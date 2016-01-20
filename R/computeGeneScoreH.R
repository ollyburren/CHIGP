## generateResourceFiles.R

## REQUIREMENTS


## FUNCTION: This script computes basic gene scores by integrating functional and chromatin conformational data
## for a given GWAS

## EXAMPLE USAGE:

# Rscript generateResourceFiles.R --prefix="test_"
# --cSNP_file="../DATA/" \\ 
# --interaction_file="../DATA/chic/misfud_et_al.pm.chr22.tab" --pchic.thresh=5 \\
# --res_frag_file='../DATA/support/Digest_Human_HindIII_chr22.bed'
# --region_bed_file='../DATA/support/0.1cM.regions.b37_chr22.bed' \\ 
# out_dir='../DATA/RDATA/'

library(GenomicRanges)
library(data.table)
library(reshape2)
library(yaml)
library(data.tree)

GRPATH<-'/Users/oliver/gitr/'

## Environmental variable setup for GIT repository location
if(!interactive())
  GRPATH<-Sys.getenv("GRPATH")
script.dir <- file.path(GRPATH,'CHIGP/R')
data.dir <- file.path(GRPATH,'CHIGP/DATA')


source(file.path(script.dir,'common.R'))


all.thresh<-0.5
BF.thresh<-3

args<-list(
  pmi_file = file.path(data.dir,'out/0.1cM_chr22.pmi'),
  out_file = file.path(data.dir,'out/0.1cM_chr22.geneScores.sets.tab'),
  csnps = file.path(data.dir,'RDATA/test_csnps.by.ld.RData'),
  int = file.path(data.dir,'RDATA/test_interactions.RData'),
  frags = file.path(data.dir,'RDATA/test_frags.by.ld.RData'),
  ## user defined sets
  sets = file.path(data.dir,'support/test_set.yaml'),
  ## if this is set to true then only coding SNPs in target gene are 
  ## included. If set to false then coding SNPs in interactions in genes
  ## other than target are counted.
  target.gene.cSNPs.only=TRUE,
  ## if set to true then we ignore coding and promoter SNPs when calculating tissue scores
  ## the special all category does include these.
  include.interactions.only = TRUE,
  ## if set to true then the this means that scores for individual tissues within a set are calculated
  ## probably desirable 
  decompose = TRUE
)

args<-list(
pmi_file = "/Users/oliver/DATA/JAVIERRE_GWAS/out/pmi/t1d.pmi",
out_file = "/Users/oliver/DATA/JAVIERRE_GWAS/out/hierarchical_geneScore//mchc_common.pmi.tab",
int = "/Users/oliver/DATA/JAVIERRE_GWAS/RDATA//javierre_interactions.RData",
frags = "/Users/oliver/DATA/JAVIERRE_GWAS/RDATA//javierre_frags.by.ld.RData",
csnps = "/Users/oliver/DATA/JAVIERRE_GWAS/RDATA//javierre_csnps.by.ld.RData",
target.gene.cSNPs.only = 1,
sets = "/Users/oliver/DATA/JAVIERRE_GWAS/support/javierre_tree.yaml",
include.interactions.only=1,
decompose=1
)

  
if(!interactive()){
  args<-getArgs(verbose=TRUE,
                numeric=c('target.gene.cSNPs.only','decompose','include.interactions.only'),
                defaults = list(include.interactions.only=0,decompose=0,sets=''))
}

# if(sum(names(args)=='sets')==0)
# 	args[['sets']]<-''
# 
# if(is.null(args[['include.interactions.only']]))
#   args[['include.interactions.only']]<- 0
# 
# if(is.null(args[['decompose']]))
#   args[['decompose']]<- 0

## FUNCTIONS 

getTreeYAML<-function(yfile){
  osList <- yaml.load_file(yfile)
  as.Node(osList)
}

GetSetsFromTree<-function(osNode){
  sets<-lapply(Traverse(osNode,filterFun=isNotLeaf),function(x) c(x$Get("name",filterFun=isLeaf),use.names = FALSE))
  names(sets)<-names(osNode$Get('levelName', filterFun = isNotLeaf))
  sets
}

## MAIN

ints<-get(load(args[['int']]))
## grab tissue names
tmp.tnames<-names(ints)[16:length(names(ints))]
frags.gr<-get(load(args[['frags']]))
cs.gr<-get(load(args[['csnps']]))

options(stringsAsFactors=FALSE)

pmi.file=args[['pmi_file']]
disease<-sub("\\.pmi$","",basename(pmi.file))

test.pmi<-fread(pmi.file,header=TRUE)
#setnames(test.pmi,c('chr','start','end','rs','r2','i.pos','maf','pval','ppi','delta'))
pmi.gr<-with(test.pmi,GRanges(seqnames=Rle(chr),ranges=IRanges(start=end,end=end),ppi=ppi,rs=rsid))


## always remove MHC as this buggers things up
mhc.gr<-GRanges(seqnames=Rle(6),ranges=IRanges(start=25e6,end=35e6))
ol<-as.matrix(findOverlaps(pmi.gr,mhc.gr))
if(length(ol)>0)
	pmi.gr<-pmi.gr[-ol[,1],]


cs.pmi.gr<-mergeByOverlaps(pmi.gr,cs.gr)

## add gene and frag info to pmi.gr

tmp<-cs.pmi.gr$pmi.gr
mcols(tmp)<-cbind(mcols(tmp),cs.pmi.gr[,c('ensg','ld.id','frag.id')])

## in this the coding SNPs that spoilt the party potentially are those
## where frag.id is !=0
cs.pmi.gr<-tmp

to.adj<-subset(cs.pmi.gr,frag.id!=0)

## here we want to make sure we don't double count.
to.adj$ppi<--to.adj$ppi
to.adj<-with(to.adj,data.table(ppi=ppi,ensg=ensg,ld.id=ld.id,frag.id=frag.id))
to.adj<-to.adj[,list(sppi=sum(ppi)),by="ensg,ld.id,frag.id"]

cs<-with(cs.pmi.gr,data.table(ppi=ppi,ensg=ensg,ld.id=ld.id,frag.id=frag.id))
cs<-cs[,list(sppi=sum(ppi)),by="ensg,ld.id"]


##### remove all coding SNPs this is an option
if(args[['target.gene.cSNPs.only']]){
	ol<-as.matrix(findOverlaps(pmi.gr,cs.gr))
	pmi.nocs.gr<-pmi.gr[-ol[,1],]
	pmi.frags<-mergeByOverlaps(pmi.nocs.gr,frags.gr)
}else{
	pmi.frags<-mergeByOverlaps(pmi.gr,frags.gr)
}

pmi.frags<-with(pmi.frags,data.table(rs=rs,ppi=ppi,ensg=ensg,ld.id,type=type,frag.id=id))

noncoding<-pmi.frags[,list(sppi=sum(ppi)),by="frag.id,ld.id,ensg,type"]    

noncoding.prom<-subset(noncoding,type=="promoter")
## this does not change between selections

noncoding.prom<-noncoding.prom[,list(sppi=sum(sppi)),by="frag.id,ensg,ld.id"]  

noncoding.interactions<-subset(noncoding,type=="interaction")
noncoding.interactions$sid<-with(noncoding.interactions,paste(ensg,frag.id,sep=":"))

ints$sid<-with(ints,paste(ensg,oeID,sep=":"))
setkey(noncoding.interactions,sid)

## we can generate one of these for each set of tissues

if(!file.exists(args[['sets']])){
	tissues<-split(tmp.tnames,tmp.tnames)
	#tissues[['all']]<-names(ints)[16:32]
	#tissues<-c(tissues,'all')
	## decompose switch allows us to compute geneScores for sets of tissues but also 
	## all indivdual tissues note that in this case if a set has one tissue its set name 
	## will be replaced with the tissue name so as to avoid duplication.
}else {
  setTree<-getTreeYAML(args[['sets']])
	## remove single tissue sets
	tissues<-GetSetsFromTree(setTree)
	## allow ease of selecting a sets 
	names(tissues)<-paste('set',names(tissues),sep=".")
	if((args[['decompose']]))
	  tissues<-c(tissues,split(tmp.tnames,tmp.tnames))
}

if(args[['include.interactions.only']])
  names(tissues)<-paste(names(tissues),'prey_only',sep="_")

## overall genescore
tissues[['overall']]<-tmp.tnames
## overall genescore ignoring coding SNPs
tissues[['noncoding']]<-tmp.tnames
## overall prey score ignoring promoter and coding SNPs
tissues[['interaction']]<-tmp.tnames

to.adj$uid<-with(to.adj,paste(ensg,frag.id,ld.id,sep=":"))
gint<-do.call("rbind",lapply(seq_along(tissues),function(i){
	t<-names(tissues)[i]
	print(paste("Processing",t))
	## grab gene and oeID back
	sids<-ints[which(rowSums(ints[,tissues[[i]],with=FALSE])>0),]$sid
	idx<-which(noncoding.interactions$sid %in% sids)
	#if(t!='all'){
	#	idx<-which(noncoding.interactions$sid %in% ints[ints[[t]],]$sid)
	#}else{
	#	idx<-1:nrow(noncoding.interactions)
	#}
	snoncoding.interactions<-noncoding.interactions[idx,list(sppi=sum(sppi)),by="frag.id,ensg,ld.id"]
	## allow a switch that allows us to promoter component so we can examine 
	## contribution of tissue specific interactions to the gene score
	if(t %in% c('overall','noncoding') | !args[['include.interactions.only']]){
	  
	  all.genes<-rbind(noncoding.prom,snoncoding.interactions)
	}else{
    all.genes<-snoncoding.interactions
  }
	all.genes$uid<-with(all.genes,paste(ensg,frag.id,ld.id,sep=":"))
	setcolorder(to.adj,names(all.genes))
	all.genes$uid<-with(all.genes,paste(ensg,frag.id,ld.id,sep=":"))
	if(!args[['target.gene.cSNPs.only']]){
		to.adj.m<-subset(to.adj,uid %in% all.genes$uid)
		all.genes<-rbind(all.genes,to.adj.m)
	}
	all.genes<-all.genes[,list(sppi=sum(sppi,na.rm=TRUE)),by="frag.id,ld.id,ensg"]
	all.genes<-all.genes[,list(sppi=sum(sppi,na.rm=TRUE)),by="ld.id,ensg"]
	setcolorder(cs,names(all.genes))
	## allow a switch that allows us to remove coding snp component so we can examine 
	## contribution of tissue specific interactions to the gene score
	if(t == 'overall' | !args[['include.interactions.only']]){
	  total<-rbind(cs,all.genes)
	}else{
	  total<-all.genes
	}
	ld.score<-total[,list(block_score=sum(sppi,na.rm=TRUE)),by="ensg,ld.id"]
	gs<-ld.score[,list(gene_score=1-prod(1-block_score)),by="ensg"]
	gs$tissue<-t
	gs
}))

## now add in scores for coding snps and promoter.snps

## promoters
noncoding.prom$uid<-with(noncoding.prom,paste(ensg,frag.id,ld.id,sep=":"))
if(!args[['target.gene.cSNPs.only']]){
	to.adj.m<-subset(to.adj,uid %in% noncoding.prom$uid)
	noncoding.prom.cor<-rbind(noncoding.prom,to.adj.m)
}else{
	noncoding.prom.cor<-noncoding.prom
}
noncoding.prom.cor<-noncoding.prom.cor[,list(sppi=sum(sppi,na.rm=TRUE)),by="frag.id,ensg,ld.id"]
noncoding.prom.cor<-noncoding.prom.cor[,list(block_score=sum(sppi,na.rm=TRUE)),by="ensg,ld.id"]
noncoding.prom.cor<-noncoding.prom.cor[,list(gene_score=1-prod(1-block_score)),by="ensg"]
noncoding.prom.cor$tissue<-'promoter'
fi<-rbind(gint,noncoding.prom.cor)

##coding snps
coding<-cs[,list(gene_score=1-prod(1-sppi)),by="ensg"]
coding$tissue<-'coding'
fi<-rbind(fi,coding)
foo<-melt(fi,id.vars=c("ensg","tissue"))
results<-data.table(dcast(foo,ensg~tissue+variable),key="ensg")

## add in details
details<-ints[,.(ensg,name,biotype,strand,baitChr)]
setkey(details,ensg)
details<-unique(details)
r<-results[details]

missing<-subset(r,is.na(overall_gene_score))
merged<-subset(r,!is.na(overall_gene_score))
merged$disease<-disease
setcolorder(merged,c('disease',names(details),names(results)[names(results)!="ensg"]))
#write.table(merged,file=args[['out_file']],sep="\t",row.nhames=FALSE,quote=FALSE)

### hierachical analysis - i.e. which is the best tissue or set of tissues to 
## select. Currently we use an average ppi selection proceedure

## eventually this will be hived off into it's own script as we can read merged
## in from the 

merged.f<-subset(merged,overall_gene_score >= all.thresh)

par_yaml<-'
name: overall
coding:
  name: coding
noncoding:
  promoter:
    name: promoter
  interaction:
    name: interaction
'

osList <- yaml.load(par_yaml)
n<-as.Node(osList)
n$noncoding$interaction$AddChildNode(setTree$Lymphoid)
n$noncoding$interaction$AddChildNode(setTree$Myeloid)

merged.f<-subset(merged,overall_gene_score>all.thresh)

convertNodeToDTName<-function(node){
  if(node$name %in% c('overall','coding','noncoding','promoter','interaction'))
    return(paste0(node$name,'_gene_score'))
  if(node$isLeaf)
    return(paste0(node$name,'_prey_only_gene_score'))
  return(paste0('set.',node$name,'_prey_only_gene_score'))
}

avPPi<-function(node){
  return(rowMeans(merged.f[,node$Get('dtName',filterFun = isLeaf),with=FALSE],na.rm=TRUE))
}

ppi<-function(node){
  merged.f[[node$dtName]]
}

BF<-function(node){
  if(!node$isLeaf){
    ch<-node$children
    return(ch[[1]]$ppi/ch[[2]]$ppi)
  }
  return(rep(-1,length(node$avPPi)))
}

## find the best path by scoring as to whether avPPi improves

bscore<-function(node){
  #if(!node$isLeaf & !node$isRoot){
  if(!node$isRoot){
    parent.avppi<-node$parent$avPPi
    current.avppi<-node$avPPi
    #we effectively terminate branch where avppi stops increasing
    term.branch.idx<-which(parent.avppi>=current.avppi)
    current.avppi[term.branch.idx]<-0
    return(current.avppi)
  }
  return(rep(-1,length(node$avPPi)))
}


## this allows us to easilly map between nodes and the datatable containing values
n$Do(function(node) node$dtName=convertNodeToDTName(node))
n$Do(function(node) node$ppi=ppi(node) )
n$Do(function(node) node$avPPi=avPPi(node) )
n$Do(function(node) node$score=bscore(node) )
n$Do(function(node) node$BF=BF(node) )


## only consider nodes where the ppi is above our threshdold
res.dt<-do.call("rbind",lapply(Traverse(n),function(node){
  dt <- data.table(ensg = merged.f$ensg, name= merged.f$name,score= node$score,ppi=merged.f[[node$dtName]],node=node$name,BF=node$BF)
}))

## next we
res.dt<-res.dt[,nBF:=ifelse(BF<1,1/BF,BF)]
h<-res.dt[res.dt$ppi >= all.thresh,.SD[which.max(.SD$score),],by=ensg]







#lapply(foo,function(x) if(is.data.table(x)) subset(x,name=="BACH2") )

## plotting fun ;)

## here we implement for a given gene a decision tree and then spit it out

#geneName<-'CD4'
geneName<- sample(h$name,1)

ppi.sg<-function(n,geneName){
  return(subset(merged.f,name==geneName)[[n$dtName]])
}

BF.sg<-function(n,geneName){
  if (!n$isLeaf)
    return(oBF<-subset(res.dt,name==geneName & node==n$name)$nBF)
  return(0)
}

getBranchPath<-function(node,geneName){
  node.name<-subset(h,name==geneName)$node
  path<-node$Get(function(x) return(x$pathString),traversal="post-order",filterFun = function(y) y$name==node.name)
  print(path)
  return(unlist(strsplit(path,split='/')))
}

selected<-function(n,node.list){
  return(n$name %in% node.list)
}

above.thresh<-function(n,thresh){
  if(is.na(n$ppi))
    return(FALSE)
  return(n$ppi>thresh)
}

gn<-Clone(n)
selected.branch.path<-getBranchPath(gn,geneName)
gn$Do(function(node) node$ppi = ppi.sg(node,geneName))
gn$Do(function(node) node$BF = BF.sg(node,geneName))
gn$Do(function(node) node$selected = selected(node,selected.branch.path))
gn$Do(function(node) node$above.thresh = above.thresh(node,all.thresh))

library(ape)
gn$Revert()
gnp <- as.phylo(gn)

Nodelabel <- function(node) {
  if(node$isLeaf)
    return(paste0( node$name,': PPi ', format(signif(node$ppi,digits=2), scientific = ifelse(node$ppi<0.001,TRUE,FALSE), big.mark = "'")))
  if(node$above.thresh)
    return(paste0( node$name,'\n BF ', ifelse(is.na(node$BF),'NA',format(signif(node$BF,digits=2), scientific = ifelse(node$BF>100,TRUE,FALSE), big.mark = "'"))))
  return(paste0(node$name,'\n'))
}

getFormat<-function(node){
    if(node$selected)
      return(list(col='red',lw=2,tcol='red'))
    if(node$above.thresh)
      return(list(col='black',lw=1,tcol='black'))
    return(list(col='grey',lw=1,tcol='white'))
}

par(mar=c(1,1,1,1))
plot(gnp, show.tip.label = TRUE, type = "cladogram",edge.color="grey",edge.width=0.5,plot=FALSE,main=geneName)

for(node in Traverse(gn)){
  f<-getFormat(node)
  print(f)
  if(!node$isRoot)
    edges(GetPhyloNr(node$parent, "node"), GetPhyloNr(node, "node"), lwd=f[['lw']],col = f[['col']])
    if(!node$isLeaf){
      nodelabels(Nodelabel(node), GetPhyloNr(node, "node"), frame = 'none', adj = c(0.5, 0.2),col = f[['col']],cex=0.8)
    }
  if(node$isLeaf)
    tiplabels(Nodelabel(node), GetPhyloNr(node, "node"), frame = "none", adj = c(0, 0),col=f[['col']],cex=0.8)
  
  
}

## another thing one might do is plot the pathweight for a complete disease for all selected genes

## get unique nodes

wc<-h[,list(ncount=nrow(.SD)),by="node"]
wn<-Clone(n)
wnp <- as.phylo(wn)
assignWeights<-function(n,dt){
  print(n$name)
  w<-subset(dt,node==n$name)$ncount
  if(length(w)==0)
    w<-0
  return(w)
}

wn$Do(function(node) node$weight = assignWeights(node,wc))

Nodelabel2<-function(node){
  if(node$isLeaf)
    return(paste0(node$name,' ',node$weight))
  return(paste0('\n',node$name))
}

cumWeight<-function(node){
  w<-sum(sum(node$Get('weight')))
  if(w==0)
    w<-return(NA)
  #return(log(w))
  return(w)
}



tw<-pmin(wn$Get(cumWeight),100)
rbPal <- colorRampPalette(c('green','red'))
tw<-split(rbPal(20)[as.numeric(cut(tw,breaks = 20))],names(tw))


par(mar=c(1,1,1,1))
plot(wnp, show.tip.label = TRUE, type = "cladogram",edge.color="grey",edge.width=0.5,plot=FALSE,main="T1D")
lw.scale.factor<-20/nrow(h)
cex.scale.factor<-2/nrow(h)
for(node in Traverse(wn)){
  pw<-sum(sum(node$Get('weight')))
  if(!node$isLeaf)
    nodelabels(Nodelabel2(node), GetPhyloNr(node, "node"), frame = 'none', adj = c(-0.3, 0.1),col = "black",cex=1)
  if(!node$isRoot)
    edges(GetPhyloNr(node$parent, "node"), GetPhyloNr(node, "node"), lwd=pw*lw.scale.factor,col = tw[[node$name]])
 
  if(node$isLeaf)
    tiplabels(Nodelabel2(node), GetPhyloNr(node, "node"), frame = "none", adj = c(0, 0),col="black",cex=0.8)
}


## took a long time but we probably don't need !
# branchTest <- function(node,leaf.names,level.no) {
#   if (!node$isLeaf) {
#     print(node$name)
#     ch <- node$children
#     #na <- names(merged.f)
#     #dt.node<-convertNodeToDTName(n)
#     dt.node<-node$dtName
#     #fs <- lapply(ch,convertNodeToDTName)
#     #print(fs)
#     if (length(ch) == 2) {
#       s1 <- ch[[1]]$dtName
#       s2 <- ch[[2]]$dtName
#       r <- merged.f[[s1]] / merged.f[[s2]]
#       dt <- data.table(ensg = merged.f$ensg, name=merged.f$name,BF = r,ppi=merged.f[[dt.node]])
#       dt$hyp <- node$name
#       dt[dt$BF < 1 / BF.thresh,]$hyp <- ch[[2]]$name
#       dt[dt$BF > BF.thresh,]$hyp <- ch[[1]]$name
#       dt[is.na(merged.f[[s1]]),]$hyp<-ch[[2]]$name
#       dt[is.na(merged.f[[s2]]),]$hyp<-ch[[1]]$name
#       dt[dt$BF < 1 / BF.thresh,]$ppi <- merged.f[dt$BF < 1 / BF.thresh,s2,with=FALSE]
#       dt[dt$BF > BF.thresh,]$ppi <- merged.f[dt$BF > BF.thresh,s1,with=FALSE]
#       dt[is.na(merged.f[[s1]]),]$ppi<-merged.f[is.na(merged.f[[s1]]),s2,with=FALSE]
#       dt[is.na(merged.f[[s2]]),]$ppi<-merged.f[is.na(merged.f[[s2]]),s1,with=FALSE]
#       dt$isLeaf <- dt$hyp %in% leaf.names
#       dt$level <- unlist(level.no[dt$hyp])
#       ## compute an average for each gene of leaf nodes
#       dt$avPPi<-rowMeans(merged.f[,node$Get('dtName',filterFun = isLeaf),with=FALSE],na.rm=TRUE)
#       dt$node<-node$name
#       return(dt)
#       #return(subset(dt,hyp==node$name))
#     }
#   }
# }
# 
# leaf.names<-sapply(n$leaves,'[[','name')
# level.no<-n$Get(function(no) no$level)
# foo<-(n$Get(branchTest,leaf.names,level.no,simplify="array"))
#h<-rbindlist(foo[!sapply(foo,is.logical)])
#h<-h[,nBF:=ifelse(BF<1,1/BF,BF)]
## we need to pick the top BF for each category as we end up
## adding BF to the wrong category sometimes
#hs<-subset(h,nBF>=BF.thresh & ppi>=all.thresh)
#res<-hs[,.SD[head(order((.SD$level^2*rank(.SD$BF))*.SD$ppi,decreasing = TRUE),n=1),],by=ensg]
#res<-hs[,.SD[head(order(rank(.SD$avPPi),decreasing = TRUE),n=1),],by=ensg]