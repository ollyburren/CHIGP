## generateResourceFiles.R

## REQUIREMENTS


## FUNCTION: This script computes support data for priortising genes in tissue contexts.

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
  pmi_file = file.path(data.dir,'out/0.1cM_chr22.ppi'),
  out_dir = file.path(data.dir,'out/hierachical_geneScore/'),
  int = file.path(data.dir,'RDATA/test_interactions.RData'),
  frags = file.path(data.dir,'RDATA/test_frags.by.ld.RData'),
  csnps = file.path(data.dir,'RDATA/test_csnps.by.ld.RData'),
  target.gene.cSNPs.only = 1,
  sets = file.path(data.dir,'support/mifsud.yaml'),
  include.interactions.only=1,
  decompose=1,
  ppi.thresh = 0.5,
  BF.thresh = 3
)


all.thresh <- args[['ppi.thresh']]
BF.thresh <- args[['BF.thresh']]
  
if(!interactive()){
  args<-getArgs(verbose=TRUE,
                numeric=c('target.gene.cSNPs.only','decompose','include.interactions.only','ppi.thresh','BF.thresh'),
                defaults = list(include.interactions.only=0,decompose=0,sets=''))
}

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

## check to see if our out dir exists if not create it

if(!dir.exists(args[['out_dir']]))
  dir.create(args[['out_dir']])

ints<-get(load(args[['int']]))
## grab tissue names
tmp.tnames<-names(ints)[16:length(names(ints))]
frags.gr<-get(load(args[['frags']]))
cs.gr<-get(load(args[['csnps']]))

options(stringsAsFactors=FALSE)

pmi.file=args[['pmi_file']]
disease<-sub("[.][^.]*$", "", basename(args[['pmi_file']]), perl=TRUE)
file.stub<-basename(args[['pmi_file']])

test.pmi<-fread(pmi.file,header=TRUE)
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
	message(paste("Processing",t))
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
## gene scores of NA mess things up set these to zero
fi[is.na(gene_score),]$gene_score<-0
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


full.results.file<-file.path(args[['out_dir']],paste0(file.stub,'_full.tab'))
write.table(merged,file=full.results.file,sep="\t",row.names=FALSE,quote=FALSE)
message(paste('written',full.results.file))
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
chil<-setTree$children
n$noncoding$interaction$AddChildNode(chil[[1]])
n$noncoding$interaction$AddChildNode(chil[[2]])
#n$noncoding$interaction$AddChildNode(setTree$Lymphoid)
#n$noncoding$interaction$AddChildNode(setTree$Myeloid)

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

wBF<-function(node){
  if(!node$isRoot){
    par<-node$parent
    ch<-par$children
    sib.name<-names(ch)[names(ch) != node$name]
    sib<-ch[[sib.name]]
    wBF<-node$ppi/sib$ppi
    # for those below threshold set to zero
    wBF[wBF<BF.thresh]<-0
    return(wBF)
  }
  return(BF.thresh)
}

cumBF<-function(node){
  if(node$isRoot)
    return(BF.thresh)
  ancwBF<-do.call("cbind",node$Get('wBF',traversal = "ancestor"))
  return(apply(ancwBF, 1, prod,na.rm=TRUE))
}




## this allows us to easilly map between nodes and the datatable containing values
n$Do(function(node) node$dtName=convertNodeToDTName(node))
n$Do(function(node) node$ppi=ppi(node) )
n$Do(function(node) node$avPPi=avPPi(node) )
n$Do(function(node) node$BF=BF(node) )
n$Do(function(node) node$wBF=wBF(node) )
n$Do(function(node) node$score=cumBF(node))
n$Do(function(node) node$nBF=ifelse(node$BF>1,node$BF,1/node$BF) )



## only consider nodes where the ppi is above our threshold
res.dt <- do.call("rbind",lapply(Traverse(n),function(node) {
  dt <-
    data.table(
      ensg = merged.f$ensg, name = merged.f$name,score = node$score,appi = node$avPPi,node =
        node$name,wBF = node$wBF,nBF = node$nBF, isLeaf = node$isLeaf
    )
}))

## next we
#res.dt<-res.dt[,nBF:=ifelse(BF<1,1/BF,BF)]
h<-res.dt[,.SD[which.max(.SD$score),],by=ensg]
h<-h[order(h$score),]
h$disease<-disease
priortised.results.file<-file.path(args[['out_dir']],paste0(file.stub,'_prioritised.tab'))
write.table(h,file=priortised.results.file,sep="\t",row.names=FALSE,quote=FALSE)
message(paste('written',priortised.results.file))






