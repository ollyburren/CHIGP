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

GRPATH<-'/Volumes/stats/GIT_REPOS/'

## Environmental variable setup for GIT repository location
if(!interactive())
  GRPATH<-Sys.getenv("GRPATH")
script.dir <- file.path(GRPATH,'CHIGP/R')
data.dir <- file.path(GRPATH,'CHIGP/DATA')


source(file.path(script.dir,'common.R'))


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
pmi_file = "/Volumes/stats/GIT_REPOS//astle/DATA//out/pmi/mchc_common.pmi",
out_file = "/Volumes/stats/GIT_REPOS//astle/DATA//out/hierarchical_geneScore//mchc_common.pmi.tab",
int = "/Volumes/stats/GIT_REPOS//astle/DATA//RDATA//astle_interactions.RData",
frags = "/Volumes/stats/GIT_REPOS//astle/DATA//RDATA//astle_frags.by.ld.RData",
csnps = "/Volumes/stats/GIT_REPOS//astle/DATA//RDATA//astle_csnps.by.ld.RData",
target.gene.cSNPs.only = 1,
sets = "/Volumes/stats/GIT_REPOS//astle/DATA//support/javierre_tree.yaml",
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
  names(tissues)<-paste(names(tissues),'interactions_only',sep="_")

## overall genescore
tissues[['overall']]<-tmp.tnames
## overall genescore ignoring coding SNPs
tissues[['noncoding']]<-tmp.tnames

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
write.table(merged,file=args[['out_file']],sep="\t",row.names=FALSE,quote=FALSE)

### hierachical analysis

all.thresh<-0.1
BF.thresh<-3
merged.f<-subset(merged,overall_gene_score >= all.thresh)

par_yaml<-'
name: overall
coding:
  name: coding
noncoding:
  promoter:
    name: promoter
  haematopoiesis:
    name: haematopoiesis
'

osList <- yaml.load(par_yaml)
n<-as.Node(osList)
n$noncoding$haematopoiesis$AddChildNode(setTree$Lymphoid)
n$noncoding$haematopoiesis$AddChildNode(setTree$Myeloid)

merged<-subset(merged,overall_gene_score>0.1)

BF <- function(node) {
  
  if(!node$isLeaf){
  ch<-node$children
  print(names(ch))
  fs<-lapply(ch,function(n){
    if(n$isLeaf){
      return(paste0(n$name,'_interactions_only_gene_score'))
    }else{
      return(paste0('set.',n$name,'_interactions_only_gene_score'))
    }
  })
  if(length(fs)==2){
    s1<-fs[[1]]
    s2<-fs[[2]]
    head(merged[[s1]],n=1)/head(merged[[s2]],n=1)
  }
  }
}

ppi <-  function(node) {
  if(node$isLeaf){
    fs<-paste0(node$name,'_interactions_only_gene_score')
  }else{
    fs<-paste0('set.',node$name,'_interactions_only_gene_score')
  }
  head(merged[[fs]],n=1)
}

#print(n,ppi=ppi)

n$Do(ppi)


## is this a coding snp or a non coding SNP (interactions or promoter)
# prev<-'none'
# s1<-'coding_gene_score'
# s2<-paste0(names(tissues[1]),'_gene_score')
# rass<-function(snode,ts,s1,s2,prev,BF.thresh=3){
#   ts$hyp<-character(length=nrow(ts))
#   ts[is.na(ts[[s1]]),]$hyp<-s2
#   ts[is.na(ts[[s2]]),]$hyp<-s1
#   ## these rows contain NA and so BF will also be NA
#   na.idx<-which(rowSums(is.na(ts[,c(s1,s2),with=FALSE]))>0)
#   ts$BF<-ts[[s1]]/ts[[s2]]
#   ts[ts$BF<=1/BF.thresh,]$hyp<-s2
#   ts[ts$BF>=BF.thresh,]$hyp<-s1
#   ts[ts$hyp=="",]$hyp<-prev
#   do.call("rbind",lapply(snode$children,function(sn){
#     rass(sn,ts,snode,)
#   }
# }


message(paste("Written",args[['out_file']]))


