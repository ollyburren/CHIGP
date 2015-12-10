## computePPi.R

## REQUIREMENTS
## tabix
## snpStats
## data.table

## FUNCTION: This script computes posterior probability of a snp being causal within an LD block
## If you have a non 1k genome imputed dataset then all of this functionality is contained within
## the PMI script.

## Note that this software uses 1KG to estimate MAF so if some SNPs are not in a given 1KG set then 
## they are likely to get dropped.

## GWAS INPUT FILES
## These are expected to be in BED5 format (chr,start,stop,name,p.value) ZERO-BASED
## that are subsequently indexed with tabix e.g.
## sed -e 's/^chr//' $filename | sort -k1,1 -k2,2n | /bin/bgzip -c > $filename.bgz
## /bin/tabix -p bed $filename.bgz

## 1KG input files - these were downloaded from -- need to check with Chris 

## EXAMPLE USAGE:

## Rscript --vanilla PMI.R \\
## --region_file=./0.1cM_shuffled_regions/0.1cMcs \\
## --out_dir=/home/uname/captureHIC_support/GWAS/PMIPV/ --gwas_tbx=./GWAS/QC_TABIX/Platelet_Volume.bed.gz \\
## --gwas_type=QUANT --n_samples=18600 --prop_cases=1 --kg_compress_dir=./1KGenome/VCF/EUR/by.chr.phase3/ALL. \\
## --kg_compress_suffix=.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz \\
## --tabix_bin=/bin/tabix --pi_1=0.0001 

## It is assumed that you would use a wrapper script such as ../perl/pmi.pl to run

## CREDITS:
## getArgs() is taken from https://github.com/chr1swallace/random-functions/blob/master/R/getArgs.R
## by Chris Wallace
## vcf2sm() was developed in collaboration with Markus Klarquist
## ppi code was adapted from coloc package https://github.com/chr1swallace/coloc


library(snpStats)
library(data.table)

## GLOBAL VARIABLES
## Environmental variable setup for GIT repository location
if(!interactive())
  GRPATH<-Sys.getenv("GRPATH")
DATA_DIR<-file.path(GRPATH,'CHIGP/DATA')
script.dir <- file.path(GRPATH,'CHIGP/R')
## read in library functions
source(file.path(script.dir,'common.R'))

## default settings for filtering 1KG genomes data
MAF.thresh<-0.01 ## default 1%
## hwe threshold
HWE.thresh<-25
## call rate threshold
CALLRATE.thresh<-0.95

source(file.path(GRPATH,'CHIGP/R/common.R'))

## required parameters

args<-list(
 	region_file=file.path(DATA_DIR,"support/0.1cM_chr22"),
 	out_dir=file.path(DATA_DIR,"out/"),
 	gwas_tbx=file.path(DATA_DIR,"gwas/OKADA_RA_22.bed.gz"),
 	gwas_type="CC",
 	n_samples=58293,
 	prop_cases=0.25,
 	kg_compress_dir=file.path(DATA_DIR,'1kgenome/VCF/EUR/by.chr.phase3/ALL.'),
 	kg_compress_suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz",
 	tabix_bin="/usr/local/bin/tabix",
 	pi_i=0.0001,
 	maf_thresh=0.01,
 	do_pmi=1
 )

numerics<-c('n_samples','prop_cases','pi_i')
if(!interactive())
args<-getArgs(verbose=TRUE,numeric=numerics)

## location of tabix bed file with pvals scores in
gwas_tbx<-args[['gwas_tbx']]
## list of regions to compute
region_file<-args[['region_file']]
## output directory
out_dir<-args[['out_dir']]
## number of samples in study
n_samples<-args[['n_samples']]
## proportion that are cases
prop_cases<-args[['prop_cases']]
## gwas type one of QUANT - quantitative trait or CC - case/control
gwas_type<-args[['gwas_type']]
## location of tabix formmated 1000 genome vcf files
kg_compress_dir <- args[['kg_compress_dir']]
## these follow a particular format and we need to know what the suffix is
kg_compress_suffix <- args[['kg_compress_suffix']]
## location of tabix binary
tabix_bin <- args[['tabix_bin']]
## prior probability that a SNP is causal usually 10e-4
pi_i<-args[['pi_i']]
## r2 threshold 
r2_thresh<-args[['r2_thresh']]
## maf thresh
maf_thresh<-ifelse(sum(grepl('maf_thresh',names(args)))!=0,args[['maf_thresh']],MAF.thresh)



## robustly sum logs
logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form)
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

## quickly convert VCF format to snpMatrix format
vcf2sm<-function(region){
	chrom <- gsub("^([^:]+).*","\\1",region)
	kg_file <- paste0(kg_compress_dir,'chr',chrom,kg_compress_suffix)
	## get the header
	my.pipe<-pipe(paste(tabix_bin,'-H',kg_file,region))
	header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
	close(my.pipe)
	cnames<-unlist(strsplit(header,"\t"))
	tmp<-as.data.frame(fread(paste(tabix_bin,kg_file,region,' | grep PASS'),sep="\t",header=FALSE,stringsAsFactors=FALSE))
	colnames(tmp)<-cnames
	## remove SNPs which have more than 2 alleles as we cannot process these
	idx<-which(nchar(tmp$ALT)>=8 | nchar(tmp$REF) >=8)
	if(length(idx)>0)
		tmp<-tmp[-idx,]
	gt<-tmp[,10:ncol(tmp)]
	if(nrow(gt)==0)
	  return(NA)
	info<-tmp[,1:9]
	sm<-apply(gt,1,function(x) sub("0\\|0","1",x))
	sm<-apply(sm,1,function(x) sub("(0\\|1)|(1\\|0)","2",x))
	sm<-apply(sm,1,function(x) sub("1\\|1","3",x))
	## set anything else to a missing value
	sm<-t(apply(sm,1,function(x) as.raw(sub("[0-9]\\|[0-9]","0",x))))
	colnames(sm)<-make.names(info$ID,unique=TRUE)
	rownames(sm)<-colnames(gt)
	sm<-new("SnpMatrix", sm)
	return(list(map=info,gt=sm))
}

## compute variance shrinkage for quantitative trait study
Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

## compute variance shrinkage for case control study
Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}

## compute approx bayes factors and resultant posterior probabilities
## based on the assumption of one causal variant in a region
approx.bf.p <- function(p,f,type, N, s,pi_i,suffix=NULL) {
  if(type=="QUANT") {
    sd.prior <- 0.15
    V <- Var.data(f, N)
  } else {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- sd.prior^2 / (sd.prior^2 + V)
  ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
  lABF = 0.5 * (log(1-r) + (r * z^2))
  sBF <- logsum(lABF + log(pi_i))
  ppi<- exp(lABF + log(pi_i))/(exp(sBF) + 1)
  ret <- data.frame(V,z,r,lABF,ppi)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)  
}



get1KGInfo<-function(r,gwas_tbx,verbose=TRUE,MAF_thresh,HWE_thresh,CALLRATE_thresh){
  if(verbose)
    message(paste("Processing region",r))
  ## get 1kg snps for this region
  smt<-vcf2sm(r)
  if(any(is.na(smt)))
    return(NA)
  sm.gt<-col.summary(smt$gt)
  idx<-which(nchar(smt$map$ALT)>=8 | nchar(smt$map$REF) >=8)
  if(!missing(HWE_thresh) & !missing(MAF_thresh) & !missing(CALLRATE_thresh))
    idx<-unique(c(idx,with(sm.gt,which(MAF< MAF_thresh | z.HWE^2>HWE_thresh | Call.rate<CALLRATE_thresh))))
   
  ## remove those variants that are med indels (i.e. > 8)
  idx<-unique(c(idx,which(duplicated(smt$map$POS))))
  if(length(idx)>0){
    sm<-list(map=smt$map[-idx,],gt=smt$gt[,-idx])
  }else{
    sm<-list(map=smt$map,gt=smt$gt)
  }
  ## read in gwas data using tabix
  tabix_pipe<-pipe(paste(tabix_bin, gwas_tbx, r))

  tmp.bed<-tryCatch(read.delim(tabix_pipe,sep="\t",header=FALSE,stringsAsFactors=FALSE),
                    error=function(e){
                      message(paste(r,'no variants found'))
                      return(NA)
                    })
  if(any(is.na(tmp.bed)))
    return(NA)
  
  
  names(tmp.bed)<-c('chr','start','end','rs','pval')
  tmp.bed<-tmp.bed[!duplicated(tmp.bed$start),]
  ## zero based !
  tmp.bed$start<-tmp.bed$start+1;
  ## next add pvals to  the info based on GWAS overlap.
  sm$map<-merge(sm$map,tmp.bed[,c('start','pval')],by.x='POS',by.y='start',all.x=TRUE)
  idx<-which(!is.na(sm$map$pval))
  ## for sparser genotyping some regions only have no variants we cannot
  ## impute these so bail !
  if(length(idx)==0)
    return(NA)
  return(sm)
}

## computePPI

computePPI<-function(r,gwas_tbx,verbose=TRUE){
  sm<-get1KGInfo(r,gwas_tbx,verbose)
  if(any(is.na(sm)))
    return(sm)
  idx<-which(!is.na(sm$map$pval))
  ## for sparser genotyping some regions only have no variants we cannot
  ## impute these so bail !
  if(length(idx)==0)
    return(NA)
  ## next filter to remove those w/o pvals
  sm$map<-sm$map[idx,]
  sm$gt<-sm$gt[,idx]
  sm$map$maf<-col.summary(sm$gt)$MAF
  ## next compute ppi using implementation of Wakfield.
  sm$map$ppi<-signif(approx.bf.p(as.numeric(sm$map$pval),sm$map$maf,gwas_type,n_samples,prop_cases,pi_i)$ppi,digits = 3)
  ret<-data.table(sm$map)
  setnames(ret,c('start','chr','rsid','ref','alt','qual','filter','info','format','pval','maf','ppi'))
  ret<-ret[,.(start,chr,rsid,pval,maf,ppi)]
  ret$end<-ret$start + 1
  setcolorder(ret,c('chr','start','end','rsid','maf','pval','ppi'))
  message(paste('Got',nrow(ret)))
  ret
}


computePMI<-function(r,gwas_tbx,verbose=TRUE,MAF.thresh,HWE.thresh,CALLRATE.thresh,R2.thresh){
  sm<-get1KGInfo(r,gwas_tbx,verbose,MAF_thresh=MAF.thresh,HWE_thresh=HWE.thresh,CALLRATE_thresh=CALLRATE.thresh)
  if(any(is.na(sm)))
    return(sm)
  idx<-which(!is.na(sm$map$pval))
  ## for sparser genotyping some regions only have no variants we cannot
  ## impute these so bail !
  if(length(idx)==0)
    return(NA)
  gt<-sm$gt
  colnames(gt)<-sm$map$POS
  tld<-ld(gt[,idx],gt[,-idx],stats="R.squared")
  ## next we select the max row for each column
  fb<-matrix(NA,nrow=ncol(tld),ncol=4)
  
  fb[,1]<-as.numeric(colnames(tld)) ## missing
  if(nrow(tld)==1){
    ## in this case they must match rownames as only 1 SNP !
    fb[,2]<-as.numeric(rownames(tld))
    fb[,3]<-as.numeric(tld)
  }else{
    fb[,2]<-as.numeric(rownames(tld)[apply(tld,2,which.max)]) # found
    fb[,3]<-as.numeric(apply(tld,2,max)) ## r2
  }
  
  ## any below cut off drop on the floor.
  
  fb<-matrix(fb[fb[,3]>=R2.thresh,],ncol=4)
  ## remove those that are not computed correctly
  fb<-fb[fb[,3]!="NaN",]
  
  ## also remove those that are NaN
  
  if(!is.matrix(fb)){
    fb<-rbind(fb,rep(NA,4))
    message("Not matrix adding NA's")
  }
  fb[,4]<-sm$map[match(fb[,2],sm$map$POS),]$pval
  
  ## finally fill in pvals
  
  
  sm$map$pval<-apply(cbind(sm$map$pval,fb[match(sm$map$POS,fb[,1]),4]),1,function(x){ 
    tmp<-x[!is.na(x)]
    ifelse(length(tmp)==1,tmp,NA)
  })	
  sm$map$imp.snp.pos<-apply(cbind(sm$map$POS,fb[match(sm$map$POS,fb[,1]),2]),1,function(x){
    tmp<-x[!is.na(x)]
    ifelse(length(tmp)>1,tmp[2],NA)
  })
  
  sm$map$imp.r2<-apply(cbind(sm$map$POS,fb[match(sm$map$POS,fb[,1]),3]),1,function(x){	
    tmp<-x[!is.na(x)]
    ifelse(length(tmp)>1,tmp[2],NA)
  })
  
  ## CURRENTLY PMI METHOD EFFECTIVELY REMOVES SNPS THAT ARE IN ORIGINAL GWAS BUT DON'T 
  ## TAG ANY OTHER SNPS (WITH THRESHOLDS). THIS CAUSES DOWN SAMPLING WHICH MIGHT BE 
  ## A REQUIRED BEHAVIOUR.
  sm$map$maf<-col.summary(sm$gt)$MAF
  ##what is the percentage imputation
  kidx<-which(!is.na(sm$map$pval))
  keep<-sm$map[kidx,]
  ## message to show how well we did 
  perc<-(nrow(keep)/nrow(sm$map)) * 100
  message(paste(floor(perc),'% imputed for',nrow(keep),'/',nrow(sm$map),r))
  
  outf<-data.table(keep)
  setnames(outf,'#CHROM','chr')
  outf$ppi<-with(outf,signif(approx.bf.p(as.numeric(pval),maf,gwas_type,n_samples,prop_cases,pi_i)$ppi,digits=3))
  ## create a output file
  #f<-sm$gt[keep$ID,]$MAF
  ret<-outf[,.(POS,chr,ID,pval,maf,ppi,imp.snp.pos,imp.r2)]
  ret$start<-ret$POS-1
  setnames(ret,c('end','chr','rsid','pval','maf','ppi','imp.snp.pos','imp.r2','start'))
  setcolorder(ret,c('chr','start','end','rsid','maf','pval','ppi','imp.snp.pos','imp.r2'))
  ret
}



## MAIN ##
#regions<-c('22:1-16252600','22:16252600-17552600')
regions<-scan(file=region_file,character())
fext<-'.ppi'
if(args[['do_pmi']]==1){
  message("Running PMI")
  fext<-'.pmi'
  results<-lapply(regions,computePMI,gwas_tbx=gwas_tbx,MAF.thresh=maf_thresh,HWE.thresh=HWE.thresh,CALLRATE.thresh=CALLRATE.thresh,R2.thresh = 0.6)
}else{
  message("Running PPi")
  results<-lapply(regions,computePPI,gwas_tbx=gwas_tbx)
}
no.res<-sapply(results,is.data.table)
no.cover<-sum(!no.res)
message(paste(no.cover,ifelse(no.cover==1,"region is","regions are"),"not covered"))
if(length(no.res)>0)
  results<-results[no.res]
out.file<-paste0(out_dir,basename(region_file),fext)
results<-results[sapply(results,is.data.table)]
results<-rbindlist(results)
## remove those regions for which we have no coverage

write.table(results,file=out.file,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
message(paste("Written",out.file))

