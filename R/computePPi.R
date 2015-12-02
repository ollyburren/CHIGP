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

## minor allele threshold 
#MAF.thresh<-0.05
## hwe threshold
HWE.thresh<-25
## call rate threshold
CALLRATE.thresh<-0.95

options(warn=1)

## FUNCTIONS ##

## process arguments
getArgs <- function(verbose=FALSE, defaults=NULL, numeric=NULL) {
  myargs <- gsub("^--","",commandArgs(TRUE))
  setopts <- !grepl("=",myargs)
  if(any(setopts))
    myargs[setopts] <- paste(myargs[setopts],"=notset",sep="")
  myargs.list <- strsplit(myargs,"=")
  myargs <- lapply(myargs.list,"[[",2 )
  names(myargs) <- lapply(myargs.list, "[[", 1)

  ## logicals
  if(any(setopts))
    myargs[setopts] <- TRUE

  ## defaults
  if(!is.null(defaults)) {
    defs.needed <- setdiff(names(defaults), names(myargs))
    if(length(defs.needed)) {
      myargs[ defs.needed ] <- defaults[ defs.needed ]
    }
  }

  ## numerics
  if(!is.null(numeric)) {
    numeric <- intersect(numeric, names(myargs))
    if(length(numeric))
      myargs[numeric] <- lapply(myargs[numeric], as.numeric)
  }

  ## verbage
  if(verbose) {
    cat("read",length(myargs),"named args:\n")
    print(myargs)
  }
  myargs
}

## required parameters

# args<-list(
# 	region_file=paste0(GRPATH,"cd4chic/gwas_integration/support/0.1cM_shuffled_regions/0.1cMcs"),
# 	out_dir=paste0(CHICDATA,"/GWAS/CHIGP/PPI/SLE/"),
# 	gwas_tbx=paste0(CHICDATA,"/GWAS/QC_TABIX/SLE.bed.gz"),
# 	gwas_type="CC",
# 	n_samples=10995,
# 	prop_cases=0.37,
# 	kg_compress_dir=paste0(STATSPATH,"1KGenome/VCF/EUR/by.chr.phase3/ALL."),
# 	kg_compress_suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz",
# 	tabix_bin="/usr/local/bin/tabix",
# 	pi_i=0.0001
# )

numerics<-c('n_samples','prop_cases','pi_i')
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

## computePPI

computePPI<-function(r,gwas_tbx,verbose=TRUE){
  if(verbose)
    message(paste("Processing region",r))
  ## get 1kg snps for this region
  smt<-vcf2sm(r)
  if(any(is.na(smt)))
    return(NA)
  sm.gt<-col.summary(smt$gt)
  #idx<-with(sm.gt,which(MAF< MAF.thresh | z.HWE^2>HWE.thresh | Call.rate<CALLRATE.thresh))
  ## remove those variants that are med indels (i.e. > 8)
  idx<-which(nchar(smt$map$ALT)>=8 | nchar(smt$map$REF) >=8)
  idx<-c(idx,which(duplicated(smt$map$POS)))
  if(length(idx)>0){
    sm<-list(map=smt$map[-idx,],gt=smt$gt[,-idx])
  }else{
    sm<-list(map=smt$map,gt=smt$gt)
  }
  ## read in gwas data using tabix
  tabix_pipe<-pipe(paste(tabix_bin, gwas_tbx, r))
  tmp.bed<-tryCatch(read.delim(tabix_pipe,sep="\t",header=FALSE,stringsAsFactors=FALSE),
                    error=function(e){
                      print(paste(r,'no variants found'))
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
  ret
}


## MAIN ##

regions<-scan(file=region_file,character())
out.file<-paste0(out_dir,basename(region_file),'.imp')
results<-lapply(regions,computePPI,gwas_tbx=gwas_tbx)
results<-results[sapply(results,is.data.table)]
results<-rbindlist(results)
## remove those regions for which we have no coverage
no.res<-which(is.na(results[,1]))
print(paste(length(no.res)," regions are not covered"))
if(length(no.res)>0)
	results<-results[-no.res,]
write.table(results,file=out.file,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
print(paste("Written",out.file))

