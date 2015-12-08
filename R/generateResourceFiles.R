## generateResourceFiles.R

## REQUIREMENTS


## FUNCTION: This script generates the files neccesary to quickly compute gene scores 

## EXAMPLE USAGE:

# Rscript generateResourceFiles.R --prefix="test_"
# --cSNP_file="../DATA/support/cSNPs_w_ENSG.e75_chr22.bed" \\ 
# --interaction_file="../DATA/chic/misfud_et_al.pm.chr22.tab" --pchic.thresh=5 \\
# --res_frag_file='../DATA/support/Digest_Human_HindIII_chr22.bed'
# --region_bed_file='../DATA/support/0.1cM.regions.b37_chr22.bed' \\ 
# out_dir='../DATA/RDATA/'

library(data.table)
library(GenomicRanges)

## Environmental variable setup for GIT repository location
GRPATH<-Sys.getenv("GRPATH")
script.dir <- file.path(GRPATH,'CHIGP/R')
data.dir <- file.path(GRPATH,'CHIGP/DATA')

source(file.path(script.dir,'common.R'))

args<-list(
    prefix = 'chr22_',
    cSNP_file = file.path(data.dir,'support/cSNPs_w_ENSG.e75_chr22.bed'),
    interaction_file = file.path(data.dir,'chic/misfud_et_al.pm.chr22.tab'),
    pchic.thresh = 5,
    res_frag_file = file.path(data.dir,'support/Digest_Human_HindIII_chr22.bed'),
    region_bed_file = file.path(data.dir,'support/0.1cM.regions.b37_chr22.bed'),
    out_dir = file.path(data.dir,'RDATA/')
)

if(!interactive()){
  numerics<-c('pchic.thresh')
  args<-getArgs(verbose=TRUE,numeric=numerics)
}


cs<-fread(args[['cSNP_file']],header=FALSE)
setnames(cs,c('ensg','chr','start','end'))

cs$start<-cs$start+1

##load in the PCHIC data
int<-fread(args[['interaction_file']],header=TRUE)

##convert chicago scores to bools
for(n in names(int)[16:length(names(int))]){
	int[[n]]<-int[[n]]>args[['pchic.thresh']]
}

## load in h3 frags for doing promoter stuff

h3<-fread(args[['res_frag_file']],header=FALSE)
setnames(h3,c('chr','start','end','id'))
## get rid of y chromosome
h3<-subset(h3,chr!='Y')
h.gr<-with(h3,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),id=id))

oe.gr<-with(int,GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,end=oeEnd),id=oeID))
oe.gr<-oe.gr[!duplicated(oe.gr$id),]
proms.gr<-with(int,GRanges(seqnames=Rle(baitChr),ranges=IRanges(start=baitStart,end=baitEnd),ensg=ensg,id=paste(ensg,baitID,sep=":")))
proms.gr<-proms.gr[!duplicated(proms.gr$id),]

fm<-unique(subsetByOverlaps(h.gr,proms.gr)$id)

fm<-unique(c(fm-1,fm,fm+1))
## grab them back 
fm.gr<-h.gr[h.gr$id %in% fm,]

m<-mergeByOverlaps(proms.gr,h.gr,maxgap=1L)

rem.idx<-which(as.character(seqnames(m$proms.gr))!=as.character(seqnames(m$h.gr)))
if(length(rem.idx)>0)
	m<-m[-rem.idx,]
	
g2f.prom<-with(m,data.table(ensg=ensg,frag.id=id.1,orig.id=sub("ENSG[0-9]+\\:(.*)","\\1",id)))

pf.gr<-m$h.gr
pf.gr<-pf.gr[!duplicated(pf.gr$id),]

af.gr<-c(pf.gr,oe.gr)
af.gr<-af.gr[!duplicated(af.gr$id),]

## next we compute LD overlap


ld.blocks<-fread(args[['region_bed_file']],header=FALSE)
setnames(ld.blocks,c('chr','start','end','det'))
## bed files are zero based - need to be consistent throughout
ld.blocks$start<-ld.blocks$start+1
ld.gr<-with(ld.blocks,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),ld.id=1:nrow(ld.blocks)))

wld<-mergeByOverlaps(af.gr,ld.gr)

##frags overlapping boundary
bidx<-which(wld$id %in% wld[duplicated(wld$id),]$id)

sings<-wld[-bidx,]
reps<-wld[bidx,]

upper<-which(with(reps,end(af.gr)>=end(ld.gr)))
lower<-which(with(reps,start(af.gr)<=start(ld.gr)))

end(reps[upper,]$af.gr)<-end(reps[upper,]$ld.gr)
start(reps[lower,]$af.gr)<-start(reps[lower,]$ld.gr)

## put back together
wldf<-rbind(sings,reps)
frag.gr<-wldf$af.gr
frag.gr$ld.id<-wldf$ld.id

##flatten and turn into data.table

frag<-data.table(as.data.frame(frag.gr))

##add in promoter gene assignments
fragy<-merge(frag,g2f.prom,by.x="id",by.y="frag.id")
fragy$type<-'promoter'

##add in oe assignments
g2f.int<-int[,.(ensg,oeID)]
fragz<-merge(frag,g2f.int,by.x="id",by.y="oeID")
fragz$orig.id=fragz$id
fragz$type<-'interaction'
frags<-rbind(fragy,fragz)
frags$uid<-paste(frags$ensg,frags$id,frags$ld.id,sep=":")
setkey(frags,uid)
frags<-unique(frags)
frags.gr<-with(frags,GRanges(seqnames=Rle(seqnames),ranges=IRanges(start=start,end=end),ensg=ensg,ld.id=ld.id,type=type,id=id))


## now sort out cSNPs

cs<-subset(cs,chr!="Y")
cs$uid<-with(cs,paste(chr,start,sep=":"))
setkey(cs,uid)
cs.unique<-unique(cs)
cs.gr<-with(cs,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),cid=uid))
cs.gr<-with(cs.unique,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),cid=uid))

##assign LD blocks
ol<-as.matrix(findOverlaps(cs.gr,ld.gr))
csdf<-mergeByOverlaps(cs.gr,ld.gr)
csdf$cs.gr$ld.id<-csdf$ld.id
cs.gr<-csdf$cs.gr
lu<-data.table(cid=cs.gr$cid,ld.id=cs.gr$ld.id)
fin<-merge(cs,lu,by.x="uid",by.y="cid")
cs.gr<-with(fin,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),ensg=ensg,ld.id=ld.id))
cs.gr$frag.id<-integer(length=length(cs.gr))

## next we need to assign these to frag.id's in our chic data

ol<-as.matrix(findOverlaps(cs.gr,af.gr))
cs.gr[ol[,1],]$frag.id<-af.gr[ol[,2],]$id

## save required R objects
frag.file<-paste0(args[['out_dir']],'/',args[['prefix']],'frags.by.ld.RData')
save(frags.gr,file=frag.file)
message(paste('Written',frag.file))
cs.file<-paste0(args[['out_dir']],'/',args[['prefix']],'csnps.by.ld.RData')
save(cs.gr,file=cs.file)
message(paste('Written',cs.file))
int.file<-paste0(args[['out_dir']],'/',args[['prefix']],'interactions.RData')
save(int,file=int.file)
message(paste('Written',int.file))
