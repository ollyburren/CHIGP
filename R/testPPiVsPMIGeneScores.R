## a function to compare the output of PMI gene scores and non PMI gene scores
library(data.table)
library(ggplot2)

options(stringsAsFactors=FALSE)

if(!interactive())
  GRPATH<-Sys.getenv("GRPATH")
DATA_DIR<-file.path(GRPATH,'CHIGP/DATA')

files<-list.files(path=file.path(DATA_DIR,'/out/'),pattern="*.tab",full.names=TRUE)

r<-lapply(seq_along(files),function(i){
  tmp<-fread(files[i],header=TRUE)
  tmp<-tmp[,.(ensg,all_gene_score)]
  setkey(tmp,ensg)
  suffix<-unlist(strsplit(basename(files[i]),"\\."))[4]
  setnames(tmp,'all_gene_score',paste('gene_score',suffix,sep="."))
  tmp
})

all<-r[[1]][r[[2]]]



pdf(file=file.path(DATA_DIR,'/out/plots/pmi_vs_ppi.pdf'))
ggplot(all,aes(x=gene_score.ppi,y=gene_score.pmi)) + geom_point() + theme_bw()
dev.off()