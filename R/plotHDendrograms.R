## TO DO 

### THESE PLOTTING FUNCTIONS TAKE THE BI-PRODUCT OF computeGeneScoreH.R and make some nice plots.
## Need to parameterise them so we can easily generate plots across diseases and for specific genes.



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




#lapply(foo,function(x) if(is.data.table(x)) subset(x,name=="BACH2") )

## plotting fun ;)

## here we implement for a given gene a decision tree and then spit it out

#geneName<-'CD4'
geneName<- sample(h$name,1)
#geneName<-'WDR59'
#geneName<-'AHR'
geneName<-'CD101'

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
  return(n$ppi>=thresh)
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
plot(gnp, show.tip.label = TRUE, type = "cladogram",edge.color="grey",edge.width=0.5,plot=FALSE,main=paste(toupper(disease),geneName,sep=':'))

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
  
  cw<-node$cWeight
  w<-paste0('(',node$weight,')')
  if(!node$isLeaf)
    return(paste0(node$name,' ',ifelse(is.na(cw),0,cw),w))
  return(paste0(node$name,' ',ifelse(is.na(cw),0,cw)))
}

cumWeight<-function(node){
  w<-sum(sum(node$Get('weight')))
  if(w==0)
    w<-return(NA)
  #return(log(w))
  return(w)
}

lineWidthFactor<-function(node){
  w<-sum(sum(node$Get('weight')))
  if(w==0)
    return(0)
  if(w==1)
    w<-2
  ceiling(log(w))*1.5
}

wn$Do(function(node) node$cWeight = cumWeight(node))
wn$Do(function(node) node$lw = lineWidthFactor(node))



#tw<-pmin(wn$Get(cumWeight), mean(sapply(wn$noncoding$interaction$children,cumWeight)))
tw<-wn$Get(lineWidthFactor)
rbPal <- colorRampPalette(c('green','red'))
tw<-split(rbPal(10)[as.numeric(cut(tw,breaks = 10))],names(tw))


par(mar=c(1,1,1,1))
plot(wnp, show.tip.label = TRUE, type = "cladogram",edge.color="grey",edge.width=0.5,plot=FALSE,main=toupper(disease))
lw.scale.factor<-20/nrow(h)
cex.scale.factor<-2/nrow(h)
for(node in Traverse(wn)){
  pw<-sum(sum(node$Get('weight')))
  if(!node$isRoot)
    #edges(GetPhyloNr(node$parent, "node"), GetPhyloNr(node, "node"), lwd=pw*lw.scale.factor,col = tw[[node$name]])
    if(node$lw==0){
      edges(GetPhyloNr(node$parent, "node"), GetPhyloNr(node, "node"), lwd=0.5,col = 'grey')
    }else{
      edges(GetPhyloNr(node$parent, "node"), GetPhyloNr(node, "node"), lwd=node$lw,col = tw[[node$name]])
    }
  if(!node$isLeaf & node$lw!=0)
    nodelabels(Nodelabel2(node), GetPhyloNr(node, "node"), frame = 'none', adj = c(-0.3, 0.1),col = "black",cex=0.7)
  if(node$isLeaf){
    if(node$lw!=0){
      tiplabels(Nodelabel2(node), GetPhyloNr(node, "node"), frame = "none", adj = c(0, 0),col = tw[[node$name]],cex=0.8)
    }else{
      tiplabels(Nodelabel2(node), GetPhyloNr(node, "node"), frame = "none", adj = c(0, 0),col = 'grey',cex=0.8)
    }
  }
}