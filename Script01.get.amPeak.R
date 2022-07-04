#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# The following initializes usage of Bioc devel
#BiocManager::install(version='devel')
#BiocManager::install("sangerseqR")
library(sangerseqR)
dir_NAME="S"
path=paste("e:/cattle/A2.analysis/runspace/ab1/",dir_NAME,"/",sep="")

setwd(path)
outpath=paste("e:/cattle/A2.analysis/runspace/ampeak/",dir_NAME,"/",sep="")
ab1_files=list.files(path=".",pattern="*.ab1")

#cut 0.15, 0.2
# uncut 0.15, 0.2

for (threshold in c(0.15, 0.2)){
  for (num in 1:length(ab1_files)){
    file=ab1_files[num]
    names=strsplit(file,split="_")
    prefix=names[[1]][1]
    seq <- readsangerseq(file)
    seq <- makeBaseCalls(seq,ratio=0.2)
    peakAmp <- peakAmpMatrix(seq)
    peakAmp <- as.data.frame(peakAmp)
    colnames(peakAmp) <- c("A","C","G","T")
    peakAmp$ratio <- apply(peakAmp,1,function(x){a=sort(x,decreasing=T);a[2]/a[1]})
    peakAmp$sig <- ifelse(peakAmp$ratio>threshold,T,F)
    out_file=paste(outpath,prefix,".",threshold,".peakAmp.txt",sep="")
    write.table(peakAmp, file=out_file,row.names=F,quote=F,sep="\t")
  }
}
