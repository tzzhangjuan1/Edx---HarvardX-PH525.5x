source("https://bioconductor.org/biocLite.R")
library(BiocInstaller)
biocLite()

library(IRanges)
biocLite("GenomicRanges")
library(GenomicRanges)

biocLite("GenomicFeatures")
library(GenomicFeatures)

biocLite("devtools")
library(devtools)
install_github("genomicsclass/ERBS")
library(ERBS)

data(HepG2)
data(GM12878)

res=findOverlaps(HepG2, GM12878)
class(res)
res
HepG2[1,]
GM12878[12,]

index=queryHits(res)
erbs=HepG2[index,]
erbs

erbsR=granges(erbs)
erbsR
