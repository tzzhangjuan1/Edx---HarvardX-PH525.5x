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

## Finding Overlaps Assessment +++++++++++++++++

# Indexing individual ranges 
HepG2[17,]       # Where does the 17th HepG2 region start?

# Closest regions in distinct GRanges
#Use distanceToNearest to find the closest region in GM12878 to the 17th region in HepG2. What is the start site of this region?
d = distanceToNearest(HepG2[17],GM12878)
index= subjectHits(d)
start(GM12878[index])

# Measuring distance between closest regions
distanceToNearest(HepG2[17],GM12878)  # distance column is 2284

#Summarizing proximities of nearest regions in a pair of GRanges
d = distanceToNearest(HepG2,GM12878)
mean( mcols(d)$distance < 2000)
