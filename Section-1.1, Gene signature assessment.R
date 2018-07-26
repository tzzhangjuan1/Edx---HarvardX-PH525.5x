library(BiocInstaller)
biocLite()

biocLite("genefu")       # Install the genefu package with biocLite.
library(genefu)
data(sig.gene70)      # read sig.gene70 data.frame - information on the 70 gene signature.
dim(sig.gene70)
head(sig.gene70)[,1:6]   #diverse ways of describing the "genes" in the signature
#How many components of the signature have a missing value for the associated NCBI gene symbol? 
NAsub <-subset(sig.gene70, NCBI.gene.symbol== NA)
dim(NAsub)
sum(is.na(sig.gene70$NCBI.gene.symbol))  # the answer is 14

# How many of the members of the 70-gene signature are genes coding for kinases?
index<-grep("kinase",sig.gene70$Description)  # use the grep function
sig.gene70$Description[index]             # the answer is 4





