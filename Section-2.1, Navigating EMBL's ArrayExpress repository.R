## Navigating EMBL's ArrayExpress repository

# direct interrogation of the EMBL-EBI archive, with queryAE
library(BiocInstaller)
biocLite("ArrayExpress")
library(ArrayExpress)
sets = queryAE(keywords = "glioblastoma", species = "homo+sapiens")  # queryAE function
dim(sets)
# the DT package,  if I say data table sets, I will get a webpage
biocLite("DT")
library(DT)
datatable(sets)
sets[5:7,-c(7,8)]

# acquire the raw data labled with a PubMed ID with the getAE function
initdir = dir()
if (!file.exists("E-MTAB-5797.sdrf.txt")) nano = getAE("E-MTAB-5797")
afterget = dir()
setdiff(afterget, initdir)
# 
