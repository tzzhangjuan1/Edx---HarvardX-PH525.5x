library(BiocInstaller)

## Integrity assurances with ExpressionSet
biocLite("Biobase")
# load the genefu package again to work with breast cancer expression data
library(Biobase)
library(genefu)
data(nkis)
dim(demo.nkis)
head(demo.nkis)[,1:8]

nkes = ExpressionSet(data.nkis, phenoData=AnnotatedDataFrame(demo.nkis),
                     featureData=AnnotatedDataFrame(annot.nkis))


## Working with GEO
biocLite("GEOquery")
library(GEOquery)
# retrieve the LSC data from GEO
lstem = getGEO("GSE3725")
class(lstem) # lstem is a "list"
lstem  #read lstem, assayData: 22690 features, 28 samples 
lstem = lstem[[1]] # obtain a simple reference to the ExpressionSet instance of interest
ncol(lstem)  # 28 samples have been hybridized?

## Details of samples
pData(lstem)$title 
lstem = lstem[, -c(1:6)] # ignore the first 6 samples to generate a listing of the cell type descriptions, note position of comma!

##  Improving sample and feature labeling for a heatmap
# perform an elementary normalization
ee = exprs(lstem)
ee[ee<0] = 0
eee = log(ee+1)
# boxplot(data.frame(eee))
meds = apply(eee,2,median)
tt = t(t(eee)-meds)
# boxplot(data.frame(tt))
# assign the normalized values to ExpressionSet
exprs(lstem) = tt
# modify the feature names to be gene symbols instead of array probe names
# simplify downstream labeling with gene symbol
featureNames(lstem) = make.names(fData(lstem)$"Gene Symbol", unique=TRUE)

# reformat the naming of cell types
ct = pData(lstem)[,1]
ct = as.character(ct)
cct = gsub(".*(\\(.*\\)).*", "\\1", ct) 
cct = make.unique(cct)
cct = gsub(" enriched", "", cct)
# use the cell types as sample names
sampleNames(lstem) = cct

# select some members of the stem cell signature
inds = which(fData(lstem)$"Gene Symbol" %in% c("Stat1", "Col4a1", "Hoxa9", "Itgb5"))

# obtain a simple heatmap
heatmap(exprs(lstem[inds,]), Colv=NA)
