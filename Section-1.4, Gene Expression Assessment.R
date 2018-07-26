library(BiocInstaller)
# use devtools to install the 'tissuesGeneExpression' data 
biocLite("devtools")
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
# load the data
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e[,1:5])
table(tissue)

# Once the tissuesGeneExpression package is loaded, and data(tissuesGeneExpression) is run, you have object e, tab, and tissue in your workspace.  
library(SummarizedExperiment)
tissSE = SummarizedExperiment(list(rma=e))
colData(tissSE) = DataFrame(tab)


assay(tissSE)   # When we want the numerical values for all features and samples in a SummarizedExperiment X, we use assay(X)
# Look at the data for the feature with ID "209169_at"
mean(assay(tissSE["209169_at",])) 

## Comparing genes for tissue-specificity

IDs = c("201884_at", "209169_at", "206269_at", "207437_at", "219832_s_at", "212827_at")   # a vector of 6 IDs which index features of 'tissSE':
# Which of the following ID(s) appear to represent a gene specific to placenta?
biocLite("rafalib")
library(rafalib) 
mypar(3,2)
sapply(IDs,function(x){boxplot(e[x,]~tissue,las=2,main=x)})

## Discovery of microarray annotation in Bioconductor

## Oligo sequences on affymetrix arrays
library(BiocInstaller)
biocLite("hgu133aprobe")
library(hgu133aprobe)
head(hgu133aprobe)
# How many oligos are used to interrogate samples for gene GCM1, annotated to probe 206269_at? 
subset(hgu133aprobe, Probe.Set.Name=="206269_at")
sum(hgu133aprobe$Probe.Set.Name=="206269_at")

# a quick illustration of annotation enhancement of a SummarizedExperiment
biocLite("hgu133a.db")
library(hgu133a.db)
sym = mapIds(hgu133a.db, keys=rownames(tissSE), column="SYMBOL", keytype="PROBEID")
nm = mapIds(hgu133a.db, keys=rownames(tissSE), column="GENENAME", keytype="PROBEID")
rowData(tissSE) = DataFrame(symbol=sym, genename=nm)
# read the rowData(tissSE) dataframe
head(rowData(tissSE))
# To restrict attention to genes with 'phosphatase' in their names
tissSE[ grep("phosphatase", rowData(tissSE)$genename), ]    # 357
#How many features are annotated to genes with 'kinase' in their name? 
tissSE[ grep("kinase", rowData(tissSE)$genename), ]  # the answer is 1064



