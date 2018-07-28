##  Using a SummarizedExperiment
# acquire the tissuesGeneExpression package using biocLite("genomicsclass/tissuesGeneExpression")
library(BiocInstaller)
biocLite("tissuesGeneExpression")
data(tissuesGeneExpression)
biocLite("SummarizedExperiment")
library(SummarizedExperiment)
tissSE = SummarizedExperiment(list(rma=e))
colData(tissSE) = DataFrame(tab)
library(hgu133a.db)
sym = mapIds(hgu133a.db, keys=rownames(tissSE), column="SYMBOL", keytype="PROBEID")
nm = mapIds(hgu133a.db, keys=rownames(tissSE), column="GENENAME", keytype="PROBEID")
rowData(tissSE) = DataFrame(symbol=sym, genename=nm)
# count the number of array features that measure expression of gene GAPDH / H2AFX.
grep("GAPDH", rowData(tissSE)$symbol, value=TRUE)
grep("H2AFX", rowData(tissSE)$symbol, value=TRUE)

## Comparing expression distributions
# Verify that 205436_s_at is the affymetrix code for H2AFX and then consider the following plot

par(las=2, mar=c(10,4,2,2))
boxplot(as.numeric(assay(tissSE["205436_s_at",]))~tissSE$Tissue)