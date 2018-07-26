library(BiocInstaller)
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome.Hsapiens.UCSC.hg19

chr11seq <- BSgenome.Hsapiens.UCSC.hg19[["chr11"]]   # access chromosome 11
subseq(chr11seq,start=10^6,width=25)   #a segment of 25 bases starting  at base 1 million

## Frequencies of short sequences
# which of the following sequences is most common on chromosome 11: "ATG", "TGA", "TAA", and "TAG"
#the fuction countPattern / matchPattern
countPattern("ATG",subseq(chr11seq))
countPattern("TGA",subseq(chr11seq))
countPattern("TAA",subseq(chr11seq))
countPattern("TAG",subseq(chr11seq))
## Nucleotide frequencies
# use the alphabetFrequency function to determine what percent of chromosome 7 is T,C,G,A, and other letters
chr7seq <- BSgenome.Hsapiens.UCSC.hg19[["chr7"]]
alphabetFrequency(subseq(chr7seq),baseOnly=TRUE, as.prob=TRUE)  

## Locations of SNPs in humans
# Download and install the SNPlocs.Hsapiens.dbSNP144.GRCh37 package
if (!("SNPlocs.Hsapiens.dbSNP144.GRCh37" %in% rownames(installed.packages()))) {
        library(BiocInstaller)
        biocLite("SNPlocs.Hsapiens.dbSNP144.GRCh37")
}
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
# To see all the SNPs on, for example, chromosome 17 
snps144 = SNPlocs.Hsapiens.dbSNP144.GRCh37
s17 = snpsBySeqname(snps144, "17")
head(s17)
# What is the location on chr17 of SNP rs73971683?
subset(s17,RefSNP_id=="rs73971683")    # the answer is 135246

##  GWAS: Linking SNP genotypes to disease risk
#Install the gwascat package
#check the version of the GWAS catalog stored in GRCh37 (hg19) coordinates.
biocLite("gwascat")
library(gwascat)
data(ebicat37)
ebicat37
# sort the chromosomes based on the number of  GWAS 'verified hits' 
sort(table(ebicat37$CHR_ID),decreasing=TRUE)   # chr. 6 has the most GWAS hits in the catalog

#What is the disease/trait with the most associations?
#use the notation mcols(ebicat37)[,"DISEASE/TRAIT"] to get a vector of names of diseases with genetic associations recorded in the gwascat
sort(table(mcols(ebicat37)[,"DISEASE/TRAIT"]),decreasing=F)   # the answer is Obesity-related traits 


