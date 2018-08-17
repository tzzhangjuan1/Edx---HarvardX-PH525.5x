library(BiocInstaller)
biocLite()
biocLite("GenomicRanges")
library(GenomicRanges)

# create a set of two ranges on a made-up chromosome, chrZ
gr <- GRanges("chrZ", IRanges(start=c(5,10),end=c(35,45)),
              strand="+", seqlengths=c(chrZ=100L))
gr
genome(gr) <- "hg19"      # we will say that these ranges refer to the genome hg19
gr
# use the shift function 
shift(gr, 10)
shift(gr, 80)   # when we try to shift the range beyond the length of the chromosome, the warning
#trim the ranges -- we obtain the ranges which are left, disregarding the portion that stretched beyond the length of the chromosome
trim(shift(gr, 80))
# the mcols function (stands for metadata columns) -- add columns of information to each range 
mcols(gr)
mcols(gr)$value <- c(-1,4)   # add $value to each range, assign as c(-1,4)
gr
mcols(gr)$value <- NULL     # We can remove the columns by assigning NULL
gr

## GRangesList
# use "GRangesList" to create a list of GRanges  - useful for representing groupings, for example the exons which belong to each gene
gr2 <- GRanges("chrZ",IRanges(11:13,51:53))
grl <- GRangesList(gr, gr2)
grl

# use "elementNROWS" to get the length of each GRanges use elementNROWS
length(grl)               # The length of the GRangesList is the number of GRanges object within
elementNROWS(grl)         # get the length of each GRanges
grl[[1]]                  # index into the list using typical list indexing of two square brackets
# width
width(grl)
sum(width(grl))
# add metadata columns 
mcols(grl)$value <- c(5,7)
grl
mcols(grl)

# findOverlaps and %over%
(gr1 <- GRanges("chrZ",IRanges(c(1,11,21,31,41),width=5),strand="*"))
(gr2 <- GRanges("chrZ",IRanges(c(19,33),c(38,35)),strand="*"))
# "findOverlaps" 
fo <- findOverlaps(gr1, gr2)
fo                             # returns a Hits object which contains the information about which ranges in the query (the first argument) overlapped which ranges in the subject (the second argument). 
queryHits(fo)
subjectHits(fo)
# "%over%"
gr1 %over% gr2           # returns a logical vector of which ranges in the first argument overlapped any ranges in the second.
gr1[gr1 %over% gr2]      
# Note that both of these are strand-specific, although findOverlaps has an ignore.strand option.
gr1 <- GRanges("chrZ",IRanges(1,10),strand="+")
gr2 <- GRanges("chrZ",IRanges(1,10),strand="-")
gr1 %over% gr2           # returns a logical vector "FALSE" - not same strand

## Rle and Views
# "Rle" stands for run-length encoding, which is a form of compression for repetitive data. -- Instead of storing: $[1,1,1,1]$, we would store the number 1, and the number of repeats
(r <- Rle(c(1,1,1,0,0,-2,-2,-2,rep(-1,20))))  
str(r)             # use "str" to examine the internal structure of the Rle, to show it is only storing the numeric values and the number of repeats
as.numeric(r)
# "Views" can be thought of as "windows" looking into a sequence.
(v <- Views(r, start=c(4,2), end=c(7,6))) # specify the windows of region (here, start=c(4,2), end=c(7,6)) in r
str(v)



### GRanges Assessment
## Understanding strand orientation with resize
x = GRanges("chr1", IRanges(c(1,101),c(50,150)), strand=c("+","-"))
ranges(x)
biocLite("genomicsclass/ph525x")
library(ph525x)
plotGRanges = function(x) plotRanges(ranges(x))
par(mfrow=c(2,1)) 
plotGRanges(x)
plotGRanges(resize(x,1))

## Intersecting transcripts with basic operations
# two different sets of ranges, which overlap somewhat but not entirely
x = GRanges("chr1", IRanges(c(101,201,401,501),c(150,250,450,550)), strand="+")
y = GRanges("chr1", IRanges(c(101,221,301,401,541),c(150,250,350,470,550)), strand="+")
# plot the two ranges
par(mfrow=c(2,1))
plotGRanges(x)
plotGRanges(y)
# combine the two GRanges into a GRangesList
GRangesList(x,y)
# combine the two GRanges into a single GRanges
c(x,y)
disjoined<-disjoin(c(x,y))
sum(width(x))+sum(width(y))-sum(width(disjoined))   # returns 140
#or
disjoined<-disjoin(c(x,y))
in.both<- disjoined %over% x & disjoined %over% y
sum(width(disjoined[ in.both ]))       # What is the total width which is covered by ranges in both x and y? returns 140
sum(width(disjoined[ !in.both ]))      # What is the total width which is in x or y but not in both? return 130

# Define a new genomic range, 'z', which covers range(ranges(x)) but has the opposite strand
z = GRanges("chr1", range(ranges(x)), strand="-")
sum(x %over% z)
