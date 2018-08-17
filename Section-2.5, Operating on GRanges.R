library(BiocInstaller)
biocLite()
biocLite("GenomicRanges")
library(GenomicRanges)


# A simple set of ranges
ir <- IRanges(c(3, 8, 14, 15, 19, 34, 40),
              width = c(12, 6, 6, 15, 6, 2, 7))
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...)
{
        height <- 1
        if (is(xlim, "Ranges"))
                xlim <- c(min(start(xlim)), max(end(xlim)))
        bins <- disjointBins(IRanges(start(x), end(x) + 1))
        plot.new()
        plot.window(xlim, c(0, max(bins)*(height + sep)))
        ybottom <- bins * (sep + height) - height
        rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
        title(main)
        axis(1)
}

plotGRanges = function (x, xlim = x, col = "black", sep = 0.5, xlimits = c(0, 
                                                                           60), ...) 
{
        main = deparse(substitute(x))
        ch = as.character(seqnames(x)[1])
        x = ranges(x)
        height <- 1
        if (is(xlim, "Ranges")) 
                xlim <- c(min(start(xlim)), max(end(xlim)))
        bins <- disjointBins(IRanges(start(x), end(x) + 1))
        plot.new()
        plot.window(xlim = xlimits, c(0, max(bins) * (height + sep)))
        ybottom <- bins * (sep + height) - height
        rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom + height, 
             col = col, ...)
        title(main, xlab = ch)
        axis(1)
}

# plot - visualize ir and several intra-range operations.
par(mfrow=c(4,1), mar=c(4,2,2,2))
plotRanges(ir, xlim=c(0,60))
plotRanges(reduce(ir), xlim=c(0,60))
plotRanges(disjoin(ir), xlim=c(0,60))
plotRanges(gaps(ir), xlim=c(0,60))

# Extension to GRanges
library(GenomicRanges)
gir = GRanges(seqnames="chr1", ir, strand=c(rep("+", 4), rep("-",3)))
gir
par(mfrow=c(4,1), mar=c(4,2,2,2))
plotGRanges(gir, xlim=c(0,60))
plotGRanges(resize(gir,1), xlim=c(0,60),col="green")
plotGRanges(flank(gir,3), xlim=c(0,60), col="purple")
plotGRanges(flank(gir,2,start=FALSE), xlim=c(0,60), col="brown")
