library(BiocInstaller)
biocLite()
biocLite("IRanges")
library(IRanges)

ir <- IRanges(5,10)
ir
start(ir)
end(ir)
width(ir)

ir<-IRanges(start=c(3,5,17), end=c(10,8,20))
ir
length(ir)
start(ir)

shift(ir,-2)
ir<-IRanges(5,10)
narrow(ir, start=2)
narrow(ir, end=5)
flank(ir, width=3, start=TRUE, both=FALSE)

plotir<- function(ir,i){arrows(start(ir)-.5,i,end(ir)+.5,i,code=3,angle=90,lwd=3)}
plot(0,0,xlim=c(0,15),ylim=c(0,8),type="n",xlab="",ylab="",xaxt="n")
axis(1,0:15)
abline(v=0:30 + .5, col=rgb(0,0,0, .5))

plotir(ir,1)
ir
polygon(c(start(ir)-.5, start(ir)-.5, end(ir)+.5, end(ir)+.5), c(-1,9,9,-1),col=rgb(1,0,0,.2),border=)
plotir(shift(ir,-2),2)
plotir(narrow(ir, start=2),3)
plotir(narrow(ir, end=5),3)
plotir(narrow(ir, end=5),4)
plotir(flank(ir, width=3, start=TRUE, both=FALSE),5)
plotir(flank(ir, width=3, start=FALSE, both=FALSE),6)
plotir(flank(ir, width=3, start=TRUE, both=TRUE),7)

ir<-IRanges(start=c(3,5,17), end=c(10,8,20))
ir
range(ir)
reduce(ir)
gaps(ir)
disjoin(ir)

++++++++++++++++++++++++++++++++++++++++++++++++++++++
### IRanges assessment

# Operating on ranges: zooming with *
ir <- IRanges(101,200)
ir *2   
# Narrowing
narrow(ir, start=20)
# Expanding with +
ir+25
# Range widths, vectorized
ir<-IRanges(start=c(1,11,21), end=c(3,15,27))
sum(width(ir))
# Visualizing and projecting ranges
ir<-IRanges(start=c(101,106,201,211,221,301,306,311,351,361,401,411,501), end=c(150,160,210,270,225,310,310,330,390,380,415,470,510))
biocLite("genomicsclass/ph525x")
library(ph525x)
plotRanges(ir)                  # plot ir
sum(width(gaps(reduce(ir))))    # What is the total width from 101 to 510 which is not covered by ranges in x?

# Disjoint subranges
disjoin(ir)

# Resizing
par(mfrow=c(2,1))
plotRanges(x, xlim=c(0,600))
plotRanges(resize(ir,1), xlim=c(0,600)) 


