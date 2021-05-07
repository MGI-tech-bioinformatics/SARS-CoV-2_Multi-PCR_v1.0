library(Cairo)

argv <- commandArgs(TRUE)

infile <- argv[1]
outfile <- argv[2]

CairoSVG(width = 11.5, height = 5,file=outfile)
dat <- read.table(infile,header = F)
sp1=spline(dat$V1,dat$V2,n=9000)
plot(dat$V1,dat$V2,type='l',ylim=c(0,max(dat$V2)),xlab="Windows", ylab="Log10(Depth)",col='white',lwd=1)
lines(sp1,col='blue',lwd=1.5) 
#lines(sp1,col='red2',lwd=1.5)

dev.off()
