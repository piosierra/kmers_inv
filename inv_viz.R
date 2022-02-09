#!/usr/bin/env Rscript

library(data.table)
library(Biostrings)

# Plots the exit from kinv.sh

args = commandArgs(trailingOnly=TRUE)


t1 <- fread(args[1])
colnames(t1) <- c("flag1", "pos1", "seq")
t2 <- fread(args[2])
colnames(t2) <- c("flag2", "pos2", "seq")
t <- merge(t1,t2)

t2$rc <- sapply(t2$seq, function(x) as.character(reverseComplement(DNAString(x))))
t2$seq <- t2$rc
t2 <- t2[,-("rc")]
tt <- merge(t1,t2)
if (nrow(t)>0 & nrow(tt)>0) {
  t<- rbind(t,tt)
} else if (nrow(t)==0) {
  t <- tt
} 

png(paste("plot_",args[3], "_", args[4] ,".png",sep=""),width=2000, height=2000)
plot(t$pos1,t$pos2,pch = 20, main= paste(args[3], "_", args[4],sep=""), 
  xlab = args[3], ylab= args[4],   cex.lab=2, cex.main=3, cex.axis=2)
dev.off()

# png(paste("x_plot_",args[3], "_", args[4] ,".png",sep=""),width=2000, height=2000)
# plot(t$pos1,t$pos2,pch = 20, main= paste(args[3], "_", args[4],sep=""), 
#      xlab = args[3], ylab= args[4],xlim=c(25000000,40000000), ylim=c(25000000,40000000), 
#      cex.lab=2, cex.main=3, cex.axis=2)
# dev.off()