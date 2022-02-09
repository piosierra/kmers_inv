#!/usr/bin/env Rscript

library(data.table)
library(Biostrings)
library(stringi)
library(intervals)

# args = commandArgs(trailingOnly=TRUE)

args <- c("~/rds/rds-durbin-group-8b3VcZwY7rY/projects/funestus/graphs/output_all2RL_s100000_p90_g30004000_k30/all2RL.fa.gz.a944863.f043790.61cf8a6.smooth.fix.gfa", "DA-402_03-2RL", "DA-416_04-2RL", "2RL")

cmd=paste("cat ",args[1], "| grep '^S' > segments",sep='')
system(cmd)

cmd=paste("cat ",args[1], " | grep \"^P\" | grep ",args[2], " | tr -s [:space:][,] \\\\n  > path_1" ,sep='')
system(cmd)
cmd=paste("cat ",args[1], " | grep \"^P\" | grep ",args[3], " | tr -s [:space:][,] \\\\n  > path_2" ,sep='')
system(cmd)



t1 <- fread("path_1", skip =1)

t2 <- fread("path_2", skip =1)

nrow(t1)
nrow(t2)

# removes end node (*)

t1<- t1[-nrow(t1),]
t2<- t2[-nrow(t2),]

colnames(t1)[1] <- "seg"
colnames(t2)[1] <- "seg"

# Add path pos variable
t1$pos <- as.numeric(rownames(t1))
t2$pos2 <- as.numeric(rownames(t2))

## Function to change sign of node
change_sign <- function(seg) {
  ifelse (stri_sub(seg,-1,-1)== "+",
          paste(stri_sub(seg,1,-2),"-",sep=""),
          paste(stri_sub(seg,1,-2),"+",sep=""))
}

# Get list of inverted nodes

t1[,inv:=change_sign(seg)]

# get list of nodes of t1 inverted in t2
t1_inv <- t1[t1$inv %in% t2$seg, ]

# order by position on the path
t1_inv <- t1_inv[order(pos)]

# Merge those that are consecutive (or chose another treshold)
t1_clus <- clusters(t1_inv$pos,1)

# Create table of intervals from the agrupated list of positions
t1_coord <- data.table(start=sapply(t1_clus, head, 1), end=sapply(t1_clus, tail, 1))


s <- fread("segments")

s[,l:=nchar(s$V3)]


t1[,sega:=as.numeric(stri_sub(seg,1,-2))]

t2[,sega:=as.numeric(stri_sub(seg,1,-2))]

#  get length of segment on t1
setkey(s,V2)
setkey(t1,sega)
setkey(t2,sega)
t1[s,len:=l]
t2[s,len2:=l]

setkey(t1,pos)
setkey(t2,pos2)
t1$end <- cumsum(t1$len)
t1$start <- 1+t1$end-t1$len
t2$end2 <- cumsum(t2$len2)
t2$start2 <- 1+t2$end2-t2$len2

setkey(t1_coord,start)
t1_coord[t1,s:=start]


setkey(t1_coord,end)
t1_coord[t1,e:=end]
t1_coord[,l:=1+e-s]



# Loop to merge inversions closer than x% of the size
min_inv_per <- 0.8
t1_merged <- data.table(start=integer(), end=integer())
current <- c(t1_coord[1]$s,t1_coord[1]$e)
current_l <- t1_coord[1]$l
current_i <- t1_coord[1]$l
for (i in 2:nrow(t1_coord)) {
  if ((current_i+t1_coord[i]$l)/(current_l+t1_coord[i]$e-current[2]) < min_inv_per ) {
    t1_merged <- rbindlist(list(t1_merged, data.table(start=current[1],end=current[2])))
    current <- c(t1_coord[i]$s,t1_coord[i]$e)
    current_l <- t1_coord[i]$l
    current_i <- t1_coord[i]$l  
  } 
  else {
    current_l <- current_l+t1_coord[i]$e-current[2]
    current_i <- current_i + t1_coord[i]$l  
    current[2] <- t1_coord[i]$e
  }
 
}
t1_merged <- rbindlist(list(t1_merged, data.table(start=current[1],end=current[2])))

t1_merged[,l:=1+end-start]

# Filter by min inversion size

min_inv_size <- 50
t1_merged <- t1_merged[l>min_inv_size]

setkey(t1_merged,start)
setkey(t1_coord,s)
t1_merged[t1_coord,pos_s:=start]
setkey(t1_merged,end)
setkey(t1_coord,e)
t1_merged[t1_coord,pos_e:=end]
setkey(t1,pos)
setkey(t1_merged,pos_s)
t1_merged[t1,seg_s:=inv]
setkey(t1_merged,pos_e)
t1_merged[t1,seg_e:=inv]
setkey(t2,seg)
setkey(t1_merged,seg_s)
t1_merged[t2,end2:= end2]
setkey(t1_merged,seg_e)
t1_merged[t2,start2:= start2]

t1_merged[,l2:=1+end2-start2]

b1 <- data.table(chr=args[4], start = t1_merged$start-1, end = t1_merged$end)
b1[,name:=paste("inv_",rownames(b1),sep="")]

b2<- data.table(chr=args[4], start = t1_merged$start2-1, end = t1_merged$end2)
b2[,name:=paste("inv_",rownames(b2),sep="")]

fwrite(b1, paste("inversions_",args[2],sep=""), sep="\t", col.names=FALSE )
fwrite(b2, paste("inversions_",args[3],sep=""), sep="\t", col.names=FALSE )