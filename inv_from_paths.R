#!/usr/bin/env Rscript

library(data.table)
library(Biostrings)
library(stringi)
library(intervals)

# args = commandArgs(trailingOnly=TRUE)

setwd("~/rds/rds-durbin-group-8b3VcZwY7rY/projects/funestus/graphs/output_all2RL_s100000_p90_g30004000_k30/")

args <- c("path_402","path_416")


t1 <- fread(args[1], skip =1)

t2 <- fread(args[2], skip =1)

nrow(t1)
nrow(t2)

# removes end node (*)

t1<- t1[-nrow(t1),]
t2<- t2[-nrow(t2),]

colnames(t1)[1] <- "seg"
colnames(t2)[1] <- "seg"

# Add path pos variable
# t1$pos <- paste(rownames(t1), "-1")
# t2$pos <- paste(rownames(t2), "-2")
t1$pos <- as.numeric(rownames(t1))
t2$pos <- as.numeric(rownames(t2))

## Function to change sign of node
change_sign <- function(seg) {
  ifelse (stri_sub(seg,-1,-1)== "+",
          paste(stri_sub(seg,1,-2),"-",sep=""),
          paste(stri_sub(seg,1,-2),"+",sep=""))
}

# Get list of inverted nodes

t1[,inv:=change_sign(seg)]
t1_inv <- t1[t1$inv %in% t2$seg, ]

t1_inv <- t1_inv[order(pos)]

t1_clus <- clusters(t1_inv$pos,1)

t1_coord <- data.table(start=sapply(t1_clus, head, 1), end=sapply(t1_clus, tail, 1))


s <- fread("segments")

s[,l:=nchar(s$V3)]


t1[,sega:=as.numeric(stri_sub(seg,1,-2))]

#  get length of segment on t1
setkey(t1,sega)
setkey(s,V2)
t1[s,len:=l]

setkey(t1,pos)
t1$end <- cumsum(t1$len)
t1$start <- 1+t1$end-t1$len

setkey(t1_coord,start)
setkey(t1,pos)
t1_coord[t1,s:=start]

setkey(t1_coord,end)
t1_coord[t1,e:=end]
t1_coord[,l:=e-s]

######  
OLD

# remove duplicated elements from them ??? (bad for pos later?)
t1 <- t1[!duplicated(t1),]
t2 <- t2[!duplicated(t2),]

nrow(t1)
nrow(t2)

# Add path variable
t1$path  <- colnames(t1)[1] 
t2$path  <- colnames(t2)[1] 

# Add path pos variable
# t1$pos <- paste(rownames(t1), "-1")
# t2$pos <- paste(rownames(t2), "-2")
t1$pos <- as.numeric(rownames(t1))
t2$pos <- as.numeric(rownames(t2))


#Change columns name and paste
colnames(t1)[1] <- "seg"
colnames(t2)[1] <- "seg"
tt <- rbind(t1,t2)

nrow(tt)

#remove duplicates
tt <- tt[!duplicated(tt$seg),]

nrow(tt)

#Get sega
tt[,sega:=stri_sub(seg,1,-2)]

#Order by path so we have later just the segments of one path.
tt[order(path)]

#remove non duplicate seg absolute
tt <- tt[duplicated(tt$sega),]

nrow(tt)

#get only the pos, regardless of path
# tt[,posn:=as.numeric(stri_sub(pos,1,-4))]

tt[order(pos)]

## These are the inversions.

s <- fread("segments")


#############################################


t1[,sense1:=!(stri_sub(`DA-402_03-2RL`,-1)=="+")]

t1[,pos:=stri_sub(`DA-402_03-2RL`,1,-2)]

t2[,sense2:=(stri_sub(`DA-416_04-2RL`,-1)=="+")]

t2[,pos:=stri_sub(`DA-416_04-2RL`,1,-2)]
setkey(t1,pos)
setkey(t2,pos)

tt <- merge(t1,t2, sort=FALSE, allow.cartesian=TRUE)

tm <- tt[sense1 != sense2]
tm <- tm[!duplicated(tm),]




t1[,sega:=stri_sub(seg,1,-2)]
t2[,sega:=stri_sub(seg,1,-2)]
t1[,sense:=stri_sub(seg,-1,-1)== "+"]
t2[,sense:=stri_sub(seg,-1,-1)]

change_sign <- function(seg) {
  ifelse (stri_sub(seg,-1,-1)== "+",
              paste(stri_sub(seg,1,-2),"-",sep=""),
              paste(stri_sub(seg,1,-2),"+",sep=""))
}

t1[,inv:=change_sign(seg)]
t1_inv <- t1[t1$inv %in% t2$seg, ]
