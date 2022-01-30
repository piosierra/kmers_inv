#!/usr/bin/env Rscript

library(data.table)
library(Biostrings)
library(stringi)

# args = commandArgs(trailingOnly=TRUE)

setwd("~/rds/rds-durbin-group-8b3VcZwY7rY/projects/funestus/graphs/output_all2RL_s100000_p90_g30004000_k30/")

args <- c("path_402","path_416")


t1 <- fread(args[1], skip =1)

t2 <- fread(args[2], skip =1)

# Keep only only element from each
t1 <- t1[!duplicated(t1),]
t2 <- t2[!duplicated(t2),]

# Add path variable
t1$path  <- colnames(t1)[1] 
t2$path  <- colnames(t2)[1] 

# Add path pos variable
t1$pos <- paste(rownames(t1), "-1")
t2$pos <- paste(rownames(t2), "-2")

#Change columns name and paste
colnames(t1)[1] <- "seg"
colnames(t2)[1] <- "seg"
tt <- rbind(t1,t2)

#remove duplicates
tt <- tt[!duplicated(tt$seg),]

#Get sega
tt[,sega:=stri_sub(seg,1,-2)]

#remove non duplicate seg absolute
tt <- tt[duplicated(tt$sega),]

#get only the pos, regardless of path
tt[,posn:=as.numeric(stri_sub(pos,1,-4))]

tt[order(pos)]

## These are the inversions.

s <- fread("segments")





t1[,sense1:=!(stri_sub(`DA-402_03-2RL`,-1)=="+")]

t1[,pos:=stri_sub(`DA-402_03-2RL`,1,-2)]

t2[,sense2:=(stri_sub(`DA-416_04-2RL`,-1)=="+")]

t2[,pos:=stri_sub(`DA-416_04-2RL`,1,-2)]
setkey(t1,pos)
setkey(t2,pos)

tt <- merge(t1,t2, sort=FALSE, allow.cartesian=TRUE)

tm <- tt[sense1 != sense2]
tm <- tm[!duplicated(tm),]



