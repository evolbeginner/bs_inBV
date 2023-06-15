#! /bin/env Rscript


############################################
library(MASS)


############################################
args = commandArgs(trailingOnly=TRUE)

infile = args[1]
outfile = args[2]


############################################
m <- as.matrix(read.table(infile))

#m <- m * rnorm(ncol(m), 1, 8e-3)

# this is to ensure that those columns where all have exactly the same value will not as they will cause the covariance matrix of
# that line all equal to zero, an irreversible matrix.
for(col in 1:ncol(m)){
	unique_count <- length(unique(m[,col]))
	#if(unique_count <= 50){
	if(unique_count == 1 & m[1,col] <= 1e-3){
	#if(unique_count <= 2){
		print(col)
		m[,col] <- m[,col] * rnorm(nrow(m), 1, 8e-3)
		#print(m[m[,col]<0, col])
		#m[,col] <- rexp(length(m[,col]), 1)
	}
}

c <- cov(m)

h <- -solve(c)

write.table(h, file=outfile, col.names=FALSE, row.names=F)


