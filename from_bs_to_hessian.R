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

m <- m + diag(ncol(m))*1e-6

c <- cov(m)

h <- -solve(c)

write.table(h, file=outfile, col.names=FALSE, row.names=F)


