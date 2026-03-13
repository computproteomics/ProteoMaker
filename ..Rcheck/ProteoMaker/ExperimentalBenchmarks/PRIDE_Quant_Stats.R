library(moments)load("FullProtQuant.RData")

# setting all NAs to zero values (these are the empty ones in a sparse matrix)
fullProtTable@x[is.na(fullProtTable@x)] <- 0
fullProtTable <- drop0(fullProtTable)

# calculate means without the zeroes
RowInd <- fullProtTable@i + 1
means <- sapply(split(fullProtTable@x, RowInd), mean)
lengths <- sapply(split(fullProtTable@x, RowInd), length)
sds <- sapply(split(fullProtTable@x, RowInd), sd)
skewnesses <- sapply(split(fullProtTable@x, RowInd), skewness)
kurtoses <- sapply(split(fullProtTable@x, RowInd), kurtosis)
prots <- rownames(fullProtTable)


# check protein with more coverage
submatr <- fullProtTable[lengths>1000,]
