################################################################################
#                     TO RUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################


#####################
## Set path:
#####################

wd <- getwd()
pathToInput <- "/Output"

#--------------------

pathToInput <- paste0(wd, pathToInput)

#####################

#####################
## Load results
#####################

fnames <- list.files(pathToInput, full.names = T, recursive = T)
fnames <- fnames[!grepl(".txt$", fnames)] # Remove the reports
lf <- vector(mode = "list", length = length(fnames))
for (i in seq_along(fnames)) {
  load(fnames[i])
  lf[[i]] <- output
  rm(output)
}

names(lf) <- gsub(pathToInput, "", fnames)
names(lf) <- gsub("/output", "", names(lf), fixed = T)
names(lf) <- gsub(".RData", "", names(lf), fixed = T)
names(lf) <- gsub("/", "", names(lf), fixed = T)

#####################

#####################
# Identifications:
#####################

li <- lf[grepl("IDs", names(lf), fixed = T)]

npar <- length(li[[1]]$Param)
cat("Column names for results of identification:\n")
print(names(li[[1]]))
df <- matrix(nrow = length(li), ncol = length(li[[1]]) + npar - 1)
row.names(df) <- names(li)
colnames(df) <- c(names(li[[1]])[names(li[[1]]) != "NumAccPerMinPepNum"], names(li[[1]]$Param))
colnames(df)[colnames(df) == "Param"] <- "outputName"
df <- as.data.frame(df)
for (r in seq_along(li)) {
  df$Fasta[r] <- li[[r]]$Fasta
  df$outputName[r] <- names(li)[r]
  df[r,3:ncol(df)] <- c(unlist(li[[r]][!(names(li[[r]]) %in% c("Fasta", "Param", "NumAccPerMinPepNum"))]), unlist(li[[r]]$Param))
}
df$NumAccPerMinPepNum_1 <- sapply(li, function(x) {
  x$NumAccPerMinPepNum[1]
})
df$NumAccPerMinPepNum_2 <- sapply(li, function(x) {
  x$NumAccPerMinPepNum[2]
})
df$NumAccPerMinPepNum_3 <- sapply(li, function(x) {
  x$NumAccPerMinPepNum[3]
})

#--------------------

pairs(as.matrix(df[,3:ncol(df)]))

#####################

#####################
# Quantification:
#####################

li <- lf[grepl("Quan", names(lf), fixed = T)]

npar <- length(li[[1]]$Param)
cat("Column names for results of identification:\n")
print(names(li[[1]]))
df <- matrix(nrow = length(li), ncol = length(li[[1]]) + npar - 1)
row.names(df) <- names(li)
colnames(df) <- c(names(li[[1]])[names(li[[1]]) != "NumberRegPerAmplitude"], names(li[[1]]$Param))
colnames(df)[colnames(df) == "Param"] <- "outputName"
df <- as.data.frame(df)
for (r in seq_along(li)) {
  df$Fasta[r] <- li[[r]]$Fasta
  df$outputName[r] <- names(li)[r]
  df[r,3:ncol(df)] <- c(unlist(li[[r]][!(names(li[[r]]) %in% c("Fasta", "Param", "NumberRegPerAmplitude"))]), unlist(li[[r]]$Param))
}
for (i in seq_len(nrow(li[[1]]$NumberRegPerAmplitude))) {
  df <- cbind(df, sapply(li, function(x) {
    x$NumberRegPerAmplitude[i,2]
  }))
  names(df)[ncol(df)] <- li[[1]]$NumberRegPerAmplitude$Regulation_Amplitude[i]
}


#--------------------

pairs(as.matrix(df[,3:ncol(df)]), col = df$NumReps)

library(ggplot2)
ggplot(data = df, aes(x = NumberTotRegulated, y = NumberTrueRegulated)) +
  geom_point() +
  facet_wrap(~NumReps)

#####################

sessionInfo()
