################################################################################
#                     TO RUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

library(ggplot2)
library(reshape2)

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

parOfInterest <- c( "NumReps", "QuantNoise", "ThreshNAQuantileProt", "AbsoluteQuanSD", "Threshqval" )
npar <- length(parOfInterest)
cat("Column names for results of identification:\n")
print(names(li[[1]]))
df <- matrix(nrow = length(li), ncol = length(li[[1]]) + npar - 1)
row.names(df) <- names(li)
colnames(df) <- c(names(li[[1]])[names(li[[1]]) != "NumberRegPerAmplitude"], parOfInterest)
colnames(df)[colnames(df) == "Param"] <- "outputName"
df <- as.data.frame(df)
for (r in seq_along(li)) {
  df$Fasta[r] <- li[[r]]$Fasta
  df$outputName[r] <- names(li)[r]
  cn <- c(unlist(li[[r]][!(names(li[[r]]) %in% c("Fasta", "Param", "NumberRegPerAmplitude"))]), unlist(li[[r]]$Param)[names(unlist(li[[r]]$Param)) %in% parOfInterest])
  df[r,3:ncol(df)] <- cn[match(colnames(df)[3:ncol(df)], names(cn))]
}
for (i in seq_len(nrow(li[[1]]$NumberRegPerAmplitude))) {
  df <- cbind(df, sapply(li, function(x) {
    x$NumberRegPerAmplitude[i,2]
  }))
  names(df)[ncol(df)] <- li[[1]]$NumberRegPerAmplitude$Regulation_Amplitude[i]
}

df$FRR <- (as.numeric(df$NumberTotRegulated) - as.numeric(df$NumberTrueRegulated)) / as.numeric(df$NumberUniqueProteoform)
df$RelMV <- as.numeric(df$NumberMissingValues) / as.numeric(df$NumReps)

#--------------------

mat <- as.matrix(df[,3:ncol(df)])
mat <- apply(mat, 2, as.numeric)
pairs(mat, col = df$NumReps)

pcamat <- mat[,3:(ncol(mat)-1)]
pcamat <- pcamat[,colnames(pcamat) != "Threshqval"]
pcamat <- scale(pcamat)
pca <- prcomp(pcamat)
# biplot(pca)

barplot(summary(pca)$importance[2,1:5], main = "Proportion of variance")
gtab <- as.data.frame(pca$rotation)
gtab$labels <- row.names(gtab)
ggplot(gtab, aes(x = PC1, y = PC2, label = labels)) +
  geom_text() +
  theme_bw()

for (qthresh in sort(unique(df$Threshqval))) {
  gtab <- df[df$Threshqval == qthresh,]
  gtab$NumberTrueRegulated <- as.numeric(gtab$NumberTrueRegulated)
  gtab$NumberTotRegulated <- as.numeric(gtab$NumberTotRegulated)
  gtab$QuantNoise <- as.numeric(gtab$QuantNoise)
  gtab$RelMV <- as.numeric(gtab$RelMV)
  gtab$AbsoluteQuanSD <- as.numeric(gtab$AbsoluteQuanSD)
  gtab$NumberTotRegulated <- as.numeric(gtab$NumberTotRegulated)
  gtab <- gtab[order(gtab$NumberTotRegulated, decreasing = T),]
  gtab$outputName <- factor(as.character(gtab$outputName), levels = as.character(gtab$outputName))
  
  ggplot(data = gtab, aes(y = NumberTotRegulated, x = outputName, fill = NumReps)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  gtab1 <- melt(gtab[,c("RelMV", "NumberTotRegulated", "6.64385618977472", "1", "3.32192809488736", "NumReps", "ThreshNAQuantileProt", "QuantNoise", "NumberTrueRegulated", "AbsoluteQuanSD")], 
                id.vars = c("NumReps", "RelMV", "NumberTotRegulated", "ThreshNAQuantileProt", "QuantNoise", "NumberTrueRegulated", "AbsoluteQuanSD"))
  # gtab1$NumberTotRegulated[gtab1$NumberTotRegulated == 0] <- NA
  g <- ggplot(data = gtab1, aes(x = RelMV, y = NumberTotRegulated, col = QuantNoise)) +
    geom_point(alpha = 0.4) +
    facet_wrap(~NumReps + variable) +
    theme_bw() +
    labs(title = "Facets = number of replicate")
  print(g)
  g <- ggplot(data = gtab1, aes(x = QuantNoise, y = NumberTotRegulated, col = NumReps)) +
    geom_point(alpha = 0.4) +
    geom_smooth() +
    facet_wrap(~ ThreshNAQuantileProt) +
    theme_bw() +
    labs(title = "Facets = Detection threshold at the proteoform level")
  print(g)
  g <- ggplot(data = gtab1, aes(x = AbsoluteQuanSD, y = NumberTotRegulated, col = NumReps)) +
    geom_point(alpha = 0.4) +
    geom_smooth() +
    facet_wrap(~ ThreshNAQuantileProt) +
    theme_bw() +
    labs(title = "Facets = Detection threshold at the proteoform level")
  print(g)
  g <- ggplot(data = gtab1, aes(x = factor(QuantNoise), y = NumberTotRegulated, col = ThreshNAQuantileProt)) +
    geom_point(alpha = 0.6, position = "jitter") +
    # stat_summary(geom= "bar" , position = "dodge", fun.y = "mean") +
    facet_wrap(~ NumReps) +
    theme_bw() +
    labs(title = "Facets = Number of replicates")
  print(g)
  g <- ggplot(data = gtab1, aes(x = factor(QuantNoise), y = NumberTotRegulated, fill = ThreshNAQuantileProt)) +
    # geom_point(alpha = 0.6, position = "jitter") +
    # stat_summary(geom= "bar" , position = "dodge", fun.y = "mean", col = "white") +
    geom_boxplot() +
    facet_wrap(~ NumReps) +
    theme_bw() +
    labs(title = "Facets = Number of replicates")
  print(g)
}

#####################

sessionInfo()
