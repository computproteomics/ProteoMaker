################################################################################
#                     QC of the Ground truth proteoforms                       #
################################################################################

#####################
## Load parameters
#####################
setwd("/Users/rpk349/Documents/Boulot/GitRepo/PhosFake")
source("Parameter.R")
#####################


#####################
## Load ground truth data (needs to be saved after beeing generated in RunAll.R)
#####################
load("RData/GroundTruthAbs.RData")
#####################

#####################
## Statistical analysis of the proteoforms. 
#####################
conditions <- unique(gsub("_R.+", "", Param$quant_colnames))
pval <- vector()
means <- matrix(ncol = length(conditions), nrow = nrow(GroundTruth))
missval <- matrix(ncol = sum(grepl("^C_", names(GroundTruth))), nrow = nrow(GroundTruth))
colnames(means) <- paste0("Mean_", conditions)
for (i in seq_len(nrow(GroundTruth))) {
  mv <- is.na(GroundTruth[i,grepl("^C_", names(GroundTruth))])
  # Number of missing values per condition:
  numpercond <- colSums(!(sapply(conditions, function(x) mv[grepl(x, colnames(mv))])))
  if (sum(numpercond >= 2) == length(numpercond)) {
    tres <- t.test(as.numeric(GroundTruth[i,grepl(conditions[1], names(GroundTruth))]), as.numeric(GroundTruth[i,grepl(conditions[2], names(GroundTruth))]))
    pval[i] <- tres$p.value
    means[i,] <- tres$estimate
  }
  missval[i,] <- mv
}

GroundTruth <- cbind(GroundTruth, means)
GroundTruth$pvalues <- pval
library(qvalue)
GroundTruth$qvalues <- qvalue(pval)$qvalues
GroundTruth$Mean_Diff <- means[,2] - means[,1]
GroundTruth$Regulated <- !is.na(GroundTruth$Regulation_Amplitude)

if (Param$FracModProt > 0) {
  GroundTruth$modProt <- !(sapply(GroundTruth$PTMType, is.null))
}


library(ggplot2)
g <- ggplot(data = GroundTruth, aes(y = -log10(qvalues), x = Mean_Diff, col = Regulated)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_point(alpha = 0.3) +
  scale_color_manual(values = c("black", "red")) +
  theme_bw()
print(g)

passTest <- GroundTruth[GroundTruth$qvalues <= 0.05 & !is.na(GroundTruth$qvalues),]
# hist(table(passTest$Accession))
cat("Number of regulated proteoforms that are not modified and modified:", paste(table(passTest$modProt), collapse = ", "))
#####################