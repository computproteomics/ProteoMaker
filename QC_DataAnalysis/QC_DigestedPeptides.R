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
## Load peptide table before MS run (needs to be saved after beeing generated in RunAll.R)
#####################
load("RData/BeforeMS.RData")
#####################

#####################
## Peptides 
#####################
peptable <- BeforeMS$NonEnriched
peptable[,grepl("^C_", names(peptable))] <- sapply(seq_len(sum(grepl("^C_", names(peptable)))), function(x) {
  peptable[,grepl("^C_", names(peptable))][,x] <- unlist(peptable[,grepl("^C_", names(peptable))][,x])
})
boxplot(peptable[,grepl("^C_", names(peptable))])
cat("There are", sum(is.na(peptable[,grepl("^C_", names(peptable))])), "missing values.")
#####################

#####################
## Peptide -> protein
#####################
uAcc <- sort(unique(peptable$Accession))
protmat <- matrix(ncol = Param$NumCond*Param$NumReps, nrow = length(uAcc))
row.names(protmat) <- uAcc
for (acc in uAcc) {
  tmp <- peptable[peptable$Accession == acc,grepl("^C_", names(peptable))]
  tmp <- tmp[order(rowSums(tmp), decreasing = T),]
  if (nrow(tmp) >= 3) {
  protmat[row.names(protmat) == acc,] <- colSums(tmp[1:3,])
  } else {
  protmat[row.names(protmat) == acc,] <- colSums(tmp)
  }
  # protmat[row.names(protmat) == acc,] <- sapply(1:ncol(tmp), function(x) {
  #   mean(tmp[,x], na.rm = T)
  # })
}
colnames(protmat) <- names((peptable[,grepl("^C_", names(peptable))]))
boxplot(protmat, col = "cornflowerblue", main = "Protein values after summarisation")
#####################

#####################
## Statistical analysis of the proteoforms. 
#####################
conditions <- unique(gsub("_R.+", "", Param$quant_colnames))
pval <- vector()
means <- matrix(ncol = length(conditions), nrow = nrow(protmat))
missval <- matrix(ncol = ncol(protmat), nrow = nrow(protmat))
colnames(means) <- paste0("Mean_", conditions)
for (i in seq_len(nrow(protmat))) {
  mv <- is.na(protmat[i,])
  # Number of missing values per condition:
  numpercond <- colSums(!(sapply(conditions, function(x) mv[grepl(x, colnames(protmat))])))
  if (sum(numpercond >= 2) == length(numpercond)) {
    tres <- t.test(as.numeric(protmat[i,grepl(conditions[1], colnames(protmat))]), as.numeric(protmat[i,grepl(conditions[2], colnames(protmat))]))
    pval[i] <- tres$p.value
    means[i,] <- tres$estimate
  }
  missval[i,] <- mv
}

protmat <- cbind(protmat, means)
protmat <- as.data.frame(protmat)
protmat$pvalues <- pval
library(qvalue)
protmat$qvalues <- qvalue(pval)$qvalues
protmat$Mean_Diff <- means[,2] - means[,1]
protmat$Accession <- row.names(protmat)

load("RData/GroundTruthAbs.RData")

protmat$Regulated <- !is.na(GroundTruth$Regulation_Amplitude)[match(protmat$Accession, GroundTruth$Accession)]

if (Param$FracModProt > 0) {
  GroundTruth$modProt <- !(sapply(GroundTruth$PTMType, is.null))
}

protmat$colour <- 2^GroundTruth$Regulation_Amplitude[match(protmat$Accession, GroundTruth$Accession)]
protmat$colour[is.na(protmat$colour)] <- 1
protmat$colour <- factor(protmat$colour)
library(ggplot2)
g <- ggplot(data = protmat, aes(y = -log10(qvalues), x = Mean_Diff, col = colour)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("black", "red", "darkorange", "gold")) +
  theme_bw()
print(g)

passTest <- protmat[protmat$qvalues <= 0.05 & !is.na(protmat$qvalues),]
# hist(table(passTest$Accession))
cat("Number of regulated proteoforms that are not modified and modified:", paste(table(passTest$modProt), collapse = ", "))
#####################
