################################################################################
#                     Analysis of PhosFake mimicking UPS1                      #
################################################################################

#####################
## Load parameters
#####################
source("../Parameters/Parameter00.R")
name <- "UPS1_00"
pathToRes <- "../RData/"
#####################


#####################
## Load data
#####################
load(file = paste0(pathToRes, name, "/GroundTruth.RData"))
load(file = paste0(pathToRes, name, "/BeforeMS.RData"))
peptable <- BeforeMS$NonEnriched
#####################


#####################
## Ground truth data set:
#####################
cat("Number of protein accessions in the ground truth data set:\n")
print(length(unique(GroundTruth$Accession)))

cat("Number of proteoforms in the ground truth data set:\n")
print(nrow(GroundTruth))

cat("Number of missingvalues:\n")
print(sum(is.na(GroundTruth[,grepl("^C_", names(GroundTruth))])))

library(ggplot2)
library(reshape2)
# Quantification
par(mar = c(6,3,1,1))
hist(unlist(sapply(seq_along(GroundTruth$Regulation_Amplitude), function(x) GroundTruth$Regulation_Amplitude[x]*GroundTruth$Regulation_Pattern[[x]])),
     col = "#c10534",
     main = "",
     xlab = "Amplitude of regulation",
     breaks = 50)

gtab <- melt(GroundTruth[,grepl("^C_|Regulation_Amp", names(GroundTruth))], id.vars = "Regulation_Amplitude")
g <- ggplot(data = gtab, aes(y = value, x = variable)) +
  geom_violin(fill = "#c10534", alpha = 0.8) +
  theme_bw() +
  xlab("") +
  labs(title = "Quan values protein table")
print(g)
g <- ggplot(data = gtab[is.na(gtab$Regulation_Amplitude),], aes(y = value, x = variable)) +
  geom_violin(fill = "black", alpha = 0.8) +
  theme_bw() +
  xlab("") +
  labs(title = "Quan values protein table: non regulated only")
print(g)
#####################


#####################
## Peptide data set:
#####################
cat("Number of protein accessions in the peptide data set:\n")
print(length(unique(peptable$Accession)))

cat("Number of unique peptide in the peptide data set:\n")
print(length(unique(peptable$PepSequence)))

cat("Number of protein accessions unambigously assigned to a min. of one unique peptide:\n")
print(length(unique(peptable$Accession[!grepl(";", peptable$Accession, fixed = T)])))

# Count number of non-unique peptides:
peptable$AccessionCount <- sapply(peptable$Accession, function(x) {
  nchar(x) - nchar(gsub(";", "", x, fixed = T)) + 1
})
hist(peptable$AccessionCount, col = "cornflowerblue", breaks = 300, xlim = c(0,50))

cat("Number of missingvalues:\n")
print(sum(is.na(peptable[,grepl("^C_", names(peptable))])))

for (c in which(grepl("^C_|^Accession$", names(peptable)))) {
  peptable[,c] <- unlist(peptable[,c])
}

gtab <- melt(peptable[,grepl("^C_|^Accession$", names(peptable))], id.vars = "Accession")
g <- ggplot(data = gtab, aes(y = value, x = variable)) +
  geom_violin(fill = "cornflowerblue", alpha = 0.8) +
  theme_bw() +
  xlab("") +
  labs(title = "Quan values peptide table")
print(g)
#####################


#####################
## Peptide -> protein
#####################
source("../../../03_DataAnalysis.R")
cat("Number of non-ambigous accession submitted to statistical analysis:\n")
print(nrow(peptable[peptable$AccessionCount == 1,]))

cat("Number of missing proteins after in silico digestion:\n")
print(length(unique(GroundTruth$Accession)) - length(unique(peptable$Accession[peptable$AccessionCount == 1])))

protmat_minPep1 <- proteinSummarisation(peptable = peptable[peptable$AccessionCount == 1,], minUniquePep = 1, parameters = Param)
protmat_minPep2 <- proteinSummarisation(peptable = peptable[peptable$AccessionCount == 1,], minUniquePep = 2, parameters = Param)
boxplot(protmat_minPep1, col = "cornflowerblue", main = "Protein values after summarisation - min 1 unique pep.")
boxplot(protmat_minPep2, col = "cornflowerblue", main = "Protein values after summarisation - min 2 unique pep.")

cat("Number of proteins quantified with min. 1 unique peptide:\n")
print(nrow(na.omit(protmat_minPep1)))

cat("Number of proteins quantified with min. 2 unique peptide:\n")
print(nrow(na.omit(protmat_minPep2)))
#####################


# #####################
# ## Statistical analysis of the proteoforms. 
# #####################
# conditions <- unique(gsub("_R.+", "", Param$quant_colnames))
# lprot <- list(protmat_minPep1, protmat_minPep2)
# lout <- vector(mode = "list")
# library(qvalue)
# 
# for (j in 1:2) {
#   protmat <- lprot[[j]]
#   pval <- vector()
#   means <- matrix(ncol = length(conditions), nrow = nrow(protmat))
#   missval <- matrix(ncol = ncol(protmat), nrow = nrow(protmat))
#   colnames(means) <- paste0("Mean_", conditions)
#   protmat[is.infinite(protmat)] <- NA
#   for (i in seq_len(nrow(protmat))) {
#     mv <- is.na(protmat[i,])
#     # Number of missing values per condition:
#     numpercond <- colSums(!(sapply(conditions, function(x) mv[grepl(x, colnames(protmat))])))
#     if (sum(numpercond >= 2) == length(numpercond)) {
#       tres <- t.test(as.numeric(protmat[i,grepl(conditions[1], colnames(protmat))]), as.numeric(protmat[i,grepl(conditions[2], colnames(protmat))]))
#       pval[i] <- tres$p.value
#       means[i,] <- tres$estimate
#     }
#     missval[i,] <- mv
#   }
#   
#   protmat <- cbind(protmat, means)
#   protmat <- as.data.frame(protmat)
#   protmat$pvalues <- pval
#   
#   protmat$qvalues <- qvalue(pval)$qvalues
#   protmat$Mean_Diff <- means[,2] - means[,1]
#   protmat$Accession <- row.names(protmat)
#   protmat$Regulated <- !is.na(GroundTruth$Regulation_Amplitude)[match(protmat$Accession, GroundTruth$Accession)]
#   
#   protmat$colour <- 2^GroundTruth$Regulation_Amplitude[match(protmat$Accession, GroundTruth$Accession)]
#   protmat$colour[is.na(protmat$colour)] <- 1
#   protmat$colour <- factor(protmat$colour)
#   
#   lout[[j]] <- protmat
# }
# 
# g <- ggplot(data = protmat, aes(y = -log10(qvalues), x = Mean_Diff, col = colour)) +
#   geom_hline(yintercept = -log10(0.05)) +
#   geom_point(alpha = 0.8) +
#   scale_color_manual(values = c("black", "red", "darkorange", "gold")) +
#   theme_bw()
# print(g)
# #####################


sessionInfo()