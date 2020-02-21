################################################################################
#                     TO RUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

library(ggplot2)
library(reshape2)
library("wesanderson")

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
fnames <- fnames[!grepl("Output/old", fnames)] # Remove the reports
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
# My param:
#####################
col_gd <- colorRampPalette(colors = c("darkred", "red", "gold"))
col_rep <- wes_palette("FantasticFox1", 4, type = "discrete")
#####################


#####################
# Identifications:
#####################

li <- lf[grepl("IDs", names(lf), fixed = T)]
li <- li[!grepl("MCs", names(li), fixed = T)]

library(zoo)

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

pairs(as.matrix(df[,3:ncol(df)]), pch = ".")

# mat <- as.matrix(df[,3:ncol(df)])
# mat <- apply(mat, 2, as.numeric)
# 
# pcamat <- mat
# pcamat <- pcamat[,colnames(pcamat) != "PepMaxLength"]
# pcamat <- scale(pcamat)
# pca <- prcomp(pcamat)
# # biplot(pca)
# 
# barplot(summary(pca)$importance[2,1:5], main = "Proportion of variance")
# gtab <- as.data.frame(pca$rotation)
# gtab$labels <- row.names(gtab)
# ggplot(gtab, aes(x = PC1, y = PC2, label = labels)) +
#   geom_text() +
#   theme_bw()

gtab <- df

ggplot(data = gtab[gtab$MaxNumMissedCleavages != 0 & gtab$MaxNumMissedCleavages %in% c(1,2) & gtab$PropMissedCleavages < 0.9,], aes(y = NumUniquePep, x = PropMissedCleavages, col = factor(PepMinLength))) +
  geom_point(alpha = 0.8) +
  # geom_smooth(se = F) +
  stat_summary(geom = "line", fun.y = "mean", size = 1.2) +
  theme_bw() +
  facet_wrap(~MaxNumMissedCleavages) +
  labs(title = "Facets = MaxNumMissedCleavages") +
  ylab("Number of unique peptide") +
  xlab("Proportion of missed cleavages") +
  scale_color_manual(values = col_gd(length(unique(gtab$PepMinLength))))

ggplot(data = gtab[gtab$MaxNumMissedCleavages != 0 & gtab$MaxNumMissedCleavages %in% c(1,2) & gtab$PropMissedCleavages < 0.9,], aes(y = NumPepOneAcc, x = PropMissedCleavages, col = factor(PepMinLength))) +
  geom_point(alpha = 0.8) +
  # geom_smooth(se = F) +
  stat_summary(geom = "line", fun.y = "mean", size = 1.2) +
  theme_bw() +
  facet_wrap(~MaxNumMissedCleavages) +
  labs(title = "Facets = MaxNumMissedCleavages") +
  ylab("Number of unique peptide corresponding to one accession") +
  xlab("Proportion of missed cleavages") +
  scale_color_manual(values = col_gd(length(unique(gtab$PepMinLength))))

library(ggnewscale)

gtab1 <-  melt(gtab, 
               id.vars = names(gtab)[1:(ncol(gtab)-3)])

g <- ggplot(data = gtab1[gtab1$PepMinLength >= 7 & gtab$MaxNumMissedCleavages > 0,], aes(y = value, x = PropMissedCleavages)) +
  geom_point(alpha = 0.8, aes(col = PepMinLength)) +
  scale_color_gradientn(colours = c("grey80", "black")) +
  theme_bw() +
  facet_wrap(~MaxNumMissedCleavages) +
  labs(title = "Facets = MaxNumMissedCleavages") +
  ylab("Number of unique protein accession") +
  xlab("Proportion of missed cleavages") 
  
g + 
  new_scale_color() +
  geom_smooth(data = gtab1[gtab1$PepMinLength >= 7 & gtab$MaxNumMissedCleavages > 0,], se = F, aes(col = variable, group = variable)) +
  scale_color_brewer(palette = "Set1")

ggplot(data = gtab[gtab$MaxNumMissedCleavages != 0,], aes(y = NumAccPerMinPepNum_2, x = PropMissedCleavages, col = factor(PepMinLength))) +
  geom_point(alpha = 0.8) +
  # geom_smooth(se = F) +
  stat_summary(fun.y = "mean", geom = "line", size = 1.2) +
  theme_bw() +
  facet_wrap(~MaxNumMissedCleavages) +
  labs(title = "Facets = MaxNumMissedCleavages") +
  ylab("Number of unique protein accession") +
  xlab("Proportion of missed cleavages") +
  scale_color_manual(values = col_gd(length(unique(gtab$PepMinLength))))

#####################

#####################
# Quantification:
#####################

li <- lf[grepl("ProtQuan", names(lf), fixed = T)]

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


for (qthresh in sort(unique(df$Threshqval))) {
  cat("q-value threshold:", qthresh, "\n")
  gtab <- df[df$Threshqval == qthresh,]
  gtab$NumberTrueRegulated <- as.numeric(gtab$NumberTrueRegulated)
  gtab$NumberTotRegulated <- as.numeric(gtab$NumberTotRegulated)
  gtab$QuantNoise <- as.numeric(gtab$QuantNoise)
  gtab$RelMV <- as.numeric(gtab$RelMV)
  gtab$AbsoluteQuanSD <- as.numeric(gtab$AbsoluteQuanSD)
  gtab$NumberTotRegulated <- as.numeric(gtab$NumberTotRegulated)
  gtab <- gtab[order(gtab$NumberTotRegulated, decreasing = T),]
  gtab$outputName <- factor(as.character(gtab$outputName), levels = as.character(gtab$outputName))
  
  # ggplot(data = gtab, aes(y = NumberTotRegulated, x = outputName, fill = NumReps)) +
  #   geom_bar(stat = "identity") +
  #   theme_classic() +
  #   theme(axis.title.x=element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank())
  
  gtab1 <- melt(gtab[,c("FRR", "RelMV", "NumberTotRegulated", "6.64385618977472", "1", "3.32192809488736", "NumReps", "ThreshNAQuantileProt", "QuantNoise", "NumberTrueRegulated", "AbsoluteQuanSD")], 
                id.vars = c("FRR", "NumReps", "RelMV", "NumberTotRegulated", "ThreshNAQuantileProt", "QuantNoise", "NumberTrueRegulated", "AbsoluteQuanSD"))
  # gtab1$NumberTotRegulated[gtab1$NumberTotRegulated == 0] <- NA
  g <- ggplot(data = gtab1, aes(x = RelMV, y = NumberTotRegulated, col = QuantNoise)) +
    geom_point(alpha = 0.4) +
    facet_wrap(~NumReps + variable) +
    theme_bw() +
    labs(title = "Facets = number of replicate")
  print(g)
  
  
  g <- ggplot(data = gtab1, aes(x = QuantNoise, y = NumberTotRegulated, col = NumReps)) +
    geom_hline(yintercept = 300, linetype = "dashed") +
    geom_point(alpha = 0.4) +
    stat_summary(fun.y = "mean", geom = "line", size = 1.2) +
    # geom_smooth(se = F) +
    facet_wrap(~ ThreshNAQuantileProt, nrow = 1) +
    theme_bw() +
    labs(title = "Facets = Detection threshold at the proteoform level") +
    scale_color_manual(values = wes_palette("FantasticFox1", length(unique(gtab1$NumReps)), type = "discrete")) +
    xlab("Noise (standard deviation) at the protein level") +
    ylab("Number of proteins considered regulated")
    
  print(g)
 
  g <- ggplot(data = gtab1, aes(x = QuantNoise, y = NumberTotRegulated, col = ThreshNAQuantileProt)) +
    geom_hline(yintercept = 300, linetype = "dashed") +
    geom_point(alpha = 0.4) +
    stat_summary(fun.y = "mean", geom = "line", size = 1.2, aes(group = ThreshNAQuantileProt)) +
    facet_wrap(~ NumReps) +
    theme_bw() +
    labs(title = "Facets = Number of replicates") +
    # scale_color_manual(values = wes_palette("Zissou1", length(unique(gtab1$ThreshNAQuantileProt)))) +
    scale_color_manual(values = col_gd(length(unique(gtab1$ThreshNAQuantileProt)))) +
    xlab("Noise (standard deviation) at the protein level") +
    ylab("Number of proteins considered regulated")
  print(g)
  
  g <- ggplot(data = gtab[gtab$ThreshNAQuantileProt == 0 | gtab$ThreshNAQuantileProt == 0.1,], aes(x = FRR, y = NumberTotRegulated, col = factor(QuantNoise))) +
    geom_point(alpha = 0.4) +
    # stat_summary(geom = "line", fun.y = "mean") +
    geom_line(aes(y=rollmean(NumberTotRegulated, 3, na.pad=TRUE))) +
    facet_wrap(~ NumReps + ThreshNAQuantileProt) +
    theme_bw() +
    labs(title = "Facets = NumReps + ThreshNAQuantileProt")
  print(g)
}

#####################
## Plot volcano plots:
#####################
detectionThresh <- sapply(li, function(x) {x$Param$ThreshNAQuantileProt})
QuantNoise <- sapply(li, function(x) {x$Param$QuantNoise})
NumReps <- sapply(li, function(x) {x$Param$NumReps})
volcano1 <- which(detectionThresh == 0 & QuantNoise == 0.26 & NumReps == 3)
volcano2 <- which(detectionThresh == 0 & QuantNoise == 0.26 & NumReps == 5)
volcano3 <- which(detectionThresh == 0.05 & QuantNoise == 0.26 & NumReps == 3)
volcano4 <- which(detectionThresh == 0.05 & QuantNoise == 0.26 & NumReps == 5)

dat <- li[c(volcano1, volcano2, volcano3, volcano4)]
gtab <- matrix(nrow = length(dat), ncol = 10)
for (r in seq_len(length(dat))) {
  gtab[r,1] <- dat[[r]]$Param$ThreshNAQuantileProt
  gtab[r,2] <- dat[[r]]$Param$QuantNoise
  gtab[r,3] <- dat[[r]]$Param$NumReps
  gtab[r,4:6] <- dat[[r]]$NumberRegPerAmplitude[order(dat[[r]]$NumberRegPerAmplitude[,1]),2]
  gtab[r,7] <- dat[[r]]$NumberTotRegulated
  gtab[r,8] <- dat[[r]]$Fasta
  gtab[r,9] <- dat[[r]]$Param$Threshqval
  gtab[r,10] <- ((dat[[r]]$NumberTotRegulated - dat[[r]]$NumberTrueRegulated) / dat[[r]]$NumberTotRegulated)*100
}
colnames(gtab) <- c("Detection_threshold", "QuantNoise", "NumReps", "FC2", "FC10", "FC100", "NumberTotRegulated", "Fasta", "Threshqval", "FRR")

gtab1 <- melt(as.data.frame(gtab), id.vars = c("Detection_threshold", "QuantNoise", "NumReps", "Fasta", "Threshqval"))

gtab1$value <- as.numeric(gtab1$value)
gtab1$variable <- factor(as.character(gtab1$variable), levels = c("NumberTotRegulated", "FC2", "FC10", "FC100", "FRR"))
ggplot(data = gtab1[gtab1$variable != "FRR",], aes(x = Detection_threshold, y = value, fill = variable)) +
  stat_summary(geom = "bar", fun.y = "mean", position = position_dodge(width = 1), alpha = 0.8, col = "black") +
  geom_hline(yintercept = c(96, 3*96), linetype = "dashed", alpha = 0.4) +
  geom_point(position = position_dodge(width = 1)) +
  facet_wrap(~Threshqval+NumReps, ncol = 2, nrow = 3) +
  scale_fill_manual(values = c("grey30", wes_palette("GrandBudapest1", 3, type = "discrete"))) +
  theme_bw()

gtab2 <- gtab1[gtab1$variable == "FRR",]
ggplot(data = gtab2, aes(x = Detection_threshold, y = value)) +
  stat_summary(geom = "bar", fun.y = "mean", position = position_dodge(width = 1), alpha = 0.8, col = "black", fill = "cornflowerblue") +
  # geom_hline(yintercept = c(96, 3*96), linetype = "dashed", alpha = 0.4) +
  geom_point(position = position_dodge(width = 1)) +
  facet_wrap(~Threshqval+NumReps, nrow = 2) +
  theme_bw() +
  ylab("Proportion of non-regulated proteins passing the statistical thresholds (%)")

# only Human: "uniprot-proteome%3AUP000005640_20200210"
gtab <- gtab[gtab[,which(colnames(gtab) == "Fasta")] == "uniprot-proteome%3AUP000005640_20200210",]

gtab1 <- melt(as.data.frame(gtab), id.vars = c("Detection_threshold", "QuantNoise", "NumReps", "Fasta", "Threshqval"))

gtab1$value <- as.numeric(gtab1$value)
gtab1$variable <- factor(as.character(gtab1$variable), levels = c("NumberTotRegulated", "FC2", "FC10", "FC100", "FRR"))
ggplot(data = gtab1[gtab1$variable != "FRR",], aes(x = Detection_threshold, y = value, fill = variable)) +
  geom_hline(yintercept = c(96, 3*96), linetype = "dashed", alpha = 0.4) +
  stat_summary(geom = "bar", fun.y = "mean", position = position_dodge(width = 1), alpha = 0.8, col = "black") +
  geom_point(position = position_dodge(width = 1)) +
  facet_wrap(~Threshqval+NumReps, ncol = 2, nrow = 3) +
  scale_fill_manual(values = c("grey30", wes_palette("GrandBudapest1", 3, type = "discrete"))) +
  theme_bw()

gtab2 <- gtab1[gtab1$variable == "FRR",]
ggplot(data = gtab2, aes(x = Detection_threshold, y = value)) +
  stat_summary(geom = "bar", fun.y = "mean", position = position_dodge(width = 1), alpha = 0.8, col = "black", fill = "cornflowerblue") +
  # geom_hline(yintercept = c(96, 3*96), linetype = "dashed", alpha = 0.4) +
  geom_point(position = position_dodge(width = 1)) +
  facet_wrap(~Threshqval+NumReps, nrow = 2) +
  theme_bw() +
  ylab("Proportion of non-regulated proteins passing the statistical thresholds (%)")


#####################
## Peptides + summarisation

li <- lf[grepl("PepQuan", names(lf), fixed = T)]

#####################

sessionInfo()
