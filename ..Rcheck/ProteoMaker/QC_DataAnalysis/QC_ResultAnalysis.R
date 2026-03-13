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
# pathToInput <- "/QC_DataAnalysis/QC_output"
# pathToInput <- paste0(wd, pathToInput)
pathToInput <- "/Volumes/Houdini/Projects/PhosFake/QC_DataAnalysis/QC_output"

#####################

#####################
## Data description
#####################

fnames <- list.files(pathToInput, full.names = T, recursive = T, pattern = ".RData")

cat("There are", length(fnames), "output files.\n")

cat("The input parameters were the following:\n")
log <- list.files(pathToInput, full.names = T, recursive = T, pattern = ".txt")
log <- log[!grepl("TestedParam", log)]
readLines(log)

rm(log)

cat("---------------------------\n")

#####################

#####################
# Retrieve parameters:
#####################

parnames <- list.files(pathToInput, full.names = T, 
                       recursive = T, pattern = "TestedParam.txt")

param <- lapply(parnames, read.table, sep = "\t", 
                header = T, stringsAsFactors = F)

PathToFasta1 <- unique(sapply(param, function(x) {x$PathToFasta}))[,1]
PathToFasta <- PathToFasta1

library(PTXQC)

doIter <- T
while (doIter) {
  toremove <- LCSn(PathToFasta, 3)
  # print(nchar(toremove))
  if (nchar(toremove) == 0) { doIter <- FALSE }
  PathToFasta <- gsub(toremove, "", PathToFasta)
  toremove <- NULL
}

fastaMapping <- data.frame("PathToFasta" = PathToFasta1, "SpeciesName" = PathToFasta)

rm(parnames)

for (i in seq_along(param)) {
  
  cat("---------------------\n\n")
  cat("Parameter set", i, ":\n")
  
  paramlist <- sapply(seq_len(ncol(param[[i]])-1), function(x) {unique(param[[i]][,x])})
  names(paramlist) <- names(param[[i]])[seq_len(ncol(param[[i]])-1)]
  print(paramlist)
  cat("---------------------\n")
  
}

#####################

#####################
## Load results
#####################

# #--------------------
# ## Subset selection
# #--------------------
# 
# PepMinLength <- 7
# PepMaxLength <- 30
# 
# cat("I select only the outputs with a min. peptide length of", PepMinLength,
#     "and a maximum peptide length of", PepMaxLength, "\n")
# 
# outputNumbers <- param[[1]]$OutputNumber[param[[1]]$PepMinLength == PepMinLength & param[[1]]$PepMaxLength == PepMaxLength]
# 
# cat("This corresponds to", length(outputNumbers), "output files.\n")
# 
# tokeep <- sapply(outputNumbers, function(x) {
#   which(grepl(paste0("_", x, ".RData"), fnames))
# })
# 
# fnames <- fnames[tokeep]
# 
# #--------------------

lf <- vector(mode = "list", length = length(fnames))
ltab <- vector(mode = "list", length = length(fnames))
for (i in seq_along(fnames)) {
  load(fnames[i])
  cmetrics <- sapply(metrics, class)
  lmetrics <- sapply(metrics, length)
  lf[[i]] <- metrics[!(cmetrics != "list" & lmetrics == 1)]
  ltab[[i]] <- metrics[cmetrics != "list" & lmetrics == 1]
  rm(metrics)
}

library(data.table)
tab <- rbindlist(ltab)

names(lf) <- gsub(pathToInput, "", fnames)
names(lf) <- gsub(".RData", "", names(lf), fixed = T)
names(lf) <- gsub("^.+\\/", "", names(lf))

tab$AnalysisName <- sapply(names(lf), function(x) {
  strsplit(x, "_", fixed = T)[[1]][length(strsplit(x, "_", fixed = T)[[1]])]
})

names(lf) <- sapply(names(lf), function(x) {
  strsplit(x, "_", fixed = T)[[1]][length(strsplit(x, "_", fixed = T)[[1]])]
})

#####################

#####################
# Merge parameters:
#####################

tab <- merge(param, tab, by.y = "AnalysisName", by.x = "OutputNumber")

tab$SpeciesName <- fastaMapping$SpeciesName[match(tab$PathToFasta, fastaMapping$PathToFasta)]

#####################

#####################
# My param:
#####################

col_gd <- colorRampPalette(colors = c("darkred", "red", "gold"))
col_rep <- wes_palette("FantasticFox1", 4, type = "discrete")
col_mc <- wes_palette("Darjeeling2", 5, type = "discrete")[c(5,1,3,4,2)]

#####################

#####################
# General numbers based on the different sets of parameters:
#####################

col_param <- c(3:6)

ggplot(data = tab, aes(x = LeastAbundantLoss, y = NumberUniquePeptide, col = factor(PropMissedCleavages))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~SpeciesName, scale = "free_y") +
  theme_bw()

ggplot(data = tab, aes(x = PropMissedCleavages, y = NumberUniquePeptide, col = factor(LeastAbundantLoss))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~SpeciesName, scale = "free_y") +
  theme_bw()

ggplot(data = tab, aes(x = LeastAbundantLoss, y = NumberUniqueProtein1, col = factor(PropMissedCleavages))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~SpeciesName, scale = "free_y") +
  theme_bw()

ggplot(data = tab, aes(x = PropMissedCleavages, y = NumberUniqueProtein1, col = factor(LeastAbundantLoss))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~SpeciesName, scale = "free_y") +
  theme_bw()

ggplot(data = tab, aes(x = LeastAbundantLoss, y = NumberUniqueProtein2, col = factor(PropMissedCleavages))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~SpeciesName, scale = "free_y") +
  theme_bw()

ggplot(data = tab, aes(x = PropMissedCleavages, y = NumberUniqueProtein2, col = factor(LeastAbundantLoss))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~SpeciesName, scale = "free_y") +
  theme_bw()

tabMC <- tab[,grepl("NumberUniquePeptide.MC", names(tab))]
tabMC <- cbind(tabMC, tab[,c(col_param, which(names(tab) == "SpeciesName"))])
tabMC <- reshape2::melt(tabMC, id.vars = c(names(tab)[col_param], "SpeciesName"))
tabMC$NumMC <- gsub("NumberUniquePeptide", "", tabMC$variable)

ggplot(data = tabMC[tabMC$SpeciesName == "human",], aes(x = factor(LeastAbundantLoss), y = value, fill = NumMC)) +
  geom_bar(alpha = 0.8, position = "dodge", stat = "identity") +
  facet_wrap(~PropMissedCleavages) +
  theme_bw() +
  ylab("Number of peptide\nbefore MS") +
  scale_fill_manual(values = col_mc) +
  labs(title = "facets = PropMissedCleavages", subtitle = "only peptidoforms from 7 to 30 aa long")

## Peptides:

data <- reshape2::melt(tab[,c(col_param, 
                              which(grepl("NumberUniquePeptide", names(tab))))], 
                       id.vars = names(tab)[col_param])

# g <- ggplot(data = data, aes(x = variable, y = value, fill = factor(SpeciesID))) +
#   geom_bar(stat = "identity", position = position_dodge(), col = "black", alpha = 0.8) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "Peptide count",
#        fill = "Fasta ID") +
#   ylab("Number of peptidoforms") +
#   xlab("")
# print(g)

for (sp in unique(tab$PathToFasta)) {
  pairs(tab[tab$PathToFasta == sp,c(4,8:14,20:ncol(tab))], col = tab$PropMissedCleavages*10, main = sp, cex = 0.5)
}

# g <- ggplot(data = data, aes(x = variable, y = value, fill = factor(SpeciesID))) +
#   geom_bar(stat = "identity", position = position_dodge(), col = "black", alpha = 0.8) +
#   theme_bw() +
#   facet_wrap(~LeastAbundantLoss + PropMissedCleavages) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "Peptide count",
#        fill = "Fasta ID") +
#   ylab("Number of peptidoforms") +
#   xlab("")
# print(g)

#--------------------

## Proteins and proteoforms:

data <- reshape2::melt(tab[,c(col_param, 
                              which(grepl("NumberUniqueProte", names(tab))))], 
                       id.vars = names(tab)[col_param])

# g <- ggplot(data = data, aes(x = variable, y = value, fill = factor(SpeciesID))) +
#   geom_bar(stat = "identity", position = position_dodge(), col = "black", alpha = 0.8) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "",
#        subtitle = "Numbers at the end of labels are the number of unique peptides/ID",
#        fill = "Fasta ID") +
#   ylab("Number of proteins and proteoforms") +
#   xlab("")
# print(g)

# g <- ggplot(data = data, aes(x = variable, y = value, fill = factor(SpeciesID))) +
#   geom_bar(stat = "identity", position = position_dodge(), col = "black", alpha = 0.8) + 
#   theme_bw() +
#   facet_wrap(~LeastAbundantLoss + PropMissedCleavages) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "",
#        subtitle = "Numbers at the end of labels are the number of unique peptides/ID",
#        fill = "Fasta ID") +
#   ylab("Number of proteins and proteoforms") +
#   xlab("")
# print(g)

data <- lapply(seq_along(lf), function(x) {
  data.frame("NumProteinPerPep" = lf[[x]]$NumProteinPerPep, 
             "NumProteoformPerPep" = lf[[x]]$NumProteoformPerPep, 
             "outputID" = rep(names(lf)[x], length(lf[[x]]$NumProteinPerPep)))
})
data <- rbindlist(data)
data <- reshape2::melt(data, id.vars = "outputID")

# g <- ggplot(data = data, aes(x = value, fill = outputID)) +
#   geom_density(alpha = 0.4, col = "black") +
#   facet_wrap(~variable) +
#   theme_bw()
# print(g)
# g <- ggplot(data = data, aes(x = value, fill = outputID)) +
#   geom_histogram(alpha = 0.8, col = "black", position = position_dodge()) +
#   facet_wrap(~variable) +
#   theme_bw() +
#   xlim(c(0,10)) +
#   labs(title = "",
#        subtitle = "x-axis cropped after 10",
#        fill = "Output ID") +
#   # ylab("Number of proteins and proteoforms") +
#   xlab("")
# print(g)

#####################


#####################
# Impact of PropMissedCleavages and LeastAbundantLoss on intensity distributions:
#####################

# library(reshape2)

# Group plots per same PropMissedCleavages + LeastAbundantLoss:
paircond <- tab[,names(tab) %in% c("PropMissedCleavages", "LeastAbundantLoss")]
paircond <- paircond[!duplicated(paircond),]
refgroup <- lapply(seq_len(nrow(paircond)), function(x) {
  tab$OutputNumber[tab$PropMissedCleavages == paircond[x,1] & tab$LeastAbundantLoss == paircond[x,2]]
})

# Make plots:


for (iter in seq_along(refgroup)) {

  pdfname <- paste0(wd, "/QC_DataAnalysis/QC_output/MC_tests_20200408_", iter, ".pdf")
  pdf(file = pdfname, width = 9, height = 7)
  
  paramID <- refgroup[[iter]]
  
  tit <- (paste("Proportion of missed cleavage:",
                unique(tab$PropMissedCleavages[tab$OutputNumber %in% paramID]),
                "\nSubsequent signal loss:",
                unique(tab$LeastAbundantLoss[tab$OutputNumber %in% paramID])
  ))
  
  ## Peptides:
  
  data <- reshape2::melt(tab[tab$OutputNumber %in% paramID,c(col_param, 
                                which(grepl("NumberUniquePeptide", names(tab))))], 
                         id.vars = names(tab)[col_param])
  
  data$SpeciesID <- c("Human", "Yeast", "Ecoli")[match(data$SpeciesID, c(2,3,1))]
  
  g <- ggplot(data = data, aes(x = variable, y = value, fill = factor(SpeciesID))) +
    geom_bar(stat = "identity", position = position_dodge(), col = "black", alpha = 0.8) +
    theme_bw() +
    facet_wrap(~SpeciesID, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = tit,
         fill = "Fasta ID") +
    ylab("Number of peptidoforms") +
    xlab("")
  print(g)
  
  ## Protein groups/proteoforms:
  
  data <- reshape2::melt(tab[tab$OutputNumber %in% paramID,c(col_param, 
                                                             which(grepl("NumberUniqueProt", names(tab))))], 
                         id.vars = names(tab)[col_param])
  
  data$SpeciesID <- c("Human", "Yeast", "Ecoli")[match(data$SpeciesID, c(2,3,1))]
  
  g <- ggplot(data = data, aes(x = variable, y = value, fill = factor(SpeciesID))) +
    geom_bar(stat = "identity", position = position_dodge(), col = "black", alpha = 0.8) +
    theme_bw() +
    facet_wrap(~SpeciesID, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = tit,
         subtitle = "Numbers (1,2) are the minimum number of unique peptide",
         fill = "Fasta ID") +
    ylab("Number of protein groups/proteoforms") +
    xlab("")
  print(g)
  
  ldata <- lf[names(lf) %in% paramID]
  names(ldata) <- names(lf)[names(lf) %in% paramID]
  allvalues <- vector(mode = "list")
  for (i in seq_along(ldata)) {
    data <- ldata[[i]]
    values <- as.data.frame(data$PeptideIntensities$Values)
    values$numMC <-data$PeptideIntensities$numMC
    values$param <-rep(names(ldata)[i], nrow(values))
    allvalues[[length(allvalues) + 1]] <- values
  }
  
  data <- rbindlist(allvalues)
  data <- reshape2::melt(data, id.vars = c("numMC", "param"))
  # data$Proportion_MC <- tab$PropMissedCleavages[match(data$param, tab$OutputNumber)]
  # data$Proportion_Loss <- tab$LeastAbundantLoss[match(data$param, tab$OutputNumber)]
  data$Species <- tab$SpeciesID[match(data$param, tab$OutputNumber)]
  data$Species <- c("Human", "Yeast", "Ecoli")[match(data$Species, c(2,3,1))]
  data <- data[data$Species != "Ecoli",] # Remove e. coli
  
  g <- ggplot(data = data, aes(x = value, fill = factor(numMC))) +
    geom_histogram(alpha = 0.85, col = "black") +
    scale_fill_manual(values = col_mc) +
    facet_wrap(~Species) +
    # geom_density() +
    theme_bw() +
    labs(title = tit, subtitle = "facets are species (i.e. Fasta files)",
         fill = "num. MC") +
    ylab("Number of peptidoforms") +
    xlab("MS signal")
  print(g)
  dev.off()
}


#####################


#--------------------

sessionInfo()
