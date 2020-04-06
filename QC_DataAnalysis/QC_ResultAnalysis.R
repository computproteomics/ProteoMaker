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
pathToInput <- "/QC_DataAnalysis/QC_output"

#--------------------

pathToInput <- paste0(wd, pathToInput)

#####################

#####################
## Load results
#####################

fnames <- list.files(pathToInput, full.names = T, recursive = T, pattern = ".RData")
# fnames <- fnames[!grepl(".txt$", fnames)] # Remove the reports
# fnames <- fnames[!grepl("Output/old", fnames)] # Remove the reports
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
# Match parameters:
#####################
fnames <- list.files(pathToInput, full.names = T, recursive = T, pattern = "TestedParam.txt")

param <- lapply(fnames, read.table, sep = "\t", header = T, stringsAsFactors = F)

tab <- merge(param, tab, by.y = "AnalysisName", by.x = "OutputNumber")

# pairs(as.matrix(tab[,-c(1:2)]))
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

col_param <- c(3:8)

## Peptides:

data <- reshape2::melt(tab[,c(col_param, 
                              which(grepl("NumberUniquePeptide", names(tab))))], 
                       id.vars = names(tab)[col_param])

g <- ggplot(data = data, aes(x = variable, y = value, fill = factor(SpeciesID))) +
  geom_bar(stat = "identity", position = position_dodge(), col = "black", alpha = 0.8) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Peptide count",
       fill = "Fasta ID") +
  ylab("Number of peptides") +
  xlab("")
print(g)

#--------------------

## Proteins and proteoforms:

data <- reshape2::melt(tab[,c(col_param, 
                              which(grepl("NumberUniqueProte", names(tab))))], 
                       id.vars = names(tab)[col_param])

g <- ggplot(data = data, aes(x = variable, y = value, fill = factor(SpeciesID))) +
  geom_bar(stat = "identity", position = position_dodge(), col = "black", alpha = 0.8) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "",
       subtitle = "Numbers at the end of labels are the number of unique peptides/ID",
       fill = "Fasta ID") +
  ylab("Number of proteins and proteoforms") +
  xlab("")
print(g)

data <- lapply(seq_along(lf), function(x) {
  data.frame("NumProteinPerPep" = lf[[x]]$NumProteinPerPep, 
             "NumProteoformPerPep" = lf[[x]]$NumProteoformPerPep, 
             "outputID" = rep(names(lf)[x], length(lf[[x]]$NumProteinPerPep)))
})
data <- rbindlist(data)
data <- reshape2::melt(data, id.vars = "outputID")

g <- ggplot(data = data, aes(x = value, fill = outputID)) +
  geom_density(alpha = 0.4, col = "black") +
  facet_wrap(~variable) +
  theme_bw()
print(g)
g <- ggplot(data = data, aes(x = value, fill = outputID)) +
  geom_histogram(alpha = 0.8, col = "black", position = position_dodge()) +
  facet_wrap(~variable) +
  theme_bw() +
  xlim(c(0,10)) +
  labs(title = "",
       subtitle = "x-axis cropped after 10",
       fill = "Output ID") +
  # ylab("Number of proteins and proteoforms") +
  xlab("")
print(g)

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

  paramID <- refgroup[[iter]]

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
  
  tit <- (paste("Proportion of missed cleavage:",
                unique(tab$PropMissedCleavages[tab$OutputNumber == paramID]),
                "\nSubsequent signal loss:",
                unique(tab$LeastAbundantLoss[tab$OutputNumber == paramID])
  ))
  
  g <- ggplot(data = data, aes(x = value, fill = factor(numMC))) +
    geom_histogram(alpha = 0.85, col = "black") +
    scale_fill_manual(values = col_mc) +
    facet_wrap(~Species) +
    # geom_density() +
    theme_bw() +
    labs(title = tit, subtitle = "facets are species (i.e. Fasta files)",
         fill = "num. MC") +
    ylab("Number of peptides") +
    xlab("MS signal")
  print(g)
  
}

#####################


#--------------------

sessionInfo()
