################################################################################
#                     TO RUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

#####################
## Set path:
#####################

wd <- getwd()
# # For Computerome:
# wd <- "/home/projects/jensenlab/people/malopa/PhosFake/BenchmarkingNoMod"
# #
pathToRes <- paste0(wd, "/Output/BenchmarkIDs_MCs")
pathToFasta <- paste0(wd, "/input_fasta")
pathToFasta <- list.files(path = pathToFasta, full.names = T, pattern = ".fasta")
pathToFunctions <- paste0(wd, "/Functions")
if (!dir.exists(pathToRes)) {
  cat("Create result directory:", pathToRes, "\n")
  dir.create(pathToRes)
}
#####################

#####################
## Load results
#####################
sapply(list.files(pathToFunctions, full.names = T), source)
# Parameters to test:
paramToTest <- list("PathToFasta" = pathToFasta, 
                    "PropMissedCleavages" = c(0, 0.2, 0.1), # According to Trypsing MCs in Chiva et al. 2014. Journal of Proteome Research.
                    "MaxNumMissedCleavages" = c(0, 2, 1))
#####################

#####################
## Run the sample preparation simulation:
#####################
# Create the initial list of proteoforms for each fasta file tested:
lp <- vector(mode = "list")
for (i in 1:length(paramToTest$PathToFasta)) {
  Param$PathToFasta <- paramToTest$PathToFasta[i]
  lp[[i]] <- samplePreparation(parameters = Param)
  fasname <- gsub(getwd(), "", paramToTest$PathToFasta[i])
  fasname <- gsub(".fasta", "", fasname)
  names(lp)[i] <- gsub("^.+/", "", fasname)
  rm(fasname)
}

GroundTruth <- lapply(lp, addProteoformAbundance, parameters = Param)

paramToTest <- paramToTest[-1]
#####################

#####################
## Digestion and sample peptide enrichment:
#####################
library(purrr)
listtotest <- cross(paramToTest)
# Select param. of interest:
listtotest <- listtotest[c(1,5,9)]
cat("Start generation of", length(listtotest) * length(GroundTruth), "parameter sets for digestion\n")
# Digest all the proteoforms and get peptide table:
iter <- 1
for (i in seq_along(GroundTruth)) {
  f <- GroundTruth[[i]]
  # I digest the table without missed cleavages only once to gain time:
  cat("Start digestion\n")
  # Get all the peptides without missed cleavage:
  d <- lapply(seq_len(nrow(f)),function(x) {
    getDigestTablesNoMC(f[x,], parameters = c(Param[!(names(Param) %in% names(x))], x))
  })
  cat("End digestion\n")
  for (x in listtotest) {
    d1 <- d
    x <- c(Param[!(names(Param) %in% names(x))], x)
    if (x$MaxNumMissedCleavages > 0 & x$PropMissedCleavages > 0) {
      ## Generate missed cleavages:
      cat("Start generation of missed-cleavages\n")
      for (el in seq_along(d1)) {
        prot <- d1[[el]]
        # print(prot)
        if (!is.null(prot)) {
          iter_mc <- 0
          r <- 1
          while (r < nrow(prot)) {
            if (runif(1, 0, 1) <= x$PropMissedCleavages & iter_mc < x$MaxNumMissedCleavages & (nrow(prot) - r) >= iter_mc) {
              newpep <- prot[r+1,]
              prot <- prot[-(r+1),]
              newpep$peptide <- paste(c(prot$peptide[r], newpep$peptide), collapse = "")
              newpep$start <- prot$start[r]
              iter_mc <- iter_mc + 1
              newpep$mc <- iter_mc
              prot[r,] <- newpep
              r <- r + 1
            } else {
              iter_mc <- 0
              r <- r + 1
            }
          }
          d1[[el]] <- prot
        }
      }
    }  
    cat("Row-bind missed cleavages\n")
    peptable <- as.data.frame(data.table::rbindlist(d1))
    cat("Table done\n")
    
    names(peptable)[names(peptable) == "peptide"] <- "PepSequence"
    names(peptable)[names(peptable) == "start"] <- "PepStart"
    names(peptable)[names(peptable) == "stop"] <- "PepStop"
    peptable$ID <- paste(peptable$PepSequence, peptable$PTMPos, peptable$PTMType, sep="_")
    peptable <- peptable[order(paste(peptable$Accession, peptable$PepStart)),]
    cat("Start filtering\n")
    #FILTERING
    pepLength <- nchar(peptable$PepSequence)
    peptable <- peptable[(pepLength >= x$PepMinLength & pepLength <= x$PepMaxLength),]
    peptable <- peptable[!is.na(peptable$ID),]
    
    upep <- names(table(peptable$PepSequence))[table(peptable$PepSequence) == 1]
    utab <- peptable[peptable$PepSequence %in% upep,]
    output <- list("Fasta" = names(GroundTruth)[i],
                   "Param" = x, 
                   "NumUniquePep" = length(unique(peptable$PepSequence)), 
                   "NumPepOneAcc" = length(upep), 
                   "NumAccPerMinPepNum" = table(table(utab$Accession)),
                   "data" = utab)
    save(output,
         file = paste0(pathToRes, "/output", iter, ".RData"))
    cat("Save output", iter, "over", length(listtotest), "\n")
    iter <- iter + 1
  }
}
#####################

#-------------------------------------------------------------------------------

#####################
## OUTPUT
#####################

#####################
## Load results
#####################

fnames <- list.files(pathToRes, full.names = T, recursive = T)
lf <- vector(mode = "list", length = length(fnames))
for (i in seq_along(fnames)) {
  load(fnames[i])
  lf[[i]] <- output
  rm(output)
}

names(lf) <- gsub(pathToRes, "", fnames)
names(lf) <- gsub("/output", "", names(lf), fixed = T)
names(lf) <- gsub(".RData", "", names(lf), fixed = T)

#####################

#####################
# Identifications:
#####################

li <- lf

library(zoo)

npar <- length(li[[1]]$Param)
cat("Column names for results of identification:\n")
print(names(li[[1]]))
df <- matrix(nrow = length(li), ncol = length(li[[1]]))
row.names(df) <- names(li)
colnames(df) <- c(names(li[[1]])[names(li[[1]]) != "NumAccPerMinPepNum" & names(li[[1]]) != "data"], "PropMissedCleavages", "MaxNumMissedCleavages")
colnames(df)[colnames(df) == "Param"] <- "outputName"
df <- as.data.frame(df)
for (r in seq_along(li)) {
  df$Fasta[r] <- li[[r]]$Fasta
  df$outputName[r] <- names(li)[r]
  df[r,3:ncol(df)] <- c(unlist(li[[r]][!(names(li[[r]]) %in% c("Fasta", "Param", "NumAccPerMinPepNum", "data"))]), li[[r]]$Param$PropMissedCleavages, li[[r]]$Param$MaxNumMissedCleavages)
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

pairs(as.matrix(df[,3:ncol(df)]), col = c(1,2,3)[match(df$Fasta, unique(df$Fasta))])

NumMC <- sapply(li, function(x) {
  table(x$data$mc)
})
tot <- sapply(NumMC, sum)
propMC <- sapply(seq_along(tot), function(x) {
  NumMC[[x]]/tot[x]
})
gtab <- data.frame("OutputName" =  unlist(sapply(seq_along(NumMC), function(x) {
  rep(df$outputName[x], length(NumMC[[x]]))
})),
"NumberMC" = unlist(sapply(NumMC, names)),
"ProportionMC" = unlist(propMC),
"Fasta" = unlist(sapply(seq_along(NumMC), function(x) {
  rep(df$Fasta[x], length(NumMC[[x]]))
})
))
gtab$OutputName <- as.character(gtab$OutputName)
gtab$OutputName[gtab$OutputName %in% c("1","4","7")] <- "No MC"
gtab$OutputName[gtab$OutputName %in% c("3","6","9")] <- "FASP"
gtab$OutputName[gtab$OutputName %in% c("2","5","8")] <- "in solution"
gtab$OutputName <- factor(as.character(gtab$OutputName), levels = sort(unique(gtab$OutputName))[c(3, 2, 1)])

ufasta <- sort(unique(gtab$Fasta))
fname <- c("mouse", "yeast", "human")
gtab$fname <- fname[match(gtab$Fasta, ufasta)]

ggplot(data = gtab, aes(y = ProportionMC, x = NumberMC, fill = OutputName)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9, preserve = "single"), col = "grey30") +
  scale_fill_manual(values = c("grey15", "#00A1D5FF", "#DF8F44FF")) +
  facet_wrap(~fname) +
  theme_bw() +
  xlab("Number of missed cleavages") +
  ylab("Proportion of peptides")

gtab <- df


ggplot(data = gtab, aes(y = NumUniquePep, x = PropMissedCleavages, col = factor(PepMinLength))) +
  geom_point(alpha = 0.8) +
  geom_smooth() +
  theme_bw() +
  facet_wrap(~MaxNumMissedCleavages) +
  labs(title = "Facets = MaxNumMissedCleavages")

ggplot(data = gtab, aes(y = NumPepOneAcc, x = PropMissedCleavages, col = factor(PepMinLength))) +
  geom_point(alpha = 0.8) +
  geom_smooth() +
  theme_bw() +
  facet_wrap(~MaxNumMissedCleavages) +
  labs(title = "Facets = MaxNumMissedCleavages")

gtab1 <-  melt(gtab, 
               id.vars = names(gtab)[1:(ncol(gtab)-3)])

ggplot(data = gtab1, aes(y = value, x = PropMissedCleavages, col = variable)) +
  geom_point(alpha = 0.8) +
  geom_smooth() +
  theme_bw() +
  facet_wrap(~MaxNumMissedCleavages) +
  labs(title = "Facets = MaxNumMissedCleavages")

ggplot(data = gtab, aes(y = NumAccPerMinPepNum_2, x = PropMissedCleavages, col = factor(PepMinLength))) +
  geom_point(alpha = 0.8) +
  geom_smooth() +
  theme_bw() +
  facet_wrap(~MaxNumMissedCleavages) +
  labs(title = "Facets = MaxNumMissedCleavages")

#####################


sessionInfo()
