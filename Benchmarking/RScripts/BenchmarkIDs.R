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
pathToRes <- paste0(wd, "/Output/BenchmarkIDs")
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
                    "PropMissedCleavages" = seq(from = 0, to = 1, by = 0.2), 
                    "MaxNumMissedCleavages" = 0:4,
                    "PepMinLength" = seq(from = 4, to = 12, by = 1),
                    "PepMaxLength" = seq(from = 16, to = 35, by = 1))
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
cat("Start generation of", length(listtotest), "parameter sets for digestion\n")
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
    peptable <- peptable[(pepLength >= parameters$PepMinLength & pepLength <= parameters$PepMaxLength),]
    peptable <- peptable[!is.na(peptable$ID),]
    
    upep <- names(table(peptable$PepSequence))[table(peptable$PepSequence) == 1]
    utab <- peptable[peptable$PepSequence %in% upep,]
    output <- list("Fasta" = names(GroundTruth)[i],
                   "Param" = x, 
                   "NumUniquePep" = length(unique(peptable$PepSequence)), 
                   "NumPepOneAcc" = length(upep), 
                   "NumAccPerMinPepNum" = table(table(utab$Accession)))
    save(output,
         file = paste0(pathToRes, "/output", iter, ".RData"))
    cat("Save output", iter, "over", length(listtotest), "\n")
    iter <- iter + 1
  }
}
#####################

sessionInfo()
