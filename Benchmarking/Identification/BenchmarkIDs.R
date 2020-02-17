################################################################################
#                     TO RUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

#####################
## Set path:
#####################

wd <- getwd()
pathToInput <- "ProteoformsQuan/RData"

#--------------------

pathToInput <- paste0(wd, "/", pathToInput)

#####################

#####################
## Load results
#####################
sapply(list.files(pathToInput, full.names = T), load)
# Parameters to test:
paramToTest <- list("PathToFasta" = pathToFasta, 
                    "PropMissedCleavages" = seq(from = 0, to = 1, by = 0.05), 
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
  for (x in listtotest) {
    d <- DigestGroundTruth(GroundTruth = f, parameters = c(Param[!(names(Param) %in% names(x))], x))
    upep <- names(table(d$PepSequence))[table(d$PepSequence) == 1]
    utab <- d[d$PepSequence %in% upep,]
    output <- list("Fasta" = names(GroundTruth)[i],
                   "Param" = x, 
                   "NumUniquePep" = length(unique(d$PepSequence)), 
                   "NumPepOneAcc" = length(upep), 
                   "NumAccPerMinPepNum" = table(table(utab$Accession)))
    save(output,
         file = paste0(pathToRes, "/ouptut", iter, ".RData"))
    iter <- iter + 1
  }
}
#####################

sessionInfo()
