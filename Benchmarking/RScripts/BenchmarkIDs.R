################################################################################
#                     TO RUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

#####################
## Set path:
#####################

wd <- getwd()
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
         file = paste0(pathToRes, "/output", iter, ".RData"))
    cat("Save output", iter, "over", length(listtotest), "\n")
    iter <- iter + 1
  }
}
#####################

sessionInfo()
