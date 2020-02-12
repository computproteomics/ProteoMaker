################################################################################
#                     TO RUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################

#####################
## Load parameters
#####################
source("../../Parameter.R")
pathToRes <- "RData/"
# Parameters to test:
fastapath <- list.files(path = "input_data", full.names = T)
paramToTest <- list("PathToFasta" = fastapath, 
                    "PropMissedCleavages" = seq(from = 0, to = 1, by = 0.05), 
                    "MaxNumMissedCleavages" = 0:4,
                    "PepMinLength" = seq(from = 4, to = 12, by = 1),
                    "PepMaxLength" = seq(from = 16, to = 35, by = 1))
#####################

#####################
## Run the sample preparation simulation:
#####################
source("../../01_GenerateGroundTruth.R")
# Create the initial list of proteoforms for each fasta file tested:
lp <- vector(mode = "list")
for (i in 1:length(paramToTest$PathToFasta)) {
  Param$PathToFasta <- paramToTest$PathToFasta[i]
  lp[[i]] <- samplePreparation(fasta.path = Param$PathToFasta, parameters = Param)
  fasname <- gsub(getwd(), "", paramToTest$PathToFasta[i])
  fasname <- gsub(".fasta", "", fasname)
  names(lp)[i] <- gsub("^.+/", "", fasname)
  rm(fasname)
}

cat("Number of proteoforms per fasta:\n")
print(sapply(lp, dim))

GroundTruth <- lapply(lp, addProteoformAbundance, parameters = Param)

paramToTest <- paramToTest[-1]
#####################

#####################
## Digestion and sample peptide enrichment:
#####################
source("../../02_Digestion.R")
library(purrr)
listtotest <- cross(paramToTest)
cat("Start generation of", length(listtotest), "parameter sets for digestion\n")
# Digest all the proteoforms and get peptide table:
ld <- vector(mode = "list")
for (f in GroundTruth) {
  for (x in listtotest) {
    d <- DigestGroundTruth(GroundTruth = f, parameters = c(Param[!(names(Param) %in% names(x))], x))
    upep <- names(table(d$PepSequence))[table(d$PepSequence) == 1]
    utab <- d[d$PepSequence %in% upep,]
    ld[[length(ld) + 1]] <- list("Param" = x, 
                                 "NumUniquePep" = length(unique(d$PepSequence)), 
                                 "NumPepOneAcc" = length(upep), 
                                 "NumAccPerMinPepNum" = table(table(utab$Accession)))
  }
}
names(ld) <- names(GroundTruth)

save(ld, file = paste0(pathToRes, "ouptut.RData"))
#####################

sessionInfo()
