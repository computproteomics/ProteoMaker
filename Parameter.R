# Variables
ExpDesign <- rep(1:2,3)
Param <- list()
GroundTruth <- data.table()
Digested <- data.table()

### Parameter setting
## Ground truth
# Fraction of proteins getting PTMs
Param$FracModProt <- 0.3
# Fraction of modifiable proteins that get a PTM
Param$FracModPerProt <- 2
# PTM types
Param$PTMTypes <- c("ph")
# residues for PTM type and probabilities to set them
Param$ModifiableResidues <- list()
Param$ModifiableResiduesProb <- list()
for (mod in Param$PTMTypes) {
  Param$ModifiableResidues$mod <- c("S","T","Y")
  Param$ModifiableResiduesProb$mod <- c(0.7,0.2,0.1)*0.1
}
# fraction of modifiable proteins to be sampled for modifications  (might require more dedicated function 
# taking into account protein properties)
Param$ModSamplingFrac <- 2

# percentage of modifiable protein without non-modified forms
Param$RemoveNonModFormFrac <- 0.2

# Number of conditions
Param$NumCond <- 2
# Number of replicates
Param$NumReps <- 5

# General noise level of all quantitative values (standard deviation of normal distribution)
Param$QuantNoise <- 1


## Sample preparation


## MS run




