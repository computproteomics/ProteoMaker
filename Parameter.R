# Variables
ExpDesign <- rep(1:2,3)
Param <- list()
GroundTruth <- data.table()
Digested <- data.table()

### Parameter setting

## Ground truth
# Fraction of proteins selected for getting PTMs
Param$FracModProt <- 0.3
# fraction of modifiable proteins to be sampled for modifications  (might require more dedicated function 
# taking into account protein properties)
Param$FracModPerProt <- 2
# PTM types and fraction of PTMs (with respect to protein chosen to be modified)
Param$PTMTypes <- c("ph")
Param$PTMNumber <- c("2")

# Distribution of multiply modified proteins is Poisson. Setting lambda
Param$PTMMultipleLambda <- 2

# residues for PTM type and relative distribution to set them
Param$ModifiableResidues <- list()
Param$ModifiableResiduesDistr <- list()

for (mod in Param$PTMTypes) {
  Param$ModifiableResidues$mod <- c("S","T","Y")
  Param$ModifiableResiduesDistr$mod <- c(0.7,0.2,0.1)
}

# percentage of modifiable protein without non-modified forms
Param$RemoveNonModFormFrac <- 0.2

# Number of conditions
Param$NumCond <- 2
# Number of replicates
Param$NumReps <- 5

# General noise level of all quantitative values (standard deviation of normal distribution)
Param$QuantNoise <- 1

# Fraction of "differentially" regulated proteoforms
Param$DiffRegFrac <- 0.2

# max. amplitude of difference (differentially regulated proteins). 
#Will be taken from uniform distribution with randomly chosen directions
Param$DiffRegMax <- 2

# threshold to remove quantitative values (proteoform level)
Param$ThreshNAProteoform <- -2

## Sample preparation

# Higher miscleavage ratio for PTM
Param$MiscleavageRatio <- c(1.8)

# filter for min and max of peptide length
Param$PepMinLength <- 7
Param$PepMaxLength <- 30

# Loss of phosphorylated peptides during enrichment
Param$EnrichmentLoss <- 0.2

# Enrichment efficiency: how many non-modified still enter the enriched fraction
Param$EnrichmentEfficiency <- 0.8

# Noise due to enrichment protocol
Param$EnrichmentNoise <- 0.2




## MS run




