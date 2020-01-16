source("SamplePrep.R")

# Variables
ExpDesign <- rep(1:2,3)
Param <- list()
GroundTruth <- data.frame()
Digested <- data.frame()

### Parameter setting

## Ground truth
# Fraction of proteins selected for getting PTMs
Param$FracModProt <- 0.3
# fraction of modifiable proteins to be sampled for modifications  (might require more dedicated function 
# taking into account protein properties)
Param$FracModPerProt <- 2
# PTM types and fraction of PTMs (with respect to protein chosen to be modified)
Param$PTMTypes <- c("ph")
#Param$PTMNumber <- c("2")

# Distribution of multiply modified proteins is Poisson. Setting lambda
# Parameter is scaled to the number of possible PTM sites. Therefore set it to a value <1
Param$PTMMultipleLambda <- 0.1

# residues for PTM type and relative distribution to set them
Param$ModifiableResidues <- list()
Param$ModifiableResiduesDistr <- list()

for (mod in Param$PTMTypes) {
  Param$ModifiableResidues$mod <- c("S","T","Y")
  Param$ModifiableResiduesDistr$mod <- c(0.86,0.13,0.01)
}

# percentage of modifiable protein without non-modified forms
Param$RemoveNonModFormFrac <- 0.2

# Number of conditions
Param$NumCond <- 2
# Number of replicates
Param$NumReps <- 5

# vector for column names
quant_colnames <- paste0("C_",rep(1:Param$NumCond,each=Param$NumReps),"_R_", rep(1:Param$NumReps, Param$NumCond))

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
#####################
proteoforms <- samplePreparation(fasta.path = "fasta.fasta", parameters = Param)
#####################

# Proportion of missed cleavages peptides
Param$PropMissedCleavages <- 0.25

# Maximum number of missed cleavaged peptides
Param$MaxNumMissedCleavages <- 0

# filter for min and max of peptide length
Param$PepMinLength <- 7
Param$PepMaxLength <- 30

# Loss of phosphorylated peptides during enrichment
Param$EnrichmentLoss <- 0.2

# Enrichment efficiency: number of phosphorylated peptides with respect tohow many non-modified still enter the enriched fraction
Param$EnrichmentEfficiency <- 0.8

# Noise due to enrichment protocol
Param$EnrichmentNoise <- 0.2

## MS run

# Percentage of detected peptides
Param$PercDetectedPep <- 0.8

# Percentage of detected values (replicate/condition)
Param$PercDetectedVal <- 0.8

# Weights for intensity-dependence of non-detection (0 means no dependence). 
# Parameter is the power to the ranks (given by number 0 to 1)
Param$WeightDetectVal <- 0

# Wrong identifications
Param$WrongIDs <- 0.01

# Wrong localizations
Param$WrongLocalizations <- 0.01




