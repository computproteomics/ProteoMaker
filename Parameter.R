################################################################################
#                     PARAMETERS OF THE PHOSFAKE PIPELINE                      #
################################################################################

# Create list that will contain all the parameters:
Param <- list()

#####################
## Variables for generating an experimental design:
#####################
# Number of conditions
Param$NumCond <- 2
# Number of replicates
Param$NumReps <- 5
#####################

#####################
## Path to input data
#####################
Param$PathToFasta <- c("InputData/fasta_full_human.fasta")
Param$PathToProteinList <- NULL
#Param$PathToProteinList <- c("InputData/Accessions.txt")
#####################

#####################
## Parameter for generation of the ground truth proteoform-centric data
#####################
## Generation of proteoform IDs:
# Fraction of the proteins selected to be modified:
Param$FracModProt <- 0.5 # Set to 0 if no modified proteins should be generated or 1 if only modified proteins should be generated.
# fraction of modifiable proteins to be sampled for modifications >> (might require more dedicated function
# taking into account protein properties)
Param$FracModPerProt <- 3 # Here, a parameter of 2 will lead to 2 times more proteoforms than the set of selected proteins for modification
# PTM types and fraction of PTMs (with respect to protein chosen to be modified)
Param$PTMTypes <- c("ph") #Only phosphorylation
#Param$PTMNumber <- c("2") #It is not used.
#Background distribution of PTM types.
Param$PTMTypesDist <- c(1) # 100% of modifications are phosphorylation.
Param$PTMTypesMass <- c(79.9663) #Phospho
#Example for multiple modification types
# Param$PTMTypes <- c("ph", "ub") # phosphorylation and ubiquitination (Example for multiple modification types)
# Param$PTMTypesDist <- c(0.83, 0.17) # 83% of modifications are phosphorylation and 17% ubiquitination. (Example for multiple modification types)
# Param$PTMTypesMass <- c(79.996, 114.043) #Phosphorylation and ubiquitination mass shifts

# Distribution of multiply modified proteins is Poisson. Setting lambda
# Parameter is scaled to the number of possible PTM sites. Therefore set it to a value <1
Param$PTMMultipleLambda <- 0.1
# Residues for PTM type and relative distribution 
# For the moment, we only consider "usual" phosphorylation on serine, threonine and tyrosines. 
# We use the proportion of each phosphorylation type based on what is observed in the literature.
Param$ModifiableResidues <- list(c("S","T","Y"))
Param$ModifiableResiduesDistr <- list(c(0.86,0.13,0.01))
#Example for multiple modification types
# Param$ModifiableResidues <- list(c("S","T","Y"), c("K"))
# Param$ModifiableResiduesDistr <- list(c(0.86,0.13,0.01), c(1))

# percentage of modifiable protein without non-modified forms
Param$RemoveNonModFormFrac <- 0
#####################

#####################
## Generation of proteoform quantities:
#####################
# Vector for column names
Param$QuantColnames <- paste0("C_",rep(1:Param$NumCond,each=Param$NumReps),"_R_", rep(1:Param$NumReps, Param$NumCond))
# General noise level of all quantitative values (standard deviation of normal distribution)
Param$QuantNoise <- 0.5
# Fraction of "differentially" regulated proteoforms
Param$DiffRegFrac <- 0.1
# max. amplitude of difference (differentially regulated proteins) on log2 scale. 
# Will be taken from uniform distribution with randomly chosen directions
Param$DiffRegMax <- 3

# >>>> For input of custom set of regulated proteoforms:
Param$UserInputFoldChanges <- NULL
# I use the UPS1 setup (Ramus et al. 2016). KEPP NULL ID DO NOT WANT.
# Param$UserInputFoldChanges <- list("NumRegProteoforms" = rep(96, 3),
#                                    "RegulationFC" = rep(log2(c(100, 10, 2)), rep(96, 3)))

# threshold to remove quantitative values (proteoform level)
# COMMENT FROM MLP: I still think this is not right. A very abundant protein could 
# have a FC a lot higher and we'll still detect the conditions where it is the 
# lower intensity. Removing the values under the detection threshold cannot be 
# done on relative quan..
# COMMENT FROM VS:
# What is the relation of this one to the others? Seems to be an arbitrary number which is 
# difficult to put into context
Param$ThreshNAProteoform <- -4
## So if we want to add absolute quan:
Param$AbsoluteQuanMean <- 30.5
Param$AbsoluteQuanSD <- 3.6
Param$ThreshNAQuantileProt <- 0.01
## With this, we'll need to set up a strategy for proteoform distributions. Like:
# the proteoforms with the same accession will have a portion of the total signal of one
# point of the distribution? -> generate a distribution based on the parameters for
# each unique accession, and then devide it for each of the proteoforms of this accession?
#####################

#####################
## Parameters for enzymatic digestion (proteoforms -> peptides)
#####################
# Enzyme
Param$Enzyme <- "trypsin"
# Proportion of peptides with missed cleavages
Param$PropMissedCleavages <- 0.05
# Maximum number of missed cleavages per peptide
Param$MaxNumMissedCleavages <- 2
# Miss cleavage abundance proportion from parental proteoform
#Param$PropMissedCleavagesAbundance <- c(0.8, 0.15, 0.5)
# Filter for min and max of peptide length
Param$PepMinLength <- 7
Param$PepMaxLength <- 30
# For parallel computing
Param$Cores <- NA
Param$ClusterType <- NA
#Param$Cores <- 10
#Param$ClusterType <- "PSOCK" #"FORK" for linux. PSOCK works for both linux and windows
# Remove percentage of least summarized abundant peptides (0-1).
Param$LeastAbundantLoss <- 0.1
#####################

#####################
## Paramters for simulating sample preparation before MS analysis (phospho-enrichment)
#####################
# Loss of phosphorylated peptides during enrichment
Param$EnrichmentLoss <- 0.2
# Enrichment efficiency: number of phosphorylated peptides with respect tohow many non-modified still enter the enriched fraction
Param$EnrichmentEfficiency <- 1
# Loss of signal for the non-modified peptides polluting the enriched fraction
Param$EnrichmentNonModSignalLoss <- 0
# Noise due to enrichment protocol
Param$EnrichmentNoise <- 0.2
#####################

#####################
## MS run
#####################
# Percentage of detected peptides
Param$PercDetectedPep <- 0.3
# Percentage of detected values (replicate/condition)
Param$PercDetectedVal <- 0.3
# Weights for intensity-dependence of non-detection (0 means no dependence). 
# Parameter is the power to the ranks (given by number 0 to 1) (Maximum of 40)
Param$WeightDetectVal <- 10
# Add noise due to MS instrument:
Param$MSNoise <- 0.5
#####################

#####################
## MS search    
#####################
# Wrong identifications
Param$WrongIDs <- 0.01
# Wrong localizations
Param$WrongLocalizations <- 0.01
#####################

#####################
## Filter MS search results
#####################
#Removes peptides that have more NA values than a specific number.
Param$MaxNAPerPep <- 100

#####################
## Summarize to proteins
#####################
#Summarization method (until now only "sum.top3" and "medpolish")
Param$ProtSummarization <- "medpolish"
#Minimum number of available unique peptides
Param$MinUniquePep <- 1

#####################
## Statistical testing of peptides and proteins
#####################
#Paired or unpaired
Param$StatPaired <- FALSE

#####################
## Load local test parameters if available:
#####################
if (file.exists("Parameter_mytests.R")) {
  source("Parameter_mytests.R")
}
#####################

#####################
## Output description:
#####################
cat("Load Parameter:\n")
cat("Quan. noise at proteoform level =", Param$QuantNoise, "(standard deviation of mean(log2(values)) = 0)\n")
cat("Loss due to detection threshold =", Param$ThreshNAQuantileProt, "(quantile of the quan. values)\n")
cat("Max. number of missed cleavages =", Param$MaxNumMissedCleavages, "occuring on", Param$PropMissedCleavages * 100, "% of the digested peptides\n")
cat("Min. peptide length =", Param$PepMinLength, "and max. peptide length =", Param$PepMaxLength, "\n")
cat("Loss in the mass spectrometer =", (1 - Param$PercDetectedPep) * 100, "% of the peptides, and", (1 - Param$PercDetectedVal) * 100, "% of the individual quantitative values (with a weigth based on inverse signal - power =", Param$WeightDetectVal, ")\n")
cat("MS noise =", Param$MSNoise, "\n")
#####################
