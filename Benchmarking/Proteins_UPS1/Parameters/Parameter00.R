################################################################################
#                     PARAMETERS OF THE PHOSFAKE PIPELINE                      #
################################################################################

cat("Load Parameter00: No noise, no loss due to detection threshold.\n")
cat("Load Parameter00: No missed cleavages, peptide length between 7 and 30.\n")

# Create list that will contain all the parameters:
Param <- list()

#####################
## Variables for generating an experimental design:
#####################
# Number of conditions
Param$NumCond <- 2
# Number of replicates
Param$NumReps <- 4
#####################

#####################
## Path to input data
#####################
Param$PathToFasta <- c("/Users/rpk349/Documents/Boulot/GitRepo/PhosFake/InputData/fasta_full_human.fasta")
#####################

#####################
## Parameter for generation of the ground truth proteoform-centric data
#####################
## Generation of proteoform IDs:
# Fraction of the proteins selected to be modified:
Param$FracModProt <- 0 # Set to 0 if no modified proteins should be generated
# # fraction of modifiable proteins to be sampled for modifications >> (might require more dedicated function
# # taking into account protein properties)
# Param$FracModPerProt <- 2 # Here, a parameter of 2 will lead to 2 times more proteoforms than the set of selected proteins for modification
# # PTM types and fraction of PTMs (with respect to protein chosen to be modified)
# Param$PTMTypes <- c("ph")
# #Param$PTMNumber <- c("2")
# 
# # Distribution of multiply modified proteins is Poisson. Setting lambda
# # Parameter is scaled to the number of possible PTM sites. Therefore set it to a value <1
# Param$PTMMultipleLambda <- 0.1
# # residues for PTM type and relative distribution 
# Param$ModifiableResidues <- list()
# Param$ModifiableResiduesDistr <- list()
# # For the moment, we only consider "usual" phosphorylation on serine, threonine and tyrosines. 
# # We use the proportion of each phosphorylation type based on what is observed in the literature.
# for (mod in Param$PTMTypes) {
#   Param$ModifiableResidues$mod <- c("S","T","Y")
#   Param$ModifiableResiduesDistr$mod <- c(0.86,0.13,0.01)
# }
# # percentage of modifiable protein without non-modified forms
# Param$RemoveNonModFormFrac <- 0.2
#####################

#####################
## Generation of proteoform quantities:
#####################
# vector for column names
Param$quant_colnames <- paste0("C_",rep(1:Param$NumCond,each=Param$NumReps),"_R_", rep(1:Param$NumReps, Param$NumCond))
# General noise level of all quantitative values (standard deviation of normal distribution)
Param$QuantNoise <- 0
# Fraction of "differentially" regulated proteoforms
Param$DiffRegFrac <- NULL
# max. amplitude of difference (differentially regulated proteins). 
# Will be taken from uniform distribution with randomly chosen directions
Param$DiffRegMax <- NULL

# >>>> For input of custom set of regulated proteoforms:
Param$UserInputFoldChanges <- NULL
# I use the UPS1 setup (Ramus et al. 2016). KEPP NULL ID DO NOT WANT.
Param$UserInputFoldChanges <- list("NumRegProteoforms" = rep(96, 3),
                                   "RegulationFC" = rep(log2(c(100, 10, 2)), rep(96, 3)))

# threshold to remove quantitative values (proteoform level)
# COMMENT FROM MLP: I still think this is not right. A very abundant protein could 
# have a FC a lot higher and we'll still detect the conditions where it is the 
# lower intensity. Removing the values under the detection threshold cannot be 
# done on relative quan..
Param$ThreshNAProteoform <- NULL
## So if we want to add absolute quan:
Param$AbsoluteQuanMean <- 30.5
Param$AbsoluteQuanSD <- 3.6
Param$ThreshNAQuantileProt <- 0
## With this, we'll need to set up a strategy for proteoform distributions. Like:
# the proteoforms with the same accession will have a portion of the total signal of one
# point of the distribution? -> generate a distribution based on the parameters for
# each unique accession, and then devide it for each of the proteoforms of this accession?
#####################

#####################
## Parameters for enzymatic digestion (proteoforms -> peptides)
#####################
# Proportion of peptides with missed cleavages
Param$PropMissedCleavages <- 0
# Maximum number of missed cleavages per peptide
Param$MaxNumMissedCleavages <- 0
# filter for min and max of peptide length
Param$PepMinLength <- 7
Param$PepMaxLength <- 30
#####################

# #####################
# ## Paramters for simulating sample preparation before MS analysis (phospho-enrichment)
# #####################
# # Loss of phosphorylated peptides during enrichment
# Param$EnrichmentLoss <- 0.2
# # Enrichment efficiency: number of phosphorylated peptides with respect tohow many non-modified still enter the enriched fraction
# Param$EnrichmentEfficiency <- 0.8
# # Loss of signal for the non-modified peptides poluting the enriched fraction
# Param$EnrichmentNonModSignalLoss <- 0.7
# # Noise due to enrichment protocol
# Param$EnrichmentNoise <- 0.2
# #####################

#####################
## MS run
#####################
# Percentage of detected peptides
Param$PercDetectedPep <- 1
# Percentage of detected values (replicate/condition)
Param$PercDetectedVal <- 1
# Weights for intensity-dependence of non-detection (0 means no dependence). 
# Parameter is the power to the ranks (given by number 0 to 1)
Param$WeightDetectVal <- 0
# Add noise due to MS instrument:
Param$MSNoise <- 0
#####################

#####################
## MS search
#####################
# Wrong identifications
Param$WrongIDs <- 0.01
# Wrong localizations
Param$WrongLocalizations <- 0.01
#####################



