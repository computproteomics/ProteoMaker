# Function to create parameter table
param_table <- function() {
    
#     # Create vectors for each parameter and its attributes with default values
#     params_list <- list(
#         NumCond = c("Experimental Design", "paramsGroundTruth", "Number of conditions.", 2, 10, 2),
#         NumReps = c("Experimental Design", "paramsGroundTruth", "Number of replicates.", 1, 10, 3),
#         PathToFasta = c("Input Data", "paramsGroundTruth", "Path to input FASTA file.", NA, NA, "fasta_full_yeast.fasta"),
#         PathToProteinList = c("Input Data", "paramsGroundTruth", "Path to input protein list (optional).", NA, NA, NA),
#         FracModProt = c("Ground Truth Data", "paramsGroundTruth", "Fraction of proteins selected to be modified (0 for none, 1 for all).", 0, 1, 0),
#         FracModPerProt = c("Ground Truth Data", "paramsGroundTruth", "Factor for generating more proteoforms from selected proteins.", 0, 5, 0),
#         PTMTypes = c("Ground Truth Data", "paramsGroundTruth", "Types of post-translational modifications (PTMs).", NA, NA, NA),
#         PTMTypesDist = c("Ground Truth Data", "paramsGroundTruth", "Distribution of PTM types.", 0, 1, NA),
#         PTMTypesMass = c("Ground Truth Data", "paramsGroundTruth", "Mass shifts for each PTM type.", 0, 200, NA),
#         PTMMultipleLambda = c("Ground Truth Data", "paramsGroundTruth", "Poisson parameter for multiply modified proteins, scaled to the number of possible PTM sites.", 0, 5, NA),
#         ModifiableResidues = c("Ground Truth Data", "paramsGroundTruth", "Residues that can be modified for each PTM type.", NA, NA, NA),
#         ModifiableResiduesDistr = c("Ground Truth Data", "paramsGroundTruth", "Distribution of modifications across residues.", 0, 1, NA),
#         RemoveNonModFormFrac = c("Ground Truth Data", "paramsGroundTruth", "Percentage of modifiable proteins without non-modified forms.", 0, 1, 0),
#         QuantColnames = c("Proteoform Quantities", "paramsProteoformAb", "Vector of column names for quantitative data.", NA, NA, NA),
#         QuantNoise = c("Proteoform Quantities", "paramsProteoformAb", "Noise level of quantitative values (standard deviation).", 0, 5, 0.5),
#         DiffRegFrac = c("Proteoform Quantities", "paramsProteoformAb", "Fraction of differentially regulated proteoforms.", 0, 1, 0.3),
#         DiffRegMax = c("Proteoform Quantities", "paramsProteoformAb", "Maximum amplitude of difference for differentially regulated proteoforms (log2 scale).", 0, 10, 1),
#         UserInputFoldChanges = c("Proteoform Quantities", "paramsProteoformAb", "Custom set of regulated proteoforms (optional).", NA, NA, NA),
#         ThreshNAProteoform = c("Proteoform Quantities", "paramsProteoformAb", "Threshold to remove quantitative values at the proteoform level.", 0, 1000, 0),
#         AbsoluteQuanMean = c("Proteoform Quantities", "paramsProteoformAb", "Mean for absolute quantification.", 0, 100, 30.5),
#         AbsoluteQuanSD = c("Proteoform Quantities", "paramsProteoformAb", "Standard deviation for absolute quantification.", 0, 10, 3.6),
#         ThreshNAQuantileProt = c("Proteoform Quantities", "paramsProteoformAb", "Quantile threshold for removing quantitative values.", 0, 1, 0.01),
#         Enzyme = c("Enzymatic Digestion", "paramsDigest", "Enzyme used for digestion (e.g., trypsin).", NA, NA, "trypsin"),
#         PropMissedCleavages = c("Enzymatic Digestion", "paramsDigest", "Proportion of peptides with missed cleavages.", 0, 1, 0.01),
#         MaxNumMissedCleavages = c("Enzymatic Digestion", "paramsDigest", "Maximum number of missed cleavages per peptide.", 0, 10, 2),
#         PepMinLength = c("Enzymatic Digestion", "paramsDigest", "Minimum peptide length.", 1, 50, 7),
#         PepMaxLength = c("Enzymatic Digestion", "paramsDigest", "Maximum peptide length.", 1, 50, 30),
#         Cores = c("Enzymatic Digestion", "-", "Number of cores for parallel computing.", 1, 32, 1),
#         ClusterType = c("Enzymatic Digestion", "-", "Type of cluster for parallel computing (e.g., PSOCK).", NA, NA, "PSOCK"),
#         LeastAbundantLoss = c("Enzymatic Digestion", "paramsDigest", "Percentage of least abundant peptides to remove.", 0, 1, 0),
#         EnrichmentLoss = c("Sample Preparation", "paramsDigest", "Loss of phosphorylated peptides during enrichment.", 0, 1, 0.2),
#         EnrichmentEfficiency = c("Sample Preparation", "paramsDigest", "Efficiency of enrichment for phosphorylated peptides.", 0, 1, 1),
#         EnrichmentNonModSignalLoss = c("Sample Preparation", "paramsDigest", "Signal loss for non-modified peptides in the enriched fraction.", 0, 1, 0),
#         EnrichmentNoise = c("Sample Preparation", "paramsDigest", "Noise due to enrichment protocol.", 0, 1, 0.2),
#         DetectabilityThreshold = c("MS Run", "paramsMSRun", "Threshold for detectibility of peptides from PeptideRanger predictions (trained on ProteomicsDB).", 0, 1, 0.5),
# #        PercDetectedPep = c("MS Run", "paramsMSRun", "Percentage of detected peptides.", 0, 1, 0.2),
#         PercDetectedVal = c("MS Run", "paramsMSRun", "Percentage of detected values (replicate/condition).", 0, 1, 0.5),
#         WeightDetectVal = c("MS Run", "paramsMSRun", "Weights for intensity-dependence of non-detection.", 0, 1, 0.1),
#         MSNoise = c("MS Run", "paramsMSRun", "Noise due to the MS instrument.", 0, 1, 0.25),
#         WrongIDs = c("MS Search", "paramsMSRun", "Percentage of wrong identifications.", 0, 0.1, 0.01),
#         WrongLocalizations = c("MS Search", "paramsMSRun", "Percentage of wrong localizations.", 0, 0.1, 0),
#         MaxNAPerPep = c("MS Search Results Filter", "paramsMSRun", "Maximum number of NA values per peptide allowed.", 0, 100, 100),
#         ProtSummarization = c("Protein Summarization", "paramsDataAnalysis", "Method for summarizing to proteins (e.g., sum.top3, medpolish).", NA, NA, "medpolish"),
#         MinUniquePep = c("Protein Summarization", "paramsDataAnalysis", "Minimum number of unique peptides available.", 1, 10, 1),
#         StatPaired = c("Statistical Testing", "paramsDataAnalysis", "Whether the statistical testing is paired or unpaired.", FALSE, TRUE, FALSE)
#     )
    
    yaml_file <- system.file("config", "parameters.yaml", package = "PhosFake")
    
    # Read the YAML file
    params <- yaml::yaml.load_file(yaml_file)$params
    
    # Convert NA values from strings to real NA
    for (l in names(params)) {
        for (k in names(params[[l]])) {
            if (params[[l]][[k]] == "NA") {
                params[[l]][[k]] <- NA
            }
        }
    }
    
    # Convert the list to a data frames
    params_table <- do.call(rbind, lapply(params, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
    
    # Name the columns
    colnames(params_table) <- c("Category", "Group", "Explanation", "MinValue", "MaxValue", "DefaultValue")
    
    
    # Ensure MinValue and MaxValue are numeric
    # params_table$MinValue <- as.numeric(params_table$MinValue)
    # params_table$MaxValue <- as.numeric(params_table$MaxValue)
    # params_table$DefaultValue <- as.character(params_table$DefaultValue)
    # 
    params_table
}

# 
# 
# 
# # Create list that will contain all the parameters:
# Param <- list()
# #####################
# ## Variables for generating an experimental design:
# #####################
# # Number of conditions
# Param$NumCond <- 2
# # Number of replicates
# Param$NumReps <- 5
# #####################
# 
# #####################
# ## Path to input data
# #####################
# Param$PathToFasta <- c("InputData/fasta_full_human.fasta")
# Param$PathToProteinList <- NULL
# #Param$PathToProteinList <- c("InputData/Accessions.txt")
# #####################
# 
# #####################
# ## Parameter for generation of the ground truth proteoform-centric data
# #####################
# ## Generation of proteoform IDs:
# # Fraction of the proteins selected to be modified:
# Param$FracModProt <- 0.0 # Set to 0 if no modified proteins should be generated or 1 if only modified proteins should be generated.
# # fraction of modifiable proteins to be sampled for modifications >> (might require more dedicated function
# # taking into account protein properties)
# Param$FracModPerProt <- 2 # Here, a parameter of 2 will lead to 2 times more proteoforms than the set of selected proteins for modification
# # PTM types and fraction of PTMs (with respect to protein chosen to be modified)
# Param$PTMTypes <- c("ph") #Only phosphorylation
# #Param$PTMNumber <- c("2") #It is not used.
# #Background distribution of PTM types.
# Param$PTMTypesDist <- c(1) # 100% of modifications are phosphorylation.
# Param$PTMTypesMass <- c(79.9663) #Phospho
# #Example for multiple modification types
# # Param$PTMTypes <- c("ph", "ub") # phosphorylation and ubiquitination (Example for multiple modification types)
# # Param$PTMTypesDist <- c(0.83, 0.17) # 83% of modifications are phosphorylation and 17% ubiquitination. (Example for multiple modification types)
# # Param$PTMTypesMass <- c(79.996, 114.043) #Phosphorylation and ubiquitination mass shifts
# 
# # Distribution of multiply modified proteins is Poisson. Setting lambda
# # Parameter is scaled to the number of possible PTM sites. Therefore set it to a value <1
# Param$PTMMultipleLambda <- 0.5
# # Residues for PTM type and relative distribution 
# # For the moment, we only consider "usual" phosphorylation on serine, threonine and tyrosines. 
# # We use the proportion of each phosphorylation type based on what is observed in the literature.
# Param$ModifiableResidues <- list(c("S","T","Y"))
# Param$ModifiableResiduesDistr <- list(c(0.86,0.13,0.01))
# #Example for multiple modification types
# # Param$ModifiableResidues <- list(c("S","T","Y"), c("K"))
# # Param$ModifiableResiduesDistr <- list(c(0.86,0.13,0.01), c(1))
# 
# # percentage of modifiable protein without non-modified forms
# Param$RemoveNonModFormFrac <- 0
# #####################
# 
# #####################
# ## Generation of proteoform quantities:
# #####################
# # Vector for column names
# Param$QuantColnames <- paste0("C_",rep(1:Param$NumCond,each=Param$NumReps),"_R_", rep(1:Param$NumReps, Param$NumCond))
# # General noise level of all quantitative values (standard deviation of normal distribution)
# Param$QuantNoise <- 0.0
# # Fraction of "differentially" regulated proteoforms
# Param$DiffRegFrac <- 0.0
# # max. amplitude of difference (differentially regulated proteins) on log2 scale. 
# # Will be taken from uniform distribution with randomly chosen directions
# Param$DiffRegMax <- 3
# 
# # >>>> For input of custom set of regulated proteoforms:
# Param$UserInputFoldChanges <- NULL
# # I use the UPS1 setup (Ramus et al. 2016). KEEP NULL ID DO NOT WANT.
# # Param$UserInputFoldChanges_NumRegProteoforms = rep(96, 3)
# # Param$UserInputFoldChanges_NumRegProteoforms_RegulationFC = rep(log2(c(100, 10, 2)), rep(96, 3))
# 
# # threshold to remove quantitative values (proteoform level)
# # COMMENT FROM MLP: I still think this is not right. A very abundant protein could 
# # have a FC a lot higher and we'll still detect the conditions where it is the 
# # lower intensity. Removing the values under the detection threshold cannot be 
# # done on relative quan..
# # COMMENT FROM VS:
# # What is the relation of this one to the others? Seems to be an arbitrary number which is 
# # difficult to put into context
# Param$ThreshNAProteoform <- -100
# ## So if we want to add absolute quan:
# Param$AbsoluteQuanMean <- 30.5
# Param$AbsoluteQuanSD <- 3.6
# Param$ThreshNAQuantileProt <- 0.01
# ## With this, we'll need to set up a strategy for proteoform distributions. Like:
# # the proteoforms with the same accession will have a portion of the total signal of one
# # point of the distribution? -> generate a distribution based on the parameters for
# # each unique accession, and then devide it for each of the proteoforms of this accession?
# #####################
# 
# #####################
# ## Parameters for enzymatic digestion (proteoforms -> peptides)
# #####################
# # Enzyme
# Param$Enzyme <- "trypsin"
# # Proportion of peptides with missed cleavages
# Param$PropMissedCleavages <- 1
# # Maximum number of missed cleavages per peptide
# Param$MaxNumMissedCleavages <- 6
# # Miss cleavage abundance proportion from parental proteoform
# #Param$PropMissedCleavagesAbundance <- c(0.8, 0.15, 0.5)
# # Filter for min and max of peptide length
# Param$PepMinLength <- 7
# Param$PepMaxLength <- 30
# # For parallel computing
# Param$Cores <- 1
# Param$ClusterType <- "PSOCK"
# #Param$Cores <- 10
# #Param$ClusterType <- "PSOCK" #"FORK" for linux. PSOCK works for both linux and windows
# # Remove percentage of least summarized abundant peptides (0-1).
# Param$LeastAbundantLoss <- 0.0
# #####################
# 
# #####################
# ## Parameters for simulating sample preparation before MS analysis (phospho-enrichment)
# #####################
# # Loss of phosphorylated peptides during enrichment
# Param$EnrichmentLoss <- 0.2
# # Enrichment efficiency: number of phosphorylated peptides with respect tohow many non-modified still enter the enriched fraction
# Param$EnrichmentEfficiency <- 1
# # Loss of signal for the non-modified peptides polluting the enriched fraction
# Param$EnrichmentNonModSignalLoss <- 0
# # Noise due to enrichment protocol
# Param$EnrichmentNoise <- 0.2
# #####################
# 
# #####################
# ## MS run
# #####################
# # Percentage of detected peptides
# Param$PercDetectedPep <- 0.3
# # Percentage of detected values (replicate/condition)
# Param$PercDetectedVal <- 0.3
# # Weights for intensity-dependence of non-detection (0 means no dependence). 
# # Parameter is the power to the ranks (given by number 0 to 1) (Maximum of 40)
# Param$WeightDetectVal <- 10
# # Add noise due to MS instrument:
# Param$MSNoise <- 0.5
# #####################
# 
# #####################
# ## MS search    
# #####################
# # Wrong identifications
# Param$WrongIDs <- 0.01
# # Wrong localizations
# Param$WrongLocalizations <- 0.01
# #####################
# 
# #####################
# ## Filter MS search results
# #####################
# #Removes peptides that have more NA values than a specific number.
# # Veit says: might become obsolete
# Param$MaxNAPerPep <- 100
# 
# #####################
# ## Summarize to proteins
# #####################
# #Summarization method (until now only "sum.top3" and "medpolish")
# Param$ProtSummarization <- "medpolish"
# #Minimum number of available unique peptides
# Param$MinUniquePep <- 1
# 
# #####################
# ## Statistical testing of peptides and proteins
# #####################
# #Paired or unpaired
# Param$StatPaired <- FALSE
# 
# #####################
# ## Load local test parameters if available:
# #####################
# if (file.exists("Parameter_mytests.R")) {
#     source("Parameter_mytests.R")
# }
# #####################
# 
# #####################
# ## Output description:
# #####################
# cat("Load Parameter:\n")
# cat("Quan. noise at proteoform level =", Param$QuantNoise, "(standard deviation of mean(log2(values)) = 0)\n")
# cat("Loss due to detection threshold =", Param$ThreshNAQuantileProt, "(quantile of the quan. values)\n")
# cat("Max. number of missed cleavages =", Param$MaxNumMissedCleavages, "occuring on", Param$PropMissedCleavages * 100, "% of the digested peptides\n")
# cat("Min. peptide length =", Param$PepMinLength, "and max. peptide length =", Param$PepMaxLength, "\n")
# cat("Loss in the mass spectrometer =", (1 - Param$PercDetectedPep) * 100, "% of the peptides, and", (1 - Param$PercDetectedVal) * 100, "% of the individual quantitative values (with a weigth based on inverse signal - power =", Param$WeightDetectVal, ")\n")
# cat("MS noise =", Param$MSNoise, "\n")
# #####################
