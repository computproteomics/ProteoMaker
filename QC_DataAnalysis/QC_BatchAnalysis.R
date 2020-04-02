################################################################################
#                       QC batch test of the parameters                        #
################################################################################

#####################
## Load parameters
#####################
setwd("/Users/rpk349/Documents/Boulot/GitRepo/PhosFake")
source("Parameter.R")
#####################

#####################
## Load parameters
#####################

# Fasta files:
pathToFasta <- paste0(getwd(), "/Benchmarking/input_fasta")
pathToFasta <- list.files(path = pathToFasta, full.names = T, pattern = ".fasta")

#--------------------
## Generate a list of parameters to test:

paramToTest <- list("PathToFasta" = pathToFasta, 
                    "InputProportionMC" = seq(from = 0.1, to = 0.8, by = 0.1), 
                    "MaxNumMissedCleavages" = 4,
                    "PepMinLength" = c(6,7),
                    "PepMaxLength" = seq(from = 15, to = 35, by = 5),
                    "LeastAbundantLoss" = seq(from = 0, to = 0.8, by = 0.2))

#--------------------
# Name of report:

repportName <- paste0("QC_output/QC_Report_BatchAnalysis_",
                      gsub(":", "", format(Sys.time(), "%Y%b%d%X")))

#####################

#####################
## Metrics that we want to retreive:
#####################

metrics <- list()

## Peptides:

# Number of unique peptides (including modification):
metrics$NumberUniquePeptide <- vector()
# Number of unique peptides (including modification) with no miss-cleavage:
metrics$NumberUniquePeptide0MC <- vector()
# Number of unique peptides (including modification) with 1 miss-cleavage:
metrics$NumberUniquePeptide1MC <- vector()
# Number of unique peptides (including modification) with 2 miss-cleavages:
metrics$NumberUniquePeptide2MC <- vector()
# Intensities of the peptides after summarisation:
metrics$PeptideIntensities$Values <- vector(mode = "list")
# Matching number of missed-cleavages:
metrics$PeptideIntensities$numMC <- vector(mode = "list")

## Proteins:

# Number of unique proteins with min. 1 proteotypic peptide:
metrics$NumberUniqueProtein1 <- vector()
# Number of unique proteins with min. 2 proteotypic peptides:
metrics$NumberUniqueProtein2 <- vector()
# Number of proteins a peptide is matched to:
metrics$NumProteinPerPep <- vector(mode = "list")

## Proteoforms:

# Number of unique proteoforms with  min. 1 peptide specific to the proteoform:
metrics$NumberUniqueProteoforms1 <- vector()
# Number of unique proteoforms with  min. 1 peptide specific to the proteoform:
metrics$NumberUniqueProteoforms2 <- vector()
# Number of proteoforms a peptide is matched to:
metrics$NumProteoformPerPep <- vector(mode = "list")

## Quantities:

# Correlation of the protein quantities per condition:
metrics$CorProteinQuan <- vector(mode = "list")
# Correlation of the proteoform quantities per condition:
metrics$CorProteoformQuan <- vector(mode = "list")

#####################

#####################
## Run the sample preparation:
#####################

source("01_GenerateGroundTruth.R")

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
## Digestion and sample peptide enrichment according to each set of parameters:
#####################
source("02_Digestion.R")

library(purrr)
listtotest <- cross(paramToTest)

# Generate table with parameters and results:

res <- bind_rows(listtotest)
res$OutputNumber <- 1:nrow(res)

#>>> Start writing report <<<#

sink(paste0(repportName,".txt"))

cat("Parameter:\n")
cat("Quan. noise at proteoform level =", Param$QuantNoise, "(standard deviation of mean(log2(values)) = 0)\n")
cat("Loss due to detection threshold =", Param$ThreshNAQuantileProt, "(quantile of the quan. values)\n")
cat("Max. number of missed cleavages =", Param$MaxNumMissedCleavages, "occuring on", Param$PropMissedCleavages * 100, "% of the digested peptides\n")
cat("Min. peptide length =", Param$PepMinLength, "and max. peptide length =", Param$PepMaxLength, "\n")
cat("Loss in the mass spectrometer =", (1 - Param$PercDetectedPep) * 100, "% of the peptides, and", (1 - Param$PercDetectedVal) * 100, "% of the individual quantitative values (with a weigth based on inverse signal - power =", Param$WeightDetectVal, ")\n")
cat("MS noise =", Param$MSNoise, "\n")

cat("Generation of", length(listtotest), "parameter sets for digestion\n")
cat("They can be found in the following table:", paste0(repportName, "_TestedParam.txt"), "\n")
write.table(res, 
            file = paste0(repportName, "_TestedParam.txt"), 
            sep = "\t",
            row.names = F)

sink()

#<<<>>>#

cat("Start generation of", length(listtotest), "parameter sets for digestion\n")

# Digest all the proteoforms and get peptide table:
iter <- 1

for (i in seq_along(GroundTruth)) { # iter through the fasta files
  
  f <- GroundTruth[[i]]
  
  for (x in listtotest) {
    
    # Generate the proportion of abundance of missed-cleaved peptides beased on the max number of missed cleavages accepted:
    if (x$MaxNumMissedCleavages > 0) {
      x$PropMissedCleavagesAbundance <- c((1-x$InputProportionMC)^2, sapply(seq_len(x$MaxNumMissedCleavages), function(mc) {
        (1-x$InputProportionMC)^2*(x$InputProportionMC^mc)
      }))
    } else {
      x$PropMissedCleavagesAbundance <- 1
    }
    
    parametersLoop <- c(Param[!(names(Param) %in% names(x))], x)
  
    # Digest all the proteoforms and get peptide table:
    peptable <- digestGroundTruth(proteoforms = f, parameters = parametersLoop)
    peptable <- digestionProductSummarization(peptides = peptable, parameters = parametersLoop)
    
    BeforeMS <- filterDigestedProt(DigestedProt = peptable, parameters = parametersLoop)
    
    ######### Add results to res:
    
    metrics$NumberUniquePeptide[iter] <- nrow(peptable)
    metrics$NumberUniquePeptide0MC[iter] <- nrow(peptable[sapply(peptable$MC, unique) == 0,])
    metrics$NumberUniquePeptide1MC[iter] <- nrow(peptable[sapply(peptable$MC, unique) == 1,])
    metrics$NumberUniquePeptide2MC[iter] <- nrow(peptable[sapply(peptable$MC, unique) == 2,])
    metrics$PeptideIntensities$Values[iter] <- as.matrix(peptable[,grepl("^C_", names(peptable))])
    metrics$PeptideIntensities$numMC[iter] <- sapply(peptable$MC, unique)
    vec_numProtein <- sapply(peptable$Accession, function(x) {length(unique(x))})
    vec_numProteoform <- sapply(peptable$Proteoform_ID, function(x) {length(unique(x))})
    metrics$NumberUniqueProtein1[iter] <- nrow(peptable[vec_numProtein == 1,])
    metrics$NumberUniqueProtein2[iter] <- nrow(peptable[vec_numProtein <= 2,])
    metrics$NumProteinPerPep[iter] <- vec_numProtein
    metrics$NumberUniqueProteoforms1[iter] <- nrow(peptable[vec_numProteoform == 1,])
    metrics$NumberUniqueProteoforms2[iter] <- nrow(peptable[vec_numProteoform <= 2,])
    metrics$NumProteoformPerPep[iter] <- vec_numProteoform
    
    save(metrics,
         file = paste0(repportName, "_", iter, ".RData"))
    
    cat("Save output", iter, "over", length(listtotest)*length(GroundTruth), "\n")
    #########
    
    iter <- iter + 1
    
  }
  
}


#####################

#>>> End report <<<#

sink(paste0(repportName,".txt"), append = T)

print(sessionInfo())

sink()

#<<<>>>#