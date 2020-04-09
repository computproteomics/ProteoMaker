################################################################################
#                       QC batch test of the parameters                        #
################################################################################

library(purrr)
library(crayon)
library(parallel)

#####################
## Paths and directories
#####################
#Working directory should be PhosFake's main directory.

if(sum(list.files("QC_DataAnalysis") == "QC_output") != 1){
  
  dir.create("QC_DataAnalysis/QC_output")
  
}

repportName <- paste0("QC_DataAnalysis/QC_output/QC_Report_BatchAnalysis_", gsub(":|[[:space:]]", "_", format(Sys.time(), "%Y%b%d%X")))

file.create(paste0(repportName,".txt"))
#####################

#####################
## Load sources
#####################

source("Parameter.R")
source("01_GenerateGroundTruth.R")
source("02_Digestion.R")
source("03_MSRun.R")
#####################

#####################
## Create list of testing parameters
#####################

# Fasta files:
pathToFasta <- list.files(path = "QC_DataAnalysis", full.names = T, pattern = ".fasta")

# Generate a list of parameters to test:
paramToTest <- list("PathToFasta" = pathToFasta,
                    "PropMissedCleavages" = c(0.01, seq(from = 0.1, to = 0.8, by = 0.1), 0.99), 
                    "MaxNumMissedCleavages" = 4,
                    # "PepMinLength" = c(6,7),
                    # "PepMaxLength" = seq(from = 15, to = 35, by = 5),
                    "LeastAbundantLoss" = c(seq(from = 0, to = 0.6, by = 0.2), seq(from = 0.7, to = 0.9, by = 0.05)))
#####################

#####################
## Run the sample preparation:
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
#####################

#####################
## Digestion according to each set of parameters:
#####################
listtotest <- purrr::cross(paramToTest)

# Generate table with parameters and results:

res <- dplyr::bind_rows(listtotest)
listtotest <- listtotest[order(res$PathToFasta)]
res <- res[order(res$PathToFasta),]
res$SpeciesID <- unlist(lapply(1:length(pathToFasta), function(x) rep(x, length(listtotest)/length(pathToFasta))))
res <- res[, c(1, 7, 2, 3, 4, 5, 6)]
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
cat("They can be found in the following table:", paste0("QC_DataAnalysis/",repportName, "_TestedParam.txt"), "\n")
write.table(res, 
            file = paste0(repportName, "_TestedParam.txt"), 
            sep = "\t",
            row.names = F)

sink()

#<<<>>>#

cat("Start generation of", length(listtotest), "parameter sets for digestion\n\n")

#Does benchmark for each row of res and the corresponding parameters in listtotest.
RunBenchmark <- function(GroundTruth, parameters, listtotest, id){
  
  parametersLoop <- c(parameters[!(names(parameters) %in% names(listtotest))], listtotest)
  
  peptable <- digestGroundTruth(proteoforms = GroundTruth, parameters = parametersLoop)
  peptable <- digestionProductSummarization(peptides = peptable, parameters = parametersLoop)
  
  #Metrics that we want to retreive:
  metrics <- list()
  
  #Peptides:
  #Number of unique peptides (including modification):
  metrics$NumberUniquePeptide <- nrow(peptable)
  #Number of unique peptides (including modification) with no miss-cleavage:
  metrics$NumberUniquePeptide0MC <- sum(unlist(lapply(peptable$MC, function(x) x[1])) == 0)
  #Number of unique peptides (including modification) with 1 miss-cleavage:
  metrics$NumberUniquePeptide1MC <- sum(unlist(lapply(peptable$MC, function(x) x[1])) == 1)
  #Number of unique peptides (including modification) with 2 miss-cleavages:
  metrics$NumberUniquePeptide2MC <- sum(unlist(lapply(peptable$MC, function(x) x[1])) == 2)
  #Number of unique peptides (including modification) with 3 miss-cleavages:
  metrics$NumberUniquePeptide3MC <- sum(unlist(lapply(peptable$MC, function(x) x[1])) == 3)
  #Number of unique peptides (including modification) with 4 miss-cleavages:
  metrics$NumberUniquePeptide4MC <- sum(unlist(lapply(peptable$MC, function(x) x[1])) == 4)
  
  pepLowerThanMedian <- which(rowMeans(peptable[,parametersLoop$QuantColnames]) <= median(unlist(peptable[,parametersLoop$QuantColnames]), na.rm = T))
  
  #Number of unique peptides with lower abundance than the median (including modification) with no miss-cleavage:
  metrics$NumberUniquePeptideMedian0MC <- sum(unlist(lapply(peptable$MC[pepLowerThanMedian], function(x) x[1])) == 0)
  #Number of unique peptides with lower abundance than the median (including modification) with 1 miss-cleavage:
  metrics$NumberUniquePeptideMedian1MC <- sum(unlist(lapply(peptable$MC[pepLowerThanMedian], function(x) x[1])) == 1)
  #Number of unique peptides with lower abundance than the median (including modification) with 2 miss-cleavages:
  metrics$NumberUniquePeptideMedian2MC <- sum(unlist(lapply(peptable$MC[pepLowerThanMedian], function(x) x[1])) == 2)
  #Number of unique peptides with lower abundance than the median (including modification) with 3 miss-cleavages:
  metrics$NumberUniquePeptideMedian3MC <- sum(unlist(lapply(peptable$MC[pepLowerThanMedian], function(x) x[1])) == 3)
  #Number of unique peptides with lower abundance than the median (including modification) with 4 miss-cleavages:
  metrics$NumberUniquePeptideMedian4MC <- sum(unlist(lapply(peptable$MC[pepLowerThanMedian], function(x) x[1])) == 4)
  
  #Number of peptides being N or C terminus
  fasta <- protr::readFASTA(file = parametersLoop$PathToFasta, legacy.mode = TRUE, seqonly = FALSE)
  last <- nchar(unlist(fasta))
  names(last) <- sub(".*[|]([^.]+)[|].*", "\\1", names(fasta))
  rm(fasta)
  
  last.mapping <- last[unlist(sapply(peptable$Accession, function(x) x[1]))]
  metrics$CTerminusPeptides <- sum(unlist(sapply(peptable$Stop, function(x) x[1])) == last.mapping)
  metrics$NTerminusPeptides <- sum(unlist(sapply(peptable$Start, function(x) x[1])) == 1)

  # Intensities of the peptides after summarisation:
  metrics$PeptideIntensities$Values <- as.matrix(peptable[,Param$QuantColnames])
  # Matching number of missed-cleavages:
  metrics$PeptideIntensities$numMC <- sapply(peptable$MC, unique)
  
  #Proteins
  vec_numProtein <- sapply(peptable$Accession, function(x) {length(unique(x))})
  vec_numProteoform <- sapply(peptable$Proteoform_ID, function(x) {length(unique(x))})
  
  # Number of unique proteins with min. 1 proteotypic peptide:
  metrics$NumberUniqueProtein1 <- sum(vec_numProtein == 1)
  # Number of unique proteins with min. 2 proteotypic peptides:
  metrics$NumberUniqueProtein2 <- sum(vec_numProtein >= 2)
  # Number of proteins a peptide is matched to:
  metrics$NumProteinPerPep <- vec_numProtein
  
  #Proteoforms
  # Number of unique proteoforms with  min. 1 peptide specific to the proteoform:
  metrics$NumberUniqueProteoforms1 <- sum(vec_numProteoform == 1)
  # Number of unique proteoforms with  min. 1 peptide specific to the proteoform:
  metrics$NumberUniqueProteoforms2 <- sum(vec_numProteoform >= 2)
  # Number of proteoforms a peptide is matched to:
  metrics$NumProteoformPerPep <- vec_numProteoform
  
  save(metrics, file = paste0(repportName, "_", id, ".RData"))
  
  cat("\n\n")
  cat(crayon::red("Save output", id, "over", length(listtotest), "\n\n"))
  
}

#Set up parallel computing environment

#Set number of clusters, PSOCK in windows and linux, FORK for linux (faster)
cluster <- parallel::makeCluster(30, type = "FORK")
parallel::setDefaultCluster(cluster)
parallel::clusterEvalQ(cluster, library(dplyr))
parallel::clusterExport(cluster, c("repportName", "Param", "GroundTruth", "listtotest", "res", "RunBenchmark", "digestGroundTruth", "proteoformDigestion", "fastDigest", "digestionProductSummarization"), envir = environment())
peptides <- parallel::parLapply(cluster, 1:nrow(res), function(x) RunBenchmark(GroundTruth = GroundTruth[[res$SpeciesID[x] ]], listtotest = listtotest[[x]], parameters = Param, id = res$OutputNumber[x]) )
parallel::stopCluster(cluster)

#####################

#>>> End report <<<#

sink(paste0(repportName,".txt"), append = T)

print(sessionInfo())

sink()

#<<<>>>#





