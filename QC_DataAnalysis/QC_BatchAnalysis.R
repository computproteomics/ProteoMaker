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
                    "PepMaxLength" = seq(from = 15, to = 35, by = 5))

#--------------------
# Name of report:

repportName <- paste0("QC_output/QC_Report_BatchAnalysis_",
                      gsub(":", "", format(Sys.time(), "%Y%b%d%X")))

#####################

#####################
## Run the sample preparation:
#####################

source("../01_GenerateGroundTruth.R")

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

#>>> Start writing report <<<#

sink(paste0(repportName,".txt"))

cat("Generation of", length(listtotest), "parameter sets for digestion\n")
cat("They can be found in the following table:", paste0(repportName, "_TestedParam.txt"), "\n")
tmp <- bind_rows(listtotest)
tmp$OutputNumber <- 1:nrow(tmp)
write.table(tmp, 
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
    
    iter <- iter + 1
    
  }
  
}


#####################

#>>> End report <<<#

sink(paste0(repportName,".txt"), append = T)

print(sessionInfo())

sink()

#<<<>>>#