################################################################################
#                       QC batch test of the parameters                        #
################################################################################


## only once for installation
 myPaths <- .libPaths()   # get the paths
 #.libPaths("/home/rstudio/R")
 .libPaths("/work/PhosFake Program and results/PhosFake/lib/")
 
 library("BiocManager")
#BiocManager::install(c("protr","preprocessCore","matrixStats","extraDistr","fdrtool","qvalue","limma","moments","igraph","pheatmap","Biostrings", "corrr"))
 
 library(purrr)
 library(crayon)
 library(parallel)
 library(digest)
 library(tidyr)
 library(pheatmap)
 library(stringr)
 library(corrr)
 
#### YOU NEED TO BE IN THE MAIN FOLDER OF THE PHOSFAKE SCRIPTS

####################### Paths and directories
#####################
#Working directory should be PhosFake's main directory.
# File paths
Param <- NULL
fastaFilePath <- "InputData"
resultFilePath <- "KBN_output/KBN_ALBU_BOVIN"
# For parallel computing
cores <- 4
clusterType <- "FORK"
calcAllBenchmarks <- T

try(dir.create(resultFilePath))


repportName <- paste0(resultFilePath,"/ALBU_Report_BatchAnalysis_", gsub(":|[[:space:]]", "_", format(Sys.time(), "%Y%b%d%X")))

file.create(paste0(repportName,".txt"))
#####################

#####################
## Load sources
#####################

source("Parameter.R")
source("01_GenerateGroundTruth.R")
source("02_Digestion.R")
source("03_MSRun.R")
source("04_DataAnalysis.R")
source("05_Statistics.R")
source("06_Benchmarks.R")
source("viqor/Functions.R")

#####################

#####################
## Create list of testing parameters
#####################


# Fasta files:
pathToFasta <- list.files(path = fastaFilePath, full.names = T, pattern = ".fasta")

# Generate a list of parameters to test:
# Note: here we need to define ALL parameters as this will generate the necessary hash
paramGroundTruth <- list(#####   Ground truth generation
  "PathToFasta" = "QC_DataAnalysis/fasta_full_human.fasta",
  "PathToProteinList" = NULL,
  "NumReps" = c(5), #DEFAULT 3
  "NumCond" = 4,  #DEFAULT 2
  "FracModProt" = 0.2,  #DEFAULT 0.2
  "FracModPerProt" = 2,  #DEFAULT 2
  "PTMTypes" = "ph", # Phosphorylation
  "PTMTypesDist" = 1, # Frequency distribution of PTM
  "PTMTypesMass" = 79.6, # Mass difference after the PTM 
  "PTMMultipleLambda" = 0.5, #Make lower for next sim (0.05) # Lambda value for the poisson distribution used in adjusting probability wights of possible PTM locations
  "ModifiableResidues" = "S", # Which residues are modifiable 
  "ModifiableResiduesDistr" = 1, # Distribution of modifiable residues
  "RemoveNonModFormFrac" = 0 
)
paramProteoformAb <- list(
  "QuantNoise" = c(0.5),
  "DiffRegFrac" = c(0.25), 
  "DiffRegMax" = c(3), 
  "UserInputFoldChanges_NumRegProteoforms" = NULL, # rep(2, 2), #### not correctly implemented
  "UserInputFoldChanges_RegulationFC" = NULL , #rep(log2(c(100, 10)), rep(2, 2)), #### not correctly implemented
  "ThreshNAProteoform" = -100, # remove values below this threshold
  "AbsoluteQuanMean" = 30.5, #DEFAULT
  "AbsoluteQuanSD" = 3.6, #DEFAULT
  "ThreshNAQuantileProt" = 0.0 # Threshold of proteoforms to continue with
)
paramDigest <- list(
  ##### Digestion
  "Enzyme" = "trypsin",
  "PropMissedCleavages" = 0.05, 
  "MaxNumMissedCleavages" = 2,
  "PepMinLength" = 7, #DEFAULT
  "PepMaxLength" = 30, #DEFAULT
  "LeastAbundantLoss" = 0,
  "EnrichmentLoss" = 0.001,
  "EnrichmentEfficiency" = 0.2,
  "EnrichmentNonModSignalLoss" = 0.2,
  "EnrichmentNoise" = 0.5
)
paramMSRun <- list(
  ##### MSRun
  "PercDetectedPep" = 0.1,
  "PercDetectedVal" = 0.9,
  "WeightDetectVal" = 1,
  "MSNoise" = 0.2,
  "WrongIDs" = 0.05,
  "WrongLocalizations" = 0.05,
  "MaxNAPerPep" = 1000
)
paramDataAnalysis <- list(
  ##### Data analysis
  "ProtSummarization" = "medpolish",
  "MinUniquePep" = c(1),
  "StatPaired" = FALSE
)
#####################

#####################
## Run the sample preparation:
#####################

# all combinations of the parameters
#listtogroundtruth <- purrr::cross(compact(paramGroundTruth))
listtogroundtruth <- lapply(split(tidyr::expand_grid(!!!paramGroundTruth), seq(nrow(tidyr::expand_grid(!!!paramGroundTruth)))), function(x) as.list(x))


## setting up the entire thing hierarchically.
## Problem: how to parallelize

#Gather always benchmarking data
allBs <- NULL

for (hh in 1:length(listtogroundtruth)) {
  
  ## Running ground thruth generations
  # check whether file with correct parameters exists
  tParam <- listtogroundtruth[[hh]]
  Param <- "none"
  md5 <- digest(tParam, algo = "md5")
  filename <- paste0(resultFilePath,"/outputGroundTruth_",md5,".RData")
  if (file.exists(filename)) {
    load(filename)
  }
  if (all(unlist(Param) != unlist(tParam))) {
    Param <- tParam
    groundTruth <- samplePreparation(parameters = Param)
    save(Param, groundTruth, file = filename)
  }
  gtParam <- Param
  ## quantitative proteoform abundances
  #listtoproteoformab <- purrr::cross(compact(paramProteoformAb))
  listtoproteoformab <- lapply(split(tidyr::expand_grid(!!!paramProteoformAb), seq(nrow(tidyr::expand_grid(!!!paramProteoformAb)))), function(x) as.list(x))
  for (ii in 1:length(listtoproteoformab)) {
    Param <- "none"
    # create combined parameterfile
    tParam <- c(gtParam,listtoproteoformab[[ii]])
    tParam <- c(tParam, list(QuantColnames = paste0("C_",rep(1:tParam$NumCond,each=tParam$NumReps),"_R_", rep(1:tParam$NumReps, tParam$NumCond))))
    # hash code to represent parameter configuration
    md5 <- digest(tParam, algo = "md5")
    filename <- paste0(resultFilePath,"/outputProteoformAb_",md5,".RData")
    if (file.exists(filename)) {
      load(filename)
    }
    # test for exactly same parameters and run if not
    if (!all(unlist(Param) == unlist(tParam))) {
      Param <- tParam
      proteoformAb <- addProteoformAbundance(proteoforms = groundTruth, parameters = Param)
      save(Param, proteoformAb, file = filename)
    }
    pfParam <- Param
    
    ### Digestion
    #listtodigestion <- purrr::cross(compact(paramDigest))
    listtodigestion <- lapply(split(tidyr::expand_grid(!!!paramDigest), seq(nrow(tidyr::expand_grid(!!!paramDigest)))), function(x) as.list(x))
    for (jj in 1:length(listtodigestion)) {
      Param <- "none"
      tParam <- c(pfParam,listtodigestion[[jj]])
      md5 <- digest(tParam, algo = "md5")
      filename <- paste0(resultFilePath,"/outputDigestion_",md5,".RData")
      if (file.exists(filename)) {
        load(filename)
      }
      if (!all(unlist(Param) == unlist(tParam))) {
        Param <- tParam
        peptable <- digestGroundTruth(proteoforms = proteoformAb, parameters = c(Param, list(Cores=cores, ClusterType=clusterType)))
        peptable <- digestionProductSummarization(peptides = peptable, parameters = Param)
        BeforeMS <- filterDigestedProt(DigestedProt = peptable, parameters = Param)
        save(Param, BeforeMS, file = filename)
      }
      dgParam <- Param
      
      ### MS run
      #listtomsrun <- purrr::cross(compact(paramMSRun))
      listtomsrun <- lapply(split(tidyr::expand_grid(!!!paramMSRun), seq(nrow(tidyr::expand_grid(!!!paramMSRun)))), function(x) as.list(x))
      for (kk in 1:length(listtomsrun)) {
        Param <- "none"
        tParam <- c(dgParam,listtomsrun[[kk]])
        md5 <- digest(tParam, algo = "md5")
        filename <- paste0(resultFilePath,"/outputMSRun_",md5,".RData")
        if (file.exists(filename)) {
          load(filename)
        }
        if (!all(unlist(Param) == unlist(tParam))) {
          Param <- tParam
          AfterMSRun <- vector(mode = "list")
          for (i in which(sapply(BeforeMS, length) > 0)) {
            AfterMSRun[[length(AfterMSRun) + 1]] <- MSRunSim(Digested = BeforeMS[[i]], parameters = Param)
          }
          names(AfterMSRun) <- names(BeforeMS)[which(sapply(BeforeMS, length) > 0)]
          save(Param, AfterMSRun, file = filename)
        }
        msParam <- Param
        
        ### Protein abundance
        #listtodatanalysis <- purrr::cross(compact(paramDataAnalysis))
        listtodatanalysis <- lapply(split(tidyr::expand_grid(!!!paramDataAnalysis), seq(nrow(tidyr::expand_grid(!!!paramDataAnalysis)))), function(x) as.list(x))
        for (ll in 1:length(listtodatanalysis)) {
          Param <- "none"
          tParam <- c(msParam,listtodatanalysis[[ll]])
          md5 <- digest(tParam, algo = "md5")
          filename <- paste0(resultFilePath,"/outputDataAnalysis_",md5,".RData")
          if (file.exists(filename)) {
            load(filename)
          }
          if (!all(unlist(Param) == unlist(tParam))) {
            Param <- tParam
            Prots <- proteinSummarisation(peptable = AfterMSRun$NonEnriched, parameters = Param)
            # Don't accept anything below 100 proteins
            if(nrow(Prots) > 99) { 
              # Filter for having at least 1 actual value per protein group and peptide
              Prots <- Prots[rowSums(is.na(Prots[, Param$QuantColnames])) < length(Param$QuantColnames), ]
              allPeps <- as.data.frame(do.call("rbind", AfterMSRun))
              allPeps <- allPeps[rowSums(is.na(allPeps[, Param$QuantColnames])) < length(Param$QuantColnames), ]
              rownames(allPeps) <- paste0("pep", 1:nrow(allPeps))
              Stats <- runPolySTest(Prots, Param, refCond=1, onlyLIMMA=T)
              # much faster with only LIMMA tests   
              StatsPep <- runPolySTest(allPeps, Param, refCond=1, onlyLIMMA=T)
              Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
              save(Param, Stats, StatsPep, Benchmarks, file = filename)
            } else {
              print("Too few proteins!!!")
              Benchmarks <- Stats <- StatsPep <- NULL
              save(Param, Stats, StatsPep, Benchmarks, file = filename)
            }
          } else if (calcAllBenchmarks & !is.null(Stats)) {
            Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
            save(Param, Stats, StatsPep, Benchmarks, file = filename)
          }
          
          allBs[[md5]] <- list(Benchmarks, Param)
        }
      }
    }
  }
}
cat("###### Finished data set generation \n")

### This part can be used for visualizing and comparing
## Preferably calling an external script

# extracting all benchmarks (sometimes there are more or less per run)
t_allbnames <- NULL
for (i in names(allBs)) {
  t_allbnames <- c(t_allbnames, names(unlist(allBs[[i]][[1]]$globalBMs)))
}
benchNames <- unique(t_allbnames)
BenchMatrix <- data.frame(matrix(NA, ncol=length(benchNames)+length(Param)  , nrow=length(allBs)))
colnames(BenchMatrix) <- c(benchNames, names(Param))
rownames(BenchMatrix) <- names(allBs)

# writing all results and parameters into matrix
for (i in names(allBs)) {
  tglob <- unlist(allBs[[i]][[1]]$globalBMs)
  BenchMatrix[i, names(tglob)] <- tglob
  tpar <- allBs[[i]][[2]]
  BenchMatrix[i, names(tpar)] <- sapply(tpar, function(x) ifelse(length(x)>1, paste0(x,collapse="_"), x))
}

#### Simple visualizations of peptide correlation matrix analysis

afterMsTotal <- as.data.frame(rbind(AfterMSRun$NonEnriched,AfterMSRun$Enriched)) #main info

condAndRepStart <- which( colnames(afterMsTotal) == "C_1_R_1")
condAndRepEnd <- paramGroundTruth[["NumCond"]] * paramGroundTruth[["NumReps"]]+ (which( colnames(afterMsTotal) == "C_1_R_1")-1) #The number of columns to plot, taken from the number of copnditions and replicates 14 is the first column of simulated experiemnt
condTimesRep <- paramGroundTruth[["NumCond"]] * paramGroundTruth[["NumReps"]] #The total number of replicates and conditions. Used when collecting the data from afterMSTotal

AfterMSRunRegulated <- afterMsTotal %>% filter(!is.na(Regulation_Amplitude)) # AfterMSRun with the NA values filtered out
AfterMSRunRegulatedUnnest <- unnest(AfterMSRunRegulated, c(Peptide, Start, Stop, MC, Accession, Proteoform_ID, Regulation_Amplitude, Regulation_Pattern)) #using the tidyr unnest function get the data to plot every peptide, even those with multiple accession numbers to the same sequence 
AfterMSRunUnnest <- unnest(afterMsTotal, c(Peptide, Start, Stop, MC, Accession, Proteoform_ID, Regulation_Amplitude, Regulation_Pattern))

AccessionFreqTable <- as.data.frame(table(AfterMSRunUnnest$Accession))
colnames(AccessionFreqTable) <- c("Accession", "# of peptides")


mostSharedPeptide <- AccessionFreqTable[which.max(AccessionFreqTable$`# of peptides`),]

makeHeatmap = FALSE #is the data simple enough to create a heatmap of the experimental design? - less than 65k peptides

groupMaxSize <- 2

source("Visualizations.R")

proteinToAnalyzeCenteredMatplot <- VisualizeMs(afterMsTotal, groupMaxSize, condAndRepStart, condAndRepEnd, savePath = "KBN_output/", numCond = paramGroundTruth$NumCond, makeHeatmap = makeHeatmap) #MAKE A SINGLE FUNCTION THAT MAKES LINEAR PLOT + COR MATRIX AND HEATMAP FOR INPUT DATA
dev.off()




CorrelationGrouping <-  function(dfToAnalyze, peptideCor) {

  outStretch <- t(dfToAnalyze) %>% #the matrix is stretched so we can select the most correlated peptides
    correlate(quiet = TRUE) %>%
    stretch() 
  
  outAnalyze <- outStretch[outStretch$r>peptideCor,] %>% # removing NA values from the stretched corMatrix as well as only selecting the peptides with a correlation greater than peptideCor
    remove_missing() 
  
  tmpUnique <-  unique(outAnalyze$x) #get a vector of unique peptides to further group. 
  
  if(length(outAnalyze$r) != 0){
    pepToRemove <- outAnalyze[outAnalyze$r == max(outAnalyze$r),]$x[1]
    tmpUnique <- tmpUnique[!tmpUnique %in% pepToRemove]
    #outLoopAnalyze <- outAnalyze #make it so we can check for this inside the loop
    outLoopTest <-  remove_missing(outStretch) #value to test if a second selection has to be done
    counter <- 1
    
  } else {
    #corMatrixOut <- outLoop
    print("The loop has eliminated all peptides")
    #outLoopTest <- 0
    #outLoopTest$r <- 0
    #corMatrixOut <- NULL
    return(NULL)
  }
  
  while(length(outLoopTest$r ) > 1) {
    if ( length( outLoopTest$r ) == 1 )  {                #If the input already only has one pair
      
      print("The input data already only contains only one peptide pair which are correlated")
      corMatrixOut <- t(dfToAnalyze) %>%  # the input matrix is saved to the output
        correlate(quiet = TRUE)
      print("The result is saved to corMatrixOut")
      break
      
    } else if( min( outLoopTest$r ) > peptideCor && exists("outLoop")){   #IF the output of the algorithm contains only entries that are more correlated, then save it and stop running
        corMatrixOut <- outLoop
        
        #autoplot(outLoop)
        print("The result is saved to corMatrixOut")
        break
        
    } else if( min( outLoopTest$r ) < peptideCor & length(outLoopTest$r) > 1 ) {      #IF the length of the analyzed data is over 1, AND it does not only contain correlations above the threshold, then do another reduction of the correlation matrix
      
        tmpUniqueLoop <-  unique(outLoopTest$x) #get a vector of unique peptides to further group
        pepToRemoveLoop <- outLoopTest[outLoopTest$r == min((outLoopTest$r)),]$x[1] 
        tmpUniqueLoop <- tmpUniqueLoop[!tmpUniqueLoop %in% pepToRemoveLoop]
        
        outLoopSelection <-  dfToAnalyze[tmpUniqueLoop,] #select the data of the second correlation
        
        outLoop <- t(outLoopSelection) %>% #perform second cor matrix calculation
          correlate(quiet = TRUE) 
        
        
        outLoopStretch <- t(outLoopSelection) %>%  #stretch it to see if all correlations are above the starting threshold of correlation
          correlate(quiet = TRUE) %>%
          stretch() 
        
        
        outLoopTest <-  remove_missing(outLoopStretch) #value to test if a second selection has to be done
        corMatrixOut <- outLoop
        #outLoopAnalyze <- outLoopTest[outLoopTest$r>peptideCor,]# Next matrix to analyze
        
        #print(outLoopTest)
        print(paste(c("Loop", counter , "-"," The peptide which was removed was:", pepToRemoveLoop, "The current minimum correlation is:", min( outLoopTest$r ) ), collapse = " "))
        counter <- counter+1
    
    } else {
        print("The correlation calculation is done before the desired cutoff was achieved. The best pair was saved to corMatrixOut")
        corMatrixOut <- t(dfToAnalyze) %>%  # the input matrix is saved to the output
          correlate(quiet = TRUE)
        break
    }
  }
  
  return(corMatrixOut)
}
  

peptideCor <- 0.5 #How related should the peptide correlations be

##################################################
#### Peptide Correlation Matrix Analysis loop ####
##################################################

corMatrixToRemove <- proteinToAnalyzeCenteredMatplot #duplicating data to not run it again every time
iii <- 1 # iterator for correlation matrix classification
collectionOfCorrelations <- list() #initialization of list for collection of results

while ( length(corMatrixToRemove[,1]) > 2 ){ #keep looping while there is two elements or fewer to correlate
  
  
  pepRemoveFromTotalList <-  CorrelationGrouping(dfToAnalyze = corMatrixToRemove, peptideCor = peptideCor)
  if( is.null(pepRemoveFromTotalList) ){
    break
  }

  if ( length(pepRemoveFromTotalList$term) > 1 ) {
    collectionOfCorrelations[[iii]] <- pepRemoveFromTotalList
    iii <- iii+1
  }

    
  corMatrixToRemove <- corMatrixToRemove[!rownames(corMatrixToRemove) %in% pepRemoveFromTotalList$term ,]
} # loop to classify every peptide of the protein

autoplot(collectionOfCorrelations[[1]], triangular = "full", low = "navy", mid = "white", high = "red") # use to plot the different cor_df found
Reduce(intersect , afterMsTotal[collectionOfCorrelations[[1]][["term"]],]$Proteoform_ID)
afterMsTotal[collectionOfCorrelations[[1]][["term"]],]$Proteoform_ID


#### correct % of correlation matrices ####
if (length ( collectionOfCorrelations) != 0 ){ #calculating the percentage of peptides with a single group
  correctCount <- 0 #counter for correct classification percentage
  for(i in 1: length(collectionOfCorrelations)) {
    
    if (length(Reduce(intersect , afterMsTotal[collectionOfCorrelations[[i]][["term"]],]$Proteoform_ID) != 0)) {
      correctCount <- correctCount+1
    }
    
  }
  
  correctPercentage <-(correctCount / length ( collectionOfCorrelations)) *100 
  print( paste( c ("The percentage of correctly assigned peptide groups is ",correctPercentage)))
} else {
  correctPercentage <- NULL
  print("Your collection of correlations contains 0 correlations - cannot compute")
  }

