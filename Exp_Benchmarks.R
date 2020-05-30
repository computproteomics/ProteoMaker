library(moments)
library(stringr)

############# Calculating ROC curves
calcROC <- function (Stats, columnName, groundtruthColumn="min1Reg") {
  FPs <- TPs <- 0
  FNs <- sum(Stats$min1Reg)
  TNs <- nrow(Stats) - FNs
  ROC <- NULL
  FDR <- NULL
  AUC <- 0
  oldFDR <- -1
  oldFPR <- 0
  for (i in order(Stats[,columnName])) {
    if(Stats[i, groundtruthColumn]) {
      TPs <- TPs + 1
      FNs <- FNs -1
    } else {
      FPs <- FPs + 1
      TNs <- TNs - 1
    }
    currFDR <- Stats[i, columnName]  
    if (is.na(currFDR)) currFDR <- 1
    if (currFDR != oldFDR) {
      ROC <- rbind(ROC, c(FPs/(FPs+TNs), TPs/(TPs+FNs)))
      FDR <- rbind(FDR, c(currFDR, FPs/(FPs + TPs)))
      AUC <- AUC + (FPs/(FPs+TNs) - oldFPR) *  TPs/(TPs+FNs)
    } else {
      ROC[nrow(ROC), ] <- c(FPs/(FPs+TNs), TPs/(TPs+FNs))
      FDR[nrow(FDR), ] <- c(currFDR, FPs/(FPs + TPs))    
      AUC <- AUC + (FPs/(FPs+TNs) - oldFPR) *  TPs/(TPs+FNs)
    }
    
    oldFDR <- currFDR
    oldFPR <- FPs/(FPs+TNs)
  }
  colnames(ROC) <- c("FPR", "TPR")
  colnames(FDR) <- c("FDR", "tFDR")
  cbind(ROC, FDR, AUC)
}

# wrapper for calculating all metrics
calcBenchmarks <- function(Stats, StatsPep, Param)  {
  
  Benchmarks <- NULL
  
  # global means 1 number per metric
  globalBMs <- list(
    # Peptide level
    numPeptides=0, numProteins=0, propUniquePep=0, propSharedPep=0, percMissingPep=0,
    aucDiffRegPeptides=list(), tprPep0.01=list(), tprPep0.05=list(), tFDRPep0.01=list(), tFDRPep0.05=list(), 
    propMisCleavedPeps=list(),skewnessPeps=0, kurtosisPeps=0, sdPeps=0,
    # Protein level
    numQuantProtGroups=0, propUniqueProts=0, percMissingProt=0, meanPepPerProt=0, aucDiffRegProteins=list(), 
    tFDRProt0.01=list(), tFDRProt0.05=list(), tprProt0.01=list(), tprProt0.05=list(), sumSquareDiffFC=0, propMisCleavedProts=0,
    propDiffRegWrongIDProt0.01=list(),propDiffRegWrongIDProt0.05=list(),skewnessProts=0, kurtosisProts=0, sdProts=0,
    # PTM level
    numProteoforms=0, numModPeptides=0, meanProteoformsPerProt=0, propModAndUnmodPep=0, aucDiffRegAdjModPep=list(),
    tFDRAdjModPep0.01=list(), tFDRAdjModPep0.05=list(), tprAdjModPep0.01=list(), tprAdjModPep0.05=list(),
    propDiffRegPepWrong0.01=list(),propDiffRegPepWrong0.05=list(), percOverlapModPepProt=0)  
  
  
  #### Calculating peptide numbers
  globalBMs["numPeptides"] <- nrow(StatsPep)
  globalBMs["numProteins"] <- length(unique(unlist(StatsPep$Accession)))
  globalBMs["propUniquePep"] <-  sum(sapply(StatsPep$Accession, function(x) length(x) == 1)) / nrow(StatsPep)
  globalBMs["propSharedPep"] <-  sum(sapply(StatsPep$Accession, function(x) length(x) > 1)) / nrow(StatsPep)
  globalBMs["percMissingPep"] <- sum(is.na(unlist(StatsPep[,Param$QuantColnames]))) / length(unlist(StatsPep[,Param$QuantColnames])) * 100
  
  # Which tests are there?
  statCols <- grep("FDR",colnames(StatsPep), value=T)
  
  # results on basis of ground truth
  ROC <- list()
  plot(0:1, 0:1, type="n", main="Peptides")
  for (test in statCols) {
    print(test)
    testSum <-  calcROC(StatsPep, test)
    if (nrow(testSum) > 1) {
      lines(testSum[,1], testSum[,2], type="s", col=which(test==statCols))
      lines(testSum[,3], testSum[,4], type="l", col=which(test==statCols),lty=3)
      ROC[[test]] <- testSum
      at.01 <- which.min(abs(testSum[,"FDR"] - 0.01))
      at.05 <- which.min(abs(testSum[,"FDR"] - 0.05))
      
      globalBMs$aucDiffRegPeptides[[test]] <- testSum[1,"AUC"]
      globalBMs$tprPep0.01[[test]] <- testSum[at.01, "TPR"]
      globalBMs$tprPep0.05[[test]] <- testSum[at.05, "TPR"]
      globalBMs$tFDRPep0.01[[test]] <- testSum[at.01, "tFDR"]
      globalBMs$tFDRPep0.05[[test]] <- testSum[at.05, "tFDR"]
    }
  }
  legend("bottomright", lwd=1, col=1:length(statCols), legend = statCols, cex=0.6)
  Benchmarks$PepStat <- ROC
  
  # miscleavages
  globalBMs["propMisCleavedPeps"] <- list(table(sapply(StatsPep$MC, function(x) x[1])) / nrow(StatsPep))
  
  #### Calculating protein numbers
  globalBMs["numQuantProtGroups"] <- nrow(Stats)
  
  globalBMs["propUniqueProts"] <- sum(unlist(sapply(str_split(Stats$num_accs,";"), function(x) unique(as.numeric(unlist(x))) == 1))) / nrow(Stats)
  globalBMs["percMissingProt"] <- sum(is.na(as.vector(Stats[,Param$QuantColnames])))  / length(Param$QuantColnames) / nrow(Stats) * 100
  pepDistr <- sapply(str_split(Stats$Sequence,";"), function(x) length(unique(x)))
  barplot(table(pepDistr), ylab="Frequency", xlab="Peptides per protein")
  globalBMs["meanPepPerProt"] <-  mean(pepDistr)
  
  # results on basis of ground truth
  ROC <- list()
  plot(0:1, 0:1, type="n", main="Proteins", xlim=c(0,1), ylim=c(0,1))
  for (test in statCols) {
    print(test)
    testSum <-  calcROC(Stats, test)
    if (nrow(testSum) > 1) {
      lines(testSum[,1], testSum[,2], type="s", col=which(test==statCols))
      lines(testSum[,3], testSum[,4], type="l", col=which(test==statCols),lty=3)
      ROC[[test]] <- testSum
      at.01 <- which.min(abs(testSum[,"FDR"] - 0.01))
      at.05 <- which.min(abs(testSum[,"FDR"] - 0.05))
      
      globalBMs$aucDiffRegProteins[[test]] <- testSum[1,"AUC"]
      globalBMs$tprProt0.01[[test]] <- testSum[at.01, "TPR"]
      globalBMs$tprProt0.05[[test]] <- testSum[at.05, "TPR"]
      globalBMs$tFDRProt0.01[[test]] <- testSum[at.01, "tFDR"]
      globalBMs$tFDRProt0.05[[test]] <- testSum[at.05, "tFDR"]
    }
  }
  
  legend("bottomright", lwd=1, col=1:length(statCols), legend = statCols, cex=0.6)
  Benchmarks$ProtStat <- ROC
  
  ## Calculating differences between actual and "measured" fold-changes
  patterns <- lapply(Stats$Regulation_Pattern, function(x)  matrix(as.numeric(unlist(str_split(x, ";"))), byrow=T, ncol = Param$NumCond))
  amplitudes <- lapply(Stats$Regulation_Amplitude, function(x) as.numeric(unlist(str_split(x, ";"))))
  diffs <- vector("numeric",length(patterns))
  sumsquare <- 0
  for (i in 1:length(patterns)) {
    tampl <- na.omit(amplitudes[i][[1]])
    if (length(tampl)> 0) {
      tval <- patterns[i][[1]] * tampl
      if(length(tval) > 2) {
        tval <- colMeans(tval[,2:ncol(tval), drop=F] - tval[,1])
        diffs[i] <- tval
        sumsquare <- sumsquare + (tval - Stats$`log-ratios 2 vs 1`[i]) * (tval - Stats$`log-ratios 2 vs 1`[i])
      } else {
        diffs[i] <- 0
      }
    } else {
      diffs[i] <- 0 
    }
  }
  
  plot(0, xlim=range(Stats$`log-ratios 2 vs 1`, na.rm=T), ylim=range(Stats$`log-ratios 2 vs 1`, na.rm=T), type="n", xlab="Ground truth", ylab="Measured ratios")
  points(diffs, Stats$`log-ratios 2 vs 1`, pch=15, cex=0.7, col="#00000055")
  abline(0,1)
  globalBMs["sumSquareDiffFC"] <- sumsquare / sum(diffs != 0)
  
  # Counting miscleavages
  globalBMs["propMisCleavedProts"] <- sum(sapply(Stats$MC, function(x) sum(as.numeric(unlist(strsplit(x, ";"))))) > 0) / nrow(Stats)
  
  # statistics with respect to wrong identifications
  wrong_ids <- sapply(Stats$WrongID, function(x) sum(as.logical(unlist(strsplit(x, ";")))))
  for (test in statCols) {
    globalBMs$propDiffRegWrongIDProt0.01[[test]] <- sum(Stats[,test] < 0.01 & wrong_ids > 0, na.rm=T) / sum(Stats[,test] < 0.01, na.rm=T)
    globalBMs$propDiffRegWrongIDProt0.05[[test]] <- sum(Stats[,test] < 0.05 & wrong_ids > 0, na.rm=T) / sum(Stats[,test] < 0.05, na.rm=T)
  }
  
  # checking quantitative values for assymetric distribution: skewness
  globalBMs$skewnessProts <- skewness(unlist(Stats[,Param$QuantColnames]), na.rm=T)
  globalBMs$kurtosisProts <- kurtosis(unlist(Stats[,Param$QuantColnames]), na.rm=T)
  globalBMs$sdProts <- sd(unlist(Stats[,Param$QuantColnames]), na.rm=T)
  globalBMs$skewnessPeps <- skewness(unlist(StatsPep[,Param$QuantColnames]), na.rm=T)
  globalBMs$kurtosisPeps <- kurtosis(unlist(StatsPep[,Param$QuantColnames]), na.rm=T)
  globalBMs$sdPeps <- sd(unlist(StatsPep[,Param$QuantColnames]), na.rm=T)
  
  ###### metrics on PTM level
  # number of proteoforms per protein group and in total, not all identifiable
  ProteoformDistr <- sapply(Stats$Proteoform_ID, function(x) length(unique(as.numeric(unlist(strsplit(x, ";"))))))
  barplot(table(ProteoformDistr), 100)
  globalBMs$numProteoforms <- sum(ProteoformDistr)
  globalBMs$meanProteoformsPerProt <- mean(ProteoformDistr)
  
  globalBMs$numModPeptides <- sum(StatsPep$PTMType != "NULL")
  ModPeps <- StatsPep[StatsPep$PTMType != "NULL",]
  # Proportion of modified peptides with identical non-modifiedpeptide
  pepgroups <- by(StatsPep[,c("Sequence", "Accession", "PTMType", "PTMPos")], StatsPep$Sequence, function(x) x )
  pepgroups <- lapply(pepgroups, function(x) {x[x=="NULL"] <- NA; x})  
  modpepgroups <- sapply(pepgroups, function(x) {sum(as.numeric(unlist(x[,"PTMPos"])),na.rm=T) > 0})
  modpepgroups <- pepgroups[modpepgroups]
  if (length(modpepgroups) > 0) {
    modunmodgroups <- sapply(modpepgroups, function(x) sum(is.na(x[,"PTMPos"])) > 0)
    if (length(modunmodgroups) > 0) {
      modunmodgroups <- modpepgroups[modunmodgroups]
      globalBMs$propModAndUnmodPep <- sum(sapply(modunmodgroups, function(x) nrow(x)-1)) / globalBMs$numModPeptides
    }
  }
  
  # modified peptides and their proteins
  ModPeps <- cbind(ModPeps, merged_accs=sapply(ModPeps$Accession, function(x) paste(unlist(x),collapse=";")))
  ModPepsWithProt <- ModPeps[ModPeps$merged_accs  %in% rownames(Stats), ]
  globalBMs$percOverlapModPepProt <- nrow(ModPepsWithProt) / nrow(ModPeps) * 100
  
  # Adjust by protein expression (could be moved to Statistics)
  AdjModPepsWithProt <- ModPepsWithProt
  AdjModPepsWithProt[,Param$QuantColnames] <- AdjModPepsWithProt[,Param$QuantColnames] - Stats[as.character(AdjModPepsWithProt$merged_accs), Param$QuantColnames]
  StatsAdjModPep <- 0
  if (nrow(AdjModPepsWithProt) > 0) {
    StatsAdjModPep <- runPolySTest(AdjModPepsWithProt, Param, refCond=1, onlyLIMMA=F)
    
    # results on basis of ground truth
    ROC <- list()
    plot(0:1, 0:1, type="n", main="Adj. modified peptides")
    for (test in statCols) {
      print(test)
      testSum <-  calcROC(StatsAdjModPep, test)
      if (nrow(testSum) > 1) {
        lines(testSum[,1], testSum[,2], type="s", col=which(test==statCols))
        lines(testSum[,3], testSum[,4], type="l", col=which(test==statCols),lty=3)
        ROC[[test]] <- testSum
        at.01 <- which.min(abs(testSum[,"FDR"] - 0.01))
        at.05 <- which.min(abs(testSum[,"FDR"] - 0.05))
        
        globalBMs$aucDiffRegAdjModPep[[test]] <- testSum[1,"AUC"]
        globalBMs$tprAdjModPep0.01[[test]] <- testSum[at.01, "TPR"]
        globalBMs$tprAdjModPep0.05[[test]] <- testSum[at.05, "TPR"]
        globalBMs$tFDRAdjModPep0.01[[test]] <- testSum[at.01, "tFDR"]
        globalBMs$tFDRAdjModPep0.05[[test]] <- testSum[at.05, "tFDR"]
      }
    }
    legend("bottomright", lwd=1, col=1:length(statCols), legend = statCols, cex=0.6)
    Benchmarks$AdjModPepStat <- ROC
  }
  
  # Back to original modified peptides, counting the wrong differential regulations
  for (test in statCols) {
    globalBMs$propDiffRegPepWrong0.01[[test]] <- sum(ModPeps[[test]] < 0.01, na.rm=T) / nrow(ModPeps)
    globalBMs$propDiffRegPepWrong0.05[[test]] <- sum(ModPeps[[test]] < 0.05, na.rm=T) / nrow(ModPeps)
  }
  
  
  Benchmarks$globalBMs <- globalBMs
  return(Benchmarks)
  
}

# wrapper for calculating all metrics
calcBasicBenchmarks <- function(Stats, StatsPep, Param)  {
  
  Benchmarks <- NULL
  
  # global means 1 number per metric
  globalBMs <- list(
    # Peptide level
    numPeptides=0, numProteins=0, propUniquePep=0, propSharedPep=0, percMissingPep=0,
    propMisCleavedPeps=list(),skewnessPeps=0, kurtosisPeps=0, sdPeps=0,
    # Protein level
    numQuantProtGroups=0, propUniqueProts=0, percMissingProt=0, meanPepPerProt=0, 
    propMisCleavedProts=0, skewnessProts=0,kurtosisProts=0, sdProts=0,
    # PTM level
    numProteoforms=0, numModPeptides=0, meanProteoformsPerProt=0, propModAndUnmodPep=0, percOverlapModPepProt=0)  
  
  
  #### Calculating peptide numbers
  globalBMs["numPeptides"] <- nrow(StatsPep)
  globalBMs["numProteins"] <- length(unique(unlist(StatsPep$Accession)))
  globalBMs["propUniquePep"] <-  sum(sapply(StatsPep$Accession, function(x) length(x) == 1)) / nrow(StatsPep)
  globalBMs["propSharedPep"] <-  sum(sapply(StatsPep$Accession, function(x) length(x) > 1)) / nrow(StatsPep)
  globalBMs["percMissingPep"] <- sum(is.na(unlist(StatsPep[,Param$QuantColnames]))) / length(unlist(StatsPep[,Param$QuantColnames])) * 100
  
  # miscleavages
  globalBMs["propMisCleavedPeps"] <- list(table(sapply(StatsPep$MC, function(x) x[1])) / nrow(StatsPep))
  
  #### Calculating protein numbers
  globalBMs["numQuantProtGroups"] <- nrow(Stats)
  
  globalBMs["propUniqueProts"] <- sum(unlist(sapply(str_split(Stats$num_accs,";"), function(x) unique(as.numeric(unlist(x))) == 1))) / nrow(Stats)
  globalBMs["percMissingProt"] <- sum(is.na(as.vector(Stats[,Param$QuantColnames])))  / length(Param$QuantColnames) / nrow(Stats) * 100
  pepDistr <- sapply(str_split(Stats$Sequence,";"), function(x) length(unique(x)))
  barplot(table(pepDistr), ylab="Frequency", xlab="Peptides per protein")
  globalBMs["meanPepPerProt"] <-  mean(pepDistr)
  
  # Counting miscleavages
  globalBMs["propMisCleavedProts"] <- sum(sapply(Stats$MC, function(x) sum(as.numeric(unlist(strsplit(x, ";"))))) > 0, na.rm=T) / nrow(Stats)
  
  
  # checking quantitative values for assymetric distribution: skewness
  globalBMs$skewnessProts <- skewness(unlist(Stats[,Param$QuantColnames]), na.rm=T)
  globalBMs$kurtosisProts <- kurtosis(unlist(Stats[,Param$QuantColnames]), na.rm=T)
  globalBMs$sdProts <- sd(unlist(Stats[,Param$QuantColnames]), na.rm=T)
  globalBMs$skewnessPeps <- skewness(unlist(StatsPep[,Param$QuantColnames]), na.rm=T)
  globalBMs$kurtosisPeps <- kurtosis(unlist(StatsPep[,Param$QuantColnames]), na.rm=T)
  globalBMs$sdPeps <- sd(unlist(StatsPep[,Param$QuantColnames]), na.rm=T)
  
  ###### metrics on PTM level
  
  globalBMs$numModPeptides <- sum(StatsPep$PTMType != "NULL")
  ModPeps <- StatsPep[StatsPep$PTMType != "NULL",]
  # Proportion of modified peptides with identical non-modifiedpeptide
  pepgroups <- by(StatsPep[,c("Sequence", "Accession", "PTMType", "PTMPos")], StatsPep$Sequence, function(x) x )
  pepgroups <- lapply(pepgroups, function(x) {x[x=="NULL"] <- NA; x})  
  modpepgroups <- sapply(pepgroups, function(x) {sum(as.numeric(unlist(x[,"PTMPos"])),na.rm=T) > 0})
  modpepgroups <- pepgroups[modpepgroups]
  if (length(modpepgroups) > 0) {
    modunmodgroups <- sapply(modpepgroups, function(x) sum(is.na(x[,"PTMPos"])) > 0)
    if (length(modunmodgroups) > 0) {
      modunmodgroups <- modpepgroups[modunmodgroups]
      globalBMs$propModAndUnmodPep <- sum(sapply(modunmodgroups, function(x) nrow(x)-1)) / globalBMs$numModPeptides
    }
  }
  
  # modified peptides and their proteins
  ModPeps <- cbind(ModPeps, merged_accs=sapply(ModPeps$Accession, function(x) paste(unlist(x),collapse=";")))
  ModPepsWithProt <- ModPeps[ModPeps$merged_accs  %in% rownames(Stats), ]
  globalBMs$percOverlapModPepProt <- nrow(ModPepsWithProt) / nrow(ModPeps) * 100
  
  
  Benchmarks$globalBMs <- globalBMs
  return(Benchmarks)
  
}


readMaxQuant <- function(allPeps, Prots, Param=NULL) {
  allPeps$Accession <- sapply(allPeps$Proteins, function(x) strsplit(as.character(x), ";"))
  allPeps$MC <- allPeps$Missed.cleavages
  allPeps$PTMType <- as.character(allPeps$Modifications)
  allPeps$PTMType[allPeps$PTMType == "Unmodified"] <- "NULL"
  # did not find corresponding field    
  allPeps$PTMPos <- NA
  allPeps$PTMPos[allPeps$PTMType != "NULL"] <- 1
  
  Param$QuantColnames <- grep("Intensity\\.", names(allPeps), value=T)
  
  # TODO: check whether LFQ results are available 
  #protCols <- grepl("LFQ.intensity", names(Prots))
  #names(Prots)[protCols] <- Param$QuantColnames
  Prots$Accession <- Prots$Sequence <- Prots$Protein.IDs
  Prots$num_accs <- Prots$Proteins
  
  # filter for rows with no quantifications
  tquant <- allPeps[,Param$QuantColnames]
  tquant[tquant == 0] <- NA
  allPeps[, Param$QuantColnames] <- log2(tquant)
  allPeps <- allPeps[rowSums(is.na(allPeps[, Param$QuantColnames])) < length(Param$QuantColnames), ]
  tquant <- Prots[,Param$QuantColnames]
  allPeps$Sequence <- as.character(allPeps$Sequence)
  tquant[tquant == 0] <- NA
  Prots[, Param$QuantColnames] <- log2(tquant)
  Prots <- Prots[rowSums(is.na(Prots[, Param$QuantColnames])) < length(Param$QuantColnames), ]
  rownames(Prots) <- Prots$Accession
  # add column with miscleavages
  Prots$MC <- NA
  if (!is.null(allPeps$MC)) {
    mergedMCs <- unlist(by(allPeps$Missed.cleavages, as.character(allPeps$Proteins), function(x) paste(x,collapse=";")))
    Prots[names(mergedMCs), "MC"] <- mergedMCs
  } else {
    allPeps$MC <- as.character(0)
    Prots$MC <- as.character(0)
  }
  
  
  # TODO define parameters for statistics
  
  
  
  return(list(Param=Param, allPeps=allPeps, Prots=Prots))
}


############ Call the function to load and prepare data, and then to calculate the benchmarks
### It is now adopted to download files from PRIDE according to the metadata
library(jsonlite)
## General parameters
withGroundTruth <- F
tmpdir <- "/data/tmp"
allPRIDE <- readLines("../ExpBench/pxd_accessions.json")


for (pxd in 1:length(allPRIDE))  {
  
  tdat <- fromJSON(allPRIDE[[pxd]])
  modpep <- protlist <- NULL
  if (!is.null(tdat$dataProcessingProtocol) & !is.null(tdat$files)) {
    # if maxquant
    if (grepl("maxquant",tolower(tdat$dataProcessingProtocol)) ) {
      # if maxquant files available
          if (tdat$files$list$fileType == "SEARCH" & any(grepl("/modificationSpecificPeptides.txt",tdat$files$list$downloadLink)) &
                      any(grepl("/proteinGroups.txt",tdat$files$list$downloadLink))) {
            modpep <- read.csv(grep("/modificationSpecificPeptides.txt",tdat$files$list$downloadLink, value=T), sep="\t")
            protlist <- read.csv(grep("/proteinGroups.txt",tdat$files$list$downloadLink, value=T), sep="\t")
        # if maxquant files in zip file
      } else if (any(tdat$files$list$fileType == "SEARCH" & grepl(".zip",tdat$files$list$fileName) )) {
        filenames <- tdat$files$list$downloadLink[which(tdat$files$list$fileType == "SEARCH" & grepl(".zip",tdat$files$list$fileName))]
        # non-optimal: if more than one zip file, then take the first
        for (filename in filenames) {
          tempfile <- tempfile(tmpdir=tmpdir)
          print(filename)
          download.file(filename, tempfile, mode="wb")
            filelist <- unzip(tempfile, list=T)$Name
              dfile <- grep("modificationSpecificPeptides.txt",filelist, value=T)[1]
          if (any(grepl("modificationSpecificPeptides.txt",filelist)) & any(grepl("proteinGroups.txt",filelist))) {
              system(paste0("unzip ",tempfile," \"",dfile,"\""))
              modpep <- read.csv( dfile, sep="\t")
              file.remove(dfile)
              dfile <- grep("proteinGroups.txt",filelist, value=T)[1]
              system(paste0("unzip ",tempfile," \"",dfile,"\""))
              protlist <- read.csv( dfile, sep="\t")
              file.remove(dfile)
            }
                    file.remove(tempfile)
          
          # if both files are non-zero
          if (!is.null(modpep) & !is.null(protlist)) {
            try({
              MQout <- readMaxQuant(modpep, protlist)
              allPeps <- MQout$allPeps
              Prots <- MQout$Prots
              Param <- MQout$Param
              Benchmarks <- NULL
              rownames(allPeps) <- paste0("pep", 1:nrow(allPeps))
              if (withGroundTruth) {
                Stats <- runPolySTest(Prots, Param, refCond=1, onlyLIMMA=F)
                # much faster with only LIMMA tests   
                StatsPep <- runPolySTest(allPeps, Param, refCond=1, onlyLIMMA=T)
                Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
              } else {
                Benchmarks <- calcBasicBenchmarks(Prots, allPeps, Param)
              }
              save(Benchmarks, Prots, allPeps, Param, tdat, file =paste0("benchmarks_",which(filename==filenames),"_",tdat$accession, ".RData") )
            })
          }
        }
      }
    }
  }
}
  
  