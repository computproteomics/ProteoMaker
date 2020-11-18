library(moments)
library(stringr)
source ("06_Benchmarks.R")

readMaxQuant <- function(allPeps, Prots, Param=NULL) {
  allPeps$Accession <- sapply(allPeps$Proteins, function(x) strsplit(as.character(x), ";"))
  allPeps$MC <- allPeps$Missed.cleavages
  allPeps$PTMType <- as.character(allPeps$Modifications)
  allPeps$PTMType[allPeps$PTMType == "Unmodified"] <- "NULL"
  # did not find corresponding field    
  allPeps$PTMPos <- NA
  allPeps$PTMPos[allPeps$PTMType != "NULL"] <- 1
  
  # remove entries with missing protein name (should be reverse)
  allPeps <- allPeps[allPeps$Proteins != "",]
  
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
  allPeps <- allPeps[rowSums(is.na(allPeps[, Param$QuantColnames,drop=F])) < length(Param$QuantColnames), ]
  tquant <- Prots[,Param$QuantColnames]
  allPeps$Sequence <- as.character(allPeps$Sequence)
  tquant[tquant == 0] <- NA
  Prots[, Param$QuantColnames] <- log2(tquant)
  Prots <- Prots[rowSums(is.na(Prots[, Param$QuantColnames,drop=F])) < length(Param$QuantColnames), ]
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
  
  return(list(Param=Param, allPeps=allPeps, Prots=Prots))
}


############ Call the function to load and prepare data, and then to calculate the benchmarks
### It is now adopted to a set of zip files in a given folder
library(jsonlite)
## General parameters
withGroundTruth <- F
tmpdir <- "tmp"
metadatafile <- "/tmp/MNTFileParsed.txt"
metadata <- read.csv(metadatafile, sep="\t", stringsAsFactors = F)
metadata$FileName <- sapply(as.character(metadata$FileName), basename)
tmpdir <- normalizePath(tmpdir)


allzips <- grep(".zip",list.files(full.names=T),value=T)

for (zip in allzips)  {
  
  zip <- normalizePath(zip)  
  filelist <- system(paste0("jar tf ", zip), intern=T)
  # non-optimal: takes only the first of the files if multiples are available in zip-file     
  if (any(grepl("modificationSpecificPeptides.txt",filelist)) & any(grepl("proteinGroups.txt",filelist))) {
    
    # match to raw file name in metadata file
    tzip <- basename(zip)
    tzip <- sub(".zip","",tzip)
    m_ind <- grep(tzip, metadata$FileName)
    tdat <- NA
    if (length(m_ind) == 1) {
      tdat <- metadata[m_ind,]    
    } else {
      print("WARNING: no metadata available")
    }
    
    #	    print(filelist)
    dfile <- grep("modificationSpecificPeptides.txt",filelist, value=T)[1]
    system(paste0("cd ",tmpdir,";jar xvf ",zip," \"",dfile,"\""))
    tdfile <- paste0(tmpdir,"/", dfile)
    modpep <- NULL
    if (file.size(tdfile) > 0) {
    modpep <- read.csv(tdfile, sep="\t")
    }
    file.remove(tdfile)
    dfile <- grep("proteinGroups.txt",filelist, value=T)[1]
    tdfile <- paste0(tmpdir,"/", dfile)
    system(paste0("cd ",tmpdir,";jar xvf ",zip," \"",dfile,"\""))
    protlist <- NULL
    if (file.size(tdfile) > 0) {
      protlist <- read.csv(tdfile, sep="\t")
    }
    file.remove(tdfile)
    
    
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
        
        ### Get additional benchmarks only for experimental data
        ## Max. difference retention time
        Benchmarks$globalBMs$diffRetentionTime <- diff(range(MQout$allPeps$Retention.time))
        
        ## Number of accepted PSMs (count scan numbers)
        Benchmarks$globalBMs$acceptedPSMs <- sum(MQout$allPeps$MS.MS.Count)
        zip <- basename(zip)
        save(Benchmarks, Prots, allPeps, Param, tdat, file =paste0("benchmarks_",zip, ".RData") )
      })
    }
    
  }
}



benchmarks <- list.files("./","benchmarks.*RData")
AllExpBenchmarks <- NULL
for (bench in benchmarks) {
  load(bench)
  print(bench)
  Benchmarks$globalBMs <- NA
  try({
    if (withGroundTruth) {
      Stats <- runPolySTest(Prots, Param, refCond=1, onlyLIMMA=F)
      # much faster with only LIMMA tests   
      StatsPep <- runPolySTest(allPeps, Param, refCond=1, onlyLIMMA=T)
      Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
    } else {
      Benchmarks <- calcBasicBenchmarks(Prots, allPeps, Param)
    }
  })
  save(Benchmarks, Prots, allPeps, Param, tdat, file =bench)
  if (length(Benchmarks$globalBMS) == 1 & !is.na(Benchmarks$globalBMs)) {
    AllExpBenchmarks <- rbind(AllExpBenchmarks,c(unlist(Benchmarks$globalBMs), tdat))
  } else {
    AllExpBenchmarks <- rbind(AllExpBenchmarks,NA)
    
  }
}
names(AllExpBenchmarks) <- benchmarks
save(AllExpBenchmarks, file="AllExpBenchmarks.RData")
