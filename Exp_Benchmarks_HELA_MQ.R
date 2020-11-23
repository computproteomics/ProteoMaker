library(moments)
library(stringr)
source ("06_Benchmarks.R")

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
AllExpBenchmarks <- AllProts <- AllPeps <- NULL
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
    
    AllProts[[bench]] <- Prots[,c("Accession","Intensity.HeLA")]
    AllProts[[bench]] <- AllProts[[bench]][!is.na(AllProts[[bench]][,1]),]
    AllPeps[[bench]] <- cbind(SeqMod=paste0(allPeps[,"Sequence"], allPeps[,"Modifications"]),allPeps[,"Intensity.HeLA"])
    ### Get additional benchmarks only for experimental data
    ## Max. difference retention time
    Benchmarks$globalBMs$diffRetentionTime <- diff(range(allPeps$Retention.time))
    
    ## Number of accepted PSMs (count scan numbers)
    Benchmarks$globalBMs$acceptedPSMs <- sum(allPeps$MS.MS.Count)
    
    
  })
  save(Benchmarks, Prots, allPeps, Param, tdat, file =bench)
  AllExpBenchmarks[[bench]] <- c(unlist(Benchmarks$globalBMs), tdat)
  
  
}


# getting all protein names and peptide sequences
AllAccs <- unlist(sapply(AllProts, function(x) x[,1]))
AllSeqs <- unlist(sapply(AllPeps, function(x) x[,1]))
AllAccs <- unique(AllAccs)
AllSeqs <- unique(AllSeqs)


# writing full tables
fullProtTable <- matrix(NA, nrow=length(AllAccs), ncol=length(AllProts), dimnames=list(x=AllAccs, y=names(AllProts)))
fullPepTable <- matrix(NA, nrow=length(AllSeqs), ncol=length(AllPeps), dimnames=list(x=AllSeqs, y=names(AllPeps)))
for (bench in names(AllProts)) {
  fullProtTable[AllProts[[bench]][,1], bench] <- AllProts[[bench]][,2]
  fullPepTable[AllPeps[[bench]][,1], bench] <- AllPeps[[bench]][,2]
} 

names(AllExpBenchmarks) <- benchmarks
save(fullProtTable, fullPepTable, AllExpBenchmarks, file="AllExpBenchmarks.RData")



