library(moments)
library(readxl)
library(stringr)
source ("06_Benchmarks.R")

############ Call the function to load and prepare data, and then to calculate the benchmarks
### It is now adopted to a set of zip files in a given folder
library(jsonlite)
## General parameters
withGroundTruth <- F
tmpdir <- "/tmp"
metadatafile <- "../../MNTFileParsed.txt"
metadata <- read.csv(metadatafile, sep="\t", stringsAsFactors = F)
metadata$FileName <- sapply(as.character(metadata$FileName), basename)
tmpdir <- normalizePath(tmpdir)

allxlsx <- list.files("../../WOMBAT-Results/",recursive=T, pattern="tpp\\.zip",full.names=T)

for (dat in allxlsx)  {

  print(dat)
#  dat <- normalizePath(dat)  
#  file.copy(dat, "tmp/")
  unzip(dat, exdir="tmp/")
  samplename <- paste0(sub("\\.raw_tpp\\.zip","",basename(dat)), ".raw")
  protlist <- modpep <- psms <- NULL
  print("dat")
  print(dat)
  if (file.exists("tmp/stand_pep_quant_merged.csv") > 0) {
  m_ind <- grep(samplename, metadata$FileName)
  tdat <- NA
  if (length(m_ind) == 1) {
    tdat <- metadata[m_ind,]    
  } else {
    print("WARNING: no metadata available")
  }
  
  # reading data from tables
  protlist <- read.csv("tmp/stand_prot_quant_merged.csv")
  modpep <- read.csv("tmp/stand_pep_quant_merged.csv")
  }
  # Delete files in tmp folder
  unlink("tmp/*")

  # if both files are non-zero
  if (!is.null(modpep) & !is.null(protlist)) {
    try({
      MQout <- readWombat(modpep, protlist)
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
      Benchmarks$globalBMs$acceptedPSMs <- sum(MQout$allPeps$number_of_psms)
      dat <- basename(dat)
      save(Benchmarks, Prots, allPeps, Param, tdat, file =paste0("benchmarks_",dat, ".RData") )
    })
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
    
    AllProts[[bench]] <- Prots[,c("Accession",Param$QuantColnames)]
    AllProts[[bench]] <- AllProts[[bench]][!is.na(AllProts[[bench]][,1]),]
    # to avoid having NAs in the sequence endings
    #NAsremoved <- allPeps[,"Modifications"]
    #NAsremoved[is.na(NAsremoved)] <- ""
    #AllPeps[[bench]] <- cbind(SeqMod=paste0(allPeps[,"Sequence"], NAsremoved),allPeps[,Param$QuantColnames])
    AllPeps[[bench]] <- cbind(SeqMod=allPeps[,"Sequence"], allPeps[,Param$QuantColnames])
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

save(AllAccs, AllSeqs, AllProts, AllPeps, file="ttt.RData")

# writing full tables
fullProtTable <- matrix(NA, nrow=length(AllAccs), ncol=length(AllProts), dimnames=list(x=AllAccs, y=names(AllProts)))
fullPepTable <- matrix(NA, nrow=length(AllSeqs), ncol=length(AllPeps), dimnames=list(x=AllSeqs, y=names(AllPeps)))
for (bench in names(AllProts)) {
  fullProtTable[AllProts[[bench]][,1], bench] <- AllProts[[bench]][,2]
  fullPepTable[AllPeps[[bench]][,1], bench] <- AllPeps[[bench]][,2]
} 

names(AllExpBenchmarks) <- benchmarks
save(fullProtTable, fullPepTable, AllExpBenchmarks, file="AllExpBenchmarks.RData")



