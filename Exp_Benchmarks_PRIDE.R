library(bigmemory)
library(moments)
library(stringr)
source ("06_Benchmarks.R")


############ Call the function to load and prepare data, and then to calculate the benchmarks
### It is now adopted to download files from PRIDE according to the metadata
library(jsonlite)
## General parameters
withGroundTruth <- F
tmpdir <- "tmp"
tmpdir <- normalizePath(tmpdir)

#allPRIDE <- readLines("pxd_accessions.json")
allPRIDE <- read.csv(unz("pride_df.zip","pride_df.csv"), stringsAsFactors = F)

for (pxd in 1:nrow(allPRIDE))  {
  
  tdat <- allPRIDE[pxd,]
  modpep <- protlist <- NULL
  if (!is.null(tdat$dataProcessingProtocol)) {
    # if maxquant
    try({
      if (grepl("maxquant",tolower(tdat$dataProcessingProtocol)) ) {
        print(pxd)
        tfile_list <- read_json(paste0("https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession=",tdat$accession))
        file_types <- sapply(tfile_list, function(x) x$fileCategory$value)
        file_list <- sapply(tfile_list, function(x) ifelse(grepl("ftp://",x$publicFileLocations[[1]]$value), x$publicFileLocations[[1]]$value, x$publicFileLocations[[2]]$value ))
        # if maxquant files available
        if (any(grepl("/modificationSpecificPeptides.txt",file_list)) &  any(grepl("/proteinGroups.txt",file_list))) {
          modpep <- read.csv(grep("/modificationSpecificPeptides.txt",file_list, value=T), sep="\t")
          protlist <- read.csv(grep("/proteinGroups.txt",file_list, value=T), sep="\t")
          # if maxquant files in zip file
        } else if (any(file_types == "SEARCH" & grepl(".zip",file_list) )) {
          filenames <- file_list[which(file_types == "SEARCH" & grepl(".zip",file_list))]
          for (filename in filenames) {
            tempfile <- tempfile(tmpdir=tmpdir)
            print(filename)
            download.file(filename, tempfile, mode="wb")
            #              filelist <- unzip(tempfile, list=T)$Name
            filelist <- system(paste0("jar tf ", tempfile), intern=T)
            # non-optimal: takes only the first of the files if multiples are available in zip-file     
            if (any(grepl("modificationSpecificPeptides.txt",filelist)) & any(grepl("proteinGroups.txt",filelist))) {
              #	    print(filelist)
              dfile <- grep("modificationSpecificPeptides.txt",filelist, value=T)[1]
              system(paste0("cd ",tmpdir,";jar xvf ",tempfile," \"",dfile,"\""))
              tdfile <- paste0(tmpdir,"/", dfile)
              modpep <- read.csv(tdfile, sep="\t")
              file.remove(tdfile)
              dfile <- grep("proteinGroups.txt",filelist, value=T)[1]
              tdfile <- paste0(tmpdir,"/", dfile)
              system(paste0("cd ",tmpdir,";jar xvf ",tempfile," \"",dfile,"\""))
              protlist <- read.csv(tdfile, sep="\t")
              file.remove(tdfile)
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
                
                ### Get additional benchmarks only for experimental data
                ## Max. difference retention time
                Benchmarks$globalBMs$diffRetentionTime <- diff(range(MQout$allPeps$Retention.time))
                
                ## Number of accepted PSMs (count scan numbers)
                Benchmarks$globalBMs$acceptedPSMs <- sum(MQout$allPeps$MS.MS.Count)
                
                save(Benchmarks, Prots, allPeps, Param, tdat, file =paste0("benchmarks_",which(filename==filenames),"_",tdat$accession, ".RData") )
              })
            }
          }
        }
      }
    })
  }
}



benchmarks <- list.files("./","benchmarks.*RData")
AllExpBenchmarks <- AllProts <- AllPeps <- numQuantCol <- NULL
for (bench in benchmarks) {
  try({
    load(bench)
    print(bench)
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
    AllPeps[[bench]] <- cbind(SeqMod=paste0(allPeps[,"Sequence"], allPeps[,"Modifications"]),allPeps[,Param$QuantColnames])
    numQuantCol[[bench]] <- length(Param$QuantColnames)

    
    ### Get additional benchmarks only for experimental data
    ## Max. difference retention time
    Benchmarks$globalBMs$diffRetentionTime <- diff(range(allPeps$Retention.time))
    
    ## Number of accepted PSMs (count scan numbers)
    Benchmarks$globalBMs$acceptedPSMs <- sum(allPeps$MS.MS.Count)
    
    
    save(Benchmarks, Prots, allPeps, Param, tdat, file =bench)
    AllExpBenchmarks[[bench]] <- unlist(c(Benchmarks$globalBMs, tdat))
  })
}

save(AllExpBenchmarks, file="AllExpBenchmarks.RData")


# getting all protein names and peptide sequences
AllAccs <- unlist(sapply(AllProts, function(x) x[,1]))
AllSeqs <- unlist(sapply(AllPeps, function(x) x[,1]))
AllAccs <- unique(AllAccs)
AllSeqs <- unique(AllSeqs)


# save all protein names and peptide sequence
save(numQuantCol, allAccs, allSeqs, file="FullProtPepList.RData")

library(bigmemory)
load("FullProtPepList.RData")
benchmarks <- list.files("./","benchmarks.*RData")

# writing full tables
print(paste("Altogether", sum(unlist(numQuantCol)), "columns"))
cols <- NULL
for (bench in benchmarks) {
  cols <- c(cols, paste(bench,1:numQuantCol,sep="_"))
}
options(bigmemory.allow.dimnames=TRUE)
fullProtTable <- big.matrix(nrow=length(AllAccs), ncol=length(cols), type="double", init=NA, dimnames=list(x=AllAccs, y=cols), backingpath="./", backingfile="tmp.big")
fullPepTable <- big.matrix(nrow=length(AllSeqs), ncol=length(cols), type="double", init=NA, dimnames=list(x=AllSeqs, y=cols), backingpath="./", backingfile="tmp2.big")
#fullProtTable <- matrix(NA, nrow=length(AllAccs), ncol=length(benchmarks), dimnames=list(x=AllAccs, y=names(benchmarks)))
#fullPepTable <- matrix(NA, nrow=length(AllSeqs), ncol=length(benchmarks), dimnames=list(x=AllSeqs, y=names(benchmarks)))
for (bench in benchmarks) {
  try({
    load(bench)
    print(bench)
    cols <- paste(bench,1:numQuantCol,sep="_")
    AllProts <- Prots[,c("Accession",Param$QuantColnames)]
    AllProts <- AllProts[!is.na(AllProts[,1]),]
    fullProtTable[as.character(AllProts[,1]), cols] <- AllProts[,2:ncol(AllProts)]
    AllPeps <- cbind(SeqMod=paste0(allPeps[,"Sequence"], allPeps[,"Modifications"]),allPeps[,Param$QuantColnames])
    fullPepTable[as.character(AllPeps[,1]), cols] <- AllPeps[,2:ncol(AllProts)]
  })
} 

save(fullProtTable, "FullProtQuant.RData")
save(fullPepTable, "FullPepQuant.RData")

