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
  }
}



benchmarks <- list.files("./","benchmarks.*RData")
AllExpBenchmarks <- NULL
for (bench in benchmarks) {
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
  save(Benchmarks, Prots, allPeps, Param, tdat, file =bench)
  AllExpBenchmarks <- rbind(c(unlist(Benchmarks$globalBMs), tdat))
}
rownames(AllExpBenchmarks) <- benchmarks
save(AllExpBenchmarks, file="AllExpBenchmarks.RData")

