################################################################################
#                     TO RUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################


#####################
## Set path:
#####################
wd <- getwd()
# # For Computerome:
# wd <- "/home/projects/jensenlab/people/malopa/PhosFake/BenchmarkingNoMod"
# #
pathToRes <- paste0(wd, "/Output/BenchmarkPepQuan")
pathToFasta <- paste0(wd, "/input_fasta")
pathToFasta <- list.files(path = pathToFasta, full.names = T, pattern = ".fasta")
pathToFunctions <- paste0(wd, "/Functions")
if (!dir.exists(pathToRes)) {
  cat("Create result directory:", pathToRes, "\n")
  dir.create(pathToRes)
}
#####################

#####################
## Load parameters
#####################
sapply(list.files(pathToFunctions, full.names = T), source)
# Parameters to test:
paramToTest <- list("PathToFasta" = pathToFasta,
                    "NumReps" = 3:4,
                    "QuantNoise" = 0.1, # I take as max sd the sd of the proteoform quan. values.
                    "PropMissedCleavages" = seq(from = 0, to = 0.3, by = 0.1), 
                    "MaxNumMissedCleavages" = 0:2,
                    "ThreshNAQuantileProt" = seq(from = 0, to = 0.2, by = 0.1),
                    "Threshqval" = c(0.01, 0.05), 
                    "PercDetectedPep" = seq(from = 0.8, to = 1, by = 0.1),
                    "PercDetectedVal" = seq(from = 0.8, to = 1, by = 0.1),
                    "WrongIDs" = c(0,0.01),
                    "WrongLocalizations" = c(0,0.01),
                    "MSNoise" = seq(from = 0, to = 0.1, by = 0.05),
                    "minUniquePep" = 1:2) 
paramPreDig <- c("NumReps", "QuantNoise", "PropMissedCleavages", "MaxNumMissedCleavages", "ThreshNAQuantileProt")
paramPostDig <- setdiff(names(paramToTest), c(paramPreDig, "PathToFasta"))
#####################

#####################
## Run the sample preparation simulation:
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
#####################

#####################
## Add quan. values to proteoforms:
#####################
library(purrr)
listtotest_predig <- cross(paramToTest[names(paramToTest) %in% paramPreDig])
listtotest_postdig <-  cross(paramToTest[names(paramToTest) %in% paramPostDig])
cat("Start generation of", length(listtotest_predig) * length(listtotest_postdig), "parameter sets for digestion per fasta.\n")

conditions <- unique(gsub("_R.+", "", Param$quant_colnames))
library(qvalue)

iter <- 1

for (j in seq_along(lp)) {
  f <- lp[[j]]
  
  for (x in listtotest_predig) {
    x <- c(Param[!(names(Param) %in% names(x))], x)
    x <- x[!(names(x) %in% paramPostDig)]
    x$quant_colnames <- paste0("C_",rep(1:x$NumCond,each=x$NumReps),"_R_", rep(1:x$NumReps, x$NumCond))
    d <- addProteoformAbundance(proteoforms = f, parameters = x)
    for (na in c("PTMPos", "PTMType")) {
      d[,names(d) == na] <- unlist(d[,names(d) == na])
    }
    
    dprot <- d
    
    ## Digestion:
    d <- DigestGroundTruth(GroundTruth = d, parameters = x)
    npep <-nrow(d)
    
    # Count number of non-unique peptides:
    d$AccessionCount <- sapply(d$Accession, function(x) {
      nchar(x) - nchar(gsub(";", "", x, fixed = T)) + 1
    })
    d <- d[d$AccessionCount == 1,]
    
    dpep <- d
    
    for (y in listtotest_postdig) {
      y <- c(x[!(names(x) %in% names(y))], y)
      
      ## Peptide -> protein
      protmat <- proteinSummarisation(peptable = d, minUniquePep = y$minUniquePep, parameters = y)
      # head(nacc <- nrow(na.omit(protmat))
      
      ## Stat:
      pval <- vector()
      means <- matrix(ncol = length(conditions), nrow = nrow(protmat))
      missval <- matrix(ncol = ncol(protmat), nrow = nrow(protmat))
      colnames(means) <- paste0("Mean_", conditions)
      # d[is.infinite(d)] <- NA
      for (i in seq_len(nrow(protmat))) {
        mv <- is.na(protmat[i,])
        # Number of missing values per condition:
        numpercond <- colSums(!(sapply(conditions, function(x) mv[grepl(x, colnames(protmat))])))
        mytest <- try(t.test(as.numeric(d[i,grepl(conditions[1], colnames(d))]), as.numeric(d[i,grepl(conditions[2], colnames(d))])), TRUE)
        if (inherits(mytest, "try-error")) {
          pval[i] <- NA
        } else {
          # } else {
          pval[i] <- mytest$p.value
          means[i,] <- mytest$estimate
        }
        missval[i,] <- mv
      }
      
      protmat <- cbind(protmat, means)
      protmat <- as.data.frame(protmat)
      protmat$pvalues <- pval
      
      protmat$qvalues <- qvalue(pval)$qvalues
      protmat$Mean_Diff <- means[,2] - means[,1]
      protmat$Regulated <- !is.na(dprot$Regulation_Amplitude[match(row.names(protmat), dprot$Accession)])
      
      numtotprot <- nrow(protmat)
      numRegTrue <- sum(protmat$Regulated & protmat$qvalues <= protmat$Threshqval, na.rm = T)
      numRegTot <- sum(protmat$qvalues <= protmat$Threshqval, na.rm = T)
      rownames(protmat)[grepl("^C_", rownames(protmat))] <- paste0("PostDig_", rownames(protmat)[grepl("^C_", rownames(protmat))])
      outtab <- cbind(dprot, protmat[match(dprot$Accession, row.names(protmat)),])
      
      output <- list("Fasta" = names(lp)[j],
                     "Param" = y, 
                     "NumberPep" = npep,
                     "NumberTrueRegulated" = numRegTrue,
                     "NumberTotRegulated" = numRegTot, 
                     "FinalNumberOfProt" = numtotprot,
                     "data" = outtab)
      save(output,
           file = paste0(pathToRes, "/output", iter, ".RData"))
      cat("Save output", iter, "over", length(listtotest), "\n")
      iter <- iter + 1
      
      
    }
  }
}
#####################

sessionInfo()
