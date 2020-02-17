################################################################################
#                     TO RUN THE ENTIRE PHOSFAKE PIPELINE                      #
################################################################################


#####################
## Set path:
#####################
wd <- getwd()
pathToRes <- paste0(wd, "/RData")
pathToFasta <- paste0(wd, "/input_data")
pathToFasta <- list.files(path = pathToFasta, full.names = T, pattern = ".fasta")
pathToFunctions <- paste0(wd, "/Functions/")
#####################

#####################
## Load parameters
#####################
sapply(list.files(pathToFunctions, full.names = T), source)
# Parameters to test:
paramToTest <- list("PathToFasta" = pathToFasta, 
                    "NumReps" = seq(from = 3, to = 8, by = 1), 
                    "QuantNoise" = seq(from = 0.01, to = Param$AbsoluteQuanSD, by = 0.1), # I take as max sd the sd of the proteoform quan. values.
                    "ThreshNAQuantileProt" = seq(from = 0, to = 0.6, by = 0.05),
                    "AbsoluteQuanSD" = seq(from = 0.01, to = 6, by = 0.2), 
                    "Threshqval" = c(0.01, 0.05)) 
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
paramToTest <- paramToTest[-1]
listtotest <- cross(paramToTest)
cat("Start generation of", length(listtotest), "parameter sets for digestion\n")

conditions <- unique(gsub("_R.+", "", Param$quant_colnames))
library(qvalue)

iter <- 1
for (j in seq_along(lp)) {
  f <- lp[[j]]
  for (x in listtotest) {
    d <- addProteoformAbundance(proteoforms = f, parameters = c(Param[!(names(Param) %in% names(x))], x))
    for (na in c("PTMPos", "PTMType")) {
      d[,names(d) == na] <- unlist(d[,names(d) == na])
    }
    nacc <- length(unique(d$Accession))
    nprot <-nrow(d)
    nMC <-sum(is.na(d[,grepl("^C_", names(d))]))
    
    ## Stat:
    pval <- vector()
    means <- matrix(ncol = length(conditions), nrow = nrow(d))
    missval <- matrix(ncol = ncol(d), nrow = nrow(d))
    colnames(means) <- paste0("Mean_", conditions)
    # d[is.infinite(d)] <- NA
    for (i in seq_len(nrow(d))) {
      mv <- is.na(d[i,])
      # Number of missing values per condition:
      numpercond <- colSums(!(sapply(conditions, function(x) mv[grepl(x, colnames(d))])))
      if (sum(numpercond >= 2) == length(numpercond)) {
        tres <- t.test(as.numeric(d[i,grepl(conditions[1], colnames(d))]), as.numeric(d[i,grepl(conditions[2], colnames(d))]))
        pval[i] <- tres$p.value
        means[i,] <- tres$estimate
      }
      missval[i,] <- mv
    }
    
    d <- cbind(d, means)
    d <- as.data.frame(d)
    d$pvalues <- pval
    
    d$qvalues <- qvalue(pval)$qvalues
    d$Mean_Diff <- means[,2] - means[,1]
    # d$Accession <- row.names(d)
    d$Regulated <- !is.na(d$Regulation_Amplitude)
    
    numRegTrue <- sum(d$Regulated & d$qvalues <= x$Threshqval)
    numRegTot <- sum(d$qvalues <= x$Threshqval)
    numRegTruePerAmplitude <- sapply(unique(d$Regulation_Amplitude)[!is.na(unique(d$Regulation_Amplitude))], function(amp) {
      sum(d$Regulated[d$Regulation_Amplitude == amp & d$qvalues <= x$Threshqval & !is.na(d$Regulation_Amplitude)])
    })
    matRegPerAmp <- data.frame("Regulation_Amplitude" = unique(d$Regulation_Amplitude)[!is.na(unique(d$Regulation_Amplitude))],
                               "numRegTrue" = numRegTruePerAmplitude)
    
    output <- list("Fasta" = names(lp)[j],
                   "Param" = x, 
                   "numberUniqueAccessions" = nacc,
                   "NumberUniqueProteoform" = nprot,
                   "NumberMissingValues" = nMC,
                   "NumberTrueRegulated" = numRegTrue,
                   "NumberTotRegulated" = numRegTot, 
                   "NumberRegPerAmplitude" = matRegPerAmp)
    save(output,
         file = paste0(pathToRes, "/ouptut", iter, ".RData"))
    iter <- iter + 1
  }
}
#####################

sessionInfo()
