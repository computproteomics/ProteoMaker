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
pathToRes <- paste0(wd, "/Output/BenchmarkProteoformsQuan")
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
pathToFasta <- pathToFasta[grepl("5640", pathToFasta)]

paramToTest <- list("PathToFasta" = pathToFasta,
                    "QuantNoise" = c(0.01, 0.25), # I take as max sd the sd of the proteoform quan. values.
                    "ThreshNAQuantileProt" = seq(from = 0, to = 0.2, by = 0.05),
                    "Threshqval" = c(0.01, 0.05)) 
#####################

#####################
## Generate the modified forms:
#####################
pathToStoechio <- paste0(wd, "/input_Stoechio/TG_Dataset_14_PTM_Stoichiometry.txt")
st <- read.table(pathToStoechio, sep = "\t", header = T, stringsAsFactors = F)
st <- st[st$Mod == "phospho",] # I keep only the phospho for the moment
library(stringr)
st$modcound <- sapply(st$AA, str_count, ",") + 1
st$Accession <- sapply(st$Protein, function(x) {
  strsplit(x, "|", fixed = T)[[1]][2]
})
PTMType <- lapply(seq_len(nrow(st)), function(x) {
  rep("ph", st$modcound[x])
})
PTMPos <- st$Site
mats <- data.frame("Accession" = st$Accession, "PTMPos" = PTMPos, "Stoichio" = st$Adult_Retina)
mats$PTMType <- PTMType
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

# iter <- 53 # Reset to 1
iter <- 1
# for (j in seq_along(lp)) {
j = 1
f <- lp[[j]]
# Custom: # To reset
# if (j == 1) {
#   listtotest2 <- listtotest[53:length(listtotest)]
# } else (
listtotest2 <- listtotest
# )
#

mats <- mats[mats$Accession %in% f$Accession,]
mats <- mats[mats$Stoichio > 0,]

for (x in listtotest2) {
  x <- c(Param[!(names(Param) %in% names(x))], x)
  
  # I keep the info of quantnoise to use it after proteoforms generation:
  quantNoise <- x$QuantNoise
  x$QuantNoise <- 0
  
  x$quant_colnames <- paste0("C_",rep(1:x$NumCond,each=x$NumReps),"_R_", rep(1:x$NumReps, x$NumCond))
  d <- addProteoformAbundance(proteoforms = f, parameters = x)
  
  bnoise <- quantile(as.matrix(d[,grepl("^C_", names(d))]), probs = 0.05)
  
  # Generate the proteoforms:
  lnew <- vector(mode = "list")
  for (r in seq_len(nrow(mats))) {
    vec <- mats[r,]
    if (vec$Stoichio > 0) {
      newprot <- d[d$Accession == as.character(vec$Accession),]
      newprot[grepl("^C_", names(newprot))] <- log2(vec$Stoichio * 2^newprot[grepl("^C_", names(newprot))])
      if (vec$Stoichio < 1) {
        d[d$Accession == as.character(vec$Accession),grepl("^C_", names(d))] <- log2(2^newprot[grepl("^C_", names(newprot))] * (1-vec$Stoichio))
      } else {
        d[d$Accession == as.character(vec$Accession),grepl("^C_", names(d))] <- rnorm(n = sum(grepl("^C_", names(d))), mean = bnoise, sd = 0.01) # I add noise when total stoichiometry > 1
        # TODO: Correct this issue of stoichio for prot with tot. stoichio > 1
      }
      newprot$PTMPos <- vec$PTMPos
      newprot$PTMType <- vec$PTMType
      lnew[[length(lnew) + 1]] <- newprot
    }
  }
  d <- rbind(d, data.table::rbindlist(lnew))
  d <- as.data.frame(d)
  # Add noise:
  matnoise <- matrix(rnorm(n = sum(grepl("^C_", names(newprot))) * nrow(d), mean = 0, sd = quantNoise), ncol = sum(grepl("^C_", names(newprot))), nrow = nrow(d))
  for (colnum in which(grepl("^C_", names(d)))) {
    d[,colnum] <- d[,colnum] + matnoise[,(colnum - which(grepl("^C_", names(d)))[1] + 1)]
  }
  
  nacc <- length(unique(d$Accession))
  nprot <-nrow(d)
  nMV <-sum(is.na(d[,grepl("^C_", names(d))]))
  
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
    # if (sum(numpercond >= 2) == length(numpercond)) {
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
  
  
  d <- cbind(d, means)
  d <- as.data.frame(d)
  d$pvalues <- pval
  
  d$qvalues <- qvalue(pval)$qvalues
  d$Mean_Diff <- means[,2] - means[,1]
  # d$Accession <- row.names(d)
  d$Regulated <- !is.na(d$Regulation_Amplitude)
  
  numRegTrue <- sum(d$Regulated & d$qvalues <= x$Threshqval, na.rm = T)
  numRegTot <- sum(d$qvalues <= x$Threshqval, na.rm = T)
  numRegTruePerAmplitude <- sapply(unique(d$Regulation_Amplitude)[!is.na(unique(d$Regulation_Amplitude))], function(amp) {
    sum(d$Regulated[d$Regulation_Amplitude == amp & d$qvalues <= x$Threshqval & !is.na(d$Regulation_Amplitude)], na.rm = T)
  })
  matRegPerAmp <- data.frame("Regulation_Amplitude" = unique(d$Regulation_Amplitude)[!is.na(unique(d$Regulation_Amplitude))],
                             "numRegTrue" = numRegTruePerAmplitude)
  
  # library(ggplot2)
  # ggplot(d, aes(x = Mean_Diff, y = -log10(qvalues), col = Regulated)) + geom_point()
  
  output <- list("Fasta" = names(lp)[j],
                 "Param" = x, 
                 "numberUniqueAccessions" = nacc,
                 "NumberUniqueProteoform" = nprot,
                 "NumberMissingValues" = nMV,
                 "NumberTrueRegulated" = numRegTrue,
                 "NumberTotRegulated" = numRegTot, 
                 "NumberRegPerAmplitude" = matRegPerAmp,
                 "data" = d)
  save(output,
       file = paste0(pathToRes, "/output_4_", iter, ".RData"))
  cat("Save output", iter, "over", length(listtotest)*length(lp), "\n")
  iter <- iter + 1
}
# }
#####################

sessionInfo()
