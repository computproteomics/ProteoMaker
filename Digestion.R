source("Parameter.R")

library(OrgMassSpecR)

# proteoformsRow <- GroundTruth[1,]

getDigestTables <- function(proteoformsRow) {
  df <- Digest(proteoformsRow$Sequence,enzyme = "trypsin", missed = 0)
  df$Accession <- proteoformsRow$Accession
  df$PTMPos <- NA
  df$PTMType <- NA
  vecquan <- proteoformsRow[grepl( "^C_" , names( proteoformsRow ) ) ]
  quan <- matrix(nrow = nrow(df), ncol = length(vecquan), byrow = T, data = vecquan)
  df <- cbind(df, quan)
  names(df)[(ncol(df) - ncol(quan) + 1):ncol(df)] <- names(vecquan)
  if (!is.null(proteoformsRow$PTMPos)) {
    for (i in seq_len(nrow(df))) {
      sel <- unlist(proteoformsRow$PTMPos) >= df$start[i] & unlist(proteoformsRow$PTMPos) <= df$stop[i]
      if (any(sel)) {
        df$PTMType[i] <- c(unlist(proteoformsRow$PTMType)[sel])
        df$PTMPos[i] <- c(unlist(proteoformsRow$PTMPos)[sel])
      } else {
        df$PTMPos[i] <- NA
        df$PTMType[i] <- NA
      }
    }
  }
  df
}

ldig <- lapply(seq_len(nrow(GroundTruth)), function(x) {
  getDigestTables(GroundTruth[x,])
})
peptable <- ldig[[1]]
for (i in 2:length(ldig)) {
  peptable <- rbind(peptable, ldig[[i]])
}

names(peptable)[names(peptable) == "peptide"] <- "PepSequence"
names(peptable)[names(peptable) == "start"] <- "PepStart"
names(peptable)[names(peptable) == "stop"] <- "PepStop"

peptable$ID <- paste(peptable$PepSequence, peptable$PTMPos, peptable$PTMType, sep="_")
unique_pep <- names(table(peptable$ID)[table(peptable$ID) == 1])
uniquetab <- peptable[peptable$ID %in% unique_pep,]
redundanttab <- peptable[!(peptable$ID %in% unique_pep),]

iterID <- unique(redundanttab$ID)
list_quan <- vector(mode = "list")
for (i in seq_along(iterID)) {
  df <- redundanttab[redundanttab$ID == iterID[i],]
  mtx <- df[, grepl( "^C_" , names( df ) ) ]
  vec <- sapply(1:ncol(mtx),function(x) {
    log2(sum(2^as.numeric(mtx[,x]), na.rm=T))
  })
  df <- df[1,]
  df[, grepl( "^C_" , names( df ) ) ] <- vec
  list_quan[[i]] <- df
}

redundanttab <- do.call(rbind,list_quan)
peptable <- rbind(redundanttab, uniquetab)

 
#insilico_peptides <- peptable
#save(insilico_peptides, file = "data/insilicoPep.RData")

##############################################################
#FILTERING

pepLength <- nchar(peptable$PepSequence)
peptable <- peptable[(pepLength >= Param$PepMinLength & pepLength <= Param$PepMaxLength),]

peptable$Enriched <- grepl("ph", peptable$PTMType)
numRemove <- floor(sum(peptable$Enriched) * Param$EnrichmentLoss)
idx <- sample(which(peptable$Enriched == TRUE), size = numRemove)
tokeep <- setdiff(seq_len(nrow(peptable)), idx)
peptable <- peptable[tokeep,]

numRemove <- floor(sum(!(peptable$Enriched)) * (1 - Param$EnrichmentEfficiency))
idx <- sample(which(peptable$Enriched == FALSE), size = numRemove)
noise <- peptable[idx,]
noise$Enriched <- TRUE
peptable <- rbind(peptable, noise)

#unique_pep <- names(table(peptable$ID)[table(peptable$ID) == 1])
enrichedtab <- peptable[peptable$Enriched,]
nonenrichedtab <- peptable[!peptable$Enriched,]
nrowTab <- nrow(enrichedtab)
ncolTab <- length(grepl( "^C_" , names( proteoformsRow ) ))
mtx <- matrix(nrow = nrowTab, ncol= ncolTab, data=rnorm(n=ncolTab * nrowTab, mean=0, sd=Param$EnrichmentNoise))
enrichedtab[, grepl( "^C_" , names( enrichedtab ) ) ] <- apply(enrichedtab[, grepl( "^C_" , names( enrichedtab ) ) ], 2, as.numeric)
enrichedtab[, grepl( "^C_" , names( enrichedtab ) ) ] <- enrichedtab[, grepl( "^C_" , names( enrichedtab ) ) ] + mtx

#insilico_peptides_enriched <- enrichedtab
#save(insilico_peptides_enriched, file = "data/insilicoPepEnriched.RData")



