source("Parameter.R")

library(OrgMassSpecR)

# proteoformsRow <- proteoforms[58,]

getDigestTables <- function(proteoformsRow) {
  df <- Digest(proteoformsRow$Sequence,enzyme = "trypsin", missed = 0)
  df$Accession <- proteoformsRow$Accession
  df$PTMPos <- NA
  df$PTMType <- NA
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

ldig <- lapply(seq_len(nrow(proteoforms)), function(x) {
  getDigestTables(proteoforms[x,])
})
peptable <- ldig[[1]]
for (i in 2:length(ldig)) {
  peptable <- rbind(peptable, ldig[[i]])
}

names(peptable)[names(peptable) == "peptide"] <- "PepSequence"
names(peptable)[names(peptable) == "start"] <- "PepStart"
names(peptable)[names(peptable) == "stop"] <- "PepStop"

insilico_peptides <- peptable

save(insilico_peptides, file = "data/insilicoPep.RData")
