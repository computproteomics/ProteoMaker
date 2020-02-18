################################################################################
#                    PEPTIDE DIGESTION OF PROTEOFORM TABLE                     #
################################################################################

library(OrgMassSpecR)

# proteoformsRow <- GroundTruth[1,]

getDigestTablesNoMC <- function(proteoformsRow, parameters) {
  # print(proteoformsRow$Sequence)
  if (grepl("K|R", proteoformsRow$Sequence)) {
    df <- Digest(proteoformsRow$Sequence,enzyme = "trypsin", missed = 0)
    df$Accession <- proteoformsRow$Accession
    df$PTMPos <- NA
    df$PTMType <- NA
    vecquan <- proteoformsRow[grepl( "^C_" , names( proteoformsRow ) ) ]
    quan <- matrix(nrow = nrow(df), 
                   ncol = length(vecquan), 
                   byrow = T, 
                   data = rep(as.numeric(vecquan), nrow(df)))
    colnames(quan) <- names(vecquan)
    df <- cbind(df, quan)
    names(df)[(ncol(df) - ncol(quan) + 1):ncol(df)] <- names(vecquan)
    if (sum(sapply(proteoformsRow$PTMPos, is.null)) == 0) {
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
    return(df)
  }
}

DigestGroundTruth <- function(GroundTruth, parameters) {
  
  cat("Start digestion\n")
  # Get all the peptides without missed cleavage:
  d <- lapply(seq_len(nrow(GroundTruth)),function(x) {
    # print(GroundTruth[x,]$Sequence)
    getDigestTablesNoMC(GroundTruth[x,], parameters)
  })
  
  if (parameters$MaxNumMissedCleavages > 0) {
    ## Generate missed cleavages:
    cat("Start generation of missed-cleavages\n")
    for (el in seq_along(d)) {
      x <- d[[el]]
      if (nrow(x) > 1) {
        iter_mc <- 0
        r <- 1
        while (r < nrow(x)) {
          if (runif(1, 0, 1) <= parameters$PropMissedCleavages & iter_mc < parameters$MaxNumMissedCleavages & (nrow(x) - r) >= iter_mc) {
            newpep <- x[r+1,]
            x <- x[-(r+1),]
            newpep$peptide <- paste(c(x$peptide[r], newpep$peptide), collapse = "")
            newpep$start <- x$start[r]
            iter_mc <- iter_mc + 1
            newpep$mc <- iter_mc
            x[r,] <- newpep
            r <- r + 1
          } else {
            iter_mc <- 0
            r <- r + 1
          }
        }
      }
      d[[el]] <- x
    }
  }  
  cat("Row-bind missed cleavages\n")
  peptable <- as.data.frame(data.table::rbindlist(d))
  cat("Table done\n")
  
  names(peptable)[names(peptable) == "peptide"] <- "PepSequence"
  names(peptable)[names(peptable) == "start"] <- "PepStart"
  names(peptable)[names(peptable) == "stop"] <- "PepStop"
  peptable$ID <- paste(peptable$PepSequence, peptable$PTMPos, peptable$PTMType, sep="_")
  peptable <- peptable[order(paste(peptable$Accession, peptable$PepStart)),]
  cat("Start filtering\n")
  #FILTERING
  pepLength <- nchar(peptable$PepSequence)
  peptable <- peptable[(pepLength >= parameters$PepMinLength & pepLength <= parameters$PepMaxLength),]
  peptable <- peptable[!is.na(peptable$ID),]
  
  return(peptable)
  
}

mapQuanToDigestionProd <- function(DigestedProt) {
  # For each unique peptide, get the sum of intensities in all conditions
  # >> Do we want to change all I in L or reverse? << #
  unique_pep <- names(table(DigestedProt$ID)[table(DigestedProt$ID) == 1])
  uniquetab <- DigestedProt[DigestedProt$ID %in% unique_pep,]
  redundanttab <- DigestedProt[!(DigestedProt$ID %in% unique_pep),]
  
  iterID <- unique(redundanttab$ID)
  list_quan <- vector(mode = "list")
  for (i in seq_along(iterID)) {
    df <- redundanttab[redundanttab$ID == iterID[i],]
    mtx <- df[, grepl( "^C_" , names( df ) ) ]
    vec <- sapply(1:ncol(mtx),function(x) {
      log2(sum(2^as.numeric(mtx[,x]), na.rm=T))
    })
    df$Accession[1] <- paste(df$Accession, collapse = ";")
    df$ID[1] <- paste(sort(df$ID), collapse = ";")
    df <- df[1,]
    df[, grepl( "^C_" , names( df ) ) ] <- vec
    list_quan[[i]] <- df
  }
  
  redundanttab <- as.data.frame(data.table::rbindlist(list_quan))
  peptable <- as.data.frame(data.table::rbindlist(list(redundanttab, uniquetab)))
  return(peptable)
}

filterDigestedProt <- function(DigestedProt, parameters) {
  ## Enrichment loss:
  # >> Set up for phospho only << #
  ## Peptide lost at the enrichment step:
  enriched <- grepl("ph", DigestedProt$PTMType)
  if (sum(DigestedProt$Enriched) == 0) {
    return(list("NonEnriched" = DigestedProt, "Enriched" = NULL))
  } else {
    DigestedProt$Enriched <- enriched
    enrichedtab <- DigestedProt[DigestedProt$Enriched,]
    nonenrichedtab <- DigestedProt[!DigestedProt$Enriched,]
    numRemove <- floor(nrow(enrichedtab) * parameters$EnrichmentLoss)
    cat("Remove", numRemove, "modified peptides according to parameter EnrichmentLoss (", parameters$EnrichmentLoss, ")")
    idx <- sample(seq_len(nrow(enrichedtab)), size = numRemove)
    tokeep <- setdiff(seq_len(nrow(enrichedtab)), idx)
    enrichedtab <- enrichedtab[tokeep,]
    ## Add non-modified peptides to mimic pollution from non-modified peptides:
    numRemove <- floor((nrow(enrichedtab) * (1/parameters$EnrichmentEfficiency)) * (1 - parameters$EnrichmentEfficiency))
    idx <- sample(seq_len(nrow(nonenrichedtab)), size = numRemove)
    noise <- nonenrichedtab[idx,]
    noise$Enriched <- TRUE
    # Reduce signal of non-phospho:
    noise[, grepl( "^C_" , names( noise ) ) ] <- noise[, grepl( "^C_" , names( noise ) ) ] * (1 - parameters$EnrichmentNonModSignalLoss)
    enrichedtab <- rbind(enrichedtab, noise)
    ## Noise do to enrichment procedure:
    nrowTab <- nrow(enrichedtab)
    ncolTab <- length(grepl( "^C_" , names( DigestedProt ) ))
    mtx <- matrix(nrow = nrowTab, ncol= ncolTab, data=rnorm(n=ncolTab * nrowTab, mean=0, sd=parameters$EnrichmentNoise))
    enrichedtab[, grepl( "^C_" , names( enrichedtab ) ) ] <- apply(enrichedtab[, grepl( "^C_" , names( enrichedtab ) ) ], 2, as.numeric)
    enrichedtab[, grepl( "^C_" , names( enrichedtab ) ) ] <- enrichedtab[, grepl( "^C_" , names( enrichedtab ) ) ] + mtx
    return(list("NonEnriched" = nonenrichedtab, "Enriched" = enrichedtable))
  }
}




