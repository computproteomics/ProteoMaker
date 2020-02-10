################################################################################
#               Data analysis of the outputs of the peptide quan.              #
################################################################################

proteinSummarisation <- function(peptable, method = "sum.top3", minUniquePep = 2, parameters) {
  
  cat("Number of minimum unique peptide per protein:", minUniquePep, "\n")
  
  if (method == "sum.top3") {
    cat("Protein sumarisation using the", method, "approach.\n")
    
    uAcc <- sort(unique(peptable$Accession))
    protmat <- matrix(ncol = parameters$NumCond*parameters$NumReps, nrow = length(uAcc))
    row.names(protmat) <- uAcc
    
    for (acc in uAcc) {
      tmp <- peptable[peptable$Accession == acc,grepl("^C_", names(peptable))]
      if (nrow(tmp) >= minUniquePep) {
        tmp <- tmp[order(rowSums(tmp), decreasing = T),]
        if (nrow(tmp) >= 3) {
          protmat[row.names(protmat) == acc,] <- log2(colSums(2^tmp[1:3,]))
        } else {
          protmat[row.names(protmat) == acc,] <- log2(colSums(2^tmp))
        }
        # protmat[row.names(protmat) == acc,] <- sapply(1:ncol(tmp), function(x) {
        #   mean(tmp[,x], na.rm = T)
        # })
      }
    }
    colnames(protmat) <- names((peptable[,grepl("^C_", names(peptable))]))
  } else {
    
  }
  
  return(protmat)
  
}
