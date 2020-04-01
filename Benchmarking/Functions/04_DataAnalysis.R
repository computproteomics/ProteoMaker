################################################################################
#               Data analysis of the outputs of the peptide quan.              #
################################################################################

proteinSummarisation <- function(peptable, method = "sum.top3", minUniquePep = 2, parameters) {
  
  cat("Number of minimum unique peptide per protein:", minUniquePep, "\n")
  
  if (method == "sum.top3") {
    cat("Protein sumarisation using the", method, "approach.\n")
    
    uAcc <- sort(unique(unlist(peptable$Accession)))
    protmat <- matrix(ncol = length(parameters$QuantColnames), nrow = length(uAcc))
    
    for (i in 1:length(uAcc)){
      
      tmp <- peptable[which(sapply(peptable$Accession, function(x) any(x == uAcc[i]))), parameters$QuantColnames]
      
      if (nrow(tmp) >= minUniquePep) {
        
        tmp <- tmp[order(rowSums(tmp), decreasing = T),]
        
        if (nrow(tmp) >= 3) {
          
          protmat[i,] <- log2(colSums(2^tmp[1:3,], na.rm = T))
          
        } else {
          
          protmat[i,] <- log2(colSums(2^tmp, na.rm = T))
          
        }
      }
    }
    
    protmat[protmat == -Inf] <- NA
    colnames(protmat) <- parameters$QuantColnames
    protmat <- data.frame(Accession = uAcc, protmat)

  } else {
    
  }
  
  return(protmat)
  
}




