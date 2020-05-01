################################################################################
#               Data analysis of the outputs of the peptide quan.              #
################################################################################

library(preprocessCore)
  
proteinSummarisation <- function(peptable, parameters) {
  
  method <- parameters$ProtSummarization
  
minUniquePep <- parameters$MinUniquePep
  
  cat("Number of minimum unique peptide per protein:", minUniquePep, "\n")
  
  cat("Protein summarisation using the", method, "approach.\n")
  
  uAcc <- sort(unique(unlist(peptable$Accession)))
  protmat <- matrix(ncol = length(parameters$QuantColnames), nrow = length(uAcc))
  
  pb <- txtProgressBar(min=0, max=length(uAcc))
  
  for (i in 1:length(uAcc)){
    
    setTxtProgressBar(pb, i)
    tmp <- peptable %>% filter(Accession == uAcc[i]) %>% select(c("Sequence",parameters$QuantColnames))
    #tmp <- peptable[peptable$Accession == uAcc[i], parameters$QuantColnames]
    
    tmp <- as.data.frame(tmp)
    rownames(tmp) <- tmp$Sequence
    tmp <- as.matrix(tmp[,parameters$QuantColnames])
    print(tmp)

    if (nrow(tmp) >= minUniquePep) {
      
      if (method == "sum.top3") {
        
        tmp <- tmp[order(rowSums(tmp), decreasing = T),]
        
        if (nrow(tmp) >= 3) {
          
          protmat[i,] <- log2(colSums(2^tmp[1:3,], na.rm = T))
          
        } else {
          
          protmat[i,] <- log2(colSums(2^tmp, na.rm = T))
        }
        
      } else if (method == "medpolish"){
        summed <- colSummarizeMedianpolish(tmp)$Estimates
        if (length(summed) > 0)
        protmat[i,] <- summed
        
      } else {
        
      }
    }
  }
  
  close(pb)
  protmat[protmat == -Inf] <- NA
  colnames(protmat) <- parameters$QuantColnames
  protmat <- data.frame(Accession = uAcc, protmat)
  
  
  return(protmat)
  
}




