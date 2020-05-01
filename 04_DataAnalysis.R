################################################################################
#               Data analysis of the outputs of the peptide quan.              #
################################################################################

library(preprocessCore)

proteinSummarisation <- function(peptable, parameters) {
  
  method <- parameters$ProtSummarization
  
  minUniquePep <- parameters$MinUniquePep
  
  cat("Number of minimum unique peptide per protein:", minUniquePep, "\n")
  
  cat("Protein summarisation using the", method, "approach.\n")
  
  # writing new column with unlisted and merged protein names
  peptable <- cbind(peptable , merged_accs=sapply(peptable$Accession, function(x) paste(unlist(x),collapse = ";")), 
                    num_accs = sapply(peptable$Accession, length))
  
  # Sort table according to protein accession, needs to stay in this order!
  peptable <- peptable %>% arrange(merged_accs) %>% select(c("Sequence","merged_accs","num_accs",parameters$QuantColnames))
  cat(" - Sorted protein table\n")
  
  # Vector with row indices of protein groups
  all_accs <- peptable$merged_accs
  prot_ind <- 1
  names(prot_ind) <- all_accs[1]
  for (i in 2:nrow(peptable)) {
    if (all_accs[i-1] != all_accs[i]) {
      prot_ind <- c(prot_ind, i)
      names(prot_ind)[length(prot_ind)] <- all_accs[i]
    }
  }
  prot_ind <- c(prot_ind, nrow(peptable))
  cat(" - built protein index for faster summarization\n")
  pb <- txtProgressBar(min=0, max=length(prot_ind))
  
  # Initiate and fill matrix with proteins
  protmat <- matrix(ncol = length(parameters$QuantColnames), nrow = length(prot_ind))
  rownames(protmat) <- names(prot_ind)

  for (i in 1:(length(prot_ind)-1)){
    
    setTxtProgressBar(pb, i)
    
    tmp <- as.data.frame(peptable[prot_ind[i]:prot_ind[i+1],])
    rownames(tmp) <- tmp$Sequence
    
    tmp <- tmp[tmp$num_accs==1, parameters$QuantColnames]
    
    if (nrow(tmp) >= minUniquePep) {
      tmp <- as.matrix(tmp)
      
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
  protmat <- data.frame(Accession = names(prot_ind), protmat)
  protmat <- protmat[rowSums(is.na(protmat)) < ncol(protmat),]
  
  
  return(protmat)
  
}




