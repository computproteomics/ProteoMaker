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
  peptable <- peptable %>% arrange(merged_accs)
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
  other_cols <- colnames(peptable)[!colnames(peptable) %in% parameters$QuantColnames]
  cat(" - built protein index for faster summarization\n")
  pb <- txtProgressBar(min=0, max=length(prot_ind))
  
  # Initiate and fill matrix with proteins
  protmat <- as.data.frame(matrix(ncol = ncol(peptable), nrow = length(prot_ind)))
  rownames(protmat) <- names(prot_ind)
  colnames(protmat) <- colnames(peptable)

  for (i in 1:(length(prot_ind)-1)){
    
    setTxtProgressBar(pb, i)
    
    tmp <- as.data.frame(peptable[prot_ind[i]:prot_ind[i+1],])
    rownames(tmp) <- tmp$Sequence
    
    # add other information
    protmat[i,other_cols] <- apply(tmp[,other_cols], 2, paste, collapse=";")

    tmp <- tmp[tmp$num_accs==1, parameters$QuantColnames]
    if (nrow(tmp) >= minUniquePep) {
      tmp <- as.matrix(tmp)
      
      if (method == "sum.top3") {
        
        tmp <- tmp[order(rowSums(tmp), decreasing = T),]
        
        if (nrow(tmp) >= 3) {
          
          protmat[i,parameters$QuantColnames] <- log2(colSums(2^tmp[1:3,], na.rm = T))
          
        } else {
          
          protmat[i,parameters$QuantColnames] <- log2(colSums(2^tmp, na.rm = T))
        }
        
      } else if (method == "medpolish"){
        summed <- colSummarizeMedianpolish(tmp)$Estimates
        if (length(summed) > 0)
          protmat[i,parameters$QuantColnames] <- summed
        
      } else {
        
      }
      
    }
    
  }
  
  close(pb)
  protmat[protmat == -Inf] <- NA
#  for (i in parameters$QuantColnames) protmat[,i] <- as.numeric(protmat[,i]) 

  return(protmat)
  
}




