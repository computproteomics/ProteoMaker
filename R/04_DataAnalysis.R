################################################################################
#               Data analysis of the outputs of the peptide quan.              #
################################################################################

library(preprocessCore)

proteinSummarisation <- function(peptable, parameters) {
    
    method <- parameters$ProtSummarization
    
    minUniquePep <- parameters$MinUniquePep
    
    QuantColnames <- parameters$QuantColnames
    
    cat(" + Remove all modified peptides\n")
    peptable <- peptable[sapply(peptable$PTMType, function(x) length(x) == 0),]
    cat("  - Remaining number of non-modified peptides:", nrow(peptable), "\n")
    
    cat(" + Protein summarisation\n")
    
    cat("  - Number of minimum unique peptide per protein:", minUniquePep, "\n")
    cat("  - Protein summarisation using the", method, "approach.\n")
    
    # writing new column with unlisted and merged protein names
    peptable$merged_accs <- sapply(peptable$Accession, function(x) paste(unique(unlist(x)), collapse=";"))
    peptable$num_accs <- sapply(peptable$Accession, function(x) length(unique(x)))
    
    # Sort table according to protein accession, needs to stay in this order!
    peptable <- peptable[order(peptable$merged_accs), ]
    cat("  - Sorted protein table\n")
    
    # Reducing table to relevant columns
    peptable <- peptable[peptable$num_accs == 1, ]
    
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
    other_cols <- colnames(peptable)[!colnames(peptable) %in% QuantColnames]
    cat("  - built protein index for faster summarization\n")

    # Initiate and fill matrix with proteins
    protmat <- as.data.frame(matrix(ncol = ncol(peptable), nrow = length(prot_ind)))
    rownames(protmat) <- names(prot_ind)
    colnames(protmat) <- colnames(peptable)
    
    cat("  - Initiated protein matrix for", length(prot_ind), "protein groups\n")
    cat("  - Summarizing proteins, this can take a while\n")
    
    # Function to summarize protein groups
    summarizeProtein <- function(tmp) {
        out <- NULL
        if (nrow(tmp) >= minUniquePep) {
            tmp <- as.matrix(tmp)
            if (method == "sum.top3") {
                tmp <- tmp[order(rowSums(tmp), decreasing = T),]
                if (nrow(tmp) >= 3) {
                    out <- log2(colSums(2^tmp[1:3,], na.rm = T))
                } else {
                    out <- log2(colSums(2^tmp, na.rm = T))
                }
            } else if (method == "medpolish"){
                summed <- NULL
                if (nrow(tmp) == 1) {
                    summed <- tmp
                } else {
                    summed <- medpolish(tmp, na.rm=T, trace.iter=F)$col
                }
                if (length(summed) > 0)
                    out <- summed
            } else {
                # Any other method to add
            }
        }
        return(out)
    }
    
    if (!is.null(parameters$Cores)) {
        cores <- parameters$Cores
        
        if (parallel::detectCores() <= parameters$Cores) {
            cores <- parallel::detectCores() - 1
        }
        
        cluster <- parallel::makeCluster(cores, type = parameters$ClusterType)
        parallel::setDefaultCluster(cluster)
        parallel::clusterExport(cluster, c("peptable","summarizeProtein","minUniquePep","prot_ind","other_cols","QuantColnames"), envir = environment())
        proteins <- parallel::parLapply(cluster, 1:(length(prot_ind)-1), function(i) {
            tmp <- as.data.frame(peptable[prot_ind[i]:(prot_ind[i+1]-1),])
            rownames(tmp) <- tmp$Sequence
            out <- tmp[1,]
            tout <- summarizeProtein(tmp[,QuantColnames,drop=F])
            if (!is.null(tout)) {
                out[QuantColnames] <- tout
                # add other information
                out[other_cols] <- sapply(tmp[,other_cols], function(x) paste(unlist(x), collapse=";"))
            } else {
                out <- NULL
            }
            return(out)
        })
        parallel::stopCluster(cluster)
    } else {
        proteins <- lapply(1:(length(prot_ind)-1), function(i) {
            tmp <- as.data.frame(peptable[prot_ind[i]:(prot_ind[i+1]-1),])
            rownames(tmp) <- tmp$Sequence
            out <- tmp[1,]
            tout <- summarizeProtein(tmp[,QuantColnames,drop=F])
            if (!is.null(tout)) {
            out[QuantColnames] <- tout
            # add other information
            out[other_cols] <- sapply(tmp[,other_cols], function(x) paste(unlist(x), collapse=";"))
            } else {
                out <- NULL
            }
            return(out)
        })
    }
    # join all protein data
    protmat <- do.call(rbind, proteins)
    protmat[protmat == -Inf] <- NA
    protmat <- protmat[rowSums(is.na(protmat[,QuantColnames])) < length(QuantColnames), ]
    
    cat ("  - Finished summarizing into", nrow(protmat) ,"proteins\n")
    #  for (i in parameters$QuantColnames) protmat[,i] <- as.numeric(protmat[,i]) 
    
    return(protmat)
    
}




