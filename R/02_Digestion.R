################################################################################
#                    PEPTIDE DIGESTION OF PROTEOFORM TABLE                     #
################################################################################

library(dplyr)
library(stringi)
library(crayon)
library(parallel)
library(purrr)
library(scales)

#####################
## Function to perform enzymatic digestion on a single protein sequence.
## - Supports 8 different cleavage rules, for trypsin (considering proline or just cleavage after lysine and arginine), chymotrypsin (high and low specificity),
##   pepsin (pH 2 and pH 1.3), lysC and argC.
#    Enzyme rules according to: https://www.nature.com/articles/srep22286/tables/1
#                               https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#exceptions
## - Digests the sequence to the maximum possible number of miss-cleavages based on the argument missed.
## - Filters proteolytic peptides according to their length.
## - Calculates the mass for every yielded peptide for charges +1, +2 and +3.
## - Filters the single charged peptides based on a upper and lower mass threshold.
## - If there are no cleavage sites on a sequence, or if the filtering removed or peptides, the output is NULL.
#####################
fastDigest <- function(sequence, enzyme = "trypsin", missed = 0, length.max = NA, length.min = NA, mass.max = NA, mass.min = NA) {
    cleavage_rules <- dplyr::case_when(
        enzyme == "trypsin" ~ "(?!(RP|KP))(?=(K|R))(?!(K|R)$)",
        enzyme == "trypsin.strict" ~ "(?=(K|R))(?!(K|R)$)",
        enzyme == "chymotrypsin.h" ~ "(?!(FP|YP|PY|WP))(?=(F|Y|W))(?!(F|Y|W)$)",
        enzyme == "chymotrypsin.l" ~ "(?!(FP|YP|PY|WP|LP|MP))(?=(F|Y|W|L|P))(?!(F|Y|W|L|P)$)",
        enzyme == "pepsin.2" ~ "(?=(F|L|W|Y|A|E|Q))(?!(F|L|W|Y|A|E|Q)$)",
        enzyme == "pepsin.1.3" ~ "(?=(F|L))(?!(F|L)$)",
        enzyme == "lysC" ~ "(?=(K))(?!(K)$)",
        enzyme == "argC" ~ "(?!(RP))(?=(R))(?!(R)$)"
    )
    
    if (!is.na(cleavage_rules)) {
        cleave <- stringi::stri_locate_all_regex(str = sequence, pattern = cleavage_rules)[[1]][, 1]
    } else {
        cat(crayon::red(enzyme, "is invalid enzyme argument!\n\n"))
        return(NULL)
    }
    
    start <- c(1, cleave + 1)
    stop <- c(cleave, nchar(sequence))
    
    if (is.na(cleave[1])) {
        return(NULL)
    }
    
    missed <- min(length(start) - 2, missed)
    
    peptides <- lapply(0:missed, function(x) {
        data.frame(
            Peptide = substring(sequence, start[1:(length(start) - x)], stop[(x + 1):length(stop)]),
            Start = start[1:(length(start) - x)],
            Stop = stop[(x + 1):length(stop)],
            MC = x, stringsAsFactors = F
        )
    })
    peptides <- dplyr::bind_rows(peptides)
    
    if (!is.na(length.max) & !is.na(length.min)) {
        lengths <- nchar(peptides$Peptide)
        peptides <- peptides[(lengths >= length.min & lengths <= length.max), ]
    }
    
    if (nrow(peptides) > 0) {
        AAs <- strsplit(peptides$Peptide, split = "")
        
        AA.mass <- c(
            "A" = 71.03711, "R" = 156.10111, "N" = 114.04293, "D" = 115.02694, "C" = 103.00919,
            "E" = 129.04259, "Q" = 128.05858, "G" = 57.02146, "H" = 137.05891, "I" = 113.08406,
            "L" = 113.08406, "K" = 128.09496, "M" = 131.04049, "F" = 147.06841, "P" = 97.05276,
            "S" = 87.03203, "T" = 101.04768, "W" = 186.07931, "Y" = 163.06333, "V" = 99.06841
        )
        
        peptide.mass <- sapply(AAs, function(x) sum(AA.mass[x], 18.01528))
        peptides$MZ1 <- peptide.mass + 1.007276466
        peptides$MZ2 <- (peptide.mass + (1.007276466 * 2)) / 2
        peptides$MZ3 <- (peptide.mass + (1.007276466 * 3)) / 3
    } else {
        return(NULL)
    }
    
    if (!is.na(mass.max) & !is.na(mass.min)) {
        peptides <- peptides[(peptides$MZ1 >= mass.min & peptides$MZ1 <= mass.max), ]
    }
    
    if (nrow(peptides) > 0) {
        rownames(peptides) <- 1:nrow(peptides)
    } else {
        return(NULL)
    }
    
    return(peptides)
}
#####################

#####################
## Function to perform enzymatic digestion on a set of proteoforms.
## - Calls fastDigest to digest a set of proteoforms.
## - Maps the modification sites on peptide sequences.
## - Adds the mass addition per peptide based on the modifications they have.
## - Peptide abundance depends on the parental proteoform.
#####################
proteoformDigestion <- function(proteoform, parameters) {
    peptides <- fastDigest(sequence = proteoform$Sequence, enzyme = parameters$Enzyme, missed = parameters$MaxNumMissedCleavages, length.max = parameters$PepMaxLength, length.min = parameters$PepMinLength)
    
    if (!is.null(peptides)) {
        peptides$Accession <- proteoform$Accession
        peptides$Proteoform_ID <- proteoform$Proteoform_ID
        
        peptides$PTMPos <- vector(mode = "list", length = nrow(peptides))
        peptides$PTMType <- vector(mode = "list", length = nrow(peptides))
        peptides$Regulation_Amplitude <- proteoform$Regulation_Amplitude
        peptides$Regulation_Pattern <- proteoform$Regulation_Pattern
        
        # Map modification sites on peptides.
        if (!is.null(proteoform$PTMPos[[1]])) {
            proteoform.position <- unlist(proteoform$PTMPos)
            proteoform.type <- unlist(proteoform$PTMType)
            pep.indices <- lapply(proteoform.position, function(x) which(x >= peptides$Start & x <= peptides$Stop))
            
            if (sum(lengths(pep.indices)) != 0) {
                pep.position <- lapply(1:length(pep.indices), function(x) sapply(pep.indices[[x]], function(y) proteoform.position[x] - peptides$Start[y] + 1))
                pep.type <- lapply(1:length(pep.indices), function(x) rep(proteoform.type[x], length(pep.indices[[x]])))
                
                to.aggregate <- data.frame(unlist(pep.indices), unlist(pep.position), unlist(pep.type), stringsAsFactors = F)
                to.aggregate <- aggregate(to.aggregate[, 2:3], by = list(to.aggregate[, 1]), FUN = list)
                
                # Calculate and add the mass addition due to modifications per modified peptide.
                modification.mass <- parameters$PTMTypesMass
                names(modification.mass) <- parameters$PTMTypes
                to.aggregate$mass_shift <- sapply(to.aggregate[, 3], function(x) sum(modification.mass[x], na.rm = T))
                
                peptides[to.aggregate[, 1], c("PTMPos", "PTMType")] <- to.aggregate[, 2:3]
                peptides[to.aggregate[, 1], c("MZ1", "MZ2", "MZ3")] <- peptides[to.aggregate[, 1], c("MZ1", "MZ2", "MZ3")] + as.numeric(to.aggregate[, 4]) %*% t(c(1, 0.5, 1 / 3))
            }
        }
        
        # Add proteoform abundance to all peptides.
        peptides.abundance <- as.data.frame(matrix(NA, ncol = length(parameters$QuantColnames), nrow = nrow(peptides)))
        colnames(peptides.abundance) <- parameters$QuantColnames
        peptides.abundance[1:nrow(peptides.abundance), parameters$QuantColnames] <- proteoform[parameters$QuantColnames]
        
        # Bind everything.
        peptides <- dplyr::bind_cols(peptides, peptides.abundance)
    }
    
    return(peptides)
}
#####################

#####################
## Function that wraps proteoformDigestion function to perform digestion on a set of proteoforms.
## - Gives the option of parallel computing.
## - Samples peptides based on a distribution derived from PropMissedCleavages according to their MC.
#####################
digestGroundTruth <- function(proteoforms, parameters) {
    cat("#PROTEOFORM DIGESTION - Start\n\n")
    cat(" + Digestion input:\n")
    cat("  - A total number of", nrow(proteoforms), "proteoforms, is proceed for proteolytic digestion.\n")
    cat("  - Unmodified fraction contains", sum(lengths(proteoforms$PTMType) == 0), "proteoforms and modified fraction", sum(lengths(proteoforms$PTMType) != 0), "proteoforms.\n")
    cat(
        "  - Cleavage will be performed by", parameters$Enzyme, "with a maximum of", parameters$MaxNumMissedCleavages,
        "miss-cleavages, to create peptides of length", parameters$PepMinLength, "to", parameters$PepMaxLength, "amino acids.\n"
    )
    
    if (!is.null(parameters$Cores)) {
        cores <- parameters$Cores
        
        if (parallel::detectCores() <= parameters$Cores) {
            cores <- parallel::detectCores() - 1
        }
        
        cluster <- parallel::makeCluster(cores, type = parameters$ClusterType)
        #on.exit(parallel::stopCluster(cluster))
        parallel::setDefaultCluster(cluster)
        parallel::clusterExport(cluster, c("proteoforms", "parameters", "proteoformDigestion", "fastDigest"), envir = environment())
        peptides <- parallel::parLapply(cluster, 1:nrow(proteoforms), function(x) proteoformDigestion(proteoform = proteoforms[x, ], parameters = parameters))
        parallel::stopCluster(cluster)
    } else {
        peptides <- lapply(1:nrow(proteoforms), function(x) proteoformDigestion(proteoform = proteoforms[x, ], parameters = parameters))
    }
    
    cat("  - All proteoforms are digested successfully!\n\n")
    
    peptides <- dplyr::bind_rows(peptides)
    
    # Sample peptides per MC by size determined by PropMissedCleavages.
    if (parameters$MaxNumMissedCleavages > 0) {
        if (parameters$PropMissedCleavages != 0 & parameters$PropMissedCleavages != 1) {
            MC.proportions <- sapply(0:parameters$MaxNumMissedCleavages, function(x) (1 - parameters$PropMissedCleavages)^2 * parameters$PropMissedCleavages^x)
            MC.proportions <- scales::rescale(x = MC.proportions, to = c(0, 1), from = c(0, max(MC.proportions, na.rm = T)))
            peptide.indices <- lapply(1:parameters$MaxNumMissedCleavages, function(x) which(peptides$MC == x))
            peptide.indices <- unlist(lapply(1:parameters$MaxNumMissedCleavages, function(x) sample(peptide.indices[[x]], size = ceiling(sum(peptides$MC == 0) * MC.proportions[x + 1]), replace = FALSE)))
            peptide.indices <- sort(c(which(peptides$MC == 0), peptide.indices))
            peptides <- peptides[peptide.indices, ]
        } else if (parameters$PropMissedCleavages == 0) {
            peptides <- peptides[peptides$MC == 0, ]
        }
    }
    
    cat(" + Digestion output:\n")
    cat("  - A total number of", nrow(peptides), "peptides is generated.\n")
    cat("  - Unmodified fraction contains", sum(lengths(peptides$PTMType) == 0), "peptides and modified fraction", sum(lengths(peptides$PTMType) != 0), "peptides.\n")
    cat(
        "  - The amount of peptides with", paste0(0:parameters$MaxNumMissedCleavages, collapse = ", "), "miss-cleavages is",
        paste0(sapply(0:parameters$MaxNumMissedCleavages, function(x) sum(peptides$MC == x)), collapse = ", "), "respectively.\n\n"
    )
    
    cat("#PROTEOFORM DIGESTION - Finish\n\n")
    
    return(peptides)
}
#####################

#####################
## Function that groups peptides and summarizes their abundance.
## - Creates unique ids for each peptide based on the modifications type and position (pep_id column).
## - Substitutes leucine to isoleucine residues (Sequence column).
## - Create groups of peptides and aggregates them.
##   The output structure is:
##
##   Sequence                 <character>  Contains the unique peptide sequence after I to L substitution.
##   Peptide                  <list of character vectors> Contains the peptides that are grouped based on Sequence, but prior to I to L substitution
##   Start                    <list of integer vectors> Contains the starting position of the peptides in the Peptide vectors on the protein sequence.
##   Stop                     <list of integer vectors> Contains the ending position of the peptides in the Peptide vectors on the protein sequence.
##   MC                       <list of integer vectors> Contains MC of the peptides in the Peptide vectors on the protein sequence.
##   MZ1                      <integer> Peptide mass for charge +1.
##   MZ2                      <integer> Peptide mass for charge +2.
##   MZ3                      <integer> Peptide mass for charge +3.
##   Accession                <list of character vectors> Contains the parental protein Accession of the peptides in Peptide vectors.
##   Proteoform_ID            <list of integer vectors> Contains the unique protoform identifiers of the Accession vectors.
##   PTMPos                   <list of integer vectors> Contains the position of the modifications on the peptides in the Peptide vectors.
##   PTMType                  <list of character vectors> Contains the modification types of the modifications of PTMPos vectors.
##   Regulation_Amplitude     <list of numeric vectors> Contains the regulation amplitudes of the proteoforms in Accession vectors.
##   Regulation_Pattern       <list of numeric vectors> Contains the regulation patterns of the proteoforms in Accession vectors.
##   Quantitative Columns     <numeric> Contains the abundances of the peptide group of Sequence.
## - Removes a percentage of randomly selected summarized peptides.
#####################
digestionProductSummarization <- function(peptides, parameters) {
    cat("#PEPTIDE SUMMARIZATION - Start\n\n")
    cat(" + Summarization input:\n")
    cat("  - A total number of", nrow(peptides), "peptides is proceed for summarization.\n")
    
    # Create unique ID for each peptide based on the PTMType and PTMPos. No aggregation technique in any package supports lists...
    peptides$pep_id <- as.character(mapply(list, peptides$PTMType, peptides$PTMPos, SIMPLIFY = F))
    
    cat("  - Unique peptide IDs are generated.\n")
    
    # Create a Sequence column where isoleucine is substituted by leucine.
    peptides$Sequence <- gsub("[I]", "L", peptides$Peptide)
    
    cat("  - Isoleucine substitution to leucine is done.\n")
    
    # Helping functions for summarization per group for specific columns.
    log2.sum <- function(x) {
        x <- log2(sum(2^x, na.rm = T))
        if (is.finite(x)) {
            return(x)
        } else {
            return(NA)
        }
    }
    
    select.first <- function(x) {
        return(x[1])
    }
    
    # Create groups based on Sequence and pep_id
    peptides <- dplyr::group_by(.data = peptides, Sequence, pep_id)
    
    cat("  - Peptide groups are generated.\n")
    
    peptides.1 <- peptides %>% dplyr::summarise(
        Peptide = list(Peptide),
        Start = list(Start),
        Stop = list(Stop),
        MC = list(MC),
        MZ1 = select.first(MZ1),
        MZ2 = select.first(MZ2),
        MZ3 = select.first(MZ3),
        Accession = list(Accession),
        Proteoform_ID = list(Proteoform_ID),
        PTMPos = select.first(PTMPos),
        PTMType = select.first(PTMType),
        Regulation_Amplitude = list(Regulation_Amplitude),
        Regulation_Pattern = list(Regulation_Pattern)
    )
    
    
    peptides.2 <- peptides %>% dplyr::summarise_at(.vars = dplyr::vars(parameters$QuantColnames), .funs = c("log2.sum"))
    
    peptides <- dplyr::inner_join(peptides.1, peptides.2, by = c("Sequence", "pep_id"))
    peptides <- dplyr::select(peptides, -c("pep_id"))
    
    cat("  - Peptide groups summarization is done.\n")
    
    # Remove a percentage of randomly selected summarized peptides.
    remove <- sample(1:nrow(peptides), size = nrow(peptides) * parameters$LeastAbundantLoss, replace = FALSE)
    
    if (length(remove) != 0) {
        peptides <- peptides[-remove, ]
    }
    
    cat("  - Remove", parameters$LeastAbundantLoss * 100, "% of the least abundant peptides, which corresponds to", length(remove), "peptides.\n\n")
    
    cat(" + Summarization output:\n")
    cat("  - A total number of ", nrow(peptides), "summarized peptides is generated.\n\n")
    cat("#PEPTIDE SUMMARIZATION - Finish\n\n")
    
    return(peptides)
}
#####################

#####################
## Function that creates enriched and non-enriched fractions of the proteolytic peptides.
## - Separates modified peptides and reduce the set based on the EnrichmentLoss parameter.
## - Simulates enrichment efficiency as an addition of non-modified peptides to the enriched set.
##   The abundance of these spiked peptides is reduced according to EnrichmentEfficiency parameter.
## - Introduces noise due to modification enrichment process to the abundances of the enrich fraction.
## - Returns two a list that contains the two sets.
#####################
filterDigestedProt <- function(DigestedProt, parameters) {
    enriched <- lengths(DigestedProt$PTMType) != 0
    
    if (sum(enriched) == 0) {
        return(list("NonEnriched" = DigestedProt, "Enriched" = NULL))
    } else {
        cat ("#ENRICHMENT SIMULATION - Start\n\n")
        ## Exact copy of "sample"
        enrichedtab <- DigestedProt
        
        ## Removing fraction according to EnrichmentLoss parameter
        numRemove <- floor(nrow(enrichedtab) * parameters$EnrichmentLoss)
        cat(" + Enrichment loss:\n")
        cat("  - Remove", numRemove, "peptides according to parameter EnrichmentLoss (", parameters$EnrichmentLoss, ")\n")
        idx <- sample(seq_len(nrow(enrichedtab)), size = numRemove, replace = F)
        enrichedtab <- enrichedtab[-idx, ]
        
        # Calculate total sum of intensities for modified and non-modified peptides in enriched fraction
        modified <- lengths(enrichedtab$PTMType) != 0
        enrichedtab_modified <- enrichedtab[modified, parameters$QuantColnames]
        enrichedtab_nonmodified <- enrichedtab[!modified, parameters$QuantColnames]
        totalModified <- sum(enrichedtab_modified, na.rm = T)
        totalNonModified <- sum(enrichedtab_nonmodified, na.rm = T)
        
        ## Adjust the intensities of modified peptides to mimic enrichment efficiency without loss of signal
        enrichedtab[!modified, parameters$QuantColnames] <- (totalModified + totalNonModified) *
            (1-parameters$EnrichmentEfficiency) / totalNonModified
        enrichedtab[modified, parameters$QuantColnames] <- (totalModified + totalNonModified) *
            parameters$EnrichmentEfficiency / totalModified
        cat(" + Enrichment efficiency:\n")
        
        cat("Enrichment efficiency is", parameters$EnrichmentEfficiency, "leading to the modified
        peptides contributing to", parameters$EnrichmentEfficiency*100, " % of the total intensity.")
        
        ## Noise due to enrichment procedure:
        nrowTab <- nrow(enrichedtab)
        ncolTab <- length(parameters$QuantColnames)
        cat(" + Enrichment noise:\n")
        cat("  - The enrichment noise standard deviation is", parameters$EnrichmentNoise, ".\n")
        mtx <- matrix(nrow = nrowTab, ncol = ncolTab, data = rnorm(n = ncolTab * nrowTab, mean = 0, sd = parameters$EnrichmentNoise))
        enrichedtab[, parameters$QuantColnames] <- enrichedtab[, parameters$QuantColnames] + mtx
        cat("  - Noise added to all samples!\n\n")
        cat("#ENRICHMENT SIMULATION - Finish\n\n")
        
        return(list("NonEnriched" = DigestedProt, "Enriched" = enrichedtab))
    }
}
#####################
