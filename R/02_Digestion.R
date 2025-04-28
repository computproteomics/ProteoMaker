################################################################################
#                    PEPTIDE DIGESTION OF PROTEOFORM TABLE                     #
################################################################################


#' Perform enzymatic digestion on a single protein sequence
#'
#' This function simulates enzymatic digestion on a single protein sequence, supporting multiple
#' cleavage rules for enzymes such as trypsin, chymotrypsin, pepsin, lysC, and argC. The function
#' returns peptides filtered by length and mass, including calculations for charges +1, +2, and +3.
#'
#' @param sequence A character string representing the protein sequence to be digested.
#' @param enzyme A character string specifying the enzyme to use for digestion. Default is "trypsin".
#' @param missed An integer specifying the maximum number of missed cleavages. Default is 0.
#' @param length.max An integer specifying the maximum peptide length. Default is NA.
#' @param length.min An integer specifying the minimum peptide length. Default is NA.
#' @param mass.max A numeric value specifying the maximum peptide mass. Default is NA.
#' @param mass.min A numeric value specifying the minimum peptide mass. Default is NA.
#'
#' @return A data frame with the digested peptides, including their mass/charge (M/Z) ratios for charges +1, +2, and +3.
#' If no peptides meet the criteria, the function returns \code{NULL}.
#'
#' @importFrom dplyr case_when bind_rows
#' @importFrom stringi stri_locate_all_regex
#' @importFrom crayon red
#' @keywords internal
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
    message(crayon::red(enzyme, " is invalid enzyme argument!"))
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

#' Perform enzymatic digestion on a set of proteoforms
#'
#' This function performs enzymatic digestion on a set of proteoforms, mapping modification sites
#' on peptide sequences and adding mass shifts per peptide based on modifications. The peptide
#' abundance is set based on the parental proteoform.
#'
#' @param proteoform A data frame containing proteoform sequences and associated data.
#' @param parameters A list containing various parameters for the digestion process,
#' including enzyme type, maximum number of missed cleavages, peptide length limits,
#' and modification masses.
#'
#' @return A data frame containing the digested peptides with mapped modifications,
#' mass shifts, and associated abundances. If no peptides are generated, the function
#' returns \code{NULL}.
#'
#' @importFrom dplyr bind_cols
#' @importFrom stats aggregate
#' @keywords internal
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
        to.aggregate <- stats::aggregate(to.aggregate[, 2:3], by = list(to.aggregate[, 1]), FUN = list)

        # Calculate and add the mass addition due to modifications per modified peptide.
        modification.mass <- parameters$PTMTypesMass[[1]]
        names(modification.mass) <- parameters$PTMTypes[[1]]
        to.aggregate$mass_shift <- sapply(to.aggregate[, 3], function(x) sum(unlist(modification.mass[x]), na.rm = T))

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

#' Perform proteoform digestion on a set of proteoforms with optional parallel computing
#'
#' This function wraps the \code{proteoformDigestion} function to perform enzymatic digestion
#' on a set of proteoforms. It supports parallel computing and allows sampling of peptides
#' based on a distribution derived from \code{PropMissedCleavages} according to their missed
#' cleavages (MC).
#'
#' @param proteoforms A data frame containing the proteoform sequences and associated data.
#' @param parameters A list of parameters including enzyme type, number of cores for parallel
#' computing, maximum number of missed cleavages, peptide length limits, and proportion of missed cleavages.
#'
#' @return A data frame containing the digested peptides, with details on modifications,
#' mass shifts, and missed cleavages.
#'
#' @importFrom dplyr bind_rows
#' @importFrom parallel makeCluster detectCores setDefaultCluster clusterExport stopCluster parLapply
#' @importFrom scales rescale
#' @keywords internal
digestGroundTruth <- function(proteoforms, parameters) {
  message("\n#PROTEOFORM DIGESTION - Start\n")
  message(" + Digestion input:")
  message("  - A total number of ", nrow(proteoforms), " proteoforms, is proceed for proteolytic digestion.")
  message("  - Unmodified fraction contains ", sum(lengths(proteoforms$PTMType) == 0), " proteoforms and modified fraction ", sum(lengths(proteoforms$PTMType) != 0), " proteoforms.")
  message(
    "  - Cleavage will be performed by ", parameters$Enzyme, " with a maximum of ", parameters$MaxNumMissedCleavages,
    " miss-cleavages, to create peptides of length ", parameters$PepMinLength, " to ", parameters$PepMaxLength, " amino acids."
  )

  # TODO REVERT
  parameters$Cores <- NULL
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

  message("  - All proteoforms are digested successfully!")

  peptides <- dplyr::bind_rows(peptides)

  # Sample peptides per MC by size determined by PropMissedCleavages.
  if (parameters$MaxNumMissedCleavages > 0) {
    if (parameters$PropMissedCleavages != 0 & parameters$PropMissedCleavages != 1) {
      # set max number of missed cleavages for probability calculation to min 5
      max_misscleav <- ifelse(parameters$MaxNumMissedCleavages < 5, 5, parameters$MaxNumMissedCleavages)
      MC.proportions <- sapply(0:parameters$MaxNumMissedCleavages, function(x) (1 - parameters$PropMissedCleavages)^(max_misscleav-x) *
                                 parameters$PropMissedCleavages^x)
      MC.proportions <- scales::rescale(x = MC.proportions, to = c(0, 1), from = c(0, max(MC.proportions, na.rm = T)))
      peptide.indices <- lapply(1:parameters$MaxNumMissedCleavages, function(x) which(peptides$MC == x))
      peptide.indices <- unlist(lapply(1:parameters$MaxNumMissedCleavages, function(x) sample(peptide.indices[[x]], size = ceiling(sum(peptides$MC == 0) * MC.proportions[x + 1]), replace = FALSE)))
      peptide.indices <- sort(c(which(peptides$MC == 0), peptide.indices))
      peptides <- peptides[peptide.indices, ]
    } else if (parameters$PropMissedCleavages == 0) {
      peptides <- peptides[peptides$MC == 0, ]
    }
  }

  message(" + Digestion output:")
  message("  - A total number of ", nrow(peptides), " peptides is generated.")
  message("  - Unmodified fraction contains ", sum(lengths(peptides$PTMType) == 0), " peptides and modified fraction ", sum(lengths(peptides$PTMType) != 0), " peptides.")
  message(
    "  - The amount of peptides with ", paste0(0:parameters$MaxNumMissedCleavages, collapse = ", "), " miss-cleavages is ",
    paste0(sapply(0:parameters$MaxNumMissedCleavages, function(x) sum(peptides$MC == x)), collapse = ", "), " respectively.\n"
  )

  message("#PROTEOFORM DIGESTION - Finish\n")

  return(peptides)
}
#####################


#' Calculate the detectability of a peptide sequence
#'
#' This function calculates the detectability of a peptide sequence based on the amino acid
#' using the PeptideRanger package. It returns a numeric vector representing the detectability
#' of the peptides. The prediction of the detectability is based on the amino acid composition
#' and does not take into account post-translational modifications.
#'
#' @param peptides A character vector containing the peptide sequences.
#' @param parameters A list of parameters including the number of cores for parallel computing.
#'
#' @return A numeric vector representing the detectability of the peptides.
#'
#' @importFrom PeptideRanger peptide_predictions
#' @keywords internal
addDetectability <- function(peptides, parameters) {

  # get unique peptide list and be able to map back
  unique_peptides <- unique(peptides)
  peptide_map <- match(peptides, unique_peptides)

  RFScores <- NULL
  if (!is.null(parameters$Cores)) {
    cores <- parameters$Cores

    if (parallel::detectCores() <= parameters$Cores) {
      cores <- parallel::detectCores() - 1
    }

    cluster <- parallel::makeCluster(cores, type = parameters$ClusterType)
    parallel::setDefaultCluster(cluster)

    # Ensure that the necessary package is loaded on each worker
    parallel::clusterEvalQ(cluster, library(PeptideRanger))

    # Split the data into chunks of 100 peptides
    peptide_chunks <- split(unlist(unique_peptides) , ceiling(seq_along(unlist(unique_peptides))/100))

    # Run the predictions in parallel
    RFScores <- parallel::parLapply(cluster, peptide_chunks, function(subset) {
      PeptideRanger::peptide_predictions(unlist(subset), PeptideRanger::RFmodel_ProteomicsDB)
    })

    # Combine the results into a single list or data frame
    RFScores <- do.call(rbind, RFScores)
    parallel::stopCluster(cluster)
  } else {
    RFScores<- PeptideRanger::peptide_predictions(unlist(unique_peptides), PeptideRanger::RFmodel_ProteomicsDB)
  }

  # Map back to the original peptides
  RFScores <- RFScores[peptide_map,]$RF_score
  return(RFScores)
}
#####################


#' Summarize digested peptide products
#'
#' This function groups peptides by unique identifiers, summarizes their abundance,
#' and optionally removes a percentage of the least abundant peptides. It creates
#' unique peptide IDs, substitutes isoleucine with leucine, and aggregates peptides
#' based on various characteristics.
#'
#' @param peptides A data frame containing the digested peptides and associated data.
#' @param parameters A list of parameters including QuantColnames and LeastAbundantLoss.
#'
#' @return A data frame containing the summarized peptides, with the following structure:
#' \describe{
#'   \item{Sequence}{A character vector containing the unique peptide sequence after isoleucine substitution to leucine.}
#'   \item{Peptide}{A list of character vectors containing the peptides that are grouped based on Sequence, prior to isoleucine substitution.}
#'   \item{Start}{A list of integer vectors containing the starting positions of the peptides in the Peptide vectors on the protein sequence.}
#'   \item{Stop}{A list of integer vectors containing the ending positions of the peptides in the Peptide vectors on the protein sequence.}
#'   \item{MC}{A list of integer vectors containing the number of missed cleavages (MC) for the peptides in the Peptide vectors.}
#'   \item{MZ1}{A numeric vector representing the peptide mass for charge +1.}
#'   \item{MZ2}{A numeric vector representing the peptide mass for charge +2.}
#'   \item{MZ3}{A numeric vector representing the peptide mass for charge +3.}
#'   \item{Accession}{A list of character vectors containing the parental protein Accession of the peptides in the Peptide vectors.}
#'   \item{Proteoform_ID}{A list of integer vectors containing the unique proteoform identifiers of the Accession vectors.}
#'   \item{PTMPos}{A list of integer vectors containing the positions of the modifications on the peptides in the Peptide vectors.}
#'   \item{PTMType}{A list of character vectors containing the modification types of the modifications in the PTMPos vectors.}
#'   \item{Regulation_Amplitude}{A list of numeric vectors containing the regulation amplitudes of the proteoforms in the Accession vectors.}
#'   \item{Regulation_Pattern}{A list of numeric vectors containing the regulation patterns of the proteoforms in the Accession vectors.}
#'   \item{Quantitative Columns}{Numeric columns containing the abundances of the peptide group for each QuantColname. These columns are dynamically named based on the provided QuantColnames parameter.}
#' }
#'
#' @importFrom dplyr group_by summarise summarise_at vars inner_join select %>%
#' @keywords internal
digestionProductSummarization <- function(peptides, parameters) {
  message("#PEPTIDE SUMMARIZATION - Start\n")
  message(" + Summarization input:")
  message("  - A total number of ", nrow(peptides), " peptides is proceed for summarization.")

  # Create unique ID for each peptide based on the PTMType and PTMPos. No aggregation technique in any package supports lists...
  peptides$pep_id <- as.character(mapply(list, peptides$PTMType, peptides$PTMPos, SIMPLIFY = F))

  message("  - Unique peptide IDs are generated.")

  # Create a Sequence column where isoleucine is substituted by leucine.
  peptides$Sequence <- gsub("[I]", "L", peptides$Peptide)

  message("  - Isoleucine substitution to leucine is done.")

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

  message("  - Peptide groups are generated.")

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

  message("  - Peptide groups summarization is done.")

  # Remove a percentage of randomly selected summarized peptides.
  remove <- sample(1:nrow(peptides), size = nrow(peptides) * parameters$LeastAbundantLoss, replace = FALSE)

  if (length(remove) != 0) {
    peptides <- peptides[-remove, ]
  }

  message("  - Remove ", parameters$LeastAbundantLoss * 100, "% of the least abundant peptides, which corresponds to ", length(remove), " peptides.")

  # add column with peptide detectability
  message("  - Calculating/predicting peptide detectability for later filtering with PeptideRanger.")
  peptides$Detectability <- addDetectability(peptides$Sequence, parameters)


  message(" + Summarization output:")
  message("  - A total number of ", nrow(peptides), " summarized peptides is generated.\n")
  message("#PEPTIDE SUMMARIZATION - Finish\n")

  return(peptides)
}
#####################

#' Create enriched and non-enriched fractions of proteolytic peptides
#'
#' This function separates modified peptides into enriched and non-enriched fractions,
#' adjusts the peptide abundances based on enrichment efficiency, and introduces noise
#' due to the enrichment process. It returns a list containing the enriched and
#' non-enriched peptide sets.
#'
#' @param DigestedProt A data frame containing the digested proteolytic peptides and associated data.
#' @param parameters A list of parameters that includes EnrichmentLoss, ModificationLoss, EnrichmentEfficiency,
#' EnrichmentNoise, and QuantColnames.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{NonEnriched}{A data frame containing the non-enriched peptide fraction, which includes both modified and non-modified peptides.}
#'   \item{Enriched}{A data frame containing the enriched peptide fraction, where modified peptides have been enriched based on the EnrichmentEfficiency, and noise has been added to simulate the enrichment process. If no modified peptides are present, this will be \code{NULL}.}
#' }
#'
filterDigestedProt <- function(DigestedProt, parameters) {
  modified <- lengths(DigestedProt$PTMType) != 0
  if (length(modified) == 0) modified <- NA

  ## Removing fraction according to ModificationLoss parameter
  numRemove <- 0
  if (sum(modified) > 0) {
    numRemove <- floor(sum(modified) * parameters$ModificationLoss)
  }
  message("\n#ENRICHMENT SIMULATION - Start\n")
  message(" + Modification loss")
  message("  - Remove ", numRemove, " modified peptides in non-enriched fraction according to parameter ModificationLoss (",
          parameters$ModificationLoss, ")")
  idx <- sample(which(modified), size = numRemove, replace = FALSE)
  nonenrichedtab <- DigestedProt
  if (length(idx) > 0)
    nonenrichedtab <- nonenrichedtab[-idx, ]

  if (sum(modified) == 0 | is.na(parameters$EnrichPTM) | parameters$EnrichmentEfficiency == 0) {
    message("\n#ENRICHMENT SIMULATION - Finish\n")
    return(list("NonEnriched" = nonenrichedtab, "Enriched" = NULL))
  } else {
    ## Exact copy of "sample"
    enrichedtab <- data.frame(DigestedProt)

    ## Removing fraction according to EnrichmentLoss parameter
    numRemove <- floor(nrow(enrichedtab) * parameters$EnrichmentLoss)
    message(" + Enrichment loss:")
    message("  - Remove ", numRemove, " peptides according to parameter EnrichmentLoss (", parameters$EnrichmentLoss, ")")
    idx <- sample(seq_len(nrow(enrichedtab)), size = numRemove, replace = FALSE)
    enrichedtab <- enrichedtab[-idx, ]

    # Select rows with PTM to be enriched
    modified <- sapply(enrichedtab$PTMType, function(x) sum(unlist(x) == parameters$EnrichPTM) > 0)

    message(" + Enriching PTM: ", parameters$EnrichPTM, ", having ", sum(modified), " peptides with this PTM in enriched fraction.")

    # Calculate total sum and average of intensities for modified and non-modified peptides in enriched fraction
    enrichedtab_modified <- enrichedtab[modified, parameters$QuantColnames]
    enrichedtab_nonmodified <- enrichedtab[!modified, parameters$QuantColnames]
    totalModified <- sum(unlist(enrichedtab_modified), na.rm = TRUE)
    averageModified <- mean(unlist(enrichedtab_modified), na.rm = TRUE)
    totalNonModified <- sum(unlist(enrichedtab_nonmodified), na.rm = TRUE)
    averageNonModified <- mean(unlist(enrichedtab_nonmodified), na.rm = TRUE)

    ## Adjust the intensities of modified peptides to mimic the mean difference of the parameter
    enrichedtab[!modified, parameters$QuantColnames] <- enrichedtab[!modified, parameters$QuantColnames] -
      parameters$EnrichmentQuantDiff
    enrichedtab[modified, parameters$QuantColnames] <- enrichedtab[modified, parameters$QuantColnames] +
      parameters$EnrichmentQuantDiff
    # Scale the total intensity of all peptides to the one before the adjustment
    enrichedtab[parameters$QuantColnames] <- enrichedtab[parameters$QuantColnames] *
      (totalModified + totalNonModified) / sum(unlist(enrichedtab[parameters$QuantColnames]), na.rm = TRUE)

    message(" + Enrichment efficiency:")
    message("  - Adjusted quantitative values to mimic the mean difference of the parameter EnrichmentQuantDiff (",
            parameters$EnrichmentQuantDiff, ").")
    # Calculate current fraction of modified peptides
    fracMod <- sum(modified) / nrow(enrichedtab)
    message("  - Modified peptides contribute to ",  fracMod * 100, "% of all present peptides in the enriched samples.")
    if (fracMod > parameters$EnrichmentEfficiency) {
      message("  - Warning: Enrichment efficiency is lower than the current fraction of modified peptides.
                No peptides will be removed")
    } else {
      # Remove the non-modified peptides from the enriched fraction to reach given enrichment efficiency
      numRemove <- floor(length(modified) - sum(modified) / parameters$EnrichmentEfficiency)
      message("  - Remove ", numRemove, " non-modified peptides to reach enrichment efficiency of ", parameters$EnrichmentEfficiency, "")
      idx <- sample(which(!modified), size = numRemove, replace = FALSE)
      enrichedtab <- enrichedtab[-idx, ]
      message("  - Enrichment efficiency is ", parameters$EnrichmentEfficiency, " leading to modified peptides contributing to ",
              parameters$EnrichmentEfficiency*100, "% of all present peptides in the enriched samples.")
    }

    ## Noise due to enrichment procedure:
    nrowTab <- nrow(enrichedtab)
    ncolTab <- length(parameters$QuantColnames)
    message(" + Enrichment noise:")
    message("  - The enrichment noise standard deviation is ", parameters$EnrichmentNoise, ".")
    mtx <- matrix(nrow = nrowTab, ncol = ncolTab, data = rnorm(n = ncolTab * nrowTab, mean = 0, sd = parameters$EnrichmentNoise))
    enrichedtab[, parameters$QuantColnames] <- enrichedtab[, parameters$QuantColnames] + mtx
    message("  - Noise added to all samples!")
    message("#ENRICHMENT SIMULATION - Finish")

    return(list("NonEnriched" = nonenrichedtab, "Enriched" = enrichedtab))
  }
}
#####################
