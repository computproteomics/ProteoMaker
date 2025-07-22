################################################################################
#                  GENERATE PROTEOFORM-CENTRIC GROUND TRUTH                    #
################################################################################

#' Import and Process Protein Sequence Data for Simulation
#'
#' This function imports a single protein sequence FASTA file and optionally a text file containing selected protein accessions. It processes the protein sequences by filtering out invalid entries, separating them into modified and unmodified sets, and validating the input files.
#'
#' @param parameters A list containing the following elements:
#' \describe{
#'   \item{PathToFasta}{A character string specifying the path to the FASTA file containing protein sequences.}
#'   \item{PathToProteinList}{A character string specifying the path to the text file containing a list of selected protein accessions. If `NA`, a fraction of proteins will be selected randomly for modification.}
#'   \item{FracModProt}{A numeric value between 0 and 1 indicating the fraction of proteins to be modified. Only applicable if `PathToProteinList` is `NA`.}
#'   \item{ModifiableResidues}{A list of characters indicating the residues that can be modified in the protein sequences.}
#' }
#'
#' @return A list with two elements:
#' \describe{
#'   \item{to.Modify}{A data frame containing the protein sequences selected for modification, along with their corresponding metadata. If no proteins are selected for modification, this will be `NULL`.}
#'   \item{to.be.Unmodified}{A data frame containing the protein sequences that remain unmodified, along with their corresponding metadata. If no proteins remain unmodified, this will be `NULL`.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Imports and validates the FASTA file.
#'   \item Filters out protein sequences containing unusual amino acids.
#'   \item Removes duplicate protein accessions.
#'   \item Fractionates the protein set into modified and unmodified fractions, based on either a specified percentage or a list of selected protein accessions.
#'   \item If the protein set is too small to be fractionated, it returns only an unmodified set.
#'   \item If the input files are corrupted or the protein set is empty after filtering, it returns `NULL` for both modified and unmodified sets.
#' }
#'
#' @keywords internal
#'
#'
proteinInput <- function(parameters) {
  message(" + Importing data:")
  # Check is the fasta file can be loaded.
  error <- try(protr::readFASTA(file = parameters$PathToFasta, legacy.mode = TRUE, seqonly = FALSE), silent = TRUE)

  if (!inherits(error, "try-error")) {

    # Read fasta.
    fasta <- protr::readFASTA(file = parameters$PathToFasta, legacy.mode = TRUE, seqonly = FALSE)
    fasta <- data.frame(Sequence = unlist(fasta), Accession = sub(".*[|]([^.]+)[|].*", "\\1", names(fasta)), stringsAsFactors = F)
    rownames(fasta) <- 1:nrow(fasta)
    message("  - File ", parameters$PathToFasta, " imported, containing ", nrow(fasta), " protein sequences.")

    # Filter non-expressed proteins (random)
    if (parameters$PercExpressedProt < 1) {
      to.keep <- round(parameters$PercExpressedProt * nrow(fasta))
      fasta <- fasta[sample(1:nrow(fasta), size = to.keep), ]
      message("  - A total of ", nrow(fasta), " protein sequences have been randomly selected to be expressed.")
    }

    # Filter proteins carrying unusual amino acids.
    knownAA <- c("A", "L", "R", "K", "N", "M", "D", "F", "C", "P", "E", "S", "Q", "T", "G", "W", "H", "Y", "I", "V")
    unknownAA <- setdiff(LETTERS, knownAA)
    initialRows <- nrow(fasta)
    fasta <- fasta[if (initialRows > 1) {
      rowSums(sapply(unknownAA, grepl, x = fasta$Sequence)) == 0
    } else {
      sum(sapply(unknownAA, grepl, x = fasta$Sequence)) == 0
    }, ]
    message("  - A total of ", initialRows - nrow(fasta), " protein sequences have been removed due to unusual amino acids ", paste0("(", paste0(unknownAA, collapse = ","), ")"), " composition.")

    # Remove duplicated protein accessions.
    message("  - A total of ", length(which(duplicated(fasta$Accession))), " duplicated protein accessions have been removed.")
    fasta <- fasta[!duplicated(fasta$Accession), ]
    message("  - Total number of remaining protein sequences: ", nrow(fasta), "\n")

    # Proceed unless fasta df is empty after filtering.
    if (nrow(fasta) > 0) {
      message(" + Creating modified and unmodified fractions:")

      # Returns a list of the number of modifiable residues on all sequences, for each residue of ModifiableResidues.
      possible.modifiable.AAs <- as.data.frame(sapply(fasta$Sequence, function(x) {
        string <- strsplit(x, split = "")
        return(sapply(unlist(parameters$ModifiableResidues), function(y) length(which(string[[1]] == y))))
      }))

      rownames(possible.modifiable.AAs) <- 1:nrow(possible.modifiable.AAs)

      # Add an additional column, denoting the number of ModifiableResidues residue types found for each protein.
      possible.modifiable.AAs$Modifiable <- sapply(1:nrow(possible.modifiable.AAs), function(x) length(which(possible.modifiable.AAs[x, ] > 0)))

      # Proteins that contain at least one of the ModifiableResidues.
      unmodifiable.indices <- which(possible.modifiable.AAs$Modifiable == 0)

      # Proteins that contain at least a residue of ModifiableResidues, and thus are modifiable.
      modifiable.indices <- setdiff(1:nrow(fasta), unmodifiable.indices)

      message("  - A total of ", length(unmodifiable.indices), " sequences are unmodifiable.")

      # If there is not a specific protein list to be modified.
      if (is.na(parameters$PathToProteinList)) {
        # If the set of proteins can be fractionated. This covers the case when parameters$FracModProt = 0
        if ((parameters$FracModProt * length(modifiable.indices)) >= 1) {
          # Randomly select a FracModProt % of modifiable proteins (modifiable.indices).
          to.modify.indices <- sample(modifiable.indices, size = parameters$FracModProt * length(modifiable.indices))
          to.Modify <- fasta[to.modify.indices, ]
          rownames(to.Modify) <- 1:nrow(to.Modify)

          message("  - A total of ", length(modifiable.indices), " sequences are modifiable, from which ", parameters$FracModProt * 100, "% randomly selected to be modified.")

          to.be.Unmodified <- fasta[-to.modify.indices, ]

          # Add additional columns to the unmodified fraction df.
          to.be.Unmodified$PTMPos <- vector(mode = "list", length = nrow(to.be.Unmodified))
          to.be.Unmodified$PTMType <- vector(mode = "list", length = nrow(to.be.Unmodified))

          message("  - Modified fraction: ", nrow(to.Modify), " proteins.")
          message("  - Unmodified fraction: ", nrow(to.be.Unmodified), " proteins.")

          # Empty dfs replaced with NULL.This covers the case when parameters$FracModProt = 1
          if (nrow(to.be.Unmodified) > 0) {
            rownames(to.be.Unmodified) <- 1:nrow(to.be.Unmodified)
          } else {
            to.be.Unmodified <- NULL
          }

          return(list(to.Modify = to.Modify, to.be.Unmodified = to.be.Unmodified))
        } else {
          to.be.Unmodified <- fasta
          to.be.Unmodified$PTMPos <- vector(mode = "list", length = nrow(to.be.Unmodified))
          to.be.Unmodified$PTMType <- vector(mode = "list", length = nrow(to.be.Unmodified))

          message("  - The protein set cannot be fractionated, too low number of proteins or too low FracModProt %.")
          message("  - Modified fraction: 0 proteins.")
          message("  - Unmodified fraction: ", nrow(to.be.Unmodified), " proteins.\n")

          return(list(to.Modify = NULL, to.be.Unmodified = to.be.Unmodified))
        }
      } else { # If there is a specific protein list to be modified.

        # Check is the protein list file can be loaded.
        error <- try(read.csv(parameters$PathToProteinList, header = F, stringsAsFactors = F), silent = TRUE)

        if (!inherits(error, "try-error")) {
          protein.list.input <- unique(as.vector(read.csv(parameters$PathToProteinList, header = F, stringsAsFactors = F)[, 1]))
          message("  - Protein list file ", parameters$PathToProteinList, " loaded and contains ", length(protein.list.input), " unique protein accessions.")

          mapping <- unlist(lapply(protein.list.input, function(x) which(fasta$Accession == x)))
          to.modify.indices <- intersect(mapping, modifiable.indices)

          if (length(to.modify.indices) > 0) {
            message(" - A total of " , length(mapping), " protein accessions have found in fasta file from which ", length(to.modify.indices), " are modifiable.")
            to.Modify <- fasta[to.modify.indices, ]
            rownames(to.Modify) <- 1:nrow(to.Modify)
            to.be.Unmodified <- fasta[-to.modify.indices, ]

            # Add additional columns to the unmodified fraction df.
            to.be.Unmodified$PTMPos <- vector(mode = "list", length = nrow(to.be.Unmodified))
            to.be.Unmodified$PTMType <- vector(mode = "list", length = nrow(to.be.Unmodified))

            message("  - Modified fraction: ", nrow(to.Modify), " proteins.")
            message("  - Unmodified fraction: ", nrow(to.be.Unmodified), " proteins.\n")

            if (nrow(to.be.Unmodified) > 0) {
              rownames(to.be.Unmodified) <- 1:nrow(to.be.Unmodified)
            } else {
              to.be.Unmodified <- NULL
            }

            return(list(to.Modify = to.Modify, to.be.Unmodified = to.be.Unmodified))
          } else {
            to.be.Unmodified <- fasta
            to.be.Unmodified$PTMPos <- vector(mode = "list", length = nrow(to.be.Unmodified))
            to.be.Unmodified$PTMType <- vector(mode = "list", length = nrow(to.be.Unmodified))
            message("  - Proteins in ", parameters$PathToProteinList, " are not modifiable or are not found in ", parameters$PathToFasta, ".")
            message("  - Modified fraction: 0 proteins.")
            message("  - Unmodified fraction: ", nrow(to.be.Unmodified), " proteins.\n")

            return(list(to.Modify = NULL, to.be.Unmodified = to.be.Unmodified))
          }
        } else {
          message(crayon::red("  - Protein list file ", parameters$PathToProteinList, " couldn't be loaded!\n"))
          return(list(to.Modify = NULL, to.be.Unmodified = NULL))
        }
      }
    } else {
      message(crayon::red("  - There are no protein sequences left!\n"))
      return(list(to.Modify = NULL, to.be.Unmodified = NULL))
    }
  } else {
    message(crayon::red("  - Fasta file ", parameters$PathToFasta, " couldn't be loaded!\n"))
    return(list(to.Modify = NULL, to.be.Unmodified = NULL))
  }
}
#####################

#' Perform Modifications on Protein Sequences
#'
#' This function performs modifications on a fraction of protein sequences selected by the `proteinInput` function. It generates proteoforms based on the provided parameters, applies modifications, and optionally retains some sequences in their unmodified form.
#'
#' @param to.Modify A data frame containing the protein sequences that have been selected for modification. This is typically the output from the `proteinInput` function.
#' @param parameters A list containing the following elements:
#' \describe{
#'   \item{PTMTypes}{A character vector specifying the types of post-translational modifications (PTMs) to be applied.}
#'   \item{PTMTypesDistr}{A numeric vector specifying the background frequency distribution of each PTM type.}
#'   \item{ModifiableResidues}{A list of character vectors specifying the residues that can be modified for each PTM type.}
#'   \item{ModifiableResiduesDistr}{A list of numeric vectors specifying the background frequency distribution of each residue type for modification.}
#'   \item{PropModPerProt}{A numeric value indicating the proportion of proteoforms to generate per protein sequence.}
#'   \item{RemoveNonModFormFrac}{A numeric value indicating the fraction of modified sequences that should retain their unmodified counterparts.}
#' }
#'
#' @return A list with two elements:
#' \describe{
#'   \item{mod.proteoforms}{A data frame containing the modified proteoforms, including details of the modification positions and types.}
#'   \item{unmod.proteoforms}{A data frame containing the unmodified counterparts of the selected protein sequences, if any. If no sequences are retained unmodified, this will be `NULL`.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Selects a proportion of the input sequences to generate multiple proteoforms per sequence.
#'   \item Applies the specified modifications to the selected proteoforms.
#' \item Summarizes the types and residues of modifications applied to each proteoform.
#'   \item Optionally retains a fraction of the input sequences in their unmodified form, based on the specified parameters.
#' }
#'
#' @keywords internal
#'
performModification <- function(to.Modify, parameters) {

  ptmtypes <- unlist(parameters$PTMTypes)
  ptmtypes <- ptmtypes[!is.na(ptmtypes)]

  if (length(ptmtypes) > 0) {
    message(" + Performing modification:")
    message(
      "  - Selected modification type(s) ", paste0('"', paste0(ptmtypes, collapse = '", "'), '"'),
      " with background frequency distribution of ", paste0(paste0(unlist(parameters$PTMTypesDistr) * 100, "%"), collapse = ", "), " respectively."
    )


    for (i in ptmtypes) {
      message(
        "  - For modification ", paste0('"', i, '"'), " residue(s) ", paste0(parameters$ModifiableResidues[[1]][[i]], collapse = ", "),
        " can be modified with background frequency distribution of ", paste0(paste0(parameters$ModifiableResiduesDistr[[1]][[i]] * 100, "%"), collapse = ", "), " respectively."
      )
    }

  # All proteins in to.Modify set are selected to be modified at least once.
  # Then randomly a set of proteoforms from the imported to.Modify set is selected based on the fraction multiplier PropModPerProt - 1.
  # The size of mod.proteoforms set depends to PropModPerProt. (When 1 a proteoform set of size equal to to.Modify is created, when 2 the size is double etc etc)
  mod.proteoforms <- to.Modify

  if (parameters$PropModPerProt >= 2) {
    mod.proteoforms <- rbind(mod.proteoforms, to.Modify[sample(x = 1:nrow(to.Modify), size = (parameters$PropModPerProt - 1) * nrow(to.Modify), replace = T), ])
    mod.proteoforms <- mod.proteoforms[order(mod.proteoforms$Accession), ]
  }

  # Create additional columns for modification positions and modification type.
  mod.proteoforms$PTMPos <- vector(mode = "list", length = nrow(mod.proteoforms))
  mod.proteoforms$PTMType <- vector(mode = "list", length = nrow(mod.proteoforms))

  # Modify the selected proteoforms.
  selected.modifications <- modify(mod.proteoforms$Sequence, parameters)
  #selected.modifications <- as.data.frame(do.call(rbind, selected.modifications))

  # Fill the columns.
  mod.proteoforms$PTMPos <- selected.modifications$Positions
  mod.proteoforms$PTMType <- selected.modifications$Types

  # Summarize modification types and residues modified by each type.
  Type_Count <- AA_Count <- list()
  for (ptm in ptmtypes) {
    AA_Count[[ptm]] <- list()
    for (res in parameters$ModifiableResidues[[1]][[ptm]]) {
      AA_Count[[ptm]][[res]] <- sum(unlist(sapply(selected.modifications$Count, function(x) x[[ptm]][[res]])))
    }
    Type_Count[[ptm]] <- sum(unlist(AA_Count[[ptm]]))
  }
  # print(count.per.AAs)

  message("  - Sequences modified!")

  for (ptm in ptmtypes) {
    message(
      "  - For modification ", paste0('"', ptm, '"'), " and residue(s) ", paste0(parameters$ModifiableResidues[[1]][[ptm]], collapse = ", "),
      " the resulting relative frequency distribution is ", paste0(paste0(round(100 * unlist(AA_Count[[ptm]])/sum(unlist(AA_Count[[ptm]])), 3), "%"), collapse = ", "), " respectively."
    )
  }

  message(
    "  - The resulting relative frequency distribution for modification type(s) ", paste0('"', paste0(ptmtypes, collapse = '", "'), '"'),
    " is ", paste0(paste0(round(100 * unlist(Type_Count)/sum(unlist(Type_Count)), 3), "%"), collapse = ", "), " respectively."
  )

  # Select a  fraction of selected to-modify sequences, to remain unmodified too.
  unmodified.proteoforms.indices <- sample(1:nrow(to.Modify), size = (1 - parameters$RemoveNonModFormFrac) * nrow(to.Modify))
  unmod.proteoforms <- to.Modify[unmodified.proteoforms.indices, ]
  unmod.proteoforms$PTMPos <- vector(mode = "list", length = nrow(unmod.proteoforms))
  unmod.proteoforms$PTMType <- vector(mode = "list", length = nrow(unmod.proteoforms))

  message("  - A fraction of ", nrow(unmod.proteoforms), " modified protein sequences will maintain their unmodified counterpart ", paste0("(", (1 - parameters$RemoveNonModFormFrac) * 100, "%)."))

  return(list(mod.proteoforms = mod.proteoforms, unmod.proteoforms = unmod.proteoforms))

} else {
  message(" + No modifications selected.")
  return(list(mod.proteoforms = NULL, unmod.proteoforms = NULL))
}
}


####################
#' Modify single sequence
#'
#' This function modified based on the modification type and residue frequencies.
#'
#' @param seq A character vector of protein sequences to be modified.
#' @param pars A list of parameters for the modification process, including:
#'  \itemize{
#'  \item \code{ModifiableResidues}: A list of amino acid residues that can be modified for each PTM type.
#'  \item \code{PTMTypes}: A character vector of modification types.
#'  \item \code{PTMTypesDistr}: A numeric vector representing the background frequencies for each PTM type.
#'  \item \code{ModifiableResiduesDistr}: A list of numeric vectors representing the background frequencies
#'  for each modifiable residue within each PTM type.
#'  \item \code{PTMMultipleLambda}: A numeric value specifying the lambda parameter for the truncated Poisson distribution,
#'  which determines the number of modifications per sequence.
#'  }
#'
#'
#'  keywords internal
modify_seq <- function(seq, pars) {
  # Simplify parameters
  pmod_res <- pars[[1]]
  ptms <- pars[[2]]
  pmod_res_distr <- pars[[3]]
  ptms_distr <- pars[[4]]
  lambda <- pars[[5]]

  seq_string <- strsplit(seq, split = "")

  # Run over modifications
  out <- setNames(
    lapply(ptms, function(x) {
      weight <- pmod_res_distr[[x]]
      names(weight) <- pmod_res[[x]]
      # Old way, not sure whether necessary
      # weight <- length(unlist(x[[y]])) / lengths(x[[y]]) * weight
      weight[!is.finite(weight)] <- 0

      # All modifiable residues
      all_positions <- setNames(lapply(pmod_res[[x]], function(y) which(y == seq_string[[1]])), pmod_res[[x]])

      # Adjust by residue frequency in protein
      weight <- sapply(pmod_res[[x]], function(y) {
        ttt <-  weight[[y]] / length(all_positions[[y]]) * length(unlist(all_positions))
        ttt[ttt > 1] <- 1
        return(ttt)
      })

      # Run over residue type to get modified residues
      out <- lapply(pmod_res[[x]], function(y) {
        # Find positions of residues
        positions <- all_positions[[y]]
        # residues <- seq_string[[1]][positions]
        pos_len <- length(positions)
        # no modification site
        if(pos_len == 0) return(NULL)
        # getting modified sites
        out <- NULL
        if (pos_len > 1) {
          out <- sample(
            x = positions,
            size = weight[[y]] *  extraDistr::rtpois(n = 1, lambda = lambda * pos_len,
                                                     a = 1, b = pos_len)
          )
        } else if (runif(1) < weight[[y]]) {
          out <- positions
        }
        return(out)
      })
      names(out) <- pmod_res[[x]]
      nummod.per.residue <- lapply(out, function(y) {
        length(y)
      })
      return(list(modified = unlist(out), count = nummod.per.residue))
    }), ptms)
  nummod.per.ptm <- lapply(out, function(x) {
    x$count
  })
  names(nummod.per.ptm) <- ptms
  Types <- unlist(lapply(ptms, function(x) rep(x, sum(unlist(nummod.per.ptm[[x]])))))
  Positions <- unlist(lapply(ptms , function(x) unlist(out[[x]]$modified)))
  return(list(Positions = Positions, Types = Types, Count = nummod.per.ptm))

}


#####################

#' Modify Protein Sequences Based on Modification Type and Residue Frequencies
#'
#' This function modifies protein sequences by introducing post-translational modifications (PTMs) based on
#' specified modification types and residue-specific background frequencies. The function identifies all
#' possible modification sites within each sequence, adjusts the probability of modification according to
#' the global background frequencies, and then samples the modification sites based on these probabilities.
#'
#' The number of modifications per sequence is determined using a truncated Poisson distribution, ensuring
#' a realistic distribution of modifications. The function returns a report detailing the modification positions,
#' types of modifications, and statistics on the number of modified amino acids per modification type.
#'
#' @param seq A character vector of protein sequences to be modified.
#' @param param A list of parameters for the modification process, including:
#'   \itemize{
#'     \item \code{ModifiableResidues}: A list of amino acid residues that can be modified for each PTM type.
#'     \item \code{PTMTypes}: A character vector of modification types.
#'     \item \code{PTMTypesDistr}: A numeric vector representing the background frequencies for each PTM type.
#'     \item \code{ModifiableResiduesDistr}: A list of numeric vectors representing the background frequencies
#'     for each modifiable residue within each PTM type.
#'     \item \code{PTMMultipleLambda}: A numeric value specifying the lambda parameter for the truncated Poisson distribution,
#'     which determines the number of modifications per sequence.
#'   }
#'
#' @return A list of lists, where each sublist corresponds to a sequence and contains:
#'   \itemize{
#'     \item \code{Positions}: A numeric vector of positions in the sequence where modifications were made.
#'     \item \code{Types}: A character vector of the modification types applied at each position.
#'     \item \code{Count}: A list of numeric vectors representing the count of modifications per amino acid residue
#'     for each PTM type.
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Identifies candidate modification sites in each sequence based on the specified modifiable residues.
#'   \item Adjusts the probability weights for each modification type and candidate residue according to the
#'   background frequencies and the composition of each sequence.
#'   \item Samples the modification sites based on the adjusted probabilities and a truncated Poisson distribution.
#'   \item Generates a report for each sequence, detailing the modification positions, types, and the number
#'   of modified amino acids per modification type.
#' }
#'
#' @keywords internal
modify <- function(seq, param) {

  # Simplify parameters
  pmod_res <- param$ModifiableResidues[[1]]
  pmod_res_distr <- param$ModifiableResiduesDistr[[1]]
  ptms <- param$PTMTypes[[1]]
  ptms_distr <- param$PTMTypesDistr[[1]]
  # In case of multiple modification types, amino acid background frequences
  # for each type are proportionally adjusted to background frequences of modification types.
  pmod_res_distr <- setNames(lapply(ptms,
                                    function(x) unlist(pmod_res_distr[[x]]) *
                                      ptms_distr[[x]]), ptms)

  pars <- list(pmod_res, ptms, pmod_res_distr, ptms_distr, param$PTMMultipleLambda)
  # Run over all sequences
  reported.modifications <- lapply(seq, function(x) modify_seq(x, pars))

  # Combine all lists into a single data frame
  reported.modifications <- data.frame(Positions = I(lapply(reported.modifications, '[[', 'Positions')),
                                       Types = I(lapply(reported.modifications, '[[', 'Types')),
                                       Counts = I(lapply(reported.modifications, '[[', 'Count')))

  return(reported.modifications)
}
#####################

#' Prepare a Sample of Proteoforms
#'
#' This function generates a set of proteoforms by using the `proteinInput` and `performModification` functions.
#' It returns a data frame that includes unmodified proteins, modified proteoforms, and modifiable but unmodified proteins.
#'
#'
#' @param parameters A list of parameters required for sample preparation, including paths to input files and modification details.
#'
#' @return A data frame containing all proteoforms generated, including:
#' \itemize{
#'   \item Unmodified proteins,
#'   \item Modified proteoforms,
#'   \item Modifiable but unmodified proteins.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Imports the protein sequences using `proteinInput`.
#'   \item Modifies the selected proteins using `performModification`.
#'   \item Combines unmodified and modified proteoforms into a single data frame.
#' }
#'
#' @keywords internal
samplePreparation <- function(parameters) {
  message("#SAMPLE PREPARATION - Start\n")

  protein.Sets <- proteinInput(parameters = parameters)

  if (is.null(protein.Sets$to.Modify) & is.null(protein.Sets$to.be.Unmodified)) {
    message(crayon::red("#SAMPLE PREPARATION - Finish (Warning: Empty sets)"))
    proteoforms <- NULL
  } else {
    if (!is.null(protein.Sets$to.Modify)) {
      proteoforms.after.modification <- performModification(to.Modify = protein.Sets$to.Modify, parameters = parameters)

      proteoforms <- rbind(
        protein.Sets$to.be.Unmodified,
        proteoforms.after.modification$mod.proteoform,
        proteoforms.after.modification$unmod.proteoforms
      )

      proteoforms$Proteoform_ID <- make.unique(proteoforms$Accession)
      proteoforms$Proteoform_ID <- rowSums(cbind(as.numeric(sub("\\.", "", sub("^[^.]*", "", proteoforms$Proteoform_ID))), rep(1, nrow(proteoforms))), na.rm = T)
      proteoforms <- proteoforms[, c(1, 2, 5, 3, 4)]

      rownames(proteoforms) <- 1:nrow(proteoforms)

      message("#SAMPLE PREPARATION - Finish\n")
    } else {
      proteoforms <- protein.Sets$to.be.Unmodified
      rownames(proteoforms) <- 1:nrow(proteoforms)
      proteoforms$PTMPos <- vector(mode = "list", length = nrow(proteoforms))
      proteoforms$PTMType <- vector(mode = "list", length = nrow(proteoforms))
      proteoforms$Proteoform_ID <- rep(1, nrow(proteoforms))
      proteoforms <- proteoforms[, c(1, 2, 5, 3, 4)]

      message("#SAMPLE PREPARATION - Finish\n")
    }
  }

  return(proteoforms)
}
#####################

#' Create a Regulation Pattern
#'
#' This function generates a regulation pattern for a given number of conditions. The pattern
#' consists of random values representing up-regulation or down-regulation for each condition.
#'
#' @param NumCond An integer specifying the number of conditions.
#'
#' @return A numeric vector of length `NumCond` representing the regulation pattern.
#'
#'
#' @keywords internal
createRegulationPattern <- function(NumCond) {
  # select 0.5 because division by 2 is already induced this way
  # this ensures that the differentiation amplitude is as speciefied by the user
  regulation_direction <- sample(unlist(lapply(1:NumCond, function(x) c(0.5, -0.5) * x))[1:NumCond], size = NumCond)
  return(regulation_direction)
}
#####################

#' Add Abundance Values to Proteoforms
#'
#' This function adds abundance values to proteoforms based on a combination of random noise and
#' differentially regulated patterns. It supports both relative and absolute quantification.
#'
#' @param proteoforms A data frame of proteoforms to which abundance values will be added.
#' @param parameters A list of parameters controlling the abundance assignment, including the
#' number of conditions, quantification columns, noise levels, and regulation patterns.
#'
#' @return A data frame with abundance values added to the proteoforms.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Adds random noise to the abundance columns.
#'   \item If specified, selects differentially regulated proteoforms and adjusts their abundance
#'   based on regulation patterns.
#'   \item Optionally, removes values below a threshold or adjusts values to an absolute quantification scale.
#' }
#'
#' @keywords internal
addProteoformAbundance <- function(proteoforms, parameters) {
  # populate the matrix with random noise
  for (name in parameters$QuantColnames) {
    proteoforms[name] <- rnorm(n = nrow(proteoforms), mean = 0, sd = parameters$QuantNoise)
  }

  if (is.na(parameters$UserInputFoldChanges_NumRegProteoforms) | is.na(parameters$UserInputFoldChanges_RegulationFC) |
      parameters$UserInputFoldChanges_NumRegProteoforms == 0 | parameters$UserInputFoldChanges_RegulationFC == 0) {
    if (parameters$DiffRegFrac != 0) {
      # select differentially regulated proteoforms
      diff_reg_indices <- sample(1:nrow(proteoforms), size = parameters$DiffRegFrac * nrow(proteoforms))

      # print(diff_reg_indices)

      # determine amplitude of regulation for regulated proteoforms
      proteoforms[diff_reg_indices, "Regulation_Amplitude"] <- runif(min = 0, max = parameters$DiffRegMax, n = length(diff_reg_indices))

      regulationPatterns <- lapply(1:length(diff_reg_indices), function(x) createRegulationPattern(parameters$NumCond))
    } else {
      proteoforms$Regulation_Amplitude <- vector(mode = "list", length = nrow(proteoforms))
      regulationPatterns <- NULL
    }
  } else {
    # select differentially regulated proteoforms

    diff_reg_indices <- sample(1:nrow(proteoforms), size = sum(parameters$UserInputFoldChanges_NumRegProteoforms))

    # determine amplitude of regulation for regulated proteoforms
    proteoforms[diff_reg_indices, "Regulation_Amplitude"] <- parameters$UserInputFoldChanges_RegulationFC


    regulationPatterns <- lapply(1:length(diff_reg_indices), function(x) createRegulationPattern(parameters$NumCond))
  }

  proteoforms$Regulation_Pattern <- vector(mode = "list", length = nrow(proteoforms))

  # print(length(regulationPatterns))
  # print(regulationPatterns)

  if (!is.null(regulationPatterns)) {
    proteoforms$Regulation_Pattern[diff_reg_indices] <- regulationPatterns
    # [diff_reg_indices, "Regulation_Pattern"]
    proteoforms[diff_reg_indices, parameters$QuantColnames] <-
      # add regulation pattern*regulation amplitude to random noise
      proteoforms[diff_reg_indices, parameters$QuantColnames, drop = FALSE] +

      # generate regulation patterns for all regulated proteoforms
      t(sapply(1:length(diff_reg_indices), function(x) {
        rep(regulationPatterns[[x]], each = parameters$NumReps)

        # multiply regulation pattern with Regulation amplitude
      })) * (proteoforms[diff_reg_indices, "Regulation_Amplitude"])
  } else {
    proteoforms$Regulation_Pattern <- vector(mode = "list", length = nrow(proteoforms))
  }

  message(" + Add quan. distribution: Relative -> absolute")
  vec <- rnorm(n = nrow(proteoforms), mean = parameters$AbsoluteQuanMean, sd = parameters$AbsoluteQuanSD)
  for (name in parameters$QuantColnames) {
    proteoforms[name] <- proteoforms[name] + vec
  }
  if (parameters$ThreshNAQuantileProt > 0) {
    # Remove Values below the threshold set in the Parameters file
    thresh <- quantile(x = vec, probs = parameters$ThreshNAQuantileProt)
    proteoforms[, parameters$QuantColnames][proteoforms[, parameters$QuantColnames] < thresh] <- NA
  }

  return(proteoforms)
}
#####################

