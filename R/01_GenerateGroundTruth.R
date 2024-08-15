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
#' @examples
#' parameters <- list(
#'   PathToFasta = "path/to/fasta.fasta",
#'   PathToProteinList = "path/to/protein_list.txt",
#'   FracModProt = 0.5,
#'   ModifiableResidues = list("S", "T", "Y")
#' )
#' result <- proteinInput(parameters)
#' 
proteinInput <- function(parameters) {
  cat(" + Importing data:\n")
  # Check is the fasta file can be loaded.
  error <- try(protr::readFASTA(file = parameters$PathToFasta, legacy.mode = TRUE, seqonly = FALSE), silent = TRUE)

  if (class(error) != "try-error") {
    # Read fasta.
    fasta <- protr::readFASTA(file = parameters$PathToFasta, legacy.mode = TRUE, seqonly = FALSE)
    fasta <- data.frame(Sequence = unlist(fasta), Accession = sub(".*[|]([^.]+)[|].*", "\\1", names(fasta)), stringsAsFactors = F)
    rownames(fasta) <- 1:nrow(fasta)
    cat("  - File", parameters$PathToFasta, "imported, containing", nrow(fasta), "protein sequences.\n")

    # Filter proteins carrying unusual amino acids.
    knownAA <- c("A", "L", "R", "K", "N", "M", "D", "F", "C", "P", "E", "S", "Q", "T", "G", "W", "H", "Y", "I", "V")
    unknownAA <- setdiff(LETTERS, knownAA)
    initialRows <- nrow(fasta)
    fasta <- fasta[if (initialRows > 1) {
      rowSums(sapply(unknownAA, grepl, x = fasta$Sequence)) == 0
    } else {
      sum(sapply(unknownAA, grepl, x = fasta$Sequence)) == 0
    }, ]
    cat("  - A total of", initialRows - nrow(fasta), "protein sequences have been removed due to unusual amino acids", paste0("(", paste0(unknownAA, collapse = ","), ")"), "composition.\n")

    # Remove duplicated protein accessions.
    cat("  - A total of", length(which(duplicated(fasta$Accession))), "duplicated protein accessions have been removed.\n")
    fasta <- fasta[!duplicated(fasta$Accession), ]
    cat("  - Total number of remaining protein sequences:", nrow(fasta), "\n\n")

    # Proceed unless fasta df is empty after filtering.
    if (nrow(fasta) > 0) {
      cat(" + Creating modified and unmodified fractions:\n")

      # Returns a list of the number of modifiable residues on all sequences, for each residue of ModifiableResidues.
      possible.modifiable.AAs <- as.data.frame(t(sapply(fasta$Sequence, function(x) {
        string <- strsplit(x, split = "")
        return(sapply(unlist(parameters$ModifiableResidues), function(y) length(which(string[[1]] == y))))
      })))

      rownames(possible.modifiable.AAs) <- 1:nrow(possible.modifiable.AAs)

      # Add an additional column, denoting the number of ModifiableResidues residue types found for each protein.
      possible.modifiable.AAs$Modifiable <- sapply(1:nrow(possible.modifiable.AAs), function(x) length(which(possible.modifiable.AAs[x, ] > 0)))

      # Proteins that contain at least one of the ModifiableResidues.
      unmodifiable.indices <- which(possible.modifiable.AAs$Modifiable == 0)

      # Proteins that contain at least a residue of ModifiableResidues, and thus are modifiable.
      modifiable.indices <- setdiff(1:nrow(fasta), unmodifiable.indices)

      cat("  - A total of", length(unmodifiable.indices), "sequences are unmodifiable.\n")

      # If there is not a specific protein list to be modified.
      if (is.na(parameters$PathToProteinList)) {
        # If the set of proteins can be fractionated. This covers the case when parameters$FracModProt = 0
        if ((parameters$FracModProt * length(modifiable.indices)) >= 1) {
          # Randomly select a FracModProt % of modifiable proteins (modifiable.indices).
          to.modify.indices <- sample(modifiable.indices, size = parameters$FracModProt * length(modifiable.indices))
          to.Modify <- fasta[to.modify.indices, ]
          rownames(to.Modify) <- 1:nrow(to.Modify)

          cat("  - A total of", length(modifiable.indices), "sequences are modifiable, from which", parameters$FracModProt * 100, "% randomly selected to be modified.\n")

          to.be.Unmodified <- fasta[-to.modify.indices, ]

          # Add additional columns to the unmodified fraction df.
          to.be.Unmodified$PTMPos <- vector(mode = "list", length = nrow(to.be.Unmodified))
          to.be.Unmodified$PTMType <- vector(mode = "list", length = nrow(to.be.Unmodified))

          cat("  - Modified fraction:", nrow(to.Modify), "proteins.\n")
          cat("  - Unmodified fraction:", nrow(to.be.Unmodified), "proteins.\n\n")

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

          cat("  - The protein set cannot be fractionated, too low number of proteins or too low FracModProt %.\n")
          cat("  - Modified fraction: 0 proteins.\n")
          cat("  - Unmodified fraction:", nrow(to.be.Unmodified), "proteins.\n\n")

          return(list(to.Modify = NULL, to.be.Unmodified = to.be.Unmodified))
        }
      } else { # If there is a specific protein list to be modified.

        # Check is the protein list file can be loaded.
        error <- try(read.csv(parameters$PathToProteinList, header = F, stringsAsFactors = F), silent = TRUE)

        if (class(error) != "try-error") {
          protein.list.input <- unique(as.vector(read.csv(parameters$PathToProteinList, header = F, stringsAsFactors = F)[, 1]))
          cat("  - Protein list file", parameters$PathToProteinList, "loaded and contains", length(protein.list.input), "unique protein accessions.\n")

          mapping <- unlist(lapply(protein.list.input, function(x) which(fasta$Accession == x)))
          to.modify.indices <- intersect(mapping, modifiable.indices)

          if (length(to.modify.indices) > 0) {
            cat(" - A total of", length(mapping), "protein accessions have found in fasta file from which", length(to.modify.indices), "are modifiable.\n")
            to.Modify <- fasta[to.modify.indices, ]
            rownames(to.Modify) <- 1:nrow(to.Modify)
            to.be.Unmodified <- fasta[-to.modify.indices, ]

            # Add additional columns to the unmodified fraction df.
            to.be.Unmodified$PTMPos <- vector(mode = "list", length = nrow(to.be.Unmodified))
            to.be.Unmodified$PTMType <- vector(mode = "list", length = nrow(to.be.Unmodified))

            cat("  - Modified fraction:", nrow(to.Modify), "proteins.\n")
            cat("  - Unmodified fraction:", nrow(to.be.Unmodified), "proteins.\n\n")

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
            cat("  - Proteins in", parameters$PathToProteinList, "are not modifiable or are not found in", parameters$PathToFasta, ".\n")
            cat("  - Modified fraction: 0 proteins.\n")
            cat("  - Unmodified fraction:", nrow(to.be.Unmodified), "proteins.\n\n")

            return(list(to.Modify = NULL, to.be.Unmodified = to.be.Unmodified))
          }
        } else {
          cat(crayon::red("  - Protein list file", parameters$PathToProteinList, "couldn't be loaded!\n\n"))
          return(list(to.Modify = NULL, to.be.Unmodified = NULL))
        }
      }
    } else {
      cat(crayon::red("  - There are no protein sequences left!\n\n"))
      return(list(to.Modify = NULL, to.be.Unmodified = NULL))
    }
  } else {
    cat(crayon::red("  - Fasta file", parameters$PathToFasta, "couldn't be loaded!\n\n"))
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
#'   \item{PTMTypesDist}{A numeric vector specifying the background frequency distribution of each PTM type.}
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
#' @examples
#' parameters <- list(
#'   PTMTypes = c("Phosphorylation", "Acetylation"),
#'   PTMTypesDist = c(0.7, 0.3),
#'   ModifiableResidues = list(c("S", "T", "Y"), c("K")),
#'   ModifiableResiduesDistr = list(c(0.5, 0.3, 0.2), c(1)),
#'   PropModPerProt = 2,
#'   RemoveNonModFormFrac = 0.5
#' )
#' to.Modify <- data.frame(
#'   Sequence = c("MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG"),
#'   Accession = c("P01112")
#' )
#' result <- performModification(to.Modify, parameters)
performModification <- function(to.Modify, parameters) {
  cat(" + Performing modification:\n")
  cat(
    "  - Selected modification type(s)", paste0('"', paste0(parameters$PTMTypes, collapse = '", "'), '"'),
    "with background frequency distribution of", paste0(paste0(parameters$PTMTypesDist * 100, "%"), collapse = ", "), "respectively.\n"
  )

  for (i in 1:length(parameters$ModifiableResidues)) {
    cat(
      "  - For modification", paste0('"', parameters$PTMTypes[i], '"'), "residue(s)", paste0(parameters$ModifiableResidues[[i]], collapse = ", "),
      "can be modified with background frequency distribution of", paste0(paste0(parameters$ModifiableResiduesDistr[[i]] * 100, "%"), collapse = ", "), "respectively.\n"
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
  selected.modifications <- as.data.frame(do.call(rbind, selected.modifications))

  # Fill the columns.
  mod.proteoforms$PTMPos <- selected.modifications$Positions
  mod.proteoforms$PTMType <- selected.modifications$Types

  # Summarize modification types and residues modified by each type.
  count.per.AAs <- Reduce(rbind, selected.modifications$Count)
  count.per.AAs <- lapply(1:ncol(count.per.AAs), function(x) colSums(Reduce(rbind, count.per.AAs[, x])))

  AAs.percentage <- sapply(count.per.AAs, function(x) sapply(x, function(y) y / sum(x) * 100))
  if (length(parameters$PTMTypes) == 1) {
    AAs.percentage <- list(as.vector(AAs.percentage))
  }
  Type.percentage <- sapply(count.per.AAs, function(x) sum(x)) / sum(unlist(count.per.AAs)) * 100

  cat("  - Sequences modified!\n")

  for (i in 1:length(parameters$ModifiableResidues)) {
    cat(
      "  - For modification", paste0('"', parameters$PTMTypes[i], '"'), "and residue(s)", paste0(parameters$ModifiableResidues[[i]], collapse = ", "),
      "the resulted frequency distribution is", paste0(paste0(round(AAs.percentage[[i]], 3), "%"), collapse = ", "), "respectively.\n"
    )
  }

  cat(
    "  - The resulted frequency distribution for modification type(s)", paste0('"', paste0(parameters$PTMTypes, collapse = '", "'), '"'),
    "is", paste0(paste0(round(Type.percentage, 3), "%"), collapse = ", "), "respectively.\n"
  )

  # Select a  fraction of selected to-modify sequences, to remain unmodified too.
  unmodified.proteoforms.indices <- sample(1:nrow(to.Modify), size = (1 - parameters$RemoveNonModFormFrac) * nrow(to.Modify))
  unmod.proteoforms <- to.Modify[unmodified.proteoforms.indices, ]
  unmod.proteoforms$PTMPos <- vector(mode = "list", length = nrow(unmod.proteoforms))
  unmod.proteoforms$PTMType <- vector(mode = "list", length = nrow(unmod.proteoforms))

  cat("  - A fraction of", nrow(unmod.proteoforms), "modified protein sequences will maintain their unmodified counterpart", paste0("(", (1 - parameters$RemoveNonModFormFrac) * 100, "%)."), "\n\n")

  return(list(mod.proteoforms = mod.proteoforms, unmod.proteoforms = unmod.proteoforms))
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
#'     \item \code{PTMTypesDist}: A numeric vector representing the background frequencies for each PTM type.
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
#' @examples
#' # Example usage:
#' sequences <- c("MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYTPGKLAG")
#' params <- list(
#'   ModifiableResidues = list(c("K", "R"), c("S", "T", "Y")),
#'   PTMTypes = c("Acetylation", "Phosphorylation"),
#'   PTMTypesDist = c(0.6, 0.4),
#'   ModifiableResiduesDistr = list(c(0.7, 0.3), c(0.4, 0.4, 0.2)),
#'   PTMMultipleLambda = 2
#' )
#' modification_report <- modify(sequences, params)
#'
#' @keywords internal
modify <- function(seq, param) {
  # Find the positions of candidate modification sites on the sequences.
  possible.modification.sites <- lapply(seq, function(x) {
    string <- strsplit(x, split = "")
    return(lapply(param$ModifiableResidues, function(y) lapply(y, function(z) which(string[[1]] == z))))
  })

  # In case of multiple modification types, amino acid background frequences for each type are proportionally adjusted to background frequences of modification types.
  param$ModifiableResiduesDistr <- lapply(1:length(param$PTMTypesDist), function(x) param$ModifiableResiduesDistr[[x]] * param$PTMTypesDist[x])

  # Find probability weights to scale the distribution of modification per residue for every sequence to fit the background frequences globaly.
  modification.probability.weight <- lapply(possible.modification.sites, function(x) {
    lapply(1:length(x), function(y) {
      weight <- length(unlist(x)) / lengths(x[[y]]) * param$ModifiableResiduesDistr[[y]]
      weight[!is.finite(weight)] <- 0
      return(weight)
    })
  })

  # Sample possible modification sites for each sequence based on the adjusted probability weights and size based on a truncated poisson distribution.
  selected.modification.sites <- lapply(1:length(possible.modification.sites), function(x) {
    if (length(unlist(possible.modification.sites[[x]])) > 1) {
      sort(sample(
        x = unlist(possible.modification.sites[[x]]),
        prob = rep(unlist(modification.probability.weight[[x]]), unlist(lapply(possible.modification.sites[[x]], function(y) lengths(y)))),
        size = extraDistr::rtpois(n = 1, lambda = param$PTMMultipleLambda * length(unlist(possible.modification.sites[[x]])), a = 1, b = length(unlist(possible.modification.sites[[x]])))
      ))
    } else {
      return(unlist(possible.modification.sites[[x]]))
    }
  })

  # Creating a report per sequence about the modification positions, types and the total number of modified AA per modification type.
  modification.report <- lapply(1:length(selected.modification.sites), function(x) {
    all.sites <- lapply(possible.modification.sites[[x]], function(y) unlist(y))

    modifications.per.AA <- lapply(possible.modification.sites[[x]], function(y) sapply(y, function(z) length(intersect(z, selected.modification.sites[[x]]))))

    selected.modification.positions <- lapply(all.sites, function(y) intersect(y, selected.modification.sites[[x]]))

    modification.types <- unlist(lapply(1:length(selected.modification.positions), function(y) rep(param$PTMTypes[y], lengths(selected.modification.positions)[y])))

    selected.modification.positions <- unlist(selected.modification.positions)

    reported.modifications <- list(
      Positions = selected.modification.positions[order(selected.modification.positions, decreasing = F)],
      Types = modification.types[order(selected.modification.positions, decreasing = F)],
      Count = modifications.per.AA
    )

    return(reported.modifications)
  })

  return(modification.report)
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
#' @examples
#' # Example usage:
#' parameters <- list(
#'   PathToFasta = "path/to/fasta/file",
#'   FracModProt = 0.5,
#'   PTMTypes = c("Phosphorylation"),
#'   PTMTypesDist = c(0.8),
#'   ModifiableResidues = list(c("S", "T", "Y")),
#'   ModifiableResiduesDistr = list(c(0.6, 0.3, 0.1)),
#'   PropModPerProt = 2,
#'   RemoveNonModFormFrac = 0.2
#' )
#' proteoforms <- samplePreparation(parameters)
#' @importFrom crayon red
#' @keywords internal
samplePreparation <- function(parameters) {
  cat("#SAMPLE PREPARATION - Start\n\n")

  protein.Sets <- proteinInput(parameters = parameters)

  if (is.null(protein.Sets$to.Modify) & is.null(protein.Sets$to.be.Unmodified)) {
    cat(crayon::red("#SAMPLE PREPARATION - Finish (Warning: Empty sets)"))
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

      cat("#SAMPLE PREPARATION - Finish")
    } else {
      proteoforms <- protein.Sets$to.be.Unmodified
      rownames(proteoforms) <- 1:nrow(proteoforms)
      proteoforms$PTMPos <- vector(mode = "list", length = nrow(proteoforms))
      proteoforms$PTMType <- vector(mode = "list", length = nrow(proteoforms))
      proteoforms$Proteoform_ID <- rep(1, nrow(proteoforms))
      proteoforms <- proteoforms[, c(1, 2, 5, 3, 4)]

      cat("#SAMPLE PREPARATION - Finish\n\n")
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
#' @examples
#' # Example usage:
#' regulation_pattern <- createRegulationPattern(NumCond = 3)
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
#' @examples
#' # Example usage:
#' parameters <- list(
#'   QuantColnames = c("Cond1", "Cond2"),
#'   QuantNoise = 0.2,
#'   DiffRegFrac = 0.3,
#'   DiffRegMax = 2,
#'   NumCond = 2,
#'   NumReps = 3
#' )
#' proteoforms <- addProteoformAbundance(proteoforms, parameters)
#'
#' @keywords internal
addProteoformAbundance <- function(proteoforms, parameters) {
  # populate the matrix with random noise
  for (name in parameters$QuantColnames) {
    proteoforms[name] <- rnorm(n = nrow(proteoforms), mean = 0, sd = parameters$QuantNoise)
  }

  if (is.na(parameters$UserInputFoldChanges_NumRegProteoforms)) {
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
      proteoforms[diff_reg_indices, parameters$QuantColnames] +

      # generate regulation patterns for all regulated proteoforms
      t(sapply(1:length(diff_reg_indices), function(x) {
        rep(regulationPatterns[[x]], each = parameters$NumReps)

        # multiply regulation pattern with Regulation amplitude
      })) * (proteoforms[diff_reg_indices, "Regulation_Amplitude"])
  } else {
    proteoforms$Regulation_Pattern <- vector(mode = "list", length = nrow(proteoforms))
  }

  if (is.na(parameters$AbsoluteQuanMean)) {
    # Remove Values below the threshold set in the Parameters file
    proteoforms[, parameters$QuantColnames][proteoforms[, parameters$QuantColnames] < parameters$ThreshNAProteoform] <- NA
  } else {
    cat("Add quan. distribution: Relative -> absolute\n")
    vec <- rnorm(n = nrow(proteoforms), mean = parameters$AbsoluteQuanMean, sd = parameters$AbsoluteQuanSD)
    for (name in parameters$QuantColnames) {
      proteoforms[name] <- proteoforms[name] + vec
    }
    if (parameters$ThreshNAQuantileProt > 0) {
      # Remove Values below the threshold set in the Parameters file
      thresh <- quantile(x = vec, probs = parameters$ThreshNAQuantileProt)
      proteoforms[, parameters$QuantColnames][proteoforms[, parameters$QuantColnames] < thresh] <- NA
    }
  }

  return(proteoforms)
}
#####################
