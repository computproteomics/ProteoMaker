################################################################################
#                  GENERATE PROTEOFORM-CENTRIC GROUND TRUTH                    #
################################################################################

library(protr)
library(crayon)
library(extraDistr)

#####################
## Function to import a single protein sequence fasta file or along with a txt file with selected protein accessions 
## and to create protein sets to be modified or remain unmodified.
## - Validates imported files (fasta, txt). (only possible warning is "incomplete final line found on" when fasta file doesnt have final white space line)
## - Filters out sequences with unusual amino acids.
## - Filters out duplicate protein accessions.
## - Fractionates initial protein set to unmodified set and modified set, either with FracModProt % or a list of selected protein accessions.
## - Checks if the protein set is too small to be fractionated in case of FracModProt parameter is used (will return only a set of proteins to remain unmodified)
##   and additionally works for fractions of modified proteins of 0 and 100%.
## - Checks if the input of the selected protein accessions in the txt file can be modified or are included in the fasta file. If not, only a set of proteins to remain unmodified is generated.
## - When any of the files is corrupted or protein set is empty after filtering then to.Modify = NULL and to.be.Unmodified = NULL.
## - Stop error function is not used not to interupt any potential interactive interface.
#####################

##### Future development #####
# - Import modification DB data.
##############################
proteinInput <- function(parameters){
  
  cat(" + Importing data:\n")
  
  #Check is the fasta file can be loaded.
  error <- try(protr::readFASTA(file = parameters$PathToFasta, legacy.mode = TRUE, seqonly = FALSE), silent = TRUE)
  
  if(class(error) != "try-error"){
    
    #Read fasta.
    fasta <- protr::readFASTA(file = parameters$PathToFasta, legacy.mode = TRUE, seqonly = FALSE)
    fasta <- data.frame(Sequence = unlist(fasta), Accession = sub(".*[|]([^.]+)[|].*", "\\1", names(fasta)), stringsAsFactors = F)
    rownames(fasta) <- 1:nrow(fasta)
    cat("  - File", parameters$PathToFasta, "imported, containing", nrow(fasta), "protein sequences.\n")
    
    #Filter proteins carrying unusual amino acids.
    knownAA <- c("A",	"L", "R",	"K", "N",	"M", "D", "F", "C",	"P", "E",	"S", "Q", "T", "G", "W", "H", "Y", "I", "V")
    unknownAA <- setdiff(LETTERS, knownAA)
    initialRows <- nrow(fasta)
    fasta <- fasta[if(initialRows > 1){ rowSums(sapply(unknownAA, grepl, x = fasta$Sequence)) == 0}
                   else { sum(sapply(unknownAA, grepl, x = fasta$Sequence)) == 0}, ]
    cat("  - A total of", initialRows - nrow(fasta) ,"protein sequences have been removed due to unusual amino acids", paste0("(", paste0(unknownAA, collapse = ","), ")"), "composition.\n")
    
    #Remove duplicated protein accessions.
    cat("  - A total of", length(which(duplicated(fasta$Accession))) ,"duplicated protein accessions have been removed.\n")
    fasta <- fasta[!duplicated(fasta$Accession),]
    cat("  - Total number of remaining protein sequences:", nrow(fasta), "\n\n")
    
    #Proceed unless fasta df is empty after filtering.
    if(nrow(fasta) > 0){
      
      cat(" + Creating modified and unmodified fractions:\n")
      
      #Returns a list of the number of modifiable residues on all sequences, for each residue of ModifiableResidues.
      possible.modifiable.AAs <- as.data.frame(t(sapply(fasta$Sequence, function(x){ 
        string = strsplit(x, split = "")
        return(sapply(unlist(parameters$ModifiableResidues), function(y) length(which(string[[1]] == y))))
      })))
      
      rownames(possible.modifiable.AAs) <- 1:nrow(possible.modifiable.AAs)
      
      #Add an additional column, denoting the number of ModifiableResidues residue types found for each protein.
      possible.modifiable.AAs$Modifiable <- sapply(1:nrow(possible.modifiable.AAs),function(x) length(which(possible.modifiable.AAs[x,] > 0))) 
      
      #Proteins that contain at least one of the ModifiableResidues.
      unmodifiable.indices <- which(possible.modifiable.AAs$Modifiable == 0)
      
      #Proteins that contain at least a residue of ModifiableResidues, and thus are modifiable.
      modifiable.indices <- setdiff(1:nrow(fasta), unmodifiable.indices)
      
      cat("  - A total of", length(unmodifiable.indices) ,"sequences are unmodifiable.\n")
      
      #If there is not a specific protein list to be modified.
      if(is.null(parameters$PathToProteinList)){
        
        #If the set of proteins can be fractionated. This covers the case when parameters$FracModProt = 0
        if((parameters$FracModProt * length(modifiable.indices)) >= 1){
          
          # Randomly select a FracModProt % of modifiable proteins (modifiable.indices).
          to.modify.indices <- sample(modifiable.indices, size = parameters$FracModProt * length(modifiable.indices))
          to.Modify <- fasta[to.modify.indices, ]
          rownames(to.Modify) <- 1:nrow(to.Modify)
          
          cat("  - A total of", length(modifiable.indices) ,"sequences are modifiable, from which", parameters$FracModProt*100,"% randomly selected to be modified.\n")
          
          to.be.Unmodified <- fasta[-to.modify.indices,]
          
          #Add additional columns to the unmodified fraction df.  
          to.be.Unmodified$PTMPos = vector(mode = "list", length = nrow(to.be.Unmodified))
          to.be.Unmodified$PTMType = vector(mode = "list", length = nrow(to.be.Unmodified))
          
          cat("  - Modified fraction:", nrow(to.Modify) ,"proteins.\n")
          cat("  - Unmodified fraction:", nrow(to.be.Unmodified) ,"proteins.\n\n")
          
          #Empty dfs replaced with NULL.This covers the case when parameters$FracModProt = 1
          if(nrow(to.be.Unmodified) > 0){
            
            rownames(to.be.Unmodified) <- 1:nrow(to.be.Unmodified)
            
          } else {
            
            to.be.Unmodified <- NULL
            
          }
          
          return(list(to.Modify = to.Modify, to.be.Unmodified = to.be.Unmodified))
          
        } else {
          
          to.be.Unmodified <- fasta
          to.be.Unmodified$PTMPos = vector(mode = "list", length = nrow(to.be.Unmodified))
          to.be.Unmodified$PTMType = vector(mode = "list", length = nrow(to.be.Unmodified))
          
          cat("  - The protein set cannot be fractionated, too low number of proteins or too low FracModProt %.\n")
          cat("  - Modified fraction: 0 proteins.\n")
          cat("  - Unmodified fraction:", nrow(to.be.Unmodified) ,"proteins.\n\n")
          
          return(list(to.Modify = NULL, to.be.Unmodified = to.be.Unmodified))
          
        }
        
      } else { #If there is a specific protein list to be modified.
        
        #Check is the protein list file can be loaded.
        error <- try(read.csv(parameters$PathToProteinList, header = F, stringsAsFactors = F), silent = TRUE)
        
        if(class(error) != "try-error"){
          
          protein.list.input <- unique(as.vector(read.csv(parameters$PathToProteinList, header = F, stringsAsFactors = F)[,1]))
          cat("  - Protein list file", parameters$PathToProteinList, "loaded and contains", length(protein.list.input),"unique protein accessions.\n")
          
          mapping <- unlist(lapply(protein.list.input, function(x) which(fasta$Accession == x)))
          to.modify.indices <- intersect(mapping, modifiable.indices)
          
          if(length(to.modify.indices) > 0){
            
            cat(" - A total of", length(mapping), "protein accessions have found in fasta file from which", length(to.modify.indices),"are modifiable.\n")
            to.Modify <- fasta[to.modify.indices, ]
            rownames(to.Modify) <- 1:nrow(to.Modify)
            to.be.Unmodified <- fasta[-to.modify.indices,]
            
            #Add additional columns to the unmodified fraction df.  
            to.be.Unmodified$PTMPos = vector(mode = "list", length = nrow(to.be.Unmodified))
            to.be.Unmodified$PTMType = vector(mode = "list", length = nrow(to.be.Unmodified))
            
            cat("  - Modified fraction:", nrow(to.Modify) ,"proteins.\n")
            cat("  - Unmodified fraction:", nrow(to.be.Unmodified) ,"proteins.\n\n")
            
            if(nrow(to.be.Unmodified) > 0){
              
              rownames(to.be.Unmodified) <- 1:nrow(to.be.Unmodified)
              
            } else {
              
              to.be.Unmodified <- NULL
              
            }
            
            return(list(to.Modify = to.Modify, to.be.Unmodified = to.be.Unmodified))
            
          } else {
            
            to.be.Unmodified <- fasta
            to.be.Unmodified$PTMPos = vector(mode = "list", length = nrow(to.be.Unmodified))
            to.be.Unmodified$PTMType = vector(mode = "list", length = nrow(to.be.Unmodified))
            cat("  - Proteins in", parameters$PathToProteinList, "are not modifiable or are not found in", parameters$PathToFasta, ".\n")
            cat("  - Modified fraction: 0 proteins.\n")
            cat("  - Unmodified fraction:", nrow(to.be.Unmodified) ,"proteins.\n\n")
            
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

#####################
## Function to perform modification to the fraction of protein sequences selected to be modified by proteinInput function.
## - Creates a set of proteoforms from the selected protein sequences. The number of proteoforms per sequence depends on FracModPerProt.
## - Calls modify function to perform the modification for each proteoform.
## - Creates a fraction of protein sequences to be modified, that will maintain their unmodified counterpart based on RemoveNonModFormFrac.
#####################
performModification <- function(to.Modify, parameters){
  
  cat(" + Performing modification:\n")
  cat("  - Selected modification type(s)", paste0('"', paste0(parameters$PTMTypes, collapse = '", "'), '"'), 
      "with background frequency distribution of", paste0(paste0(parameters$PTMTypesDist*100, "%"), collapse = ", "), "respectively.\n")
  
  for (i in 1:length(parameters$ModifiableResidues)) {
    
    cat("  - For modification", paste0('"', parameters$PTMTypes[i], '"') , "residue(s)", paste0(parameters$ModifiableResidues[[i]], collapse = ", "), 
        "can be modified with background frequency distribution of", paste0(paste0(parameters$ModifiableResiduesDistr[[i]]*100, "%"), collapse = ", "), "respectively.\n")
    
  }
  
  # All proteins in to.Modify set are selected to be modified at least once.
  # Then randomly a set of proteoforms from the imported to.Modify set is selected based on the fraction multiplier FracModPerProt - 1.
  # The size of mod.proteoforms set depends to FracModPerProt. (When 1 a proteoform set of size equal to to.Modify is created, when 2 the size is double etc etc)
  mod.proteoforms <- to.Modify
  
  if(parameters$FracModPerProt >= 2){
    
    mod.proteoforms <- rbind(mod.proteoforms, to.Modify[sample(x = 1:nrow(to.Modify), size = (parameters$FracModPerProt - 1)*nrow(to.Modify), replace = T),])
    mod.proteoforms <- mod.proteoforms[order(mod.proteoforms$Accession),]
    
  }
  
  #Create additional columns for modification positions and modification type.
  mod.proteoforms$PTMPos <- vector(mode = "list", length = nrow(mod.proteoforms))
  mod.proteoforms$PTMType <- vector(mode = "list", length = nrow(mod.proteoforms))
  
  #Modify the selected proteoforms.
  selected.modifications <- modify(mod.proteoforms$Sequence, parameters)
  selected.modifications <- as.data.frame(do.call(rbind, selected.modifications))

  #Fill the columns.
  mod.proteoforms$PTMPos <- selected.modifications$Positions
  mod.proteoforms$PTMType <- selected.modifications$Types
  
  #Summarize modification types and residues modified by each type.
  count.per.AAs <- Reduce(rbind, selected.modifications$Count)
  count.per.AAs <- lapply(1:ncol(count.per.AAs), function(x) colSums(Reduce(rbind, count.per.AAs[,x])))
  
  AAs.percentage <- sapply(count.per.AAs, function(x) sapply(x, function(y) y/sum(x)*100))
  if(length(parameters$PTMTypes) == 1) {AAs.percentage <- list(as.vector(AAs.percentage))}
  Type.percentage <- sapply(count.per.AAs, function(x) sum(x))/sum(unlist(count.per.AAs))*100
  
  cat("  - Sequences modified!\n")
  
  for (i in 1:length(parameters$ModifiableResidues)) {
    
    cat("  - For modification", paste0('"', parameters$PTMTypes[i], '"') , "and residue(s)", paste0(parameters$ModifiableResidues[[i]], collapse = ", "), 
        "the resulted frequency distribution is", paste0(paste0(round(AAs.percentage[[i]], 3), "%"), collapse = ", "), "respectively.\n")
    
  }
  
  cat("  - The resulted frequency distribution for modification type(s)", paste0('"', paste0(parameters$PTMTypes, collapse = '", "'), '"'),
      "is", paste0(paste0(round(Type.percentage, 3), "%"), collapse = ", "), "respectively.\n")
  
  #Select a  fraction of selected to-modify sequences, to remain unmodified too.
  unmodified.proteoforms.indices <- sample(1:nrow(to.Modify), size = (1-parameters$RemoveNonModFormFrac)*nrow(to.Modify))
  unmod.proteoforms <- to.Modify[unmodified.proteoforms.indices, ]
  unmod.proteoforms$PTMPos <- vector(mode = "list", length = nrow(unmod.proteoforms))
  unmod.proteoforms$PTMType <- vector(mode = "list", length = nrow(unmod.proteoforms))
  
  cat("  - A fraction of", nrow(unmod.proteoforms), "modified protein sequences will maintain their unmodified counterpart", paste0("(", (1-parameters$RemoveNonModFormFrac)*100, "%)."), "\n\n")
  
  return(list(mod.proteoforms = mod.proteoforms, unmod.proteoforms = unmod.proteoforms))  
}
#####################

#####################
## Function to modify sequences based on modification type background frequences and modification per residue background frequences.
## - Determines all possible sites for each sequence based on ModifiableResidues list.
## - Calculates the readjusted probability weights for each modification type and candidate residues for each sequence based on the composition of the sequence and the background frequences.
## - Modifies sequences by sampling based on the calculated probability weights. Sampling size per sequence is determined by a trancated poisson distribution [1,total modfiable residues]
## - Returns all reported modification positions per sequences, along with the type of modification and statistics about amound of modification type and modified residues.  
#####################
modify <- function(seq, param){
  
  #Find the positions of candidate modification sites on the sequences.
  possible.modification.sites <- lapply(seq, function(x){ 
    string = strsplit(x, split = "")
    return(lapply(param$ModifiableResidues, function(y) lapply(y, function(z) which(string[[1]] == z))))
  })

  #In case of multiple modification types, amino acid background frequences for each type are proportionally adjusted to background frequences of modification types.
  param$ModifiableResiduesDistr <- lapply(1:length(param$PTMTypesDist ), function(x) param$ModifiableResiduesDistr[[x]] * param$PTMTypesDist[x])
  
  #Find probability weights to scale the distribution of modification per residue for every sequence to fit the background frequences globaly.
  modification.probability.weight <- lapply(possible.modification.sites, function(x){
    lapply(1:length(x), function(y){
      weight <- length(unlist(x))/lengths(x[[y]])*param$ModifiableResiduesDistr[[y]]
      weight[!is.finite(weight)] <- 0
      return(weight)
    })
  })

  #Sample possible modification sites for each sequence based on the adjusted probability weights and size based on a truncated poisson distribution.
  selected.modification.sites <- lapply(1:length(possible.modification.sites), function(x) {
    if(length(unlist(possible.modification.sites[[x]])) != 1){
      
      sort(sample(x = unlist(possible.modification.sites[[x]]),
                  prob = rep(unlist(modification.probability.weight[[x]]), unlist(lapply(possible.modification.sites[[x]], function(y) lengths(y))) ),
                  size = extraDistr::rtpois(n = 1, lambda = param$PTMMultipleLambda * length(unlist(possible.modification.sites[[x]])), a = 1, b = length(unlist(possible.modification.sites[[x]])))
      ))
      
    } else { return(unlist(possible.modification.sites[[x]]))}  
    
  })
  
  #Creating a report per sequence about the modification positions, types and the total number of modified AA per modification type.
  modification.report <- lapply(1:length(selected.modification.sites), function(x){
    
    all.sites <- lapply(possible.modification.sites[[x]], function(y) unlist(y))
    
    modifications.per.AA <- lapply(possible.modification.sites[[x]], function(y) sapply(y, function(z) length(intersect(z, selected.modification.sites[[x]]))))
    
    selected.modification.positions <- lapply(all.sites, function(y) intersect(y, selected.modification.sites[[x]]))
    
    modification.types <- unlist(lapply(1:length(selected.modification.positions), function(y) rep(param$PTMTypes[y], lengths(selected.modification.positions)[y])))
    
    selected.modification.positions <- unlist(selected.modification.positions)
    
    reported.modifications <- list(Positions = selected.modification.positions[order(selected.modification.positions, decreasing = F)], 
                                   Types = modification.types[order(selected.modification.positions, decreasing = F)],
                                   Count = modifications.per.AA)
    
    return(reported.modifications)
    
  })
  
  return(modification.report)
  
}
#####################

#####################
# Function to create a proteoform sample.
# - Uses the above functions and returns a data frame that contains all proteoforms generated (unmodified fraction, modified fraction, modifiable but not modified)
#####################
samplePreparation <- function(parameters){
  
  cat("#SAMPLE PREPARATION - Start\n\n")
  
  protein.Sets <- proteinInput(parameters = parameters)
  
  if(is.null(protein.Sets$to.Modify) & is.null(protein.Sets$to.be.Unmodified)){
    
    cat(crayon::red("#SAMPLE PREPARATION - Finish (Warning: Empty sets)"))
    proteoforms <- NULL
    
  } else {
  
    if(!is.null(protein.Sets$to.Modify)){
      
      proteoforms.after.modification <- performModification(to.Modify = protein.Sets$to.Modify, parameters = parameters)
      
      proteoforms <- rbind(protein.Sets$to.be.Unmodified, 
                           proteoforms.after.modification$mod.proteoform, 
                           proteoforms.after.modification$unmod.proteoforms)
      
      proteoforms$Proteoform_ID <- make.unique(proteoforms$Accession)
      proteoforms$Proteoform_ID  <- rowSums(cbind(as.numeric(sub("\\.", "", sub("^[^.]*", "", proteoforms$Proteoform_ID))), rep(1, nrow(proteoforms))), na.rm = T)
      proteoforms <- proteoforms[,c(1,2,5,3,4)]
      
      rownames(proteoforms) <- 1:nrow(proteoforms)
      
      cat("#SAMPLE PREPARATION - Finish")
      
    } else {
      
      proteoforms <- protein.Sets$to.be.Unmodified
      rownames(proteoforms) <- 1:nrow(proteoforms)
      proteoforms$PTMPos <- vector(mode = "list", length = nrow(proteoforms))
      proteoforms$PTMType <- vector(mode = "list", length = nrow(proteoforms))
      proteoforms$Proteoform_ID <- rep(1, nrow(proteoforms))
      proteoforms <- proteoforms[,c(1,2,5,3,4)]
      
      cat("#SAMPLE PREPARATION - Finish\n\n")
      
    }
    
  }
  
  return(proteoforms)
  
}
#####################

#Below: will be developed further.
#####################
createRegulationPattern = function(NumCond){
  # select 0.5 because division by 2 is already induced this way
  # this ensures that the differentiation amplitude is as speciefied by the user
  regulation_direction <- sample(unlist(lapply(1:NumCond, function(x) c(0.5,-0.5)*x))[1:NumCond], size = NumCond)
  return(regulation_direction)
}
#####################

#####################
addProteoformAbundance <- function(proteoforms, parameters){
  
  # populate the matrix with random noise
  for (name in parameters$QuantColnames) {
    
    proteoforms[name] = rnorm(n = nrow(proteoforms), mean = 0, sd = parameters$QuantNoise)
    
  }
  
  if (is.null(parameters$UserInputFoldChanges)) {
    
    if(parameters$DiffRegFrac != 0){
      
      # select differentially regulated proteoforms
      diff_reg_indices = sample(1:nrow(proteoforms),size = parameters$DiffRegFrac*nrow(proteoforms))
      
      #print(diff_reg_indices)
      
      # determine amplitude of regulation for regulated proteoforms
      proteoforms[diff_reg_indices, "Regulation_Amplitude"] = runif(min = 0, max = parameters$DiffRegMax, n = length(diff_reg_indices))
    
      regulationPatterns <- lapply(1:length(diff_reg_indices), function(x) createRegulationPattern(parameters$NumCond))
    
    } else {
      
      proteoforms$Regulation_Amplitude <- vector(mode = "list", length = nrow(proteoforms))
      regulationPatterns <- NULL
    
    }  
    
      
  } else {
    # select differentially regulated proteoforms
    diff_reg_indices = sample(1:nrow(proteoforms),size = sum(parameters$UserInputFoldChanges$NumRegProteoforms))
    
    # determine amplitude of regulation for regulated proteoforms
    proteoforms[diff_reg_indices, "Regulation_Amplitude"] = parameters$UserInputFoldChanges$RegulationFC
    
    
    regulationPatterns <- lapply(1:length(diff_reg_indices), function(x) createRegulationPattern(parameters$NumCond))
  }
  
  proteoforms$Regulation_Pattern <- vector(mode = "list", length = nrow(proteoforms))
  
  #print(length(regulationPatterns))
  #print(regulationPatterns)
  
  if(!is.null(regulationPatterns)){
  
    proteoforms$Regulation_Pattern[diff_reg_indices] = regulationPatterns
    #[diff_reg_indices, "Regulation_Pattern"]
    proteoforms[diff_reg_indices, parameters$QuantColnames] = 
      # add regulation pattern*regulation amplitude to random noise
      proteoforms[diff_reg_indices, parameters$QuantColnames] +
      
      #generate regulation patterns for all regulated proteoforms
      t(sapply(1:length(diff_reg_indices), function(x) {
        rep(regulationPatterns[[x]], each = parameters$NumReps)
        
        # multiply regulation pattern with Regulation amplitude
      })) * (proteoforms[diff_reg_indices, "Regulation_Amplitude"])
    
  } else {
    
    proteoforms$Regulation_Pattern <- vector(mode = "list", length = nrow(proteoforms))
    
  }
  
  if (is.null(parameters$AbsoluteQuanMean)) {
    # Remove Values below the threshold set in the Parameters file
    proteoforms[,parameters$QuantColnames][proteoforms[,parameters$QuantColnames] < parameters$ThreshNAProteoform]  = NA
  } else {
    cat("Add quan. distribution: Relative -> absolute\n")
    vec <- rnorm(n = nrow(proteoforms), mean = parameters$AbsoluteQuanMean, sd = parameters$AbsoluteQuanSD)
    for (name in parameters$QuantColnames) {
      proteoforms[name] = proteoforms[name] + vec
    }
    if (parameters$ThreshNAQuantileProt > 0) {
      # Remove Values below the threshold set in the Parameters file
      thresh <- quantile(x = vec, probs = parameters$ThreshNAQuantileProt)
      proteoforms[,parameters$QuantColnames][proteoforms[,parameters$QuantColnames] < thresh]  = NA
    }
  }
  
  return(proteoforms)
}
#####################