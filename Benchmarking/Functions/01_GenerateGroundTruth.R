################################################################################
#                  GENERATE PROTEOFORM-CENTRIC GROUND TRUTH                    #
################################################################################

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

############FUTURE############
# - Make this work for multiple modification types.
# - Import modification DB data.
##############################

proteinInput <- function(parameters){
  
  cat("+ Importing data:\n")
  
  #Check is the fasta file can be loaded.
  error <- try(protr::readFASTA(file = parameters$PathToFasta, legacy.mode = TRUE, seqonly = FALSE), silent = TRUE)
  
  if(class(error) != "try-error"){
    
    #Read fasta.
    fasta <- protr::readFASTA(file = parameters$PathToFasta, legacy.mode = TRUE, seqonly = FALSE)
    fasta <- data.frame(Sequence = unlist(fasta), Accession = sub(".*[|]([^.]+)[|].*", "\\1", names(fasta)), stringsAsFactors = F)
    rownames(fasta) <- 1:nrow(fasta)
    cat(" - File", parameters$PathToFasta, "imported, containing", nrow(fasta), "protein sequences.\n")
    
    #Filter proteins carrying unusual amino acids.
    knownAA <- c("A",	"L", "R",	"K", "N",	"M", "D", "F", "C",	"P", "E",	"S", "Q", "T", "G", "W", "H", "Y", "I", "V")
    unknownAA <- setdiff(LETTERS, knownAA)
    initialRows <- nrow(fasta)
    fasta <- fasta[if(initialRows > 1){ rowSums(sapply(unknownAA, grepl, x = fasta$Sequence)) == 0}
                   else { sum(sapply(unknownAA, grepl, x = fasta$Sequence)) == 0}, ]
    cat(" - A total of", initialRows - nrow(fasta) ,"protein sequences have been removed due to unusual amino acids", paste0("(", paste0(unknownAA, collapse = ","), ")"), "composition.\n")
    
    #Remove duplicated protein accessions.
    cat(" - A total of", length(which(duplicated(fasta$Accession))) ,"duplicated protein accessions have been removed.\n")
    fasta <- fasta[!duplicated(fasta$Accession),]
    cat(" - Total number of remaining protein sequences:", nrow(fasta), "\n\n")
    
    #Proceed unless fasta df is empty after filtering.
    if(nrow(fasta) > 0){
      
      cat("+ Creating modified and unmodified fractions:\n")
      
      #Returns a list of the number of modifiable residues on all sequences, for each residue of ModifiableResidues$mod.
      possible.modifiable.AAs <- as.data.frame(t(sapply(fasta$Sequence, function(x){ 
        string = strsplit(x, split = "")
        return(sapply(parameters$ModifiableResidues$mod, function(y) length(which(string[[1]] == y))))
      })))
      
      rownames(possible.modifiable.AAs) <- 1:nrow(possible.modifiable.AAs)
      
      #Add an additional column, denoting the number of ModifiableResidues$mod residue types found for each protein.
      possible.modifiable.AAs$Modifiable <- sapply(1:nrow(possible.modifiable.AAs),function(x) length(which(possible.modifiable.AAs[x,] > 0))) 
      
      #Proteins that do not contain any of the ModifiableResidues$mod.
      unmodifiable.indices <- which(possible.modifiable.AAs$Modifiable < length(parameters$ModifiableResidues$mod))
      
      #Proteins that contain all residues in ModifiableResidues$mod, and thus are modifiable. (Should fix this)
      modifiable.indices <- setdiff(1:nrow(fasta), unmodifiable.indices)
      
      cat(" - A total of", length(unmodifiable.indices) ,"sequences are unmodifiable.\n")
      
      #If there is not a specific protein list to be modified.
      if(is.null(parameters$PathToProteinList)){
        
        #If the set of proteins can be fractionated. This covers the case when parameters$FracModProt = 0
        if((parameters$FracModProt * length(modifiable.indices)) >= 1){
          
          # Randomly select a FracModProt % of modifiable proteins (modifiable.indices).
          to.modify.indices <- sample(modifiable.indices, size = parameters$FracModProt * length(modifiable.indices))
          to.Modify <- fasta[to.modify.indices, ]
          rownames(to.Modify) <- 1:nrow(to.Modify)
          
          cat(" - A total of", length(modifiable.indices) ,"sequences are modifiable, from which", parameters$FracModProt*100,"% randomly selected to be modified.\n")
          
          to.be.Unmodified <- fasta[-to.modify.indices,]
          
          #Add additional columns to the unmodified fraction df.  
          to.be.Unmodified$PTMPos = vector(mode = "list", length = nrow(to.be.Unmodified))
          to.be.Unmodified$PTMType = vector(mode = "list", length = nrow(to.be.Unmodified))
          
          cat(" - Modified fraction:", nrow(to.Modify) ,"proteins.\n")
          cat(" - Unmodified fraction:", nrow(to.be.Unmodified) ,"proteins.\n\n")
          
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
          
          cat(" - The protein set cannot be fractionated, too low number of proteins or too low FracModProt %.\n")
          cat(" - Modified fraction: 0 proteins.\n")
          cat(" - Unmodified fraction:", nrow(to.be.Unmodified) ,"proteins.\n\n")
          
          return(list(to.Modify = NULL, to.be.Unmodified = to.be.Unmodified))
          
        }
        
      } else { #If there is a specific protein list to be modified.
        
        #Check is the protein list file can be loaded.
        error <- try(read.csv(parameters$PathToProteinList, header = F, stringsAsFactors = F), silent = TRUE)
        
        if(class(error) != "try-error"){
          
          protein.list.input <- unique(as.vector(read.csv(parameters$PathToProteinList, header = F, stringsAsFactors = F)[,1]))
          cat(" - Protein list file", parameters$PathToProteinList, "loaded and contains", length(protein.list.input),"unique protein accessions.\n")
          
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
            
            cat(" - Modified fraction:", nrow(to.Modify) ,"proteins.\n")
            cat(" - Unmodified fraction:", nrow(to.be.Unmodified) ,"proteins.\n\n")
            
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
            cat(" - Proteins in", parameters$PathToProteinList, "are not modifiable or are not found in", parameters$PathToFasta, ".\n")
            cat(" - Modified fraction: 0 proteins.\n")
            cat(" - Unmodified fraction:", nrow(to.be.Unmodified) ,"proteins.\n\n")
            
            return(list(to.Modify = NULL, to.be.Unmodified = NULL))
            
          }
          
        } else {
          
          cat(crayon::red(" - Protein list file", parameters$PathToProteinList, "couldn't be loaded!\n\n"))
          return(list(to.Modify = NULL, to.be.Unmodified = NULL))
          
        }
        
      }
      
    } else {
      
      cat(crayon::red(" - There are no protein sequences left!\n\n"))
      return(list(to.Modify = NULL, to.be.Unmodified = NULL))
      
    }
    
  } else {
    
    cat(crayon::red(" - Fasta file", parameters$PathToFasta, "couldn't be loaded!\n\n"))
    return(list(to.Modify = NULL, to.be.Unmodified = NULL))
    
  }
  
}
#####################

#####################
## Performs phosphorylation to selected protein sequences.
## Determines number of modified proteoforms for each protein.
#####################
performModification <- function(to.Modify, parameters){
  
  # Randomly selects a set of proteoforms from the imported to.Modify set. The size of mod.proteoforms set
  # depends to FracModPerProt. (When 1 a proteoform set of size equal to to.Modify is created, when 2 the size is double etc etc)
  mod.proteoforms <- to.Modify[sample(x = 1:nrow(to.Modify), size = parameters$FracModPerProt*nrow(to.Modify), replace = T),]
  mod.proteoforms <- mod.proteoforms[order(mod.proteoforms$Accession),]

  mod.proteoforms$PTMPos <- vector(mode = "list", length = nrow(mod.proteoforms))
  mod.proteoforms$PTMType <- vector(mode = "list", length = nrow(mod.proteoforms))
  
  selected.modif <- modify(mod.proteoforms$Sequence, parameters)
  mod.proteoforms$PTMPos <- selected.modif$site
  mod.proteoforms$PTMType <- lapply(selected.modif$site, function(x) sapply(1:length(x), function(y) parameters$PTMTypes))
  
  AAcounts <- sapply(1:length(parameters$ModifiableResidues$mod), function(x) 100*length(which(unlist(selected.modif$count) == x))/length(unlist(selected.modif$count)) )
  # cat("Percentages of the in silico modified", paste(parameters$ModifiableResidues$mod, collapse = ", "), "are:", paste(AAcounts, collapse = ", "), "\n")
  
  unmodified.proteoforms.indices = sample(1:nrow(to.Modify), size = (1-parameters$RemoveNonModFormFrac)*nrow(to.Modify))
  unmod.proteoforms <- to.Modify[unmodified.proteoforms.indices,]
  
  
  unmod.proteoforms$PTMPos <- vector(mode = "list", length = nrow(unmod.proteoforms))
  unmod.proteoforms$PTMType <- vector(mode = "list", length = nrow(unmod.proteoforms))
  
  return(list(mod.proteoform = mod.proteoforms, counts = AAcounts, unmod.proteoforms = unmod.proteoforms))  
}
#####################

#####################
## Determines the total number of phosphorylations per proteoform based on the parameter ModifiableResidues.
## Determines the phosphosites per AA type, so the distribution of appearences in the resulted data follows the natural observed frequences.
#####################
modify <- function(seq, param){
  
  possible.phospho.sites <- lapply(seq, function(x){ 
    string = strsplit(x, split = "")
    return(lapply(param$ModifiableResidues$mod, function(y) which(string[[1]] == y)))
  }) # -> position of modifiable sites in the sequence.
  
  adjusted.phospho.probability <- lapply(possible.phospho.sites, function(x) length(unlist(x))/lengths(x)*param$ModifiableResiduesDistr$mod)
  
  selected.phospho.sites <- lapply(1:length(seq), function(x) {
    sort(sample(
      x = unlist(possible.phospho.sites[[x]]),
      prob = unlist(sapply(1:length(param$ModifiableResidues$mod), function(y) rep(adjusted.phospho.probability[[x]][y],lengths(possible.phospho.sites[[x]][y])))),
      size = min(c(
        rpois(1, lambda = param$PTMMultipleLambda * length(unlist(
          possible.phospho.sites[[x]]
        ))) + 1,
        length(unlist(possible.phospho.sites[[x]]))
      ))
    ))
  })
  
  
  selected.phospho.count <- lapply(1:length(seq), function(selected){
    
    sapply(unlist(selected.phospho.sites[[selected]]), function(y)
      
      which(lengths(sapply(possible.phospho.sites[[selected]], function(x) which(x == y))) == 1)
      
    )
    
  })
  
  return(list(site = selected.phospho.sites, count = selected.phospho.count))
  
}
#####################

#####################
# Uses the above functions and returns all proteins contained in the sample (unmodified, modified, modifiable but not modified)
#####################
samplePreparation <- function(parameters){
  
  protein.Sets <- proteinInput(parameters = parameters)
  
  if (Param$FracModProt > 0) {
    proteoforms.after.phosphorylation <- performModification(to.Modify = protein.Sets$to.Modify, parameters = parameters)
    
    proteoforms <- rbind(protein.Sets$to.be.Unmodified, 
                         proteoforms.after.phosphorylation$mod.proteoform, 
                         proteoforms.after.phosphorylation$unmod.proteoforms)
    rownames(proteoforms) <- 1:nrow(proteoforms)

    cat("Modified AA frequencies for: ", parameters$ModifiableResidues$mod, "\n")
    cat(proteoforms.after.phosphorylation$counts, "\n")

  } else {
    proteoforms <- protein.Sets$to.be.Unmodified
    rownames(proteoforms) <- 1:nrow(proteoforms)
    proteoforms$PTMPos <- vector(mode = "list", length = nrow(proteoforms))
    proteoforms$PTMType <- vector(mode = "list", length = nrow(proteoforms))
    
    cat("No modified proteins in these data\n")
  }
  
  return(proteoforms)
  
}
#####################

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
  for (name in parameters$quant_colnames) {
    
    proteoforms[name] = rnorm(n = nrow(proteoforms), mean = 0, sd = parameters$QuantNoise)
    
  }
  
  if (is.null(parameters$UserInputFoldChanges)) {
    # select differentially regulated proteoforms
    diff_reg_indices = sample(1:nrow(proteoforms),size = parameters$DiffRegFrac*nrow(proteoforms))
    
    # determine amplitude of regulation for regulated proteoforms
    proteoforms[diff_reg_indices, "Regulation_Amplitude"] = runif(min =1, max = parameters$DiffRegMax, n = length(diff_reg_indices))
    
    regulationPatterns <- lapply(1:length(diff_reg_indices), function(x) createRegulationPattern(parameters$NumCond))
  } else {
    # select differentially regulated proteoforms
    diff_reg_indices = sample(1:nrow(proteoforms),size = sum(parameters$UserInputFoldChanges$NumRegProteoforms))
    
    # determine amplitude of regulation for regulated proteoforms
    proteoforms[diff_reg_indices, "Regulation_Amplitude"] = parameters$UserInputFoldChanges$RegulationFC
    
    
    regulationPatterns <- lapply(1:length(diff_reg_indices), function(x) createRegulationPattern(parameters$NumCond))
  }
  
  proteoforms$Regulation_Pattern <- vector(mode = "list", length = nrow(proteoforms))
  proteoforms$Regulation_Pattern[diff_reg_indices] = regulationPatterns
  #[diff_reg_indices, "Regulation_Pattern"]
  proteoforms[diff_reg_indices, parameters$quant_colnames] = 
    # add regulation pattern*regulation amplitude to random noise
    proteoforms[diff_reg_indices, parameters$quant_colnames] +
    
    #generate regulation patterns for all regulated proteoforms
    t(sapply(1:length(diff_reg_indices), function(x) {
      rep(regulationPatterns[[x]], each = parameters$NumReps)
      
      # multiply regulation pattern with Regulation amplitude
    })) * proteoforms[diff_reg_indices, "Regulation_Amplitude"]
  
  if (is.null(parameters$AbsoluteQuanMean)) {
    # Remove Values below the threshold set in the Parameters file
    proteoforms[,parameters$quant_colnames][proteoforms[,parameters$quant_colnames] < parameters$ThreshNAProteoform]  = NA
  } else {
    cat("Add quan. distribution: Relative -> absolute\n")
    vec <- rnorm(n = nrow(proteoforms), mean = parameters$AbsoluteQuanMean, sd = parameters$AbsoluteQuanSD)
    for (name in parameters$quant_colnames) {
      proteoforms[name] = proteoforms[name] + vec
    }
    if (parameters$ThreshNAQuantileProt > 0) {
      # Remove Values below the threshold set in the Parameters file
      thresh <- quantile(x = vec, probs = parameters$ThreshNAQuantileProt)
      proteoforms[,parameters$quant_colnames][proteoforms[,parameters$quant_colnames] < thresh]  = NA
    }
  }
  
  return(proteoforms)
}
#####################