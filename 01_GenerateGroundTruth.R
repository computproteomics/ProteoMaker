################################################################################
#                  GENERATE PROTEOFORM-CENTRIC GROUND TRUTH                    #
################################################################################


#####################
## Import protein sequences and split proteins to two sets (one set will be modified and the other remain unmodified)
#####################
proteinInput <- function(fasta.path, parameters){
  
  fasta <- protr::readFASTA(file = fasta.path, legacy.mode = TRUE, seqonly = FALSE)
  fasta <- data.frame(Sequence = unlist(fasta), Accession = sub(".*[|]([^.]+)[|].*", "\\1", names(fasta)), stringsAsFactors = F)
  rownames(fasta) <- 1:nrow(fasta)
  
  knownAA <- c("A",	"L", "R",	"K", "N",	"M", "D", "F", "C",	"P", "E",	"S", "Q", "T", "G", "W", "H", "Y", "I", "V")
  unknownAA <- setdiff(LETTERS, knownAA)
  cat("Remove from fasta the protein sequences with non-usual amino-acids:", unknownAA, "\n")
  
  fasta <- fasta[rowSums(sapply(unknownAA, grepl, x = fasta$Sequence)) == 0,]
  
  if (parameters$FracModProt > 0) {
    # >> Should be properly done in future. Not considering proteins that cannot be modified
    # at all parameters$ModifiableResidues$mod AAs, for modification process.
    possible.modif <- as.data.frame(t(sapply(fasta$Sequence, function(x){ 
      string = strsplit(x, split = "")
      return(sapply(parameters$ModifiableResidues$mod, function(y) length(which(string[[1]] == y))))
    }))) # -> returns a list with the numper of all possible modified residues per modifiable aminoacid.
    possible.modif$zeros <- sapply(1:nrow(possible.modif),function(x) length(which(possible.modif[x,] == 0))) # -> lists when there is no modification
    zero.indices <- which(possible.modif$zeros > 0) # -> proteins with no modifiable aminoacid.
    # Randomly select a set of proteins to modify (only in the proteins that can be modified) according to the parameter FracModProt (fraction of proteins to modify)
    # >> CAN BE CHANGED TO A LIST OF PROTEINS OF INTEREST IN THE FUTURE << #
    to.modify.indices <- sample(setdiff(1:nrow(fasta), zero.indices), size = parameters$FracModProt * nrow(fasta))
    to.Modify <- fasta[to.modify.indices,]
    rownames(to.Modify) <- 1:nrow(to.Modify)
    to.be.Unmodified <- fasta[-to.modify.indices,]
    rownames(to.be.Unmodified) <- 1:nrow(to.be.Unmodified)
    
    to.be.Unmodified$PTMPos = vector(mode = "list", length = nrow(to.be.Unmodified))
    to.be.Unmodified$PTMType = vector(mode = "list", length = nrow(to.be.Unmodified))
    
    return(list(to.Modify = to.Modify, to.be.Unmodified = to.be.Unmodified))
  } else {
    to.be.Unmodified <- fasta
    rownames(to.be.Unmodified) <- 1:nrow(to.be.Unmodified)
    return(list(to.Modify = NULL, to.be.Unmodified = to.be.Unmodified))
  }
}
#####################

#####################
## Performs phosphorylation to selected protein sequences.
## Determines number of modified proteoforms for each protein.
#####################
perform.modification <- function(to.Modify, parameters){
  
  # Randomly select a set of residue to modify according to the parameter FracModPerProt (number of modifications)
  # >> CAN BE CHANGED TO A LIST OF PROTEINS OF INTEREST IN THE FUTURE << #
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
samplePreparation <- function(fasta.path, parameters){
  
  protein.Sets <- proteinInput(fasta.path = fasta.path, parameters = parameters)
  
  if (Param$FracModProt > 0) {
    proteoforms.after.phosphorylation <- perform.modification(to.Modify = protein.Sets$to.Modify, parameters = parameters)
    
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



