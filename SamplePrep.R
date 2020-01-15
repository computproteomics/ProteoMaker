
#Import protein sequences and split proteins to two sets (to get modified and to remain unmodified)
proteinInput <- function(fasta.path, parameters){
  
  
  fasta <- protr::readFASTA(file = fasta.path, legacy.mode = TRUE, seqonly = FALSE)
  fasta <- data.frame(Sequence = unlist(fasta), Accession = sub(".*[|]([^.]+)[|].*", "\\1", names(fasta)), stringsAsFactors = F)
  rownames(fasta) <- 1:nrow(fasta)
  
  to.modify.indices <- sample(1:nrow(fasta), size = parameters$FracModProt * nrow(fasta))
  
  to.Modify <- fasta[to.modify.indices,]
  rownames(to.Modify) <- 1:nrow(to.Modify)
  to.be.Unmodified <- fasta[-to.modify.indices,]
  rownames(to.be.Unmodified) <- 1:nrow(to.be.Unmodified)
  
  to.be.Unmodified$PTMPos = vector(mode = "list", length = nrow(to.be.Unmodified))
  to.be.Unmodified$PTMType = vector(mode = "list", length = nrow(to.be.Unmodified))
  
  
  return(list(to.Modify = to.Modify, to.be.Unmodified = to.be.Unmodified))
}

#Performs phosphorylation to selected protein sequences.
#Determines number of modified proteoforms for each protein.
perform.phosphorylation <- function(to.Modify, parameters){
  
  
  mod.proteoforms <- to.Modify[sample(x = 1:nrow(to.Modify), size = parameters$FracModPerProt*nrow(to.Modify), replace = T),]
  mod.proteoforms <- mod.proteoforms[order(mod.proteoforms$Accession),]
  
  mod.proteoforms$PTMPos <- vector(mode = "list", length = nrow(mod.proteoforms))
  mod.proteoforms$PTMType <- vector(mode = "list", length = nrow(mod.proteoforms))
  
  selected.phospho <- phosphorylate(mod.proteoforms$Sequence, parameters)
  mod.proteoforms$PTMPos <- selected.phospho$site
  mod.proteoforms$PTMType <- lapply(selected.phospho$site, function(x) sapply(1:length(x), function(y) parameters$PTMTypes))
  
  AAcounts <- sapply(1:length(parameters$ModifiableResidues$mod), function(x) 100*length(which(unlist(selected.phospho$count) == x))/length(unlist(selected.phospho$count)) )
  
  unmodified.proteoforms.indices = sample(1:nrow(to.Modify), size = (1-parameters$RemoveNonModFormFrac)*nrow(to.Modify))
  unmod.proteoforms <- to.Modify[unmodified.proteoforms.indices,]
  
  
  unmod.proteoforms$PTMPos <- vector(mode = "list", length = nrow(unmod.proteoforms))
  unmod.proteoforms$PTMType <- vector(mode = "list", length = nrow(unmod.proteoforms))
  
  
  return(list(mod.proteoform = mod.proteoforms, counts = AAcounts, unmod.proteoforms = unmod.proteoforms))  
  
}

#Determines the total number of phosphorylations per proteoform.
#Determines the phosphosites per AA type, so the distribution of appearences in the resulted data follows the natural observed frequences.
phosphorylate <- function(seq, param){
  
  possible.phospho.sites <- lapply(seq, function(x){ string = strsplit(x, split = "")
  return(lapply(param$ModifiableResidues$mod, function(y) which(string[[1]] == y)))})
  
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

#Uses the above functions and returns all proteins contained in the sample (unmodified, modified, modifiable but not modified)
samplePreparation <- function(fasta.path, parameters){
  
  protein.Sets <- proteinInput(fasta.path = fasta.path, parameters = parameters)
  
  proteoforms.after.phosphorylation <- perform.phosphorylation(to.Modify = protein.Sets$to.Modify, parameters = parameters)
  
  proteoforms <- rbind(protein.Sets$to.be.Unmodified, proteoforms.after.phosphorylation$mod.proteoform, proteoforms.after.phosphorylation$unmod.proteoforms)
  rownames(proteoforms) <- 1:nrow(proteoforms)
  
  cat("Modified AA frequencies for: ", parameters$ModifiableResidues$mod, "\n")
  cat(proteoforms.after.phosphorylation$counts)
  
  return(proteoforms)
  
}
