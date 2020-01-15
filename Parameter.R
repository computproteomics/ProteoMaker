# Variables
ExpDesign <- rep(1:2,3)
Param <- list()
GroundTruth <- data.table()
Digested <- data.table()

### Parameter setting
## Ground truth
# Fraction of proteins getting PTMs
Param$FracModProt <- 1
# Fraction of modifiable proteins that get a PTM
Param$FracModPerProt <- 2
# PTM types
Param$PTMTypes <- c("ph")
# residues for PTM type
Param$ModifiableResidues <- list()
for (mod in Param$PTMTypes) 
  Param$ModifiableResidues$mod <- c("S","T","Y")

Param$NumProteoforms <- 0

#

## Sample preparation


## MS run




