source("Parameter.R")

library(OrgMassSpecR)
fragments <- Digest(proteoforms$Sequence[36],enzyme = "trypsin", missed = 0)
PTMs <- data.frame(key=proteoforms$PTMPos[36], value=proteoforms$PTMType[36])
names(PTMs) <- c("position", "modification")

# if ptm in range of start to stop for fragment add ptm
fragments$start
fragments$stop

# row wise comparison with PTMs table
fragments[(PTMs$position[1]>=fragments$start) & (PTMs$position[1]<=fragments$stop), ]$PTMType <- toString(PTMs$modification[1])

#Param$PepMinLength
#Param$PepMaxLength