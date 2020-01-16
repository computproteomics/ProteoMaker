
#setwd("/Users/Eva/Desktop/PhD/EuBIC")
 
# Set filepath:
#filepath <- "SharmaTxtFiles/pYSharma/txt/evidence.txt" # pTyr
#filepath <- "SharmaTxtFiles/proteomeSharma/txt/evidence.txt" # unmodified
#filepath <- "SharmaTxtFiles/TiO2Sharma/txt/evidence.txt" # pSerThr

protSelectionPath <- "phosphoSTY_subset.txt"
protSelection <- read.csv(protSelectionPath, allowEscapes = TRUE, check.names = FALSE,sep = "\t")
protSelectionVec <- unique(sub(".*\\|","", sub('\\|([^\\|]*)$', '', protSelection$Leading.proteins)))

ex <- strsplit(basename(filepath), split="\\.")[[1]][2]
filepaths <- c("SharmaTxtFiles/pYSharma/txt/", # pTyr
                "SharmaTxtFiles/proteomeSharma/txt/", # unmodified
                "SharmaTxtFiles/TiO2Sharma/txt/") # pSerThr

library(data.table)

getCleanTable <- function(filepath){

        filepath_evidence <- paste0(filepath, "evidence.txt")
        filepath_peptides <- paste0(filepath, "peptides.txt")

        # Read file
        mtx <- fread(filepath_evidence, select = c("Sequence", "Leading razor protein", "Modified sequence", "Experiment", "Intensity", "Modifications"))
        pepStart <- fread(filepath_peptides, select = c("Sequence", "Start position"))
        mtx <- merge(mtx, pepStart, by="Sequence")
        
        if (grepl("proteomeSharma", filepath)){
            # remove unmodified peptides
            mtx <- mtx[!grepl("Phospho", mtx$Modifications),]
            mtx$prep <- "proteome"
        } else if (grepl("TiO2Sharma", filepath)){
            mtx <- mtx[grepl("Phospho", mtx$Modifications),]
            mtx$prep <- "TiO2"
        } else if (grepl("pYSharma", filepath)){
            mtx <- mtx[grepl("Y(ph)", mtx$`Modified sequence`, fixed = T),]
            mtx$prep <- "pY"
        }
        
        mtx$`Leading razor protein` <- sub(".*\\|","", sub('\\|([^\\|]*)$', '', mtx$`Leading razor protein`))
        mtx <- subset(mtx, `Leading razor protein` %in% protSelectionVec)
        
        # Summarize rows with same modification and same experiment id by taking average of intensities
        mtx.aggr <- aggregate(Intensity ~ Sequence + `Start position` + `Leading razor protein` +`Modified sequence`+ Experiment  + prep, data=mtx, sum, na.rm=TRUE)
}

# combine the data frames corresonding to the three datasets
formattedDF <- do.call(rbind,lapply(filepaths, getCleanTable))

# order prep increasing (phTyr before pThrSer)
formattedDF <- formattedDF[order(formattedDF$prep),]
# keep only non-duplicated
formattedDF <- formattedDF[!duplicated(formattedDF[,-which(names(formattedDF) =="prep")]),]

exPlanPath <- "ExperimentalPlanSharma.txt"
exPlan <- read.csv(exPlanPath, allowEscapes = TRUE, check.names = FALSE,sep = "\t")
mapping <- data.frame(key=exPlan$label,value=exPlan$condition)
replicate <- gsub(".*(\\d+)", "\\1", formattedDF$Experiment)
mapping <- setNames(exPlan$condition, exPlan$label)[formattedDF$Experiment]
formattedDF$expGroup <- paste(mapping, replicate, sep="_")

# keep only 4 replicates
formattedDF <- formattedDF[!(as.numeric(str_sub(formattedDF$expGroup, start= -1)) > 4),]

formattedDF <- formattedDF[grepl("control|EGF15|mitosis",formattedDF$expGroup),]

# control -> C_1_R_; EGF15 -> C_2_R_; mitosis -> C_3_R_
formattedDF$expGroup <- gsub("control", "C_1_R", formattedDF$expGroup)
formattedDF$expGroup <- gsub("EGF15", "C_2_R", formattedDF$expGroup)
formattedDF$expGroup <- gsub("mitosis", "C_3_R", formattedDF$expGroup)



library(reshape2)
formattedDF$Intensity <- as.numeric(formattedDF$Intensity)
formattedDFcast <- dcast(formattedDF, Sequence + `Start position` + `Leading razor protein` + `Modified sequence` ~ expGroup, value.var = c("Intensity"))

library(stringr)
parse_modseq <- function(modified_sequences) {
    PTM_pos <- vector(mode = "list")
    PTM_type <- vector(mode = "list")
    for(i in seq_along(modified_sequences)) {
        seq <- modified_sequences[i]
        if (grepl(")", seq, fixed = T)) {
            modif <- strsplit(seq, ")", fixed = T)[[1]]
            modif <- modif[grepl("(", modif, fixed = T)]
            PTM_type[[i]] <- gsub("^.+\\(", "", modif, perl = T)
            PTM_pos[[i]] <- as.numeric(str_locate(gsub("_", "", modif, fixed = T), "\\(..$")[,1] - 1)
        } else {
            PTM_pos[[i]] <- NA
            PTM_type[[i]] <- NA
        }
    }
    return(list("PTM_type" = PTM_type, "PTM_pos" = PTM_pos))
}

PTMinfo <- parse_modseq(formattedDFcast$`Modified sequence`)
formattedDFcast$PTMtype <- PTMinfo$PTM_type

formattedDFcast$PTMpos <- lapply(seq_along(PTMinfo$PTM_pos), function(x) {
    PTMinfo$PTM_pos[[x]] + formattedDFcast$`Start position`[x] - 1
})

formattedDFcast$`Modified sequence` <- gsub("_","", formattedDFcast$`Modified sequence`)

library(dplyr)
finalDF <- formattedDFcast %>%
    select("Sequence", "Leading razor protein", "PTMpos", "PTMtype", everything())
finalDF$`Modified sequence` <- NULL

names(finalDF)[names(finalDF) == 'Sequence'] <- 'PepSequence'
names(finalDF)[names(finalDF) == 'Leading razor protein'] <- 'Accession'
names(finalDF)[names(finalDF) == 'PTMpos'] <- 'PTMPos'
names(finalDF)[names(finalDF) == 'PTMtype'] <- 'PTMType'
names(finalDF)[names(finalDF) == 'Start position'] <- 'PepStart'

# save(finalDF, file="expDataFrame.RData")
