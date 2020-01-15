

#setwd("/Users/Eva/Desktop/PhD/EuBIC")
 
# Set filepath:
#filepath <- "SharmaTxtFiles/pYSharma/txt/evidence.txt" # pTyr
#filepath <- "SharmaTxtFiles/proteomeSharma/txt/evidence.txt" # unmodified
#filepath <- "SharmaTxtFiles/TiO2Sharma/txt/evidence.txt" # pSerThr

#filepath <- "Phospho_data_TMT_PXD007871/HBT_PARK2_nonmod_PSMs.xlsx"
#filepath <- "Phospho_data_TMT_PXD007871/HBT\ PARK2\ phospho_PSMs.xlsx"

protSelectionPath <- "phosphoSTY_subset.txt"
protSelection <- read.csv(protSelectionPath, allowEscapes = TRUE, check.names = FALSE,sep = "\t")
protSelectionVec <- unique(sub(".*\\|","", sub('\\|([^\\|]*)$', '', protSelection$Leading.proteins)))


ex <- strsplit(basename(filepath), split="\\.")[[1]][2]

filepaths <- c("SharmaTxtFiles/pYSharma/txt/evidence.txt", # pTyr
                "SharmaTxtFiles/proteomeSharma/txt/evidence.txt", # unmodified
                "SharmaTxtFiles/TiO2Sharma/txt/evidence.txt") # pSerThr

library(data.table)

getCleanTable <- function(filepath){
    if (ex == "txt"){
        # Read file
        mtx <- fread(filepath, select = c("Sequence", "Leading razor protein", "Modified sequence", "Experiment", "Intensity", "Modifications"))

        #dat <- read.csv(filepath, allowEscapes = TRUE, check.names = FALSE,sep = "\t")
        
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
        
        #mtx <- as.matrix(dat[, grepl("^LFQ", names(dat))])
        #mtx <- dat[, grepl("Sequence|Leading\ razor\ protein$|Modified\ sequence$|Experiment|Intensity", names(dat))]
        mtx$`Leading razor protein` <- sub(".*\\|","", sub('\\|([^\\|]*)$', '', mtx$`Leading razor protein`))
        
        mtx <- subset(mtx, `Leading razor protein` %in% protSelectionVec)
        
        #mtx$`Modified sequence`[duplicated(mtx)]
        #temp <- dat[dat$`Modified sequence` == "_EEDEEPES(ph)PPEK_",]
        
        # Summarize rows with same modification and same experiment id by taking average of intensities
        mtx.aggr <- aggregate(Intensity ~ Sequence + `Leading razor protein` +`Modified sequence`+ Experiment  + prep, data=mtx, sum, na.rm=TRUE)
    
        # Transform to table with columns pepseq, PTMs, PTM type, accs, quant1, quant2, ...
        
    } else if (ex == "xlsx"){
        library("readxl")
        dat <- read_excel(filepath)
    }
}

formattedDF <- do.call(rbind,lapply(filepaths, getCleanTable))
#write.csv(formattedDF, file=)
#save(formattedDF, file="formattedDF.RData")

# order prep increasing (phTyr before pThrSer)
formattedDF <- formattedDF[order(formattedDF$prep),]
# keep only non-duplicated
formattedDF <- formattedDF[!duplicated(formattedDF[,-which(names(formattedDF) =="prep")]),]

#experimentMapping <- formattedDF[,c("Raw file", "Experiment")]
#experimentMapping <- experimentMapping[!duplicated(experimentMapping),]


exPlanPath <- "ExperimentalPlanSharma.txt"
exPlan <- read.csv(exPlanPath, allowEscapes = TRUE, check.names = FALSE,sep = "\t")
mapping <- data.frame(key=exPlan$label,value=exPlan$condition)
replicate <- gsub(".*(\\d+)", "\\1", formattedDF$Experiment)
mapping <- setNames(exPlan$condition, exPlan$label)[formattedDF$Experiment]

formattedDF$expGroup <- paste(mapping, replicate, sep="_")
formattedDF <- formattedDF[formattedDF$prep == "TiO2",]

library(reshape2)
formattedDFcast <- dcast(formattedDF, Sequence + `Leading razor protein` + `Modified sequence` ~ expGroup, value.var = c("Intensity"))

