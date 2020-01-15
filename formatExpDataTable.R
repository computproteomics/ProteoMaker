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


getCleanTable <- function(filepath){
    if (ex == "txt"){
        # Read file
        dat <- read.csv(filepath, allowEscapes = TRUE, check.names = FALSE,sep = "\t")
        
        if (!grepl("proteomeSharma", filepath)){
            # remove unmodified peptides
            dat<-dat[!(dat$Modifications == "Unmodified"),]
        }
        
        #mtx <- as.matrix(dat[, grepl("^LFQ", names(dat))])
        mtx <- dat[, grepl("Sequence|Leading\ razor\ protein$|Modified\ sequence$|Experiment|Intensity", names(dat))]
        mtx$`Leading razor protein` <- sub(".*\\|","", sub('\\|([^\\|]*)$', '', mtx$`Leading razor protein`))
        
        mtx <- subset(mtx, `Leading razor protein` %in% protSelectionVec)
        
        #mtx$`Modified sequence`[duplicated(mtx)]
        #temp <- dat[dat$`Modified sequence` == "_EEDEEPES(ph)PPEK_",]
        
        # Summarize rows with same modification and same experiment id by taking average of intensities
        mtx.aggr <- aggregate(Intensity ~ Sequence + `Leading razor protein` +`Modified sequence` + Experiment, data=mtx, sum, na.rm=TRUE)
    
        # Transform to table with columns pepseq, PTMs, PTM type, accs, quant1, quant2, ...
        
    } else if (ex == "xlsx"){
        library("readxl")
        dat <- read_excel(filepath)
    }
}

# mtx[mtx == 0] <- NA
# mtx[mtx == "NaN"] <- NA

