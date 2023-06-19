source("viqor/Functions.R")
### Method which takes a dataframe of simulated peptide abundance values and visualizes it using linear plotting, correlation maps and heatmaps

VisualizeMs <- function(afterMsTotal, groupMaxSize, condAndRepStart, condAndRepEnd, filename="", savePath, numCond, makeHeatmap = F){
  proteinGroupData <- as.data.frame(select(afterMsTotal, c(1, condAndRepStart:condAndRepEnd)))
  
  groups <- protein.Grouping(peptide.Data.Frame = proteinGroupData, fasta = paramGroundTruth$PathToFasta, parsimony = "no")
  member_count <- count_numbers(groups[["cc"]][["membership"]][1: length(groups[["proteins"]])])# Count the number of times each protein group value appears in the grouping

  AccessionsOfInterest <- member_count[member_count$`Protein group size` <= groupMaxSize, ]
  AccessionsOfInterestSeperated <- tidyr::separate_rows(AccessionsOfInterest, Accession, sep = ",")
  
  AfterMSRunUnnest <- unnest(afterMsTotal, c(Peptide, Start, Stop, MC, Accession, Proteoform_ID, Regulation_Amplitude, Regulation_Pattern))
  AfterMSRunRegulatedUnnest <- unnest(AfterMSRunRegulated, c(Peptide, Start, Stop, MC, Accession, Proteoform_ID, Regulation_Amplitude, Regulation_Pattern))
  
  AccessionFreqTable <- as.data.frame(table(AfterMSRunUnnest$Accession))
  colnames(AccessionFreqTable) <- c("Accession", "# of peptides")

  joinDataset <- merge(AccessionFreqTable, AccessionsOfInterestSeperated, by = "Accession")

  topProteinVal <- joinDataset[which.max(joinDataset$`# of peptides`),]
  proteinToAnalyze <- afterMsTotal[grep(topProteinVal$Accession, afterMsTotal$Accession),]
  proteinToAnalyzeCentered <- proteinToAnalyze
  proteinToAnalyzeCentered[,condAndRepStart:condAndRepEnd] <- proteinToAnalyzeCentered[,condAndRepStart:condAndRepEnd] - rowMeans(proteinToAnalyzeCentered[,condAndRepStart:condAndRepEnd], na.rm = T)
  
  proteinToAnalyzeCenteredOutput <- proteinToAnalyzeCentered[,condAndRepStart:condAndRepEnd]
  cor_matr <- cor(t(proteinToAnalyzeCenteredOutput), use="pairwise.complete.obs")
  
  #Correlation Matrix
  pheatmap(cor_matr,
           main = "Correlation of proteoform group",
           color = colorRampPalette(c("navy", "white", "red"))(50),
           #cutree_rows = length(unique(proteinToAnalyzeCentered$Proteoform_ID)),
           #cutree_cols = length(unique(proteinToAnalyzeCentered$Proteoform_ID)),
           symm = T,
           show_rownames = T, 
           show_colnames = T,
           filename = paste(savePath,format(Sys.time(), "%c"),filename,"Cor_matrix_heatmap.pdf"))    
  
  proteinToAnalyzeCenteredGgplot <- melt(t(proteinToAnalyzeCentered[,c(condAndRepStart:condAndRepEnd)])) #makes the data change from wide format to long format
  proteinToAnalyzeCenteredGgplot$Proteoform_ID <-  rep(as.character(proteinToAnalyzeCentered$Proteoform_ID), each = condTimesRep) #the list of proteoform ID's are made to be characters and then used as a factor when plotting
  proteinToAnalyzeCenteredGgplot$Sequence <-  rep(proteinToAnalyzeCentered$Sequence, each = condTimesRep)
  
  #Linear plot
  ggplot(proteinToAnalyzeCenteredGgplot, aes(x = Var1,
                                             y = value,
                                             group = Var2,
                                             col = Proteoform_ID)) +
    ggtitle( paste("Linear plot of centered abundances of protein", topProteinVal$Accession)) +
    geom_line(data = proteinToAnalyzeCenteredGgplot[!is.na(proteinToAnalyzeCenteredGgplot$value),]) + 
    geom_point() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(col="Proteoform group", x = "Experimental Condition and Replicate", y = "Centered abundances")
  
    ggsave(filename = paste(format(Sys.time(), "%c"),filename,"Linear_plot.pdf"), path = savePath, width = 12, height = 8,  units = "in")

  if (makeHeatmap) {  
    NewAfterFiltering <- AfterMSRunRegulatedUnnest[rowSums(is.na(AfterMSRunRegulatedUnnest[,condAndRepStart:condAndRepEnd])) < 3 , ] #CURRENTLY ONLY WITH 3 OR LESS NA
    NewAfterFilteringAll <- AfterMSRunUnnest[rowSums(is.na(AfterMSRunUnnest[,condAndRepStart:condAndRepEnd])) < 3 , ] #getting the data we want CURRENTLY ONLY WITH 3 OR LESS NA
    NewAfterFilteringAllRowmeans <- rowMeans(NewAfterFilteringAll[,condAndRepStart:condAndRepEnd], na.rm = T) #calculating the rowmeans
    dataForPheatmap <- NewAfterFilteringAll[,condAndRepStart:condAndRepEnd] #saving the information of interest in a separate dataframe
    dataForPheatmap <- dataForPheatmap[,1:condTimesRep] - NewAfterFilteringAllRowmeans # subtracting rowmeans
    rownames(dataForPheatmap) <- rownames(dataForPheatmap)
    isDiffReg <- as.data.frame(NewAfterFilteringAll$Accession %in% NewAfterFiltering$Accession) #saving the information of wether a peptide was diffregulated
    isDiffReg[isDiffReg == T] <- "Regulated"
    isDiffReg[isDiffReg == F] <- "Not regulated"
    colnames(isDiffReg) <- "Regulation" 
    
    #Heatmap of experiment
    pheatmap(dataForPheatmap[,1:condTimesRep],
             main = "All peptides identified in simulated MS run",
             color = colorRampPalette(c("navy", "white", "red"))(50),
             cutree_rows = length(unique(proteinToAnalyzeCentered$Proteoform_ID)),
             cutree_cols = numCond,
             annotation_row = isDiffReg,
             show_rownames = F,
             legend = T,
             kmeans_k = NA,
             filename = paste(savePath,format(Sys.time(), "%c"),"All_pep_heatmap.pdf"))
  }
  return(proteinToAnalyzeCenteredOutput)
}


count_numbers <- function(vec){
  df <- data.frame(num=vec, names=names(vec), stringsAsFactors = FALSE) #make the vector a data frame
  df <- aggregate(df$names, by = list(df$num), paste, collapse = ", ") #use the aggregate function to subset the data easily
  colnames(df) <- c("num", "name")
  df$name_count <- sapply(strsplit(as.character(df$name), ", "), length)
  result <- merge(df, data.frame(num=unique(vec)), all = TRUE)
  result[is.na(result)] <- ""
  colnames(result) <- c("protein group","Accession","Protein group size")
  return(result)
}