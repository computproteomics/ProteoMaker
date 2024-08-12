

# Wrapper to call statistical tests (set to only RefCond condition as reference)
runPolySTest <- function(fullData, Param, refCond, onlyLIMMA=F, cores=1) {
    library(SummarizedExperiment)
    Data <- fullData[,Param$QuantColnames]
    NumCond <- Param$NumCond
    NumReps <- Param$NumReps
    Reps <- rep(seq_len(NumCond), NumReps)
    isPaired <- Param$StatPaired
    
    # Set number of threads
    Sys.setenv(SHINY_THREADS=cores)
    
    cat(" + Running statistical tests\n")
    
    
    # Rearrange order of colums to grouped replicates
    act_cols <- rep(0:(NumCond-1),NumReps)*NumReps+rep(1:(NumReps), each=NumCond)
    Data <- Data[,act_cols]
    
    # Make SummarizedExperiment object
    sampleMetadata <- data.frame(Condition = rep(paste("Condition", 1:NumCond), each=NumReps),
                                 Replicate = rep(1:NumReps, each=NumCond))
    
    fulldata <- SummarizedExperiment(assays = list(quant = Data), 
                                     colData = sampleMetadata)
    rowData(fulldata) <- rownames(Data)
    metadata(fulldata) <- list(NumReps = NumReps, NumCond = NumCond)
    
    assay(fulldata, "quant") <- Data
    
    # Generate experimental design
    conditions <- unique(colData(fulldata)$Condition)
    allComps <- suppressMessages(PolySTest::create_pairwise_comparisons(conditions, 1))
    
    
    # Run tests
    testNames <- NULL
    if (onlyLIMMA) {
        testNames <- "limma"
    } else {
        testNames <- c("limma", "Miss_Test", "t_test", "rank_products", "permutation_test")
    }
    results <- NULL
    if (isPaired) {
        results <- PolySTest::PolySTest_paired(fulldata, allComps, testNames)
    } else {
        results <- PolySTest::PolySTest_unpaired(fulldata, allComps, testNames)
    }
    
    # Define comparisons to visualize from available ones
    compNames <- metadata(results)$compNames
    
    
    # Preparing data
    cnames <- colnames(rowData(results))
    rdata <- as.data.frame(rowData(results))
    LogRatios <- rdata[, grep("^log_ratios_", cnames), drop=F]
    names(LogRatios) <- gsub("^log_ratios", "log-ratios", names(LogRatios))
    names(LogRatios) <- gsub("_Condition\\.", " ", names(LogRatios))
    names(LogRatios) <- gsub("_vs", " vs", names(LogRatios))
    Pvalue <- rdata[, grep("^p_values_", cnames), drop=F]
    if (onlyLIMMA) {
        Qvalue <- rdata[, grep("^FDR_limma_", cnames), drop=F]
    } else {
        Qvalue <- rdata[, grep("^FDR_PolySTest_", cnames), drop=F]
    }
    names(Qvalue) <- gsub("_Condition\\.", " ", names(Qvalue))
    names(Qvalue) <- gsub("_vs", " vs", names(Qvalue))
    
    
    fullData$Regulation_Amplitude <- sapply(fullData$Regulation_Amplitude, function(x) paste(x, collapse=";"))
    fullData$Regulation_Pattern <- sapply(fullData$Regulation_Pattern, function(x) paste(x, collapse=";"))
    FullReg <- cbind(LogRatios, Qvalue, as.data.frame(fullData))#, WhereReg)
    
    FullReg <- cbind(FullReg, 
                     min1Reg=sapply(str_split(fullData$Regulation_Amplitude, ";"),
                                    function(x) {
                                        y <- as.numeric(ifelse(x == "NA", NA, x))
                                        sum(as.numeric(y),na.rm=T)
                                    }) != 0, 
                     allReg=sapply(str_split(fullData$Regulation_Amplitude, ";"),
                                   function(x) {
                                       y <- as.numeric(ifelse(x == "NA", NA, x))
                                       !is.na(sum(as.numeric(y),na.rm=T))
                                   }))

    return(FullReg)
}
