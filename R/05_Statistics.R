#' Run PolySTest Statistical Tests on Quantitative Data
#'
#' This function runs a set of statistical tests on quantitative data using the `PolySTest` package. The function handles paired and unpaired tests, reorganizes the data according to experimental conditions, and returns a comprehensive table with log ratios, p-values, and FDR values for the comparisons. It allows the user to specify the number of cores for parallel computation and to limit the tests to only the `limma` method if desired.
#'
#' @param fullData A data frame containing quantitative data, with rows representing features (e.g., proteins or peptides) and columns representing samples. The data frame should also include metadata columns such as `Regulation_Amplitude` and `Regulation_Pattern`.
#' @param Param A list of parameters, including:
#' \describe{
#'   \item{QuantColnames}{The names of the columns containing the quantitative data to be analyzed.}
#'   \item{NumCond}{The number of conditions in the experiment.}
#'   \item{NumReps}{The number of replicates per condition.}
#'   \item{StatPaired}{A logical value indicating whether the experiment is paired.}
#' }
#' @param refCond A numeric or character value indicating the reference condition to be used in the statistical comparisons.
#' @param onlyLIMMA A logical value indicating whether to restrict the analysis to only the `limma` test. Default is `FALSE`, which means all available tests in `PolySTest` are run.
#' @param cores An integer specifying the number of cores to use for parallel processing. Default is 1.
#'
#' @return A data frame with the original data and additional columns for log ratios, p-values, FDR values, and ground truth indicators for differential regulation.
#'
#' @import SummarizedExperiment 
#' @importFrom PolySTest create_pairwise_comparisons PolySTest_paired PolySTest_unpaired
#' @importFrom stringr str_split
#' @importFrom parallel detectCores makeCluster setDefaultCluster clusterExport parLapply stopCluster
#'
#' @keywords internal
#' 
runPolySTest <- function(fullData, Param, refCond, onlyLIMMA=F, cores=1) {
    # library(SummarizedExperiment)
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
    
    fulldata <- SummarizedExperiment::SummarizedExperiment(assays = list(quant = Data), 
                                     colData = sampleMetadata)
    SummarizedExperiment::rowData(fulldata) <- rownames(Data)
    metadata(fulldata) <- list(NumReps = NumReps, NumCond = NumCond)
    
    SummarizedExperiment::assay(fulldata, "quant") <- Data
    
    # Generate experimental design
    conditions <- unique(SummarizedExperiment::colData(fulldata)$Condition)
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
        results <- suppressWarnings(PolySTest::PolySTest_paired(fulldata, allComps, testNames))
    } else {
        results <- suppressWarnings(PolySTest::PolySTest_unpaired(fulldata, allComps, testNames))
    }
    
    # Define comparisons to visualize from available ones
    compNames <- metadata(results)$compNames
    
    
    # Preparing data
    cnames <- colnames(SummarizedExperiment::rowData(results))
    rdata <- as.data.frame(SummarizedExperiment::rowData(results))
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
                     min1Reg=sapply(stringr::str_split(fullData$Regulation_Amplitude, ";"),
                                    function(x) {
                                        y <- as.numeric(ifelse(x == "NA", NA, x))
                                        sum(as.numeric(y),na.rm=T)
                                    }) != 0, 
                     allReg=sapply(stringr::str_split(fullData$Regulation_Amplitude, ";"),
                                   function(x) {
                                       y <- as.numeric(ifelse(x == "NA", NA, x))
                                       !is.na(sum(as.numeric(y),na.rm=T))
                                   }))

    return(FullReg)
}
