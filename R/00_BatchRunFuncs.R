################################################################################
#                       QC batch test of the parameters                        #
################################################################################

# Function to install all necessary packages
install_phosfake <- function(pkgs = c(
    "purrr", "crayon", "parallel", "digest", "protr", "preprocessCore",
    "matrixStats", "extraDistr", "fdrtool", "qvalue", "limma", "moments",
    "Hmisc", "gplots"
), libpath = NULL) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    if (!is.null(libpath)) {
        .libPaths(libpath)
    }
    BiocManager::install(pkgs)
}

# Function to load all PhosFake scripts from given path and load libraries
load_phosfake <- function(path = "./") {
    source(file.path(path, "R/Parameter.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "R/01_GenerateGroundTruth.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "R/02_Digestion.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "R/03_MSRun.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "R/04_DataAnalysis.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "R/05_Statistics.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "R/06_Benchmarks.R"), echo = T, print.eval = TRUE, )
    library(crayon)
    library(parallel)
    library(digest)
    library(Hmisc)
    library(gplots)
}

# function to set the paths and general configs
set_phosfake <- function(fastaFilePath = "Proteomes", resultFilePath = "SimulatedDatasets", cores = 8, clusterType = "FORK", calcAllBenchmarks = T) {
    try(dir.create(resultFilePath))
    return(list(fastaFilePath = fastaFilePath, resultFilePath = resultFilePath, cores = cores, clusterType = clusterType, calcAllBenchmarks = calcAllBenchmarks))
}

# function to generate default parameters
def_param <- function() {
    paramsGroundTruth <- list(
        "PathToFasta" = "fasta_full_yeast.fasta",
        "PathToProteinList" = NA,
        "NumReps" = c(3),
        "NumCond" = 2,
        "FracModProt" = 0,
        "FracModPerProt" = 0,
        "PTMTypes" = NA,
        "PTMTypesDist" = NA,
        "PTMTypesMass" = NA,
        "PTMMultipleLambda" = NA,
        "ModifiableResidues" = NA,
        "ModifiableResiduesDistr" = NA,
        "RemoveNonModFormFrac" = 0
    )
    paramsProteoformAb <- list(
        "QuantNoise" = seq(0.5),
        "DiffRegFrac" = c(0.3),
        "DiffRegMax" = c(1),
        "UserInputFoldChanges" = NA,
        "UserInputFoldChanges_NumRegProteoforms" = NA,
        "UserInputFoldChanges_FoldChange" = NA,
        "ThreshNAProteoform" = -100,
        "AbsoluteQuanMean" = 30.5,
        "AbsoluteQuanSD" = 3.6,
        "ThreshNAQuantileProt" = 0.01
    )
    paramsDigest <- list(
        "Enzyme" = "trypsin",
        "PropMissedCleavages" = 0.01,
        "MaxNumMissedCleavages" = 2,
        "PepMinLength" = 7,
        "PepMaxLength" = 30,
        "LeastAbundantLoss" = 0,
        "EnrichmentLoss" = 0.2,
        "EnrichmentEfficiency" = 1,
        "EnrichmentNonModSignalLoss" = 0,
        "EnrichmentNoise" = 0.2
    )
    
    paramsMSRun <- list(
        "PercDetectedPep" = c(0.2),
        "PercDetectedVal" = c(0.5),
        "WeightDetectVal" = 0.1,
        "MSNoise" = c(0.25),
        "WrongIDs" = c(0.01),
        "WrongLocalizations" = 0.0,
        "MaxNAPerPep" = 1000
    )
    
    paramsDataAnalysis <- list(
        "ProtSummarization" = "medpolish",
        "MinUniquePep" = c(1),
        "StatPaired" = FALSE
    )
    
    ## Provide nice printout
    cat("--------------------\nGround truth generation parameters:\n")
    list.tree(paramsGroundTruth)
    cat("--------------------\nProteoform abundance parameters:\n")
    list.tree(paramsProteoformAb)
    cat("--------------------\nDigestion parameters:\n")
    list.tree(paramsDigest)
    cat("--------------------\nMSRun parameters:\n")
    list.tree(paramsMSRun)
    cat("--------------------\nData analysis parameters:\n")
    list.tree(paramsDataAnalysis)
    cat("--------------------\n")
    
    Param <- list(
        paramGroundTruth = paramsGroundTruth,
        paramProteoformAb = paramsProteoformAb,
        paramDigest = paramsDigest,
        paramMSRun = paramsMSRun,
        paramDataAnalysis = paramsDataAnalysis
    )
    
    return(Param)
}


# Function to generate combinations and convert rows to lists
generate_combinations <- function(params) {
    combinations <- expand.grid(params, stringsAsFactors = F)
    # seems to be the only way to get to a list of lists
    out_comb <- list()
    for (i in 1:nrow(combinations)) {
        out_comb[[i]] <- combinations[i, ]
    }
    out_comb
}



# Function to run all simulations
run_sims <- function(Parameters, Config) {
    # Generate combinations for each parameter set
    listtogroundtruth <- generate_combinations(Parameters$paramGroundTruth)
    listtoproteoformab <- generate_combinations(Parameters$paramProteoformAb)
    listtodigestion <- generate_combinations(Parameters$paramDigest)
    listtomsrun <- generate_combinations(Parameters$paramMSRun)
    listtodatanalysis <- generate_combinations(Parameters$paramDataAnalysis)
    
    all_params <- c(
        Parameters$paramGroundTruth,
        Parameters$paramProteoformAb,
        Parameters$paramDataAnalysis,
        Parameters$paramDigest,
        Parameters$paramMSRun
    )
    
    # Generate combinations for all parameters
    listall <- generate_combinations(all_params)
    totalbench <- length(listall)
    cat("Total number of simulations to run: ", totalbench, "\n")
    benchcounter <- 0
    
    # Gather always benchmarking data
    allBs <- NULL
    for (hh in 1:length(listtogroundtruth)) {
        ## Running ground thruth generations
        # check whether file with correct parameters exists
        tParam <- listtogroundtruth[[hh]]
        Param <- "none"
        groundTruth <- NULL
        md5 <- digest(tParam, algo = "md5")
        filename <- paste0(Config$resultFilePath, "/outputGroundTruth_", md5, ".RData")
        if (file.exists(filename)) {
            load(filename)
        } else {
            Param <- tParam
            # Add path to fasta file
            ttParam <- Param
            ttParam$PathToFasta <- paste0(Config$fastaFilePath, "/", Param$PathToFasta)
            groundTruth <- samplePreparation(parameters = ttParam)
            if (!is.null(groundTruth)) {
                save(Param, groundTruth, file = filename)
            }
        }
        gtParam <- Param
        ## quantitative proteoform abundances
        for (ii in 1:length(listtoproteoformab)) {
            Param <- "none"
            proteoformAb <- NULL
            # create combined parameterfile
            tParam <- c(gtParam, listtoproteoformab[[ii]])
            tParam <- c(tParam, list(QuantColnames = paste0("C_", rep(1:tParam$NumCond, each = tParam$NumReps), "_R_", rep(1:tParam$NumReps, tParam$NumCond))))
            # hash code to represent parameter configuration
            md5 <- digest(tParam, algo = "md5")
            filename <- paste0(Config$resultFilePath, "/outputProteoformAb_", md5, ".RData")
            if (file.exists(filename)) {
                load(filename)
            } else {
                Param <- tParam
                proteoformAb <- addProteoformAbundance(proteoforms = groundTruth, parameters = Param)
                save(Param, proteoformAb, file = filename)
            }
            pfParam <- Param
            
            ### Digestion
            for (jj in 1:length(listtodigestion)) {
                Param <- "none"
                BeforeMS <- NULL
                tParam <- c(pfParam, listtodigestion[[jj]])
                md5 <- digest(tParam, algo = "md5")
                filename <- paste0(Config$resultFilePath, "/outputDigestion_", md5, ".RData")
                if (file.exists(filename)) {
                    load(filename)
                } else {
                    Param <- tParam
                    peptable <- digestGroundTruth(proteoforms = proteoformAb, parameters = c(Param, list(Cores = Config$cores, ClusterType = Config$clusterType)))
                    peptable <- digestionProductSummarization(peptides = peptable, parameters = Param)
                    BeforeMS <- filterDigestedProt(DigestedProt = peptable, parameters = Param)
                    save(Param, BeforeMS, file = filename)
                }
                dgParam <- Param
                
                ### MS run
                for (kk in 1:length(listtomsrun)) {
                    Param <- "none"
                    AfterMSRun <- NULL
                    tParam <- c(dgParam, listtomsrun[[kk]])
                    md5 <- digest(tParam, algo = "md5")
                    filename <- paste0(Config$resultFilePath, "/outputMSRun_", md5, ".RData")
                    if (file.exists(filename)) {
                        load(filename)
                    } else {
                        Param <- tParam
                        AfterMSRun <- vector(mode = "list")
                        for (i in which(sapply(BeforeMS, length) > 0)) {
                            AfterMSRun[[length(AfterMSRun) + 1]] <- MSRunSim(Digested = BeforeMS[[i]], parameters = Param)
                        }
                        names(AfterMSRun) <- names(BeforeMS)[which(sapply(BeforeMS, length) > 0)]
                        save(Param, AfterMSRun, file = filename)
                    }
                    msParam <- Param
                    
                    ### Protein abundance
                    for (ll in 1:length(listtodatanalysis)) {
                        Param <- "none"
                        tParam <- c(msParam, listtodatanalysis[[ll]])
                        Benchmarks <- NULL
                        md5 <- digest(tParam, algo = "md5")
                        filename <- paste0(Config$resultFilePath, "/outputDataAnalysis_", md5, ".RData")
                        if (file.exists(filename)) {
                            load(filename)
                        } else {
                            # counter
                            Param <- tParam
                            Prots <- proteinSummarisation(peptable = AfterMSRun$NonEnriched, parameters = Param)
                            # Don't accept anything below 100 proteins
                            if (nrow(Prots) > 99) {
                                # Filter for having at least 1 actual value per protein group and peptide
                                Prots <- Prots[rowSums(is.na(Prots[, Param$QuantColnames])) < length(Param$QuantColnames), ]
                                allPeps <- as.data.frame(do.call("rbind", AfterMSRun))
                                allPeps <- allPeps[rowSums(is.na(allPeps[, Param$QuantColnames])) < length(Param$QuantColnames), ]
                                rownames(allPeps) <- paste0("pep", 1:nrow(allPeps))
                                Stats <- runPolySTest(Prots, Param, refCond = 1, onlyLIMMA = F)
                                # much faster with only LIMMA tests
                                StatsPep <- runPolySTest(allPeps, Param, refCond = 1, onlyLIMMA = T)
                                Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
                                save(Param, Stats, StatsPep, Benchmarks, file = filename)
                            } else {
                                print("Too few proteins!!!")
                                Benchmarks <- Stats <- StatsPep <- NULL
                                save(Param, Stats, StatsPep, Benchmarks, file = filename)
                            }
                            
                            if (Config$calcAllBenchmarks & !is.null(Stats)) {
                                Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
                                save(Param, Stats, StatsPep, Benchmarks, file = filename)
                            }
                        }
                        allBs[[md5]] <- list(Benchmarks = Benchmarks, Param = Param)
                    }
                }
            }
        }
    }
    return(allBs)
    cat("###### Finished data set generation \n")
}

# Retrieve intermediate outputs for a simulation
get_simulation <- function(Param, Config, stage="DataAnalysis") {
    # check whether file with correct parameters exists
    md5 <- digest(as.list(Param), algo = "md5")
    filename <- paste0(Config$resultFilePath, "/output", stage, "_", md5, ".RData")
    if (file.exists(filename)) {
        load(filename)
    } else {
        print("No simulation found with these parameters.")
    }
    return(list(Stats, StatsPep, Benchmarks))
}


# make matrix of benchmarks
matrix_benchmarks <- function(allBs, Config) {
    # extracting all benchmarks (sometimes there are more or less per run)
    t_allbnames <- NULL
    for (i in names(allBs)) {
        t_allbnames <- c(t_allbnames, names(unlist(allBs[[i]][[1]]$globalBMs)))
    }
    benchNames <- unique(t_allbnames)
    BenchMatrix <- data.frame(matrix(NA, ncol = length(benchNames) + length(Param), nrow = length(allBs)))
    colnames(BenchMatrix) <- c(benchNames, names(Param))
    rownames(BenchMatrix) <- names(allBs)
    
    # writing all results and parameters into matrix
    for (i in names(allBs)) {
        tglob <- unlist(allBs[[i]][[1]]$globalBMs)
        BenchMatrix[i, names(tglob)] <- tglob
        tpar <- allBs[[i]][[2]]
        BenchMatrix[i, names(tpar)] <- sapply(tpar, function(x) ifelse(length(x) > 1, paste0(x, collapse = "_"), x))
    }
    BenchMatrix
}

# Function to visualize benchmarks TODO libary(plotly)
visualize_benchmarks <- function(BenchMatrix) {
    # get parameter names
    param_table <- phosfake_params()
    param_names <- rownames(params_table)
    # Visualize roughly
    par(mfrow = c(nrow(BenchMatrix), 1), mar = c(2, 5, 2, 1), xpd = T, oma = c(5, 1, 1, 1))
    # Separate parameters and filter for actual values
    reds <- colnames(BenchMatrix) %in% param_names
    param_values <- BenchMatrix[, reds]
    BenchMatrix <- BenchMatrix[, -reds]
    to_del <- c()
    for (i in 1:ncol(BenchMatrix)) {
        tt <- unlist(BenchMatrix[, i])
        if (all(is.na(tt))) {
            to_del <- c(to_del, i)
        }
    }
    if(length(to_del) > 0)
        BenchMatrix <- BenchMatrix[, -to_del]
    # remove column QuantColnames if it exists
    if ("QuantColnames" %in% colnames(BenchMatrix)) {
        BenchMatrix <- BenchMatrix[, -which(colnames(BenchMatrix) == "QuantColnames")]
    }
    tBenchMatrix <- BenchMatrix
    # convert characters to factors
    for (i in 1:ncol(BenchMatrix)) {
        tt <- BenchMatrix[, i]
        if (is.numeric(tt)) {
            BenchMatrix[, i] <- round(tt, 2)
        }
        if (is.character(tt)) {
            tt <- as.factor(tt)
        }
        tt <- tt / max(as.numeric(tt), na.rm = T)
        tt[is.na(tt)] <- 0
        tBenchMatrix[, i] <- unlist(tt)
    }
    # define reference for x-axis
    for (sim in rownames(BenchMatrix)) {
        dat <- tBenchMatrix[sim, ]
        dat2 <- BenchMatrix[sim, ]
        if (is.numeric(dat2)) {
            dat2 <- round(dat2, 2)
        }
        
        midpoints <- barplot(unlist(dat),
                             las = 2, cex.names = 0.5, col = colorpanel(100, "#AA3333", "#3333AA")[as.numeric(dat) * 99 + 1],
                             main = sim, ylab = "Normalized values", xlab = "", ylim = c(0, 1),
                             xaxt = "none"
        )
        if (sim == rownames(BenchMatrix)[nrow(BenchMatrix)]) {
            axis(1, at = midpoints, labels = colnames(BenchMatrix), las = 2, cex.axis = 0.7)
        }
        # Add grid
        abline(v = midpoints, col = "lightgray", lty = 2)
        # Add vertical real number on top of each bar
        text(midpoints, dat, labels = as.character(dat2), pos = 3, cex = 0.7, col = 1, xpd = NA, srt = 45, offset = 0.5, adj = 0.5)
    }
    par(mfrow = c(1, 1))
}


plot_params <- function(BenchMatrix, current_frow) {
    
    # get parameter names
    param_table <- phosfake_params()
    param_names <- rownames(params_table)
    
    #filter out every NA or NULL parameters
    to_remove <- c()
    for (i in which(colnames(BenchMatrix) %in% param_names)) {
        val <- as.numeric(BenchMatrix[, i])
        if (!(length(val) ==0)) {
        if (all(is.na(val))) {
            to_remove <- append(to_remove, i)
        }
        }
    }
    if (length(to_remove) > 0)
        BenchMatrix <- BenchMatrix[, -to_remove]
    param_names <- param_names[param_names %in% colnames(BenchMatrix)]
    param_table <- param_table[param_names, ]

    # Set the number of columns
    ncols <- 2
    # Calculate number of rows based on number of plots and columns
    nrows <- ceiling(length(param_names) / ncols)
    
    # Create an empty list to store meters
    meters <- list()
    # Adjust spacing factors for left/right and top/bottom margins
    hspacing_factor <- 0.2
    vspacing_factor <- 0.5
    
    # Loop to create each meter plot
    for (i in 1:length(param_names)) {
        curr_par <- param_table[i, ]
        val <- as.numeric(BenchMatrix[current_row, param_names[i]])
        row_index <- ceiling(i / ncols)
        col_index <- (i - 1) %% ncols + 1
        x_domain <- c((col_index - 1 + hspacing_factor) / ncols, (col_index - hspacing_factor) / ncols)
        y_domain <- c(1 - (row_index - vspacing_factor) / nrows, 1 - ((row_index - 1 + vspacing_factor) / nrows))
        print(x_domain)
        print(y_domain)
        fig1 <- plot_ly(
            domain = list(x = x_domain, y = y_domain),
            value = val,
            title = list(text = param_names[i], align="center", font = list(size = 12)),
            type = "indicator",
            mode = "number+gauge",
            gauge = list(
                shape = "bullet",
                axis = list(range = list(curr_par$MinValue, curr_par$MaxValue),
                            tickmode = "auto",
                            ticks = "inside"  # Position ticks inside
                            #tickfont = list(size = 10, color = "black")  # Adjust tickfont size and color for better visibility
                ),
                bar = list(
                    color = colorpanel(100, "#33CC33", "#999966", "#CC3333")[(val - curr_par$MinValue) / (curr_par$MaxValue - curr_par$MinValue) * 100 + 1],
                    thickness = 0.6  # Adjust this value to make the gauge bar thicker
                )
            )
        ) %>%
            layout(margin = list(l = 20, r = 20))
        meters[[i]] <- fig1
    }
    # Combine all meter plots
    subplot(meters, nrows = nrows)
}
