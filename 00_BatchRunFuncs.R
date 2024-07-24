################################################################################
#                       QC batch test of the parameters                        #
################################################################################

# Function to install all necessary packages
install_phosfake <- function(pkgs = c(
                                 "purrr", "crayon", "parallel", "digest", "protr", "preprocessCore",
                                 "matrixStats", "extraDistr", "fdrtool", "qvalue", "limma", "moments",
                                 "Hmisc","gplots"
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
    source(file.path(path, "Parameter.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "01_GenerateGroundTruth.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "02_Digestion.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "03_MSRun.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "04_DataAnalysis.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "05_Statistics.R"), echo = T, print.eval = TRUE)
    source(file.path(path, "06_Benchmarks.R"), echo = T, print.eval = TRUE, )
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
get_simulation <- function(Param, Config) {
    # check whether file with correct parameters exists
    md5 <- digest(as.list(Param), algo = "md5")
    filename <- paste0(Config$resultFilePath, "/outputDataAnalysis_", md5, ".RData")
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

# Function to visualize benchmarks
visualize_benchmarks <- function(BenchMatrix) {
    # Visualize roughly
    par(mfrow = c(nrow(BenchMatrix), 1))
    tBenchMatrix <- BenchMatrix
    # convert characters to factors
    for (i in 1:ncol(BenchMatrix)) {
        tt <- BenchMatrix[, i]
        if (is.numeric(tt)) {
            BenchMatrix[,i] <- round(tt, 2)
        }
        if (is.character(tt)) {
            tt <- as.factor(tt)
        }
        tt <- tt / max(as.numeric(tt), na.rm = T)
        tt[is.na(tt)] <- 0
        tBenchMatrix[, i] <- tt
    }
    # define reference for x-axis
    for (sim in rownames(BenchMatrix)) {
        dat <- tBenchMatrix[sim, ]
        dat2 <- BenchMatrix[sim, ]
        if (is.numeric(dat2)) {
            dat2 <- round(dat2, 2)
        }

        print(as.vector(dat))
        midpoints <- barplot(as.vector(dat),
            las = 2, cex.names = 0.5, col = colorpanel(100, "#AA3333", "#3333AA")[as.numeric(dat)*99+1],
            main = sim, ylab = "Normalized values", xlab = "", ylim = c(0, 1),
            xaxt="none"
        )
        if(sim == rownames(BenchMatrix)[nrow(BenchMatrix)]){
            axis(1, at = midpoints, labels = colnames(BenchMatrix), las = 2, cex.axis = 0.5)
        }
        # Add grid
        abline(h = midpoints, col = "lightgray", lty = 2)
        # Add vertical real number on top of each bar
        text(midpoints, dat, labels = as.character(dat2), pos = 3, cex = 0.5, col = 1, xpd = NA, srt = 90)
    }
    par(mfrow = c(1, 1))
}
