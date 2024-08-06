################################################################################
#                       QC batch test of the parameters                        #
################################################################################
# 
# # Function to install all necessary packages
# install_phosfake <- function(pkgs = c(
#     "purrr", "crayon", "parallel", "digest", "protr", "preprocessCore",
#     "matrixStats", "extraDistr", "fdrtool", "qvalue", "limma", "moments",
#     "Hmisc", "gplots", "plotly", "yaml"
# ), libpath = NULL) {
#     if (!requireNamespace("BiocManager", quietly = TRUE)) {
#         install.packages("BiocManager")
#     }
#     if (!is.null(libpath)) {
#         .libPaths(libpath)
#     }
#     BiocManager::install(pkgs)
# }
# 
# # Function to load all PhosFake scripts from given path and load libraries
# load_phosfake <- function(path = "./") {
#     source(file.path(path, "R/Parameter.R"), echo = T, print.eval = TRUE)
#     source(file.path(path, "R/01_GenerateGroundTruth.R"), echo = T, print.eval = TRUE)
#     source(file.path(path, "R/02_Digestion.R"), echo = T, print.eval = TRUE)
#     source(file.path(path, "R/03_MSRun.R"), echo = T, print.eval = TRUE)
#     source(file.path(path, "R/04_DataAnalysis.R"), echo = T, print.eval = TRUE)
#     source(file.path(path, "R/05_Statistics.R"), echo = T, print.eval = TRUE)
#     source(file.path(path, "R/06_Benchmarks.R"), echo = T, print.eval = TRUE, )
#     library(crayon)
#     library(parallel)
#     library(digest)
#     library(Hmisc)
#     library(gplots)
# }

#' Set paths and general configuration for PhosFake
#'
#' This function sets the paths and general configuration parameters needed for
#' running simulations with the PhosFake package. It creates the necessary output
#' directories and returns a list of configuration settings.
#'
#' @param fastaFilePath A character string specifying the path to the directory
#' containing the FASTA files. Default is `"Proteomes"`.
#' @param resultFilePath A character string specifying the path to the directory
#' where the simulation results will be saved. Default is `"SimulatedDatasets"`.
#' @param cores An integer specifying the number of cores to use for parallel
#' processing. Default is 2.
#' @param clusterType A character string specifying the type of cluster to use
#' for parallel computing (e.g., `"FORK"`, `"PSOCK"`). Default is `"FORK"`.
#' @param calcAllBenchmarks A logical value indicating whether to calculate all
#' benchmarks. Default is `TRUE`.
#' 
#' @return A list containing the configuration settings:
#' \describe{
#'   \item{fastaFilePath}{The path to the FASTA files.}
#'   \item{resultFilePath}{The path to the result files.}
#'   \item{cores}{The number of cores used.}
#'   \item{clusterType}{The type of cluster used.}
#'   \item{calcAllBenchmarks}{Whether to calculate all benchmarks.}
#' }
#' 
#' @export
#'
#' @examples
#' config <- set_phosfake()
#' config <- set_phosfake(fastaFilePath = "CustomProteomes", resultFilePath = "Results", cores = 4, clusterType = "PSOCK")
set_phosfake <- function(fastaFilePath = "Proteomes", resultFilePath = "SimulatedDatasets", cores = 2, clusterType = "FORK", calcAllBenchmarks = T) {
    try(dir.create(resultFilePath))
    return(list(fastaFilePath = fastaFilePath, resultFilePath = resultFilePath, cores = cores, clusterType = clusterType, calcAllBenchmarks = calcAllBenchmarks))
}

#' #' Generate default parameters for PhosFake simulations
#' #'
#' #' This function generates and returns a list of default parameters for running
#' #' simulations with the PhosFake package. The parameters are organized into
#' #' several categories, including ground truth generation, proteoform abundance,
#' #' digestion, MS run, and data analysis.
#' #'
#' #' @return A list of default parameters categorized as follows:
#' #' \describe{
#' #'   \item{paramGroundTruth}{Parameters related to ground truth data generation.}
#' #'   \item{paramProteoformAb}{Parameters related to proteoform abundance.}
#' #'   \item{paramDigest}{Parameters related to enzymatic digestion.}
#' #'   \item{paramMSRun}{Parameters related to the mass spectrometry run.}
#' #'   \item{paramDataAnalysis}{Parameters related to data analysis.}
#' #' }
#' #' 
#' #' @export
#' #'
#' #' @examples
#' #' params <- def_param()
#' #' params$paramGroundTruth$NumReps <- c(2, 4, 6)
#' def_param <- function() {
#'     paramsGroundTruth <- list(
#'         "PathToFasta" = "fasta_full_yeast.fasta",
#'         "PathToProteinList" = NA,
#'         "NumReps" = c(3),
#'         "NumCond" = 2,
#'         "FracModProt" = 0,
#'         "FracModPerProt" = 0,
#'         "PTMTypes" = NA,
#'         "PTMTypesDist" = NA,
#'         "PTMTypesMass" = NA,
#'         "PTMMultipleLambda" = NA,
#'         "ModifiableResidues" = NA,
#'         "ModifiableResiduesDistr" = NA,
#'         "RemoveNonModFormFrac" = 0
#'     )
#'     paramsProteoformAb <- list(
#'         "QuantNoise" = seq(0.5),
#'         "DiffRegFrac" = c(0.3),
#'         "DiffRegMax" = c(1),
#'         "UserInputFoldChanges" = NA,
#'         "UserInputFoldChanges_NumRegProteoforms" = NA,
#'         "UserInputFoldChanges_FoldChange" = NA,
#'         "ThreshNAProteoform" = -100,
#'         "AbsoluteQuanMean" = 30.5,
#'         "AbsoluteQuanSD" = 3.6,
#'         "ThreshNAQuantileProt" = 0.01
#'     )
#'     paramsDigest <- list(
#'         "Enzyme" = "trypsin",
#'         "PropMissedCleavages" = 0.01,
#'         "MaxNumMissedCleavages" = 2,
#'         "PepMinLength" = 7,
#'         "PepMaxLength" = 30,
#'         "LeastAbundantLoss" = 0,
#'         "EnrichmentLoss" = 0.2,
#'         "EnrichmentEfficiency" = 1,
#'         "EnrichmentNonModSignalLoss" = 0,
#'         "EnrichmentNoise" = 0.2
#'     )
#'     
#'     paramsMSRun <- list(
#'         "PercDetectedPep" = c(0.2),
#'         "PercDetectedVal" = c(0.5),
#'         "WeightDetectVal" = 0.1,
#'         "MSNoise" = c(0.25),
#'         "WrongIDs" = c(0.01),
#'         "WrongLocalizations" = 0.0,
#'         "MaxNAPerPep" = 1000
#'     )
#'     
#'     paramsDataAnalysis <- list(
#'         "ProtSummarization" = "medpolish",
#'         "MinUniquePep" = c(1),
#'         "StatPaired" = FALSE
#'     )
#'     
#'     ## Provide nice printout
#'     cat("--------------------\nGround truth generation parameters:\n")
#'     list.tree(paramsGroundTruth)
#'     cat("--------------------\nProteoform abundance parameters:\n")
#'     list.tree(paramsProteoformAb)
#'     cat("--------------------\nDigestion parameters:\n")
#'     list.tree(paramsDigest)
#'     cat("--------------------\nMSRun parameters:\n")
#'     list.tree(paramsMSRun)
#'     cat("--------------------\nData analysis parameters:\n")
#'     list.tree(paramsDataAnalysis)
#'     cat("--------------------\n")
#'     
#'     Param <- list(
#'         paramGroundTruth = paramsGroundTruth,
#'         paramProteoformAb = paramsProteoformAb,
#'         paramDigest = paramsDigest,
#'         paramMSRun = paramsMSRun,
#'         paramDataAnalysis = paramsDataAnalysis
#'     )
#'     
#'     return(Param)
#' }


#' Generate default parameters or read from yaml file
#'
#' This function reads a YAML file containing parameter settings and extracts
#' them into a structured list categorized by their types. The YAML file should
#' have a specific structure, with parameters organized under a top-level 
#' `params` key, and each parameter containing metadata including its type 
#' and default value.
#'
#' @param yaml_file A character string specifying the path to the YAML file. 
#' NULL means setting default values.
#' 
#' @return A list containing categorized default parameters:
#' \describe{
#'   \item{paramGroundTruth}{Parameters related to ground truth data generation.}
#'   \item{paramProteoformAb}{Parameters related to proteoform abundance.}
#'   \item{paramDigest}{Parameters related to enzymatic digestion.}
#'   \item{paramMSRun}{Parameters related to the mass spectrometry run.}
#'   \item{paramDataAnalysis}{Parameters related to data analysis.}
#' }
#' 
#' @export
#'
#' @examples
#' # Read YAML file from inst folder
#' yaml_path <- system.file("config", "params.yaml", package = "PhosFake")
#' params <- def_param_from_yaml(yaml_path)
def_param <- function(yaml_file=NULL) {
    
    if (is.null(yaml_file)) {
        yaml_file <- system.file("config", "parameters.yaml", package = "PhosFake")
    }
    # Read the YAML file
    params <- yaml::yaml.load_file(yaml_file)$params
    
    # Convert NA values from strings to real NA
    for (l in names(params)) {
        print(params[[l]])
        for (k in names(params[[l]])) {
            if (params[[l]][[k]] == "NA") {
                params[[l]][[k]] <- NA
            }
        }
    }
    
    
    # Initialize empty lists for each category
    Param <- list(
        paramGroundTruth = list(),
        paramProteoformAb = list(),
        paramDigest = list(),
        paramMSRun = list(),
        paramDataAnalysis = list()
    )
    
    # Iterate over each parameter and place it into the correct category based on "type"
    for (param_name in names(params)) {
        param_info <- params[[param_name]]
        category <- param_info$type
        default_value <- param_info$default
        
        if (!is.null(category)) {
            Param[[category]][[param_name]] <- default_value
        }
    }
    
    # Provide a nice printout for each category
    cat("--------------------\nGround truth generation parameters:\n")
    list.tree(Param$paramGroundTruth)
    cat("--------------------\nProteoform abundance parameters:\n")
    list.tree(Param$paramProteoformAb)
    cat("--------------------\nDigestion parameters:\n")
    list.tree(Param$paramDigest)
    cat("--------------------\nMSRun parameters:\n")
    list.tree(Param$paramMSRun)
    cat("--------------------\nData analysis parameters:\n")
    list.tree(Param$paramDataAnalysis)
    cat("--------------------\n")
    
    return(Param)
}



#' Generate all combinations of parameter sets
#'
#' This function takes a list of parameter vectors and generates all possible 
#' combinations. The result is converted into a list of lists, where each 
#' sublist represents a unique combination of the provided parameters.
#'
#' @param params A list of parameter vectors, where each vector represents 
#' different possible values for a parameter.
#' 
#' @return A list of lists, each sublist containing a unique combination of 
#' parameters.
#' 
#' @export
#'
#' @examples
#' params <- def_param()
#' params$paramGroundTruth$NumReps <- c(2, 4, 6)
#' combinations <- generate_combinations(params)
generate_combinations <- function(params) {
    combinations <- expand.grid(params, stringsAsFactors = F)
    # seems to be the only way to get to a list of lists
    out_comb <- list()
    for (i in 1:nrow(combinations)) {
        out_comb[[i]] <- combinations[i, ]
    }
    out_comb
}



#' Run all simulations for given parameters and configurations
#'
#' This function runs simulations based on the provided parameter sets and 
#' configuration settings. It runs the given parameter settings, processes 
#' ground truth data, simulates proteoform abundance, performs digestion, 
#' conducts MS runs, and carries out data analysis. Results and benchmarks are 
#' saved to the output directory specified in \code{def_param}.
#'
#' @param Parameters A list containing categorized parameter sets, typically 
#' generated using \code{def_param} or \code{def_param_from_yaml}.
#' @param Config A list containing configuration settings, such as file paths 
#' and computational settings, typically generated using \code{set_phosfake}.
#' 
#' @return A list of results, where each element contains simulation results 
#' and benchmarks for a specific parameter combination.
#' 
#' @importFrom digest digest
#' @export
#'
#' @examples
#' params <- def_param()
#' config <- set_phosfake()
#' results <- run_sims(params, config)
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



#' Retrieve intermediate outputs for a simulation
#'
#' This function retrieves intermediate outputs for a simulation based on the
#' provided parameters and configuration. It checks for the existence of the 
#' output file corresponding to a specific simulation stage and loads the 
#' data if available.
#'
#' @param Param A list of parameters used in the simulation, typically a subset 
#' of those generated by functions like \code{def_param} or \code{def_param_from_yaml}.
#' @param Config A list of configuration settings, including file paths and 
#' computational settings, typically generated using \code{set_phosfake}.
#' @param stage A character string indicating the stage of the simulation for 
#' which the output is requested. Default is \code{"DataAnalysis"}. Values are \code{"GroundTruth"},
#' \code{"ProteoformAb"}, \code{"Digestion"}, \code{"MSRun"}, and \code{"DataAnalysis"}.
#' 
#' @return A list containing the loaded data objects, including \code{Stats}, 
#' \code{StatsPep}, and \code{Benchmarks}, if they exist. If no file is found, 
#' the function returns \code{NULL}.
#' 
#' @importFrom digest digest
#' @export
#'
#' @examples
#' config <- set_phosfake()
#' param <- def_param()$paramGroundTruth
#' results <- run_sims(param, config)
#' output <- get_simulation(param, config, stage = "MSRun")
get_simulation <- function(Param, Config, stage="DataAnalysis") {
    # check whether file with correct parameters exists
    md5 <- digest(as.list(Param), algo = "md5")
    filename <- paste0(Config$resultFilePath, "/output", stage, "_", md5, ".RData")
    if (file.exists(filename)) {
        # Create a temporary environment to load the objects
        temp_env <- new.env()
        # Load the objects into the temporary environment
        load(filename, envir = temp_env)
        
        # Get the names of the objects loaded into the temporary environment
        object_names <- ls(temp_env)
        
        # Create a list to store the objects
        objects_list <- lapply(object_names, function(x) get(x, envir = temp_env))
        
        # Set the names of the list elements
        names(objects_list) <- object_names
        
        # Return the list of objects
        return(objects_list)        
        
    } else {
        print("No simulation found with these parameters.")
    }
}


#' Create a matrix of benchmarks from simulation results
#'
#' This function compiles a matrix of benchmark results from the outputs of multiple 
#' simulations. It extracts benchmark metrics and associated parameter values, 
#' creating a comprehensive data frame for analysis.
#'
#' @param allBs A list of results from the \code{run_sims} function, where each 
#' element contains simulation results and benchmarks for a specific parameter 
#' combination.
#' @param Config A list of configuration settings, including file paths and 
#' computational settings, typically generated using \code{set_phosfake}.
#' 
#' @return A data frame with rows representing each simulation and columns 
#' representing benchmark values and parameter settings.
#' 
#' @export
#'
#' @examples
#' results <- run_sims(def_param(), set_phosfake())
#' benchmark_matrix <- matrix_benchmarks(results, set_phosfake())
#' 
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

#' Visualize benchmarks for a specific simulation
#'
#' This function creates visualizations of the benchmark results for a specific
#' simulation. It uses various plotting methods to display the normalized values
#' of benchmarks and parameters, aiding in the analysis of simulation outcomes.
#'
#' @param BenchMatrix A data frame containing the benchmark results for multiple 
#' simulations, typically generated by \code{matrix_benchmarks}.
#' @param current_row An integer or a character string specifying the row number 
#' or the simulation identifier for which the benchmarks are to be visualized. 
#' Default is 1.
#' 
#' @return A plot object showing the visualized benchmarks and parameters for 
#' the specified simulation.
#' 
#' @importFrom plotly plot_ly subplot
#' @importFrom gplots colorpanel
#' @export
#'
#' @examples
#' benchmarks <- matrix_benchmarks(run_sims(def_param(), set_phosfake()), set_phosfake())
#' visualize_benchmarks(benchmarks, current_row = 1)
visualize_benchmarks <- function(BenchMatrix, current_row = 1) {
    # get parameter names
    param_table <- param_table()
    param_names <- rownames(param_table)
    
    # Visualize roughly
    par(mfrow = c(nrow(BenchMatrix), 1), mar = c(2, 5, 2, 1), xpd = T, oma = c(5, 1, 1, 1))
    # Separate parameters and filter for actual values
    reds <- colnames(BenchMatrix) %in% param_names
    param_values <- BenchMatrix[, reds]
    BenchMatrix <- BenchMatrix[, !reds]
    to_del <- c()
    for (i in 1:ncol(BenchMatrix)) {
        tt <- unlist(BenchMatrix[, i])
        if (all(is.na(tt))) {
            to_del <- c(to_del, i)
        }
    }
    if(length(to_del) > 0)
        BenchMatrix <- BenchMatrix[, -to_del]
    BenchMatrix[is.na(BenchMatrix)] <- 0
    # remove column QuantColnames if it exists
    if ("QuantColnames" %in% colnames(BenchMatrix)) {
        BenchMatrix <- BenchMatrix[, -which(colnames(BenchMatrix) == "QuantColnames")]
    }
    tBenchMatrix <- BenchMatrix
    nr <- 2# nrow(BenchMatrix)
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
    
    # Define color palette
    color_palette <- colorpanel(100, "#AA3333", "#3333AA")
    
    # Create plots
    plots <- list()
    sim <- current_row
    dat <- tBenchMatrix[sim, ]
    if (is.numeric(dat)) {
        dat <- round(dat, 2)
    }
    if (!is.numeric(dat)) {
        dat <- 0
    }
    dat2 <- BenchMatrix[sim, ]
    
    if (is.numeric(dat2)) {
        dat2 <- round(dat2, 2)
    }
    
    plot <- plot_ly(
        x = names(dat),
        y = as.numeric(dat),
        type = 'bar',
        marker = list(
            color = color_palette[as.numeric(dat) * 99 + 1]
        ),
        text = as.character(dat2),
        textposition = 'auto',
        hoverinfo = 'text'
        
    ) %>%
        layout(
            title = paste("hash:", sim),
            yaxis = list(title = 'Normalized values', range = c(0, 1), tickfont = list(size = 18/nr), titlefont = list(size = 20/nr)),
            xaxis = list(title = '', tickangle = -45, tickfont = list(size = 16/nr)),
            margin = list(t= 100, b = 100, l = 100, r = 100),
            showlegend = FALSE
        )
    
    param_plot <- plot_params(param_values, sim)
    
    
    # Combine all plots into a subplot
    subplot(param_plot, plot, nrows = 2, shareX = TRUE, shareY = TRUE,
            margin = 0.01, titleX = TRUE, titleY = TRUE) %>%
        layout(showlegend = FALSE)
    
}

#' Plot parameters for a specific simulation
#'
#' This function generates visualizations of the parameter values for a specific
#' simulation. It uses a gauge-like representation to show the range and actual 
#' values of parameters, providing an intuitive understanding of the simulation 
#' setup.
#'
#' @param BenchMatrix A data frame containing the parameter values for multiple 
#' simulations.
#' @param current_row An integer or a character string specifying the row number 
#' or the simulation identifier for which the parameters are to be visualized.
#' 
#' @return A plot object visualizing the parameter values for the specified 
#' simulation.
#' 
#' @importFrom plotly plot_ly subplot
#' @export
#'
#' @examples
#' \dontrun{
#' benchmarks <- matrix_benchmarks(run_sims(def_param(), set_phosfake()), set_phosfake())
#' plot_params(benchmarks, current_row = 1)
#' }
plot_params <- function(BenchMatrix, current_row = 1) {
    
    # get parameter names
    param_table <- param_table()
    param_names <- rownames(param_table)
    
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
    val[is.na(val)] <- 0
    param_names <- param_names[param_names %in% colnames(BenchMatrix)]
    param_table <- param_table[param_names, ]
    
    # Set the number of columns
    ncols <- 4
    # Calculate number of rows based on number of plots and columns
    nrows <- ceiling(length(param_names) / ncols)
    
    # Create an empty list to store meters
    meters <- list()
    # Adjust spacing factors for left/right and top/bottom margins
    hspacing_factor <- 0.3
    vspacing_factor <- 0.3
    
    # get position of the current row
    if (!is.numeric(current_row)) {
        current_row <- which(rownames(BenchMatrix) == current_row)
    }
    nr <- 1
    
    # Loop to create each meter plot
    for (i in 1:length(param_names)) {
        curr_par <- param_table[i, ]
        val <- as.numeric(BenchMatrix[current_row, param_names[i]])
        val[is.na(val)] <- 0
        row_index <- ceiling(i / ncols)
        col_index <- (i - 1) %% ncols + 1
        x_domain <- c((col_index - 1 + hspacing_factor) / ncols, (col_index - hspacing_factor) / ncols)
        y_domain <- c(1 - (row_index - vspacing_factor) / nrows, 1 - ((row_index - 1 + vspacing_factor) / nrows))
        # Scale to row position
        # x_domain <- x_domain/2 + 0.5
        y_domain <- y_domain/2 + 0.5
        #y_domain <- y_domain/(nrow(BenchMatrix)*1.1) + (current_row-1)*1.05 / nrow(BenchMatrix)
        fig1 <- plot_ly(
            domain = list(x = x_domain, y = y_domain),
            value = val,
            title = list(text = param_names[i], align="center", font = list(size = 10/nr)),
            type = "indicator",
            mode = "gauge",
            gauge = list(
                shape = "bullet",
                axis = list(range = list(curr_par$MinValue, curr_par$MaxValue),
                            tickmode = "auto",
                            ticks = "inside",  # Position ticks inside
                            tickfont = list(size = 10/nr, color = "black")  # Adjust tickfont size and color for better visibility
                ),
                bar = list(
                    color = colorpanel(100, "#33CC33", "#999966", "#CC3333")[(val - curr_par$MinValue) / (curr_par$MaxValue - curr_par$MinValue) * 100 + 1],
                    thickness = 0.8  # Adjust this value to make the gauge bar thicker
                )
            )
        ) 
        meters[[i]] <- fig1
    }
    # Combine all meter plots
    subplot(meters, nrows = nrows, margin = 0.01)
}


