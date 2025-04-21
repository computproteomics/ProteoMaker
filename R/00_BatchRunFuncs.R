################################################################################
#                       QC batch test of the parameters                        #
################################################################################
#
#' @keywords internal
#' @importFrom stats rnorm runif quantile sd na.omit medpolish setNames
#' @importFrom graphics barplot legend lines par
#' @importFrom utils read.csv
#' @importFrom stringr str_count
NULL


#' Set paths and general configuration for ProteoMaker
#'
#' This function sets the paths and general configuration parameters needed for
#' running simulations with the ProteoMaker package. It creates the necessary output
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
#' @param runStatTests A logical value indicating whether to run statistical tests. Default is `TRUE`.
#'
#'
#' @return A list containing the configuration settings:
#' \describe{
#'   \item{fastaFilePath}{The path to the FASTA files.}
#'   \item{resultFilePath}{The path to the result files.}
#'   \item{cores}{The number of cores used.}
#'   \item{clusterType}{The type of cluster used.}
#'   \item{calcAllBenchmarks}{Whether to calculate all benchmarks.}
#'   \item{runStatTests}{Whether to run statistical tests.}
#' }
#'
#' @export
#'
#' @examples
#' config <- set_proteomaker()
#' config <- set_proteomaker(fastaFilePath = "CustomProteomes",
#' resultFilePath = "Results", cores = 4, clusterType = "PSOCK")
set_proteomaker <- function(fastaFilePath = system.file("Proteomes", package = "ProteoMaker"),
                            resultFilePath = "SimulatedDatasets",
                            cores = 2, clusterType = "FORK",
                            runStatTests = TRUE, calcAllBenchmarks = TRUE) {
  dir.create(resultFilePath, showWarnings = FALSE)
  return(list(fastaFilePath = fastaFilePath, resultFilePath = resultFilePath,
              cores = cores, clusterType = clusterType,
              runStatTests = runStatTests, calcAllBenchmarks = calcAllBenchmarks))
}


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
#' @importFrom yaml yaml.load_file
#' @importFrom Hmisc list.tree
#'
#' @export
#'
#' @examples
#' # Read YAML file from inst folder
#' yaml_path <- system.file("config", "params.yaml", package = "ProteoMaker")
#' params <- def_param(yaml_path)
def_param <- function(yaml_file=NULL) {
  if (is.null(yaml_file)) {
    yaml_file <- system.file("config", "parameters.yaml", package = "ProteoMaker")
  }
  # Read the YAML file
  params <- yaml::yaml.load_file(yaml_file)$params

  # Convert NA values from strings to real NA
  for (l in names(params)) {
    params[[l]]$class <- params[[l]]$choices <- NULL
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
    value <- param_info$value

    if (!is.null(category)) {
      if (!is.null(value)) {
        Param[[category]][[param_name]] <- value
      } else {
        Param[[category]][[param_name]] <- default_value
      }
    }
  }

  # Provide a nice printout for each category
  base::message("--------------------\nGround truth generation parameters:")
  Hmisc::list.tree(Param$paramGroundTruth, maxcomp = 100, maxlen = 100)
  base::message("--------------------\nProteoform abundance parameters:")
  Hmisc::list.tree(Param$paramProteoformAb, maxcomp = 100, maxlen = 100)
  base::message("--------------------\nDigestion parameters:")
  Hmisc::list.tree(Param$paramDigest, maxcomp = 100, maxlen = 100)
  base::message("--------------------\nMSRun parameters:")
  Hmisc::list.tree(Param$paramMSRun, maxcomp = 100, maxlen = 100)
  base::message("--------------------\nData analysis parameters:")
  Hmisc::list.tree(Param$paramDataAnalysis, maxcomp = 100, maxlen = 100)
  base::message("--------------------")
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
#' generated using \code{def_param}.
#' @param Config A list containing configuration settings, such as file paths
#' and computational settings, typically generated using \code{set_proteomaker}.
#' @param overwrite A logical value indicating whether to overwrite existing
#' simulation results. Default is \code{FALSE}.
#'
#' @return A list of results, where each element contains simulation results
#' and benchmarks for a specific parameter combination.
#'
#' @importFrom digest digest
#' @export
#'
#' @examples
#' params <- def_param()
#' config <- set_proteomaker()
#' results <- run_sims(params, config)
run_sims <- function(Parameters, Config, overwrite = FALSE) {
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
  message("Total number of simulations to run: ", totalbench)
  benchcounter <- 0

  # Gather always benchmarking data
  allBs <- list()
  for (hh in 1:length(listtogroundtruth)) {
    ## Running ground thruth generations
    # check whether file with correct parameters exists
    tParam <- listtogroundtruth[[hh]]
    # turn into list
    tParam <- lapply(tParam, function(x) x)
    Param <- "none"
    groundTruth <- NULL
    md5 <- digest::digest(tParam, algo = "md5")
    filename <- paste0(Config$resultFilePath, "/outputGroundTruth_", md5, ".RData")
    if (!overwrite & file.exists(filename)) {
      load(filename)
    } else {
      Param <- tParam
      # Add path to fasta file
      ttParam <- Param
      ttParam$PathToFasta <- paste0(Config$fastaFilePath, ifelse(Config$fastaFilePath == "", "", "/"), Param$PathToFasta)
      groundTruth <- samplePreparation(parameters = ttParam)
      if (!is.null(groundTruth)) {
        save(groundTruth, Param, file = filename)
      }
    }
    gtParam <- Param
    ## quantitative proteoform abundances
    for (ii in 1:length(listtoproteoformab)) {
      Param <- "none"
      proteoformAb <- NULL
      # create combined parameterfile
      tParam <- c(gtParam, listtoproteoformab[[ii]])
      # hash code to represent parameter configuration
      md5 <- digest::digest(tParam, algo = "md5")
      tParam <- c(tParam,
                  list(QuantColnames = paste0("C_",
                                              rep(1:tParam$NumCond, each = tParam$NumReps),
                                              "_R_",
                                              rep(1:tParam$NumReps, tParam$NumCond))))
      filename <- paste0(Config$resultFilePath, "/outputProteoformAb_", md5, ".RData")
      if (!overwrite & file.exists(filename)) {
        load(filename)
      } else {
        Param <- tParam
        proteoformAb <- addProteoformAbundance(proteoforms = groundTruth, parameters = Param)
        save(proteoformAb, Param, file = filename)
      }
      pfParam <- Param
      gc()

      ### Digestion
      for (jj in 1:length(listtodigestion)) {
        Param <- "none"
        BeforeMS <- NULL
        tParam <- c(pfParam, listtodigestion[[jj]])
        md5 <- digest::digest(tParam, algo = "md5")
        filename <- paste0(Config$resultFilePath, "/outputDigest_", md5, ".RData")
        if (!overwrite & file.exists(filename)) {
          load(filename)
        } else {
          Param <- tParam
          peptable <- digestGroundTruth(
            proteoforms = proteoformAb,
            parameters = c(
              Param, list(Cores = Config$cores,
                          ClusterType = Config$clusterType)))
          peptable <- digestionProductSummarization(peptides = peptable,
                                                    parameters = c(Param, list(Cores = Config$cores, ClusterType = Config$clusterType)))
          BeforeMS <- filterDigestedProt(DigestedProt = peptable,
                                         parameters = Param)
          save(Param, BeforeMS, file = filename)
        }
        dgParam <- Param
        gc()

        ### MS run
        for (kk in 1:length(listtomsrun)) {
          Param <- "none"
          AfterMSRun <- NULL
          tParam <- c(dgParam, listtomsrun[[kk]])
          md5 <- digest::digest(tParam, algo = "md5")
          filename <- paste0(Config$resultFilePath, "/outputMSRun_", md5, ".RData")
          if (!overwrite & file.exists(filename)) {
            load(filename)
          } else {
            Param <- tParam
            AfterMSRun <- vector(mode = "list")
            for (i in which(sapply(BeforeMS, length) > 0)) {
              AfterMSRun[[length(AfterMSRun) + 1]] <- MSRunSim(
                Digested = BeforeMS[[i]],
                parameters = c(Param,
                               list(Cores = Config$cores,
                                    ClusterType = Config$clusterType)))
            }
            names(AfterMSRun) <- names(BeforeMS)[which(sapply(BeforeMS, length) > 0)]
            save(Param, AfterMSRun, file = filename)
          }
          msParam <- Param
          gc()

          ### Protein abundance
          for (ll in 1:length(listtodatanalysis)) {
            Param <- "none"
            tParam <- c(msParam, listtodatanalysis[[ll]])
            Benchmarks <- NULL
            md5 <- digest::digest(tParam, algo = "md5")
            filename <- paste0(Config$resultFilePath, "/outputDataAnalysis_", md5, ".RData")
            if (!overwrite & file.exists(filename)) {
              load(filename)
            } else if (Config$runStatTests) {
              # counter
              Param <- tParam
              Prots <- proteinSummarisation(peptable = AfterMSRun$NonEnriched,
                                            parameters = c(Param, list(Cores = Config$cores, ClusterType = Config$clusterType)))
              # Don't accept anything below 100 proteins
              if (nrow(Prots) > 99) {
                # Filter for having at least 1 actual value per protein group and peptide
                Prots <- Prots[rowSums(is.na(Prots[, Param$QuantColnames])) < length(Param$QuantColnames), ]
                allPeps <- as.data.frame(do.call("rbind", AfterMSRun))
                allPeps <- allPeps[rowSums(is.na(allPeps[, Param$QuantColnames])) < length(Param$QuantColnames), ]
                rownames(allPeps) <- paste0("pep", 1:nrow(allPeps))
                Stats <- runPolySTest(Prots, Param, refCond = 1, onlyLIMMA = F, cores = Config$cores)
                # much faster with only LIMMA tests
                StatsPep <- runPolySTest(allPeps, Param, refCond = 1, onlyLIMMA = T, cores = Config$cores)
                #Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
                save(Param, Stats, StatsPep, Benchmarks, file = filename)
              } else {
                message("Too few proteins or no statistical tests requested.
                                        Skipping this protein summarization, statistical testing and
                                        benchmark calculation.")
                Benchmarks <- Stats <- StatsPep <- NULL
              }

              if (Config$calcAllBenchmarks & !is.null(Stats)) {
                Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
                save(Param, Stats, StatsPep, Benchmarks, file = filename)
              }
            } else {
              Param <- tParam
            }
            allBs[[md5]] <- list(Benchmarks = Benchmarks, Param = Param)
          }
        }
        gc()
      }
    }
  }
  return(allBs)
  message("###### Finished data set generation")
}



#' Retrieve intermediate outputs for a simulation
#'
#' This function retrieves intermediate outputs for a simulation based on the
#' provided parameters and configuration. It checks for the existence of the
#' output file corresponding to a specific simulation stage and loads the
#' data if available.
#'
#' @param Param A list of parameters used in the simulation, typically a subset
#' of those generated by functions like \code{def_param}.
#' @param Config A list of configuration settings, including file paths and
#' computational settings, typically generated using \code{set_proteomaker}.
#' @param stage A character string indicating the stage of the simulation for
#' which the output is requested. Default is \code{"DataAnalysis"}. Values are \code{"GroundTruth"},
#' \code{"ProteoformAb"}, \code{"Digest"}, \code{"MSRun"}, and \code{"DataAnalysis"}.
#'
#' @return A list containing the loaded data objects, including \code{Stats},
#' \code{StatsPep}, and \code{Benchmarks}, if they exist. If no file is found,
#' the function returns \code{NULL}.
#'
#' @importFrom digest digest
#'
#' @export
#'
#' @examples
#' config <- set_proteomaker(resultFilePath = tempdir())
#' Param <- def_param()
#' Param$paramGroundTruth$NumReps <- 5
#' benchmarks <- run_sims(Param, config)
#' result <- get_simulation(benchmarks[[1]]$Param, config, stage="MSRun")
get_simulation <- function(Param, Config, stage="DataAnalysis") {

  # Check for valid stage name
  if (!(stage %in% c("GroundTruth", "ProteoformAb", "Digest", "MSRun", "DataAnalysis"))) {
    stop("Invalid stage name. Please provide a valid stage name.")
  }

  # Get all parameter names
  param_names <- param_table()

  # Reduce the parameter set to the relevant stage
  param_names <- param_names[1:max(which(param_names$Group == paste0("param", stage))), ]
  # max of which is the last element of the vector
  max_pname <- max(which(names(Param) %in% rownames(param_names)))
  tParam  <- Param[1:max_pname]
  # print(tParam)

  # check whether file with correct parameters exists
  md5 <- digest::digest(as.list(tParam), algo = "md5")

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
    message("No simulation found with these parameters.")
  }
}

#' Gather parameters and benchmarks from all available runs
#'
#' This function retrieves the parameters and benchmarks (if available) from all available runs
#' and compiles them into a list for further analysis.
#'
#' @param Config A list of configuration settings, including file paths and
#' computational settings, typically generated using \code{set_proteomaker}.
#' @param stage A character string indicating the stage of the simulation for
#' which the output is requested. Default is \code{"DataAnalysis"}. Values are \code{"GroundTruth"},
#' \code{"ProteoformAb"}, \code{"Digest"}, \code{"MSRun"}, and \code{"DataAnalysis"}. Default is \code{"DataAnalysis"}.
#'
#' @return A list containing the parameters and benchmarks from all available runs.
#'
#' @export
#'
#' @examples
#' config <- set_proteomaker(resultFilePath = tempdir())
#' all_results <- gather_all_sims(config)
#'
gather_all_sims <- function(Config, stage = "DataAnalysis") {
  # Get all files in the result directory
  all_files <- list.files(Config$resultFilePath, full.names = TRUE)

  # Filter for RData files
  rdata_files <- all_files[grep(".RData", all_files)]
  rdata_files <- rdata_files[grep(paste0("output", stage), rdata_files)]

  # Initialize an empty list to store the results
  all_results <- list()

  # Iterate over each RData file
  for (file in rdata_files) {
    # Load the RData file
    load(file)

    # get hash
    hash <- gsub(paste0(Config$resultFilePath, "/output", stage, "_"), "", file)
    hash <- gsub(".RData", "", hash)

    # Get the parameter and benchmark objects
    param <- get("Param")

    if(exists("Benchmarks")) {
      benchmarks <- get("Benchmarks")

      # Create a list of the objects
      result_list <- list(Param = param, Benchmarks = benchmarks)
    } else {
      result_list <- list(Param = param)
    }
    # Append the list to the results list
    all_results[[hash]] <- result_list
  }

  # Return the list of results
  return(all_results)
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
#' computational settings, typically generated using \code{set_proteomaker}.
#'
#' @return A data frame with rows representing each simulation and columns
#' representing benchmark values and parameter settings.
#'
#' @export
#'
#' @examples
#' conf <- set_proteomaker(resultFilePath = tempdir())
#' results <- run_sims(def_param(), conf)
#' benchmark_matrix <- matrix_benchmarks(results, conf)
#'
matrix_benchmarks <- function(allBs, Config) {
  # extracting all benchmarks (sometimes there are more or less per run)
  t_allbnames <- NULL
  for (i in names(allBs)) {
    t_allbnames <- c(t_allbnames, names(unlist(allBs[[i]][[1]]$globalBMs)))
  }
  benchNames <- unique(t_allbnames)
  parNames <- unique(names(allBs[[1]]$Param))
  BenchMatrix <- data.frame(matrix(NA, ncol = length(benchNames) + length(parNames), nrow = length(allBs)))
  colnames(BenchMatrix) <- c(benchNames, parNames)
  rownames(BenchMatrix) <- names(allBs)

  # writing all results and parameters into matrix
  for (i in names(allBs)) {
    tglob <- unlist(allBs[[i]]$Benchmarks$globalBMs)
    BenchMatrix[i, names(tglob)] <- tglob
    tpar <- allBs[[i]]$Param
    tpar <- lapply(tpar, function(x) paste0(unlist(x), collapse = ";"))
    BenchMatrix[i, names(tpar)] <- sapply(tpar, function(x) ifelse(length(x) > 1, paste0(x, collapse = "_"), x))
  }
  BenchMatrix
}

#' Visualize benchmarks vs 1 or two parameters
#'
#' This function visualizes the benchmarks from a matrix of benchmark results. It creates
#' bar plots for the specified benchmarks and parameters, allowing for the
#' comparison of benchmark values across different parameter settings.
#'
#' @param benchmatrix A data frame containing benchmark results, typically
#' generated by the \code{matrix_benchmarks} function.
#' @param benchmarks A character vector specifying the names of the benchmarks
#' to be visualized. If \code{NULL}, all benchmarks will be plotted.
#' @param ref_par A character vector specifying the names of the parameters to
#' be used for visualization. It can contain one or two parameters.
#' @param cols A character vector specifying the colors to be used for the
#' visualization. If \code{NULL}, default colors will be used.
#'
#' @return A bar plot or a colored map visualizing the specified benchmarks against the
#' reference parameters.
#'
#' @examples
#' conf <- set_proteomaker(resultFilePath = tempdir())
#' results <- run_sims(def_param(), conf)
#' benchmark_matrix <- matrix_benchmarks(results, conf)
#' visualize_benchmarks(benchmark_matrix, ref_par = "NumReps")
#'
#' @importFrom colorspace qualitative_hcl
#' @importFrom lattice levelplot
#'
#' @export
visualize_benchmarks <- function(benchmatrix,
                                 benchmarks = NULL,
                                 ref_par = "WrongIDs",
                                 cols = NULL) {
  # Require at least one row in benchmatrix
  if (nrow(benchmatrix) < 1) {
    stop("No data available for visualization.")
  }

  ref_par <- make.names(ref_par)
  colnames(benchmatrix) <- make.names(colnames(benchmatrix))

  titles <- c(
    numPeptides = "I: #Total peptides (mod+unmod)",
    numProteins = "II: #Proteins (from peptides)",
    propUniquePep = "III: % Unique peptides (prot-specific)",
    uniqueStrippedPep = "IV: #Stripped sequences",
    percMissingPep = "V: % Missing peptide values",
    aucDiffRegPeptides.FDR_limma.2.vs.1.AUC = "VI: Peptide AUC (truth vs FDR)",
    tprPep0.01.FDR_limma.2.vs.1.TPR = "VII: TPR peptide (FDR < 0.01)",
    tprPep0.05.FDR_limma.2.vs.1.TPR = "VIII: TPR peptide (FDR < 0.05)",
    tFDRPep0.01.FDR_limma.2.vs.1.tFDR = "IX: True FDR (peptides, 0.01)",
    tFDRPep0.05.FDR_limma.2.vs.1.tFDR = "X: True FDR (peptides, 0.05)",
    propMisCleavedPeps.0 = "XI: No miscleavage fraction",
    dynRangePep = "XII: Dynamic range (peptides)",
    sumSquareDiffFCPep = "XIII: Fold-change error (peptides)",
    sdWithinRepsPep = "XIV: Replicate SD (peptides)",
    skewnessPeps = "XV: Skewness (peptides)",
    kurtosisPeps = "XVI: Kurtosis (peptides)",
    sdPeps = "XVII: Overall SD (peptides)",
    numQuantProtGroups = "XVIII: #Quantified proteins",
    propUniqueProts = "XIX: % Unique proteins",
    percMissingProt = "XX: % Missing protein values",
    meanPepPerProt = "XXI: Peptides per protein",
    aucDiffRegProteins.FDR_PolySTest.2.vs.1.AUC = "XXII: Protein AUC (truth vs FDR)",
    tprProt0.01.FDR_PolySTest.2.vs.1.TPR = "XXIII: TPR protein (FDR < 0.01)",
    tprProt0.05.FDR_PolySTest.2.vs.1.TPR = "XXIV: TPR protein (FDR < 0.05)",
    tFDRProt0.01.FDR_PolySTest.2.vs.1.tFDR = "XXV: True FDR (proteins, 0.01)",
    tFDRProt0.05.FDR_PolySTest.2.vs.1.tFDR = "XXVI: True FDR (proteins, 0.05)",
    sumSquareDiffFCProt = "XXVII: Fold-change error (proteins)",
    sdWithinRepsProt = "XXVIII: Replicate SD (proteins)",
    propMisCleavedProts = "XXIX: % Miscleaved proteins",
    propDiffRegWrongIDProt0.01.FDR_PolySTest.2.vs.1 = "XXX: % Wrong-ID (proteins, 0.01)",
    propDiffRegWrongIDProt0.05.FDR_PolySTest.2.vs.1 = "XXXI: % Wrong-ID (proteins, 0.05)",
    skewnessProts = "XXXII: Skewness (proteins)",
    kurtosisProts = "XXXIII: Kurtosis (proteins)",
    sdProts = "XXXIV: Overall SD (proteins)",
    numProteoforms = "XXXV: #Proteoforms",
    meanProteoformsPerProt = "XXXVI: Proteoforms per protein",
    numModPeptides = "XXXVII: #Modified peptides",
    propModAndUnmodPep = "XXXVIII: % Modified with unmodified match",
    aucDiffRegAdjModPep = "XXXIX: AUC (adj. mod. peptides)",
    tprAdjModPep0.01 = "XL: TPR (adj. mod. peptides, 0.01)",
    tprAdjModPep0.05 = "XLI: TPR (adj. mod. peptides, 0.05)",
    tFDRAdjModPep0.01 = "XLII: True FDR (mod. peptides, 0.01)",
    tFDRAdjModPep0.05 = "XLIII: True FDR (mod. peptides, 0.05)",
    propDiffRegPepWrong0.01.FDR_PolySTest.2.vs.1 = "XLIV: % Wrongly significant (mod. peptides, 0.01)",
    propDiffRegPepWrong0.05.FDR_PolySTest.2.vs.1 = "XLV: % Wrongly significant (mod. peptides, 0.05)",
    percOverlapModPepProt = "XLVI: % Mod peptides with protein quant",
    sumSquareDiffFCModPep = "XLVII: Fold-change error (mod. peptides)"
  )

  titles_params <- c(
    NumCond = "Number of conditions",
    NumReps = "Number of replicates",
    PathToFasta = "Path to FASTA file",
    PathToProteinList = "Path to protein list (optional)",
    PercExpressedProt = "Fraction of expressed proteins",
    FracModProt = "Fraction of proteins to be modified",
    PropModPerProt = "Max #Modifications per protein",
    PTMMultipleLambda = "Poisson lambda for multiple PTMs",
    RemoveNonModFormFrac = "Fraction without non-modified form",
    PTMTypes = "PTM types",
    PTMTypesDistr = "PTM type distribution",
    PTMTypesMass = "PTM mass shifts",
    ModifiableResidues = "Modifiable residues",
    ModifiableResiduesDistr = "Modifiable residue distribution",
    QuantNoise = "Quantification noise (SD)",
    DiffRegFrac = "Fraction of differentially regulated",
    DiffRegMax = "Max fold-change (log2)",
    UserInputFoldChanges_NumRegProteoforms = "#User-specified regulated proteoforms",
    UserInputFoldChanges_RegulationFC = "User-defined fold change",
    AbsoluteQuanMean = "Mean absolute quantity",
    AbsoluteQuanSD = "SD of absolute quantity",
    ThreshNAQuantileProt = "NA removal quantile threshold",
    Enzyme = "Digestion enzyme",
    PropMissedCleavages = "Proportion with missed cleavages",
    MaxNumMissedCleavages = "Max missed cleavages",
    PepMinLength = "Min peptide length",
    PepMaxLength = "Max peptide length",
    LeastAbundantLoss = "Fraction of least abundant removed",
    EnrichPTM = "Enriched PTM type",
    EnrichmentLoss = "Enrichment loss",
    EnrichmentEfficiency = "Enrichment efficiency",
    EnrichmentQuantDiff = "Quant difference due to enrichment",
    EnrichmentNoise = "Enrichment noise",
    ModificationLoss = "Loss of modified peptides",
    PercDetectability = "% Peptides detected (model)",
    PercDetectedVal = "% Detected intensities",
    WeightDetectVal = "Intensity-dependence of detection",
    MSNoise = "Instrumental noise",
    WrongIDs = "Wrong identification rate",
    WrongLocalizations = "Wrong localization rate",
    MaxNAPerPep = "Max NAs per peptide",
    ProtSummarization = "Protein summarization method",
    MinUniquePep = "Min unique peptides per protein",
    StatPaired = "Paired testing enabled"
  )

  # Setting fixed color code
  dark_colors <- colorspace::qualitative_hcl(n = length(titles), palette = "Dark 3")
  M <- matrix(1:100, nrow = 10, byrow = TRUE)
  M <- M[M <= length(dark_colors) & M > 0]
  dark_colors <- dark_colors[as.vector(M)]
  names(dark_colors) <- names(titles)
  # Setting fixed character for pch
  pch.vals <- rep(c(1, 16, 17, 15, 3, 4, 8), length.out = length(titles))
  names(pch.vals) <- names(titles)

  # Allow at most two parameters
  if (length(ref_par) > 2) {
    stop("Please provide at most two parameters for visualization.")
  }
  # Assign titles to parameters
  for (i in ref_par) {
    if (!(i %in% names(titles_params))) {
      stop(paste("Parameter name", i, "is not correct."))
    }
  }
  params <- titles_params[ref_par]

  # Setting benchmarking metrics to be plotted
  if(length(benchmarks) == 0) {
    titles <- titles[1:12]
  } else if (is.character(benchmarks)) {
    # Check if the provided benchmarks are valid
    benchmarks <- make.names(benchmarks)
    for (i in benchmarks) {
      if (!(i %in% names(titles))) {
        stop(paste("Benchmark name", i, "is not correct."))
      }
    }
    titles <- titles[benchmarks]
  } else  if (is.numeric(benchmarks)) {
    # check for range of indices
    if (max(benchmarks) > length(titles)) {
      stop(paste("Too many benchmarks selected. Please select a maximum of", length(titles)))
    }
    titles <- titles[benchmarks]

    # check whether in benchmatrix
    if (length(titles) > 0) {
      for(i in names(titles)) {
        if (!(i %in% names(benchmatrix))) {
          warning(paste("Benchmark name", i, "is not available as benchmark. Will be removed"))
          titles <- titles[-which(names(titles) == i)]
        }
      }
    }
  } else {
    stop("Please provide a character vector of benchmark names or numeric vector of indices.")
  }
  benchmarks <- names(titles)
  dark_colors <- dark_colors[benchmarks]
  pch.vals <- pch.vals[benchmarks]

  n_plots <- length(titles)
  # Adapt to plot number with 5 columns max and max. 25 plots per page
  if (n_plots > 16) {
    plotmfrow <- c(4, 4)
  } else {
    plotmfrow <- c(4, ceiling(n_plots / 4))
  }
  par(mfrow=plotmfrow, cex.main=0.9, cex.lab = 0.75, cex.axis = 0.7,
      mgp = c(1.5, 0.3, 0), mar=c(3, 3, 1.5, 1.5), xpd = TRUE, font.main = 2)

  # Generate a dark qualitative color palette
  if (!is.null(cols)) {
    dark_colors <- rep(cols, length.out = n_plots)
  }

  # Loop through all metrics
  for (i in names(titles)) {
    col <- dark_colors[i]
    pch.use <- pch.vals[i]
    # Single plots when having one parameter
    if (length(ref_par) == 1) {
      if (all(is.finite(range(benchmatrix[, i])))) {
        plot(benchmatrix[, ref_par], benchmatrix[, i],
             xlab = params, ylab = i,
             pch = pch.use, col = col, cex = 1.5, cex.lab = 1, cex.axis = 1)

        abline(h = pretty(benchmatrix[, i]), col = "gray90", lty = "dotted")
        abline(v = pretty(benchmatrix[, ref_par]), col = "gray90", lty = "dotted")

        title(main = paste0(titles[i]),
              col.main = col, font.main = 2)
      } else {
        plot(benchmatrix[, ref_par], rep(0, length(benchmatrix[, i])), type = "n",
             xlab = params, ylab = i,
             main = paste0(titles[i]),
             col.main = col, font.main = 2)
      }
    } else {
      # colormaps showing the benchmarks in a 2-dim plot as colors
      # Create matrix from benchmatrix ref_par columns with values from i
      # and plot it

      if (length(ref_par) == 2) {
        x_vals <- sort(unique(benchmatrix[, ref_par[1]]))
        y_vals <- sort(unique(benchmatrix[, ref_par[2]]))

        z_mat <- matrix(NA, nrow = length(x_vals), ncol = length(y_vals))

        for (j in 1:nrow(benchmatrix)) {
          x_idx <- which(x_vals == benchmatrix[j, ref_par[1]])
          y_idx <- which(y_vals == benchmatrix[j, ref_par[2]])
          z_mat[x_idx, y_idx] <- benchmatrix[j, i]
        }

        # Define a color ramp
        n_colors <- 100
        col_ramp <- colorRampPalette(c("white", col))
        image_colors <- col_ramp(n_colors)
        zlim <- range(z_mat, na.rm = TRUE)
        if(any(!is.finite(zlim))) {
          zlim <- c(0, 0)
        }

        # Plot image
        image(x = x_vals, y = y_vals, z = z_mat,
              col = image_colors, zlim = zlim,
              xlab = params[1], ylab = params[2], axes = TRUE)


        # Colorized title
        title(main = paste0(titles[i]),
              col.main = col, font.main = 2)

        # Expand singleton axes to avoid collapsed plots
        if (length(x_vals) == 1) {
          x_vals <- x_vals + c(-0.5, 0.5)
        }
        if (length(y_vals) == 1) {
          y_vals <- y_vals + c(-0.5, 0.5)
        }


        # Add colorbar along right side
        bar_x <- max(x_vals) -  diff(range(x_vals)) * 0.5
        bar_w <- diff(range(x_vals)) * 0.03
        rect_x <- c(bar_x, bar_x + bar_w)

        # Color bar segments
        color_steps <- length(image_colors)
        bar_y_vals <- seq(min(y_vals, na.rm=T), max(y_vals, na.rm=T), length.out = color_steps + 1)

        for (k in 1:color_steps) {
          rect(rect_x[1], bar_y_vals[k], rect_x[2], bar_y_vals[k+1],
               col = image_colors[k], border = NA)
        }
        # rectangle around all rectangles with white border
        rect(rect_x[1], min(y_vals), rect_x[2], max(y_vals),
             col = NA, border = "#333333", lwd = 0.5)

        # Add legend labels on the colorbar
        legend_labels <- pretty(zlim, n = 3)
        if (diff(range(zlim)) > 0) {
          legend_pos <- approx(zlim, range(bar_y_vals), xout = legend_labels)$y
          text(x = rect_x[2] + 0.02 * diff(range(x_vals)),
               y = legend_pos,
               labels = round(legend_labels, 2),
               cex = 0.6, adj = 0)
        }
      }
    }
  }
  # Reset plotting layout
  par(mfrow = c(1, 1), cex.main = 1.2, cex.lab = 1, cex.axis = 1, xpd = FALSE,
      font.main = 1)

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
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' benchmarks <- matrix_benchmarks(run_sims(def_param(), set_proteomaker()), set_proteomaker())
#' visualize_one_sim(benchmarks, current_row = 1)
visualize_one_sim <- function(BenchMatrix, current_row = 1) {
  # get parameter names
  param_t <- param_table()
  param_names <- rownames(param_t)

  # Visualize roughly
  par(mfrow = c(nrow(BenchMatrix), 1), mar = c(2, 5, 2, 1), xpd = T, oma = c(5, 1, 1, 1))
  # Separate parameters and filter for actual values
  reds <- colnames(BenchMatrix) %in% param_names
  param_values <- BenchMatrix[, reds]
  BenchMatrix <- BenchMatrix[, !reds]
  to_del <- c()
  if (is.null(ncol(BenchMatrix))) {
    message("No benchmarks have been calculated for this simulation.")
    return()
  }
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
  if (is.null(ncol(BenchMatrix))) {
    message("No benchmarks available for this simulation.")
    return()
  }
  for (i in 1:ncol(BenchMatrix)) {
    tt <- BenchMatrix[, i]
    if (is.numeric(tt)) {
      BenchMatrix[, i] <- round(tt, 2)
    }
    if (is.character(tt) || is.factor(tt)) {
      tt <- as.numeric(as.factor(tt))
    }
    tt <- tt - min(tt, 0)
    tt <- tt / max(tt, na.rm = T)
    tt[is.na(tt)] <- 0
    tBenchMatrix[, i] <- unlist(tt)
  }

  # Define color palette
  color_palette <- colorpanel(100, "#AA3333", "#3333AA")

  # Create plots
  plots <- list()
  sim <- current_row
  if (is.numeric(sim))
    sim <- rownames(BenchMatrix)[sim]
  dat <- tBenchMatrix[sim, ]
  dat <- sapply(dat, function(x) {
    if (is.numeric(x)) {
      dat <- round(x, 2)
    }
    if (!is.numeric(x)) {
      x <- 0
    }
    x
  })
  dat2 <- BenchMatrix[sim, ]
  dat2 <- sapply(dat2, function(x) {    if (is.numeric(x)) {
    x <- round(x, 2)
  }
    x
  })
  plot <- plotly::plot_ly(
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
    plotly::layout(
      title = paste("hash:", sim),
      yaxis = list(title = 'Normalized values', range = c(0, 1), tickfont = list(size = 18/nr), titlefont = list(size = 20/nr)),
      xaxis = list(title = '', tickangle = -45, tickfont = list(size = 16/nr)),
      margin = list(t = 100, b = 100, l = 100, r = 100),
      showlegend = FALSE
    )
  param_plot <- plot_params(param_values, sim)


  # Combine all plots into a subplot
  subplot(param_plot, plot, nrows = 2, shareX = TRUE, shareY = TRUE,
          margin = 0.01, titleX = TRUE, titleY = TRUE) %>%
    plotly::layout(showlegend = FALSE)

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
#' benchmarks <- matrix_benchmarks(run_sims(def_param(), set_proteomaker()), set_proteomaker())
#' plot_params(benchmarks, current_row = 1)
#' }
plot_params <- function(BenchMatrix, current_row = 1) {

  # get parameter names
  param_t <- param_table()
  param_names <- rownames(param_t)

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
  param_t <- param_t[param_names, ]

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
    curr_par <- param_t[i, ]
    curr_par$MaxValue <- unlist(curr_par$MaxValue)
    curr_par$MinValue <- unlist(curr_par$MinValue)
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


