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
#' config <- set_proteomaker(
#'   fastaFilePath = "CustomProteomes",
#'   resultFilePath = "Results", cores = 4, clusterType = "PSOCK"
#' )
set_proteomaker <- function(fastaFilePath = system.file("Proteomes", package = "ProteoMaker"),
                            resultFilePath = "SimulatedDatasets",
                            cores = 2, clusterType = "PSOCK",
                            runStatTests = TRUE, calcAllBenchmarks = TRUE) {
  dir.create(resultFilePath, showWarnings = FALSE)
  return(list(
    fastaFilePath = fastaFilePath, resultFilePath = resultFilePath,
    cores = cores, clusterType = clusterType,
    runStatTests = runStatTests, calcAllBenchmarks = calcAllBenchmarks
  ))
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
def_param <- function(yaml_file = NULL) {
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
  combinations <- expand.grid(params, stringsAsFactors = FALSE)
  # seems to be the only way to get to a list of lists
  out_comb <- list()
  for (i in seq_len(nrow(combinations))) {
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
  # Ensuring ProteoMaker version in parameters (for hashing)
  Parameters$paramGroundTruth$ProteoMakerVersion <- packageVersion("ProteoMaker")

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

  # Gather always benchmarking data
  allBs <- list()
  for (hh in seq_len(length(listtogroundtruth))) {
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
      tParam <- c(
        tParam,
        list(QuantColnames = paste0(
          "C_",
          rep(1:tParam$NumCond, each = tParam$NumReps),
          "_R_",
          rep(1:tParam$NumReps, tParam$NumCond)
        ))
      )
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
        SearchIndex <- NULL
        tParam <- c(pfParam, listtodigestion[[jj]])
        md5 <- digest::digest(tParam, algo = "md5")
        filename <- paste0(Config$resultFilePath, "/outputDigest_", md5, ".RData")
        if (!overwrite & file.exists(filename)) {
          load(filename)
          # Backward compatibility: older digest files may not contain SearchIndex.
          if (!exists("SearchIndex") || is.null(SearchIndex)) {
            idxParam <- Param
            idxParam$PathToFasta <- paste0(Config$fastaFilePath, ifelse(Config$fastaFilePath == "", "", "/"), Param$PathToFasta)
            idxParam$Cores <- Config$cores
            idxParam$ClusterType <- Config$clusterType
            SearchIndex <- buildSearchIndexFromFasta(parameters = idxParam)
          }
        } else {
          Param <- tParam
          # Build full FASTA search index for later reuse (per digestion parameter set)
          idxParam <- Param
          idxParam$PathToFasta <- paste0(Config$fastaFilePath, ifelse(Config$fastaFilePath == "", "", "/"), Param$PathToFasta)
          idxParam$Cores <- Config$cores
          idxParam$ClusterType <- Config$clusterType
          SearchIndex <- buildSearchIndexFromFasta(parameters = idxParam)
          peptable <- digestGroundTruth(
            proteoforms = proteoformAb,
            parameters = c(
              Param, list(
                Cores = Config$cores,
                ClusterType = Config$clusterType
              )
            ),
            searchIndex = SearchIndex
          )
          peptable <- digestionProductSummarization(
            peptides = peptable,
            parameters = c(Param, list(Cores = Config$cores, ClusterType = Config$clusterType))
          )
          BeforeMS <- filterDigestedProt(
            DigestedProt = peptable,
            parameters = Param
          )
          save(Param, BeforeMS, SearchIndex, file = filename)
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
                parameters = c(
                  Param,
                  list(
                    Cores = Config$cores,
                    ClusterType = Config$clusterType
                  )
                ),
                searchIndex = SearchIndex
              )
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
            Benchmarks <- Occupancies <- NULL
            md5 <- digest::digest(tParam, algo = "md5")
            filename <- paste0(Config$resultFilePath, "/outputDataAnalysis_", md5, ".RData")
            if (!overwrite & file.exists(filename)) {
              load(filename)
            } else if (Config$runStatTests) {
              # counter
              Param <- tParam
              Prots <- proteinSummarisation(
                peptable = AfterMSRun$NonEnriched,
                parameters = c(Param, list(Cores = Config$cores, ClusterType = Config$clusterType))
              )
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
                # get occupancies for all PTM types (for later benchmarking)
                Occupancies <- calcPTMOccupancy(allPeps, Param)
                save(Param, Stats, StatsPep, Occupancies, Benchmarks, file = filename)
              } else {
                message("Too few proteins or no statistical tests requested.
                                        Skipping this protein summarization, statistical testing and
                                        benchmark calculation.")
                Benchmarks <- Stats <- StatsPep <- NULL
              }

              if (Config$calcAllBenchmarks & !is.null(Stats)) {
                Benchmarks <- calcBenchmarks(Stats, StatsPep, Param)
                save(Param, Stats, StatsPep, Occupancies, Benchmarks, file = filename)
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
#' result <- get_simulation(benchmarks[[1]]$Param, config, stage = "MSRun")
get_simulation <- function(Param, Config, stage = "DataAnalysis") {
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
  tParam <- Param[1:max_pname]
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

#' Retrieve hashes for intermediate simulation stages
#'
#' This function computes the MD5 hashes used to identify intermediate stage
#' outputs for a given parameter set. The hashes match the filenames created
#' by \code{run_sims} and used by \code{get_simulation}.
#'
#' @param Param A list of parameters used in the simulation, typically a subset
#' of those generated by functions like \code{def_param}.
#' @param stages A character vector of stages to compute hashes for. Default is
#' all stages: \code{"GroundTruth"}, \code{"ProteoformAb"}, \code{"Digest"},
#' \code{"MSRun"}, and \code{"DataAnalysis"}.
#'
#' @return A named character vector of MD5 hashes, one per stage.
#'
#' @importFrom digest digest
#'
#' @export
#'
#' @examples
#' Param <- def_param()
#' hashes <- get_stage_hashes(Param)
#' hashes["MSRun"]
get_stage_hashes <- function(Param, stages = c("GroundTruth", "ProteoformAb", "Digest", "MSRun", "DataAnalysis")) {
  if (any(!stages %in% c("GroundTruth", "ProteoformAb", "Digest", "MSRun", "DataAnalysis"))) {
    stop("Invalid stage name. Please provide a valid stage name.")
  }

  param_names <- param_table()

  # Accept both nested (def_param) and flat (run_sims) parameter lists
  if (all(c("paramGroundTruth", "paramProteoformAb", "paramDigest", "paramMSRun", "paramDataAnalysis") %in% names(Param))) {
    param_flat <- c(
      Param$paramGroundTruth,
      Param$paramProteoformAb,
      Param$paramDigest,
      Param$paramMSRun,
      Param$paramDataAnalysis
    )
  } else {
    param_flat <- Param
  }

  hashes <- vapply(stages, function(stage) {
    stage_names <- param_names[1:max(which(param_names$Group == paste0("param", stage))), ]
    stage_order <- rownames(stage_names)
    stage_order <- stage_order[stage_order %in% names(param_flat)]
    if (length(stage_order) == 0) {
      stop("No parameters found for stage: ", stage)
    }
    tParam <- param_flat[stage_order]
    digest::digest(as.list(tParam), algo = "md5")
  }, character(1))

  hashes
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

    if (exists("Benchmarks")) {
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
#' to be visualized. If \code{NULL}, all grouped benchmarks will be plotted.
#' @param benchmark_level A character vector specifying benchmark level groups
#' to include. Supported values are \code{"peptidoform"},
#' \code{"protein_group"}, and \code{"ptm_proteoform"}.
#' @param benchmark_category A character vector specifying benchmark purpose
#' groups to include. Supported values are \code{"coverage"},
#' \code{"completeness"}, \code{"quantification_quality"},
#' \code{"differential_performance"}, \code{"artifact_sensitivity"}, and
#' \code{"ptm_adjusted_performance"}.
#' @param ref_par A character vector specifying the names of the parameters to
#' be used for visualization. It can contain one or two parameters.
#' @param fullrange A logical value indicating whether to use the maximal range of
#' set of benchmarks (works only for one `ref_par`)
#' @param cols A character vector specifying the colors to be used for the
#' visualization. If \code{NULL}, default colors will be used.
#' @param errorbar A logical value indicating whether to show mean +/- standard
#' deviations in the
#' visualization. Default is \code{FALSE}. This parameter is not applied to the color maps.
#' @param errorstyle A character value specifying how uncertainty should be shown
#' when \code{errorbar = TRUE}. Use \code{"bar"} for vertical error bars or
#' \code{"area"} for a shaded mean +/- SD band.
#' @param compare_par A character vector specifying the names of an additional parameter for subsetting and then
#' comparison of the benchmark values. This parameter is only applied when `ref_par` contains one parameter.
#' If defulat `NULL`, the full set will be shown.
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
#' @importFrom gplots plotCI
#'
#' @export
visualize_benchmarks <- function(benchmatrix,
                                 benchmarks = NULL,
                                 benchmark_level = NULL,
                                 benchmark_category = NULL,
                                 ref_par = "WrongIDs",
                                 fullrange = FALSE,
                                 cols = NULL,
                                 errorbar = FALSE,
                                 errorstyle = "bar",
                                 compare_par = NULL) {
  if (nrow(benchmatrix) < 1) {
    stop("No data available for visualization.")
  }

  ref_par <- make.names(ref_par)
  if (length(ref_par) > 2) {
    stop("Please provide at most two parameters for visualization.")
  }
  colnames(benchmatrix) <- make.names(colnames(benchmatrix))
  errorstyle <- match.arg(errorstyle, c("bar", "area"))

  bm_meta <- get_bmmeta()
  titles <- stats::setNames(bm_meta$title, bm_meta$benchmark)
  axis_labels <- stats::setNames(bm_meta$axis_label, bm_meta$benchmark)
  ranges <- set_bmranges(titles)
  titles_params <- get_paramtitles()
  level_order <- c("peptidoform", "protein_group", "ptm_proteoform")
  category_order <- c(
    "coverage", "completeness", "quantification_quality",
    "differential_performance", "artifact_sensitivity",
    "ptm_adjusted_performance"
  )
  category_colors <- c(
    coverage = "#1B6CA8",
    completeness = "#2A9D8F",
    quantification_quality = "#E9C46A",
    differential_performance = "#F4A261",
    artifact_sensitivity = "#D95D39",
    ptm_adjusted_performance = "#7A8E3A"
  )

  metric_colors <- colorspace::qualitative_hcl(n = length(titles), palette = "Dark 3")
  color_idx <- matrix(seq_along(metric_colors), nrow = 10, byrow = TRUE)
  color_idx <- as.vector(color_idx[color_idx <= length(metric_colors)])
  metric_colors <- metric_colors[color_idx]
  names(metric_colors) <- names(titles)
  pch.vals <- rep(c(1, 16, 17, 15, 3, 4, 8), length.out = length(titles))
  names(pch.vals) <- names(titles)

  bad_ref <- setdiff(ref_par, names(titles_params))
  if (length(bad_ref) > 0) {
    stop(paste("Parameter name", bad_ref[1], "is not correct."))
  }
  params <- titles_params[ref_par]

  if (!is.null(compare_par)) {
    compare_par <- make.names(compare_par)
    if (!(compare_par %in% names(titles_params))) {
      stop(paste("Compare parameter name", compare_par, "is not correct."))
    }
    if (!(compare_par %in% colnames(benchmatrix))) {
      stop(paste("Compare parameter", compare_par, "not found in benchmark matrix."))
    }
  }

  titles <- resolve_benchmark_titles(
    benchmarks = benchmarks,
    benchmark_level = benchmark_level,
    benchmark_category = benchmark_category,
    available_names = names(benchmatrix)
  )
  benchmarks <- names(titles)
  if (length(benchmarks) == 0) {
    stop("No benchmarks selected for visualization.")
  }
  group_meta <- bm_meta[bm_meta$benchmark %in% benchmarks, c("benchmark", "level", "category"), drop = FALSE]
  group_meta$level <- factor(group_meta$level, levels = level_order)
  group_meta$category <- factor(group_meta$category, levels = category_order)
  group_meta <- group_meta[order(group_meta$category, group_meta$level, group_meta$benchmark), , drop = FALSE]
  benchmarks <- group_meta$benchmark
  titles <- titles[benchmarks]
  axis_labels <- axis_labels[benchmarks]
  metric_colors <- metric_colors[benchmarks]
  pch.vals <- pch.vals[benchmarks]
  pretty_label <- function(x) {
    gsub("\\b([a-z])", "\\U\\1", gsub("_", " ", x), perl = TRUE)
  }

  # Draw a plain text header strip without panel boxes.
  draw_header <- function(label, col, cex = 1, right_label = NULL,
                          legend_title = NULL, legend_labels = NULL, legend_cols = NULL, legend_pch = NULL) {
    par(mar = c(0.2, 0.8, 0.2, 0.8))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    segments(0.04, 0.18, 0.96, 0.18, col = grDevices::adjustcolor(col, alpha.f = 0.45), lwd = 1)
    text(
      x = 0.07, y = 0.58,
      labels = label,
      adj = c(0, 0.5), cex = cex, font = 2, col = col
    )
    if (!is.null(right_label)) {
      text(x = 0.93, y = 0.58, labels = right_label, adj = c(1, 0.5), cex = 0.76, col = "#333333")
    }
    if (!is.null(legend_labels) && length(legend_labels) > 0) {
      legend(
        "bottom",
        title = legend_title,
        legend = legend_labels,
        col = legend_cols,
        pch = legend_pch,
        horiz = TRUE,
        xpd = NA,
        cex = 0.72,
        pt.cex = 0.72,
        x.intersp = 0.6,
        y.intersp = 0.7,
        bty = "n",
        inset = c(0, 0.02)
      )
    }
    par(mar = c(3, 3, 1.5, 1.5))
  }
  # Split one category into pages while keeping level blocks intact.
  build_pages <- function(category_meta, n_cols = 4, max_rows = 3) {
    # Keep category pages compact by carrying whole level blocks together.
    level_rows <- setNames(
      vapply(unique(as.character(category_meta$level)),
             function(level_id) ceiling(sum(category_meta$level == level_id) / n_cols),
             numeric(1)),
      unique(as.character(category_meta$level))
    )
    pages <- list()
    keep_levels <- character()
    row_count <- 0
    for (level_id in names(level_rows)) {
      if (length(keep_levels) > 0 && row_count + level_rows[[level_id]] > max_rows) {
        pages[[length(pages) + 1]] <- category_meta[category_meta$level %in% keep_levels, , drop = FALSE]
        keep_levels <- character()
        row_count <- 0
      }
      keep_levels <- c(keep_levels, level_id)
      row_count <- row_count + level_rows[[level_id]]
    }
    if (length(keep_levels) > 0) {
      pages[[length(pages) + 1]] <- category_meta[category_meta$level %in% keep_levels, , drop = FALSE]
    }
    pages
  }
  # Build the layout matrix for one page of grouped benchmark panels.
  build_layout <- function(page_meta, category_key, n_cols = 4, max_rows = 3) {
    # Each page starts with a category banner, then inserts one level strip
    # before the corresponding 4-column metric rows. Blank rows pad short
    # pages so metric panels keep the same size across pages.
    level_ids <- unique(as.character(page_meta$level))
    mat <- heights <- NULL
    items <- character()
    id <- 1L
    metric_rows_used <- 0L
    level_headers_used <- 0L
    mat <- rbind(mat, rep(id, n_cols))
    heights <- c(heights, 0.22)
    items <- c(items, paste0("page::", category_key))
    id <- id + 1L
    for (level_id in level_ids) {
      level_meta <- page_meta[page_meta$level == level_id, , drop = FALSE]
      mat <- rbind(mat, rep(id, n_cols))
      heights <- c(heights, 0.18)
      items <- c(items, paste0("level::", level_id))
      id <- id + 1L
      level_headers_used <- level_headers_used + 1L
      for (start_idx in seq(1, nrow(level_meta), by = n_cols)) {
        row_items <- level_meta$benchmark[start_idx:min(start_idx + n_cols - 1L, nrow(level_meta))]
        row <- rep(0L, n_cols)
        for (col_idx in seq_along(row_items)) {
          row[col_idx] <- id
          items <- c(items, row_items[col_idx])
          id <- id + 1L
        }
        mat <- rbind(mat, row)
        heights <- c(heights, 1)
        metric_rows_used <- metric_rows_used + 1L
      }
    }
    if (level_headers_used < max_rows) {
      for (i in seq_len(max_rows - level_headers_used)) {
        mat <- rbind(mat, rep(0L, n_cols))
        heights <- c(heights, 0.18)
      }
    }
    if (metric_rows_used < max_rows) {
      for (i in seq_len(max_rows - metric_rows_used)) {
        mat <- rbind(mat, rep(0L, n_cols))
        heights <- c(heights, 1)
      }
    }
    list(layout = matrix(as.integer(mat), nrow = nrow(mat)), heights = heights, items = items)
  }
  compare_values <- compare_colors <- NULL
  if (!is.null(compare_par)) {
    compare_values <- unique(benchmatrix[, compare_par])
    compare_values <- if (is.numeric(compare_values)) sort(compare_values) else sort(as.character(compare_values))
    compare_colors <- colorspace::qualitative_hcl(n = length(compare_values), palette = "Dark 3")
    if (!is.null(cols)) {
      compare_colors <- rep(cols, length.out = length(compare_values))
    }
  } else if (!is.null(cols)) {
    metric_colors <- rep(cols, length.out = length(metric_colors))
  }

  category_pages <- split(group_meta, factor(as.character(group_meta$category), levels = category_order))
  category_pages <- lapply(category_pages[vapply(category_pages, nrow, integer(1)) > 0], build_pages)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    layout(matrix(1L))
    par(old_par)
  }, add = TRUE)

  # Draw points, error bars, or uncertainty bands on 1D benchmark panels.
  draw_uncertainty <- function(x, y, sd, col, pch, lty = 1) {
    ok <- is.finite(x) & is.finite(y)
    ok_sd <- is.finite(sd)
    if (errorbar && errorstyle == "area" && sum(ok & ok_sd) >= 2) {
      ord <- order(x)
      polygon(
        c(x[ord][ok_sd[ord] & ok[ord]], rev(x[ord][ok_sd[ord] & ok[ord]])),
        c((y - sd)[ord][ok_sd[ord] & ok[ord]], rev((y + sd)[ord][ok_sd[ord] & ok[ord]])),
        col = grDevices::adjustcolor(col, alpha.f = 0.2),
        border = col, lty = lty, lwd = 0.6
      )
    }
    if (errorbar && errorstyle == "bar" && any(ok & ok_sd)) {
      segments(
        x[ok & ok_sd], y[ok & ok_sd] - sd[ok & ok_sd],
        x[ok & ok_sd], y[ok & ok_sd] + sd[ok & ok_sd],
        col = grDevices::adjustcolor(col, alpha.f = 0.6),
        lwd = 1
      )
    }
    if (sum(ok) >= 2) {
      ord <- order(x)
      lines(x[ord][ok[ord]], y[ord][ok[ord]], col = col, lty = lty)
    }
    points(x[ok], y[ok], pch = pch, col = col)
  }
  # Set the plotting range for one metric, optionally using predefined bounds.
  get_range <- function(metric) {
    out <- range(benchmatrix[, metric], na.rm = TRUE)
    if (fullrange) {
      if (!is.na(ranges[[metric]][1])) out[1] <- ranges[[metric]][1]
      if (!is.na(ranges[[metric]][2])) out[2] <- ranges[[metric]][2]
    }
    out
  }
  # Draw a single benchmark against one reference parameter.
  draw_metric_1d <- function(metric, col, pch) {
    x <- benchmatrix[, ref_par]
    y <- benchmatrix[, metric]
    yr <- get_range(metric)
    if (!all(is.finite(range(y, na.rm = TRUE)))) {
      plot.new()
      title(main = titles[metric], col.main = "black", font.main = 2)
      return(invisible())
    }
    if (!is.null(compare_par)) {
      # Compare mode overlays subsets for a second parameter on the same panel.
      x_base <- unique(x)
      x_pos <- if (is.numeric(x_base)) as.numeric(x_base) else seq_along(x_base)
      pch_comp <- rep(pch.vals, length.out = length(compare_values))
      plot(x_pos, rep(NA, length(x_pos)), type = "n", xaxt = "n", xlab = params, ylab = axis_labels[metric], ylim = yr)
      axis(1, at = x_pos, labels = x_base, las = 1)
      title(main = titles[metric], col.main = "black", font.main = 2)
      for (j in seq_along(compare_values)) {
        sel <- as.character(benchmatrix[, compare_par]) == compare_values[j]
        x_raw <- benchmatrix[sel, ref_par]
        y_raw <- benchmatrix[sel, metric]
        x_plot <- x_pos[match(as.character(x_raw), as.character(x_base))]
        if (errorbar) {
          y_mean <- tapply(y_raw, x_raw, mean, na.rm = TRUE)
          y_sd <- tapply(y_raw, x_raw, sd, na.rm = TRUE)
          y_sd[!is.finite(y_sd)] <- 0
          draw_uncertainty(x_pos[match(names(y_mean), as.character(x_base))], y_mean, y_sd, compare_colors[j], pch_comp[j], j)
        } else {
          points(x_plot, y_raw, pch = pch_comp[j], col = compare_colors[j])
        }
      }
      return(invisible())
    }
    if (errorbar) {
      # Error bars summarize repeated runs at each x value.
      y_sd <- tapply(y, x, sd, na.rm = TRUE)
      y_mean <- tapply(y, x, mean, na.rm = TRUE)
      x_ord <- names(y_mean)
      x_plot <- if (is.numeric(x_ord)) as.numeric(x_ord) else seq_along(x_ord)
      plot(x_plot, y_mean, type = "n", xaxt = "n", xlab = params, ylab = axis_labels[metric], ylim = yr)
      draw_uncertainty(x_plot, y_mean, y_sd, col, pch)
      axis(1, at = x_plot, labels = x_ord, las = 1)
    } else {
      x_plot <- x
      axis_at <- x
      if (any(is.character(x_plot))) {
        xf <- factor(x, levels = unique(x))
        names(xf) <- x
        x_plot <- xf[x_plot]
        axis_at <- x_plot
      }
      gplots::plotCI(x_plot, y, uiw = NA, gap = 0, xaxt = "n", sfrac = 0.02,
                     xlab = params, ylab = axis_labels[metric], pch = pch, col = col,
                     cex = 1.2, cex.lab = 0.9, cex.axis = 0.85, ylim = yr)
      axis(1, at = axis_at, labels = unique(x), las = 1)
    }
    title(main = titles[metric], col.main = "black", font.main = 2)
  }
  # Draw a single benchmark against two reference parameters as a heatmap.
  draw_metric_2d <- function(metric, col) {
    # Two reference parameters are shown as a heatmap plus a compact side scale.
    x_vals <- sort(unique(benchmatrix[, ref_par[1]]))
    y_vals <- sort(unique(benchmatrix[, ref_par[2]]))
    z_mat <- matrix(NA, nrow = length(x_vals), ncol = length(y_vals))
    for (j in seq_len(nrow(benchmatrix))) {
      z_mat[which(x_vals == benchmatrix[j, ref_par[1]]), which(y_vals == benchmatrix[j, ref_par[2]])] <- benchmatrix[j, metric]
    }
    image_colors <- colorRampPalette(c("white", col))(100)
    zlim <- range(z_mat, na.rm = TRUE)
    if (any(!is.finite(zlim))) zlim <- c(0, 0)
    image(x_vals, y_vals, z_mat, col = image_colors, zlim = zlim, xlab = params[1], ylab = params[2], axes = TRUE)
    title(main = titles[metric], col.main = "black", font.main = 2)
    if (length(x_vals) == 1) x_vals <- x_vals + c(-0.5, 0.5)
    if (length(y_vals) == 1) y_vals <- y_vals + c(-0.5, 0.5)
    bar_x <- max(x_vals) - diff(range(x_vals)) * 0.5
    bar_w <- diff(range(x_vals)) * 0.03
    bar_y <- seq(min(y_vals), max(y_vals), length.out = length(image_colors) + 1)
    for (k in seq_along(image_colors)) {
      rect(bar_x, bar_y[k], bar_x + bar_w, bar_y[k + 1], col = image_colors[k], border = NA)
    }
    rect(bar_x, min(y_vals), bar_x + bar_w, max(y_vals), col = NA, border = "#333333", lwd = 0.5)
    if (diff(zlim) > 0) {
      labels <- pretty(zlim, n = 3)
      label_y <- approx(zlim, range(bar_y), xout = labels)$y
      text(bar_x + bar_w + 0.02 * diff(range(x_vals)), label_y, labels = round(labels, 2), cex = 0.6, adj = 0)
    }
  }

  for (category_key in names(category_pages)) {
    pages <- category_pages[[category_key]]
    for (page_idx in seq_along(pages)) {
      page_meta <- pages[[page_idx]]
      lay <- build_layout(page_meta, category_key)
      # layout() is needed here because the category and level headers span full rows.
      layout(matrix(1L))
      layout(lay$layout, heights = lay$heights)
      par(cex.main = 0.9, cex.lab = 0.75, cex.axis = 0.7,
          mgp = c(1.5, 0.3, 0), mar = c(2.2, 2.2, 1.2, 0.8), xpd = FALSE, font.main = 2)
      for (item in lay$items) {
        if (startsWith(item, "page::")) {
          cc <- category_colors[sub("^page::", "", item)]
          draw_header(pretty_label(sub("^page::", "", item)),
                      cc,
                      cex = 1.02,
                      right_label = paste0("Page ", page_idx, "/", length(pages)),
                      legend_title = if (!is.null(compare_par)) compare_par else NULL,
                      legend_labels = if (!is.null(compare_par)) compare_values else NULL,
                      legend_cols = if (!is.null(compare_par)) compare_colors else NULL,
                      legend_pch = if (!is.null(compare_par)) rep(pch.vals, length.out = length(compare_values)) else NULL)
        } else if (startsWith(item, "level::")) {
          cc <- category_colors[category_key]
          draw_header(pretty_label(sub("^level::", "", item)),
                      cc,
                      cex = 0.82)
        } else if (length(ref_par) == 1) {
          draw_metric_1d(item, "black", pch.vals[[item]])
        } else {
          draw_metric_2d(item, metric_colors[[item]])
        }
      }
    }
  }
  invisible(NULL)
}


#' Render a Benchmark Summary Table
#'
#' This function summarizes benchmark results from a matrix of metrics
#' (typically produced by proteomics simulations or evaluation pipelines).
#' It returns a human-readable table with either raw values (if one row)
#' or value ranges (if multiple rows). It also prints a markdown version
#' of the table in the console
#'
#' @param benchmatrix A data frame or matrix containing benchmark metrics.
#' @param benchmarks A character vector of benchmark variable names, or a numeric vector of column indices.
#'        If NULL, all grouped benchmarks will be used.
#' @param benchmark_level A character vector specifying benchmark level groups
#'   to include.
#' @param benchmark_category A character vector specifying benchmark purpose
#'   groups to include.
#'
#' @return A data frame summarizing the selected benchmark metrics.
render_benchmark_table <- function(benchmatrix,
                                   benchmarks = NULL,
                                   benchmark_level = NULL,
                                   benchmark_category = NULL) {
  colnames(benchmatrix) <- make.names(colnames(benchmatrix))

  titles <- resolve_benchmark_titles(
    benchmarks = benchmarks,
    benchmark_level = benchmark_level,
    benchmark_category = benchmark_category,
    available_names = colnames(benchmatrix)
  )

  selected_cols <- intersect(names(titles), colnames(benchmatrix))
  if (length(selected_cols) == 0) {
    stop("No benchmarks selected for the summary table.")
  }
  out_data <- benchmatrix[, selected_cols, drop = FALSE]

  # Create summary output
  summarize_column <- function(x) {
    if (length(x) == 1 || nrow(benchmatrix) == 1) {
      round(x[1], 4)
    } else {
      rng <- range(x, na.rm = TRUE)
      if (isTRUE(all.equal(rng[1], rng[2], tolerance = 1e-6))) {
        round(rng[1], 4)
      } else {
        paste0(round(rng[1], 4), "-", round(rng[2], 4))
      }
    }
  }

  summary_df <- as.data.frame(lapply(out_data, summarize_column))
  colnames(summary_df) <- titles[selected_cols]

  # Print markdown table if requested
  if (requireNamespace("knitr", quietly = TRUE)) {
    print(knitr::kable(summary_df, format = "markdown", align = "c"))
  } else {
    warning("knitr package not installed. Install it or set print_table = FALSE.")
  }

  return(summary_df)
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
  bm_titles <- get_bmtitles()

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
  if (length(to_del) > 0) {
    BenchMatrix <- BenchMatrix[, -to_del]
  }
  BenchMatrix[is.na(BenchMatrix)] <- 0
  # remove column QuantColnames if it exists
  if ("QuantColnames" %in% colnames(BenchMatrix)) {
    BenchMatrix <- BenchMatrix[, -which(colnames(BenchMatrix) == "QuantColnames")]
  }
  tBenchMatrix <- BenchMatrix
  nr <- 2 # nrow(BenchMatrix)
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
  if (is.numeric(sim)) {
    sim <- rownames(BenchMatrix)[sim]
  }
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
  dat2 <- sapply(dat2, function(x) {
    if (is.numeric(x)) {
      x <- round(x, 2)
    }
    x
  })
  x_labels <- names(dat)
  matching_titles <- bm_titles[x_labels]
  x_labels[!is.na(matching_titles)] <- matching_titles[!is.na(matching_titles)]
  plot <- plotly::plot_ly(
    x = x_labels,
    y = as.numeric(dat),
    type = "bar",
    marker = list(
      color = color_palette[as.numeric(dat) * 99 + 1]
    ),
    text = as.character(dat2),
    textposition = "auto",
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      title = paste("hash:", sim),
      yaxis = list(title = "Normalized values", range = c(0, 1), tickfont = list(size = 18 / nr), titlefont = list(size = 20 / nr)),
      xaxis = list(title = "", tickangle = -45, tickfont = list(size = 16 / nr)),
      margin = list(t = 100, b = 100, l = 100, r = 100),
      showlegend = FALSE
    )
  param_plot <- plot_params(param_values, sim)


  # Combine all plots into a subplot
  subplot(param_plot, plot,
          nrows = 2, shareX = TRUE, shareY = TRUE,
          margin = 0.01, titleX = TRUE, titleY = TRUE
  ) %>%
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

  # filter out every NA or NULL parameters
  to_remove <- c()
  for (i in which(colnames(BenchMatrix) %in% param_names)) {
    val <- as.numeric(BenchMatrix[, i])
    if (!(length(val) == 0)) {
      if (all(is.na(val))) {
        to_remove <- append(to_remove, i)
      }
    }
  }
  if (length(to_remove) > 0) {
    BenchMatrix <- BenchMatrix[, -to_remove]
  }
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
    y_domain <- y_domain / 2 + 0.5
    # y_domain <- y_domain/(nrow(BenchMatrix)*1.1) + (current_row-1)*1.05 / nrow(BenchMatrix)
    fig1 <- plot_ly(
      domain = list(x = x_domain, y = y_domain),
      value = val,
      title = list(text = param_names[i], align = "center", font = list(size = 10 / nr)),
      type = "indicator",
      mode = "gauge",
      gauge = list(
        shape = "bullet",
        axis = list(
          range = list(curr_par$MinValue, curr_par$MaxValue),
          tickmode = "auto",
          ticks = "inside", # Position ticks inside
          tickfont = list(size = 10 / nr, color = "black") # Adjust tickfont size and color for better visibility
        ),
        bar = list(
          color = colorpanel(100, "#33CC33", "#999966", "#CC3333")[(val - curr_par$MinValue) / (curr_par$MaxValue - curr_par$MinValue) * 100 + 1],
          thickness = 0.8 # Adjust this value to make the gauge bar thicker
        )
      )
    )
    meters[[i]] <- fig1
  }
  # Combine all meter plots
  subplot(meters, nrows = nrows, margin = 0.01)
}



#' Function to get human-readable titles for the different parameters
#'
get_bmmeta <- function() {
  meta <- data.frame(
    benchmark = c(
      "numPeptides", "numProteins", "propUniquePep", "uniqueStrippedPep", "percMissingPep",
      "aucDiffRegPeptides.FDR_limma.2.vs.1.AUC", "tprPep0.01.FDR_limma.2.vs.1.TPR", "tprPep0.05.FDR_limma.2.vs.1.TPR",
      "tFDRPep0.01.FDR_limma.2.vs.1.tFDR", "tFDRPep0.05.FDR_limma.2.vs.1.tFDR", "propMisCleavedPeps.0", "dynRangePep",
      "meanSquareDiffFCPep", "sdWithinRepsPep", "skewnessPeps", "kurtosisPeps", "sdPeps",
      "numQuantProtGroups", "propUniqueProts", "percMissingProt", "meanPepPerProt",
      "aucDiffRegProteins.FDR_PolySTest.2.vs.1.AUC", "tprProt0.01.FDR_PolySTest.2.vs.1.TPR", "tprProt0.05.FDR_PolySTest.2.vs.1.TPR",
      "tFDRProt0.01.FDR_PolySTest.2.vs.1.tFDR", "tFDRProt0.05.FDR_PolySTest.2.vs.1.tFDR", "meanSquareDiffFCProt",
      "sdWithinRepsProt", "propMisCleavedProts", "propDiffRegWrongIDProt0.01.FDR_PolySTest.2.vs.1",
      "propDiffRegWrongIDProt0.05.FDR_PolySTest.2.vs.1", "skewnessProts", "kurtosisProts", "sdProts",
      "numProteoforms", "meanProteoformsPerProt", "numModPeptides", "propModAndUnmodPep",
      "aucDiffRegAdjModPep", "tprAdjModPep0.01", "tprAdjModPep0.05", "tFDRAdjModPep0.01", "tFDRAdjModPep0.05",
      "propDiffRegPepWrong0.01.FDR_PolySTest.2.vs.1", "propDiffRegPepWrong0.05.FDR_PolySTest.2.vs.1",
      "percOverlapModPepProt", "meanSquareDiffFCModPep"
    ),
    title = c(
      "#Total peptidoforms (mod+unmod)", "#Protein accessions (all peptidoforms)", "Fraction unique peptidoforms (single protein)",
      "#Unique peptide sequences", "% Missing peptidoform values", "Peptidoform AUC (truth vs FDR)",
      "TPR peptidoforms (FDR < 0.01)", "TPR peptidoforms (FDR < 0.05)", "True FDR (peptidoforms, 0.01)",
      "True FDR (peptidoforms, 0.05)", "No miscleavage fraction", "Dynamic range (peptidoforms)",
      "Fold-change error (peptidoforms)", "Replicate SD (regulated peptidoforms)", "Skewness (peptidoforms)",
      "Kurtosis (peptidoforms)", "Overall SD (peptidoforms)", "#Quantified protein groups", "Fraction single-protein groups",
      "Fraction missing protein-group values", "Peptide sequences per protein group", "Protein-group AUC (truth vs FDR)",
      "TPR protein groups (FDR < 0.01)", "TPR protein groups (FDR < 0.05)", "True FDR (protein groups, 0.01)",
      "True FDR (protein groups, 0.05)", "Fold-change error (protein groups)", "Replicate SD (regulated protein groups)",
      "Fraction miscleaved protein groups", "Wrong-ID fraction (protein groups, FDR<0.01)", "Wrong-ID fraction (protein groups, FDR<0.05)",
      "Skewness (protein groups)", "Kurtosis (protein groups)", "Overall SD (protein groups)", "#Proteoforms",
      "Proteoforms per protein", "#Modified peptidoforms", "Fraction modified with unmodified match",
      "AUC (adj. mod. peptidoforms)", "TPR (adj. mod. peptidoforms, 0.01)", "TPR (adj. mod. peptidoforms, 0.05)",
      "True FDR (mod. peptidoforms, 0.01)", "True FDR (mod. peptidoforms, 0.05)",
      "Wrong-signif. fraction (mod. peptidoforms, FDR<0.01)", "Wrong-signif. fraction (mod. peptidoforms, FDR<0.05)",
      "Fraction mod peptidoforms with protein quant", "Fold-change error (mod. peptidoforms)"
    ),
    axis_label = c(
      "Count", "Count", "Fraction", "Count", "% Missing", "AUC", "TPR", "TPR", "True FDR", "True FDR",
      "Fraction", "log2 Range", "FC Error", "SD", "Skewness", "Kurtosis", "SD",
      "Count", "Fraction", "% Missing", "Count", "AUC", "TPR", "TPR", "True FDR", "True FDR",
      "FC Error", "SD", "Fraction", "Fraction", "Fraction", "Skewness", "Kurtosis", "SD",
      "Count", "Count", "Count", "Fraction", "AUC", "TPR", "TPR", "True FDR", "True FDR",
      "Fraction", "Fraction", "Fraction", "FC Error"
    ),
    level = c(rep("peptidoform", 17), rep("protein_group", 17), rep("ptm_proteoform", 13)),
    category = c(
      "coverage", "coverage", "coverage", "coverage", "completeness",
      "differential_performance", "differential_performance", "differential_performance",
      "differential_performance", "differential_performance", "coverage", "quantification_quality",
      "differential_performance", "quantification_quality", "quantification_quality", "quantification_quality", "quantification_quality",
      "coverage", "coverage", "completeness", "coverage", "differential_performance", "differential_performance",
      "differential_performance", "differential_performance", "differential_performance", "differential_performance",
      "quantification_quality", "coverage", "artifact_sensitivity", "artifact_sensitivity", "quantification_quality",
      "quantification_quality", "quantification_quality", "coverage", "coverage", "coverage", "coverage",
      "ptm_adjusted_performance", "ptm_adjusted_performance", "ptm_adjusted_performance", "ptm_adjusted_performance",
      "ptm_adjusted_performance", "artifact_sensitivity", "artifact_sensitivity", "coverage", "ptm_adjusted_performance"
    ),
    direction = c(
      "higher_better", "higher_better", "higher_better", "higher_better", "lower_better",
      "higher_better", "higher_better", "higher_better", "lower_better", "lower_better", "depends", "higher_better",
      "lower_better", "lower_better", "depends", "depends", "depends",
      "higher_better", "higher_better", "lower_better", "higher_better", "higher_better", "higher_better", "higher_better",
      "lower_better", "lower_better", "lower_better", "lower_better", "lower_better", "lower_better", "lower_better",
      "depends", "depends", "depends", "higher_better", "higher_better", "higher_better", "higher_better",
      "higher_better", "higher_better", "higher_better", "lower_better", "lower_better", "lower_better", "lower_better",
      "higher_better", "lower_better"
    ),
    range_min = c(
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, NA, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, NA, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ),
    range_max = c(
      NA, NA, 1, NA, 100, 1, 1, 1, 1, 1, 1, NA, NA, NA, NA, NA, NA,
      NA, 1, NA, NA, 1, 1, 1, 1, 1, NA, NA, 1, 1, 1, NA, NA, NA,
      NA, NA, NA, 1, NA, NA, NA, NA, NA, 1, 1, NA, NA
    ),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$benchmark
  meta
}

get_bmtitles <- function() {
  meta <- get_bmmeta()
  stats::setNames(meta$title, meta$benchmark)
}

get_bmgroups <- function() {
  meta <- get_bmmeta()
  meta[, c("benchmark", "level", "category", "direction", "title", "axis_label"), drop = FALSE]
}

get_bmaxislabels <- function() {
  meta <- get_bmmeta()
  stats::setNames(meta$axis_label, meta$benchmark)
}

resolve_benchmark_titles <- function(benchmarks = NULL,
                                     benchmark_level = NULL,
                                     benchmark_category = NULL,
                                     available_names = NULL,
                                     default_first = NULL) {
  meta <- get_bmmeta()
  titles <- stats::setNames(meta$title, meta$benchmark)
  groups <- meta[, c("benchmark", "level", "category"), drop = FALSE]

  if (!is.null(benchmark_level)) {
    benchmark_level <- unique(as.character(benchmark_level))
    invalid_levels <- setdiff(benchmark_level, unique(groups$level))
    if (length(invalid_levels) > 0) {
      stop(paste("Unknown benchmark level(s):", paste(invalid_levels, collapse = ", ")))
    }
    titles <- titles[groups$benchmark[groups$level %in% benchmark_level]]
  }

  if (!is.null(benchmark_category)) {
    benchmark_category <- unique(as.character(benchmark_category))
    invalid_categories <- setdiff(benchmark_category, unique(groups$category))
    if (length(invalid_categories) > 0) {
      stop(paste("Unknown benchmark category(s):", paste(invalid_categories, collapse = ", ")))
    }
    groups <- groups[groups$benchmark %in% names(titles), , drop = FALSE]
    titles <- titles[groups$benchmark[groups$category %in% benchmark_category]]
  }

  if (is.null(benchmarks) || length(benchmarks) == 0) {
    groups <- groups[groups$benchmark %in% names(titles), , drop = FALSE]
    titles <- titles[groups$benchmark]
  } else if (is.character(benchmarks)) {
    benchmarks <- make.names(benchmarks)
    invalid_benchmarks <- setdiff(benchmarks, names(titles))
    if (length(invalid_benchmarks) > 0) {
      stop(paste("Benchmark name(s)", paste(invalid_benchmarks, collapse = ", "), "are not correct."))
    }
    titles <- titles[benchmarks]
  } else if (is.numeric(benchmarks)) {
    if (length(titles) == 0 || max(benchmarks) > length(titles)) {
      stop(paste("Too many benchmarks selected. Please select a maximum of", length(titles)))
    }
    titles <- titles[benchmarks]
  } else {
    stop("Please provide a character vector of benchmark names or numeric vector of indices.")
  }

  if (!is.null(available_names)) {
    missing_benchmarks <- setdiff(names(titles), available_names)
    if (length(missing_benchmarks) > 0) {
      warning(paste(
        "Benchmark name(s)", paste(missing_benchmarks, collapse = ", "),
        "are not available in the benchmark matrix. They will be removed."
      ))
      titles <- titles[intersect(names(titles), available_names)]
    }
  }

  titles
}

#' setting maximal ranges for benchmarking metrics
#'
#' @param titles A named character vector of benchmark metric titles, as
#'   returned by \code{get_bmtitles()}.
#' @return A named list of numeric vectors specifying the \code{c(min, max)}
#'   display range for each benchmark metric.
set_bmranges <- function(titles) {
  meta <- get_bmmeta()
  selected <- meta[names(titles), c("range_min", "range_max"), drop = FALSE]
  stats::setNames(lapply(seq_len(nrow(selected)), function(i) as.numeric(selected[i, ])), rownames(selected))
}

get_paramtitles <- function() {
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
    ThreshRemoveProteoforms = "Proteoform removal threshold",
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
    PercDetectability = "Fraction of peptides detected (model)",
    PercDetectedVal = "Fraction of detected intensities",
    WeightDetectVal = "Intensity-dependence of detection",
    MSNoise = "Instrumental noise",
    WrongIDs = "Wrong identification rate",
    WrongLocalizations = "Wrong localization rate",
    MaxNAPerPep = "Max NAs per peptide",
    ProtSummarization = "Protein summarization method",
    MinUniquePep = "Min unique peptides per protein",
    IncludeModPep = "Include modified peptidoforms in protein quantification",
    SharedPep = "Include shared peptides in protein inference",
    StatPaired = "Paired testing enabled"
  )
  titles_params
}
