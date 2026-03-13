pkgname <- "ProteoMaker"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ProteoMaker')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("def_param")
### * def_param

flush(stderr()); flush(stdout())

### Name: def_param
### Title: Generate default parameters or read from yaml file
### Aliases: def_param

### ** Examples

# Read YAML file from inst folder
yaml_path <- system.file("config", "params.yaml", package = "ProteoMaker")
params <- def_param(yaml_path)



cleanEx()
nameEx("gather_all_sims")
### * gather_all_sims

flush(stderr()); flush(stdout())

### Name: gather_all_sims
### Title: Gather parameters and benchmarks from all available runs
### Aliases: gather_all_sims

### ** Examples

config <- set_proteomaker(resultFilePath = tempdir())
all_results <- gather_all_sims(config)




cleanEx()
nameEx("generate_combinations")
### * generate_combinations

flush(stderr()); flush(stdout())

### Name: generate_combinations
### Title: Generate all combinations of parameter sets
### Aliases: generate_combinations

### ** Examples

params <- def_param()
params$paramGroundTruth$NumReps <- c(2, 4, 6)
combinations <- generate_combinations(params)



cleanEx()
nameEx("get_simulation")
### * get_simulation

flush(stderr()); flush(stdout())

### Name: get_simulation
### Title: Retrieve intermediate outputs for a simulation
### Aliases: get_simulation

### ** Examples

config <- set_proteomaker(resultFilePath = tempdir(), cores = NULL)
Param <- def_param()
Param$paramGroundTruth$NumReps <- 5
benchmarks <- run_sims(Param, config)
result <- get_simulation(benchmarks[[1]]$Param, config, stage = "MSRun")



cleanEx()
nameEx("get_stage_hashes")
### * get_stage_hashes

flush(stderr()); flush(stdout())

### Name: get_stage_hashes
### Title: Retrieve hashes for intermediate simulation stages
### Aliases: get_stage_hashes

### ** Examples

Param <- def_param()
hashes <- get_stage_hashes(Param)
hashes["MSRun"]



cleanEx()
nameEx("matrix_benchmarks")
### * matrix_benchmarks

flush(stderr()); flush(stdout())

### Name: matrix_benchmarks
### Title: Create a matrix of benchmarks from simulation results
### Aliases: matrix_benchmarks

### ** Examples

conf <- set_proteomaker(resultFilePath = tempdir(), cores = NULL)
results <- run_sims(def_param(), conf)
benchmark_matrix <- matrix_benchmarks(results, conf)




cleanEx()
nameEx("plot_params")
### * plot_params

flush(stderr()); flush(stdout())

### Name: plot_params
### Title: Plot parameters for a specific simulation
### Aliases: plot_params

### ** Examples

## Not run: 
##D benchmarks <- matrix_benchmarks(run_sims(def_param(), set_proteomaker(cores = NULL)), set_proteomaker(cores = NULL))
##D plot_params(benchmarks, current_row = 1)
## End(Not run)



cleanEx()
nameEx("run_sims")
### * run_sims

flush(stderr()); flush(stdout())

### Name: run_sims
### Title: Run all simulations for given parameters and configurations
### Aliases: run_sims

### ** Examples

params <- def_param()
config <- set_proteomaker(cores = NULL)
results <- run_sims(params, config)



cleanEx()
nameEx("set_proteomaker")
### * set_proteomaker

flush(stderr()); flush(stdout())

### Name: set_proteomaker
### Title: Set paths and general configuration for ProteoMaker
### Aliases: set_proteomaker

### ** Examples

config <- set_proteomaker()
config <- set_proteomaker(
  fastaFilePath = "CustomProteomes",
  resultFilePath = "Results", cores = 4, clusterType = "PSOCK"
)



cleanEx()
nameEx("visualize_benchmarks")
### * visualize_benchmarks

flush(stderr()); flush(stdout())

### Name: visualize_benchmarks
### Title: Visualize benchmarks vs 1 or two parameters
### Aliases: visualize_benchmarks

### ** Examples

conf <- set_proteomaker(resultFilePath = tempdir(), cores = NULL)
results <- run_sims(def_param(), conf)
benchmark_matrix <- matrix_benchmarks(results, conf)
visualize_benchmarks(benchmark_matrix, ref_par = "NumReps")




cleanEx()
nameEx("visualize_one_sim")
### * visualize_one_sim

flush(stderr()); flush(stdout())

### Name: visualize_one_sim
### Title: Visualize benchmarks for a specific simulation
### Aliases: visualize_one_sim

### ** Examples

benchmarks <- matrix_benchmarks(run_sims(def_param(), set_proteomaker(cores = NULL)), set_proteomaker(cores = NULL))
visualize_one_sim(benchmarks, current_row = 1)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
