library(PhosFake)
# Load the YAML file
yaml_file <- system.file("config", "parameters.yaml", package = "PhosFake")
if (!file.exists(yaml_file) || file.size(yaml_file) == 0) {
  stop("YAML file could not be loaded or is empty.")
}
def_params <- yaml::yaml.load_file(yaml_file)$params
if (is.null(def_params) || length(def_params) == 0) {
  stop("The params object is empty or invalid.")
}

phosfake_config <- set_phosfake(fastaFilePath = system.file("Proteomes", package = "PhosFake"),
                                resultFilePath = "SimulatedDataSets",
                                cores = 4, clusterType = "PSOCK",
                                runStatTests = F,
                                calcAllBenchmarks = T
)


print("Loaded parameter defaults")
options(shiny.maxRequestSize=1*1024^3)
