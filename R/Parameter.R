# Function to create parameter table
param_table <- function() {
    
    yaml_file <- system.file("config", "parameters.yaml", package = "PhosFake")
    
    # Read the YAML file
    params <- yaml::yaml.load_file(yaml_file)$params
    
    # Convert NA values from strings to real NA
    for (l in names(params)) {
        for (k in names(params[[l]])) {
            if (params[[l]][[k]] == "NA") {
                params[[l]][[k]] <- NA
            }
        }
    }
    
    # Convert the list to a data frames
    params_table <- do.call(rbind, lapply(params, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
    
    # Name the columns
    colnames(params_table) <- c("Category", "Group", "Explanation", "MinValue", "MaxValue", "DefaultValue")
    
    
    # Ensure MinValue and MaxValue are numeric
    params_table
}

