#' Create Parameter Table from YAML Configuration
#'
#' This function reads a YAML file containing parameter settings, processes them
#' to handle `NA` values, and converts the parameters into a data frame. The resulting
#' table includes categories, groups, explanations, and min/max/default values for each parameter.
#'
#' @return A data frame containing the parameters from the YAML file with the following columns:
#' \describe{
#'   \item{Category}{The category to which the parameter belongs.}
#'   \item{Group}{The group within the category.}
#'   \item{Explanation}{A description of the parameter.}
#'   \item{MinValue}{The minimum allowable value for the parameter.}
#'   \item{MaxValue}{The maximum allowable value for the parameter.}
#'   \item{DefaultValue}{The default value for the parameter.}
#' }
#'
#' @importFrom yaml yaml.load_file
#' @keywords internal
param_table <- function() {

    yaml_file <- system.file("config", "parameters.yaml", package = "ProteoMaker")

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

    # Convert the list to a data frames
    params_table <- do.call(rbind, lapply(params, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))

    # Name the columns
    colnames(params_table) <- c("Category", "Group", "Explanation", "MinValue", "MaxValue", "DefaultValue")


    # Ensure MinValue and MaxValue are numeric
    params_table
}

