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

#' Render Simulation Parameters as a Table
#'
#' Creates a compact table of parameter values used in a simulation and joins
#' them to the descriptions from the package parameter metadata.
#'
#' @param Param A nested parameter list returned by \code{def_param()} or a flat
#'   simulation parameter list.
#' @param columns Character vector of columns to return. Available columns are
#'   \code{"Parameter"}, \code{"Value"}, \code{"Description"}, \code{"Default"},
#'   \code{"Range"}, \code{"List"}, and \code{"Section"}.
#' @param print_table Logical; if \code{TRUE}, print a markdown table when
#'   \pkg{knitr} is available.
#'
#' @return A data frame with the selected columns.
#' @export
#'
#' @examples
#' params <- def_param()
#' render_parameter_table(params, print_table = FALSE)
render_parameter_table <- function(Param,
                                   columns = c("Parameter", "Value", "Description"),
                                   print_table = TRUE) {
  if (is.null(Param) || !is.list(Param)) {
    stop("Param must be a parameter list.")
  }

  param_groups <- c("paramGroundTruth", "paramProteoformAb", "paramDigest", "paramMSRun",
                    "paramDataAnalysis")
  if (any(names(Param) %in% param_groups)) {
    Param <- do.call(c, unname(Param[names(Param) %in% param_groups]))
  }

  format_value <- function(x) {
    if (is.null(x) || length(x) == 0) return("NULL")
    if (is.list(x) && !is.data.frame(x)) {
      if (is.null(names(x))) return(paste(unlist(x), collapse = ";"))
      return(paste(vapply(names(x), function(nm) {
        paste0(nm, "=", paste(unlist(x[[nm]]), collapse = ","))
      }, character(1)), collapse = ";"))
    }
    paste(as.character(x), collapse = ";")
  }

  values <- data.frame(
    Parameter = names(Param),
    Value = vapply(Param, format_value, character(1)),
    stringsAsFactors = FALSE
  )
  meta <- param_table()
  meta$Parameter <- rownames(meta)
  meta <- meta[match(values$Parameter, meta$Parameter), , drop = FALSE]
  missing_meta <- is.na(meta$Parameter)
  if (any(missing_meta)) {
    meta$Parameter[missing_meta] <- values$Parameter[missing_meta]
  }
  meta$Value <- values$Value

  meta$Default <- if ("DefaultValue" %in% colnames(meta)) {
    vapply(meta$DefaultValue, format_value, character(1))
  } else {
    NA_character_
  }
  meta$Range <- if (all(c("MinValue", "MaxValue") %in% colnames(meta))) {
    paste(vapply(meta$MinValue, format_value, character(1)),
          vapply(meta$MaxValue, format_value, character(1)), sep = "-")
  } else {
    NA_character_
  }
  meta$Description <- if ("Explanation" %in% colnames(meta)) {
    vapply(meta$Explanation, format_value, character(1))
  } else {
    NA_character_
  }
  meta$List <- if ("Group" %in% colnames(meta)) {
    vapply(meta$Group, format_value, character(1))
  } else {
    NA_character_
  }
  meta$Section <- if ("Category" %in% colnames(meta)) {
    vapply(meta$Category, format_value, character(1))
  } else {
    NA_character_
  }

  available <- c("Parameter", "Value", "Description", "Default", "Range", "List", "Section")
  bad_columns <- setdiff(columns, available)
  if (length(bad_columns) > 0) {
    stop("Unknown column(s): ", paste(bad_columns, collapse = ", "))
  }

  out <- meta[, columns, drop = FALSE]
  rownames(out) <- NULL
  if (print_table && requireNamespace("knitr", quietly = TRUE)) {
    print(knitr::kable(out, format = "markdown", row.names = FALSE))
  }
  out
}
