library(shiny)
library(shinyWidgets)
library(bslib)
library(shinyalert)
library(yaml)
library(shinyjs)
library(DT)
library(ProteoMaker)
library(dplyr)

# Function to create UI elements based on parameter type
create_param_ui <- function(param_name, param_info) {

  group <- param_info$group
  type <- param_info$type
  class <- param_info$class
  description <- param_info$description
  min_val <- param_info$min
  max_val <- param_info$max
  default_val <- param_info$default
  ui_element <- NULL


  if (is.null(class)) {
    # error alert in shiny app
    ui_element <- tags$div(
      class = "alert alert-danger",
      "Error: Parameter class not defined for parameter ", param_name

    )
    print("No class given!!")
  }

  ## type specific UI elements
  if (class == "file") {
    # File input
    ui_element <- tagList(fileInput(param_name, label = description), textOutput(paste0(param_name, "_file")))
  } else if (class == "numeric") {
    # Numeric input
    ui_element <- numericInput(param_name, label = description, value = default_val, min = min_val, max = max_val)
  } else if (class == "boolean") {
    # Checkbox input
    ui_element <- checkboxInput(param_name, label = description, value = default_val)
  } else if (class == "ptm") {
    # Provide summary of list
    # ui_element <- tags$div(
    #   class = "alert alert-info",
    #   paste("PTM list: ", paste(param_info$default, collapse = ", "))
    # )
    # ui_element <- selectInput(param_name, label = description, choices = param_info$choices, multiple = TRUE,
    # selected = default_val)
    ui_element <- NULL
  } else if (class == "string") {
    # Select input
    ui_element <- selectInput(param_name, label = description, choices = param_info$choices, selected = default_val)
  } else {
    # State a problem
    ui_element <- tags$div(
      class = "alert alert-danger",
      "Something is not OK for this paramter: ", param_name
    )
  }

  return(ui_element)
}

# Module for setting PTMs
create_ptm_ui <- function(def_params) {
  ptm_types <- def_params$PTMTypes$choices
  ptm_dist <- def_params$PTMTypes$default
  ptm_ui <- card(tagList(
    h4("Post-Translational Modifications"),
    selectInput("ptm_types", "Select PTM", choices = ptm_types, selected = "ph"),
    # How much of this PTM?
    numericInput("ptm_fraction", "Fraction of this PTMs within all PTMs, use '1 for equal contributions.\nWill be
                 normalized to a total of one.", value = 0, min = 0, max = 1),
    # Input for distribution over all amino acids
    h5("Distribution PTMs over amino acids, will be normalized to a total sum of 1"),
    do.call(layout_column_wrap, c(
      width = 1/4,  # This will place 4 numeric inputs per row 25% width each)
      lapply(def_params$ModifiableResidues$choices, function(aa) {
        numericInput(inputId = paste0("ptm_aa_", aa), label = aa, value = 0, min = 0, max = 1)
      })
    )),
    # button to submit PTM
    layout_column_wrap(
      width = 1/2,
      actionButton("add_ptm", "Add selected PTM"),
      actionButton("del_ptm", "Delete selected PTM")
    ),
    h5("Current PTMs"),
    # Table to show current PTMs and their configuration
    htmlOutput("ptm_table")
  ))

  return(ptm_ui)

}

# Summarize ui elements for each parameter group
create_param_panels <- function(group, def_params) {
  # Filter parameters by group
  group_params <- def_params[sapply(def_params, function(x) x$group) == group]
  # Create UI elements for each parameter
  param_ui <- do.call(layout_column_wrap,
                      c(width = 1/2, heights_equal = "row",
                        lapply(names(group_params), function(param_name) {
                          create_param_ui(param_name, group_params[[param_name]])
                        })
                      )

  )
  if (group == "Ground Truth Data") {
    # Add PTM UI
    ptm_ui <- create_ptm_ui(def_params)
    param_ui <- tagList(param_ui, ptm_ui)
  }

  return(param_ui)
}

# UI Section
ui <- page_sidebar(
  # Load the audio file and play it automatically
  title = tagList(
  tags$img(src = "ProteoMakerLogo.png", height = "70px"),
  h2("ProteoMaker Simulator")
  ),
  theme = bs_theme(bootswatch = "vapor"),
  shinyjs::useShinyjs(),
  #bootstrapLib(bs_theme()),
  sidebar = sidebar(
    # title in bold
    layout_column_wrap(
      width = 1 / 2,
      h5("Parameter Groups"),
      actionButton(
        inputId = "help_button",
        label = "Details",
        icon = icon("question-circle"),  # Font Awesome question mark icon
        onclick = "window.open('Parameters.html', '_blank')",  # JavaScript to open in a new tab
        style = "border: none; background: none; color: #007bff; font-size: 16px;"
      )
    ),
    layout_column_wrap(
      width = 1 / 2,
      # Upload button to provide parameters through yaml file
      fileInput("yaml_file", "Upload YAML File"),
      # Download button to download parameters as yaml file
      downloadButton("download_yaml", "Download YAML File")
    ),
    # Capsible Panels for each group
    do.call(accordion, c(
      id = " Parameter groups",
      # Question mark button to open documentation in a new tab
      open = "Experimental Design",
      lapply(unique(sapply(def_params, function(x) x$group)),
             function(x) {
               accordion_panel(title = x, style = "info",
                               create_param_panels(x, def_params)
               )  # Adjust padding and margin to control spacing
             })
    )),
    width = 600
  ),
  card(
    h3("Run Simulation"),
    layout_column_wrap(
      width = 1 / 3,
      # Title for simulation section
      # Button to start simulation
      actionButton("start_sim", "Start Simulation"),
      checkboxInput("run_stat", "Run statistical tests on results", value = FALSE),
      tags$audio(src = "ProteoMaker Dreams.mp3", type = "audio/mp3", controls = NA,
                 style = "width: 150px; height: 30px;"), #, autoplay = NA),
      checkboxInput("overwrite", "Overwrite existing files", value = FALSE),
      checkboxInput("run_benchmarks", "Calculate benchmarking metrics", value = FALSE),
      textOutput("parameters_valid"),
    ),
    # Log window about progress of simulation
    # A text area to display the output
    card(
    verbatimTextOutput("sim_log"),
    tags$head(tags$style("#sim_log{overflow-y:scroll; max-height: 200px;}")),
    height="400px"
  )),
  card(id="results",
       h3("Get Results"),
       # List to select type of results
       layout_column_wrap(
         width = 1/4,
         heights_equal = "row",
         selectInput("result_type", "Select Simulation Stage", choices = c("GroundTruth", "ProteoformAb", "Digest", "MSRun", "DataAnalysis"),
                     selected = "MSRun"),
         selectInput("result_subtable", "Table", choices = c()),
         # Disabled button to download results
         downloadButton("download_results", "Download Results"),
         downloadButton("download_all_results", "Download All Results (zip file)"),
       ),br(),
       dataTableOutput("results_table")
  )
)
