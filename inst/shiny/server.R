options(shiny.fullstacktrace=TRUE)
# Server Section
server <- function(input, output, session) {

  toggle_dark_mode(mode = "dark", session = session)

  # General whether all parameters are valid
  all_valid <- reactiveVal(TRUE)

  # General whether simulations is available
  sim_available <- reactiveVal(FALSE)

  # Simulation results
  sim_results <- reactiveVal(NULL)

  # Table with (intermediate) results
  sim_table <- reactiveVal(NULL)

  # Collected information about PTMs
  ptm_info <- reactiveVal(NULL)

  # Currently active simulation stage
  current_stage <- reactiveVal(NULL)

  parameters <- reactiveVal(def_params)


  output$PathToFasta_file <- renderText("fasta_example.fasta")

  output$parameters_valid <- renderText(ifelse(all_valid(), "All parameters are valid", "Please correct parameter values"))

  # Function to validate parameters based on their type and constraints
  validate_param_values <- function(param_name, param_info, input_value) {
    #print(paste("Validating parameter: ", param_name, param_info$class, input_value))
    class <- param_info$class
    min_val <- param_info$min
    max_val <- param_info$max
    choices <- param_info$choices
    valid <- TRUE  # Default state to be valid
    error_message <- NULL

    # move on if now input_value
    if (is.null(input_value)) {
      return(list(valid = TRUE, message = paste(param_name, "has not been given.")))
    }

    # Check if class is defined
    if (is.null(class)) {
      return(list(valid = FALSE, message = paste("Error: Parameter not defined: ", param_name)))
    }
    # Validation based on class/type of the parameter
    if (class == "numeric") {
      # Check if numeric value is within range
      if (!is.numeric(input_value) & !is.na(input_value)) {
        valid <- FALSE
        error_message <- paste("Error: Value for", param_name, "should be numeric:", input_value)
      } else if (!is.na(input_value)) {
        if (!is.null(min_val) && input_value < min_val) {
          valid <- FALSE
          error_message <- paste("Error: Value for", param_name, "should be greater than or equal to", min_val)
        } else if (!is.null(max_val) && input_value > max_val) {
          valid <- FALSE
          error_message <- paste("Error: Value for", param_name, "should be less than or equal to", max_val)
        }
      }
    } else if (class == "file") {
      # File check: Ensure the file is not NULL (if required)

      if (!all(is.null(unlist(input_value))) & !all(is.na(unlist(input_value)))) {
        # Check if file exists
        if ("datapath" %in% names(input_value)) {
          if (file.exists(input_value$datapath)) {
            valid <- TRUE
            proteomaker_config$fastaFilePath <<- ""
            # calculate sha256 of file
            filehash <- digest::digest(file = input_value$datapath, algo = "crc32", serialize = FALSE)
            # New file name
            new_file_name <- paste0(proteomaker_config$resultFilePath, "/",
                                    basename(input_value$name), "_", filehash,
                                    ".", tools::file_ext(input_value$datapath))
            # Copy file to original name and hash, keep extension at the end
            file.copy(input_value$datapath, new_file_name, overwrite = T)
            input_value <- new_file_name
          } else {
            valid <- FALSE
            error_message <- paste("Error: File not available, please upload existing", param_name)
          }
          # Check whether file is in the inst/Proteomes folder
        } else if (file.exists(input_value)) {
          valid <- TRUE
          proteomaker_config$fastaFilePath <<- ""
        } else if (any(grepl(input_value, list.files(path = system.file("Proteomes", package = "ProteoMaker"))))) {
          valid <- TRUE
          proteomaker_config$fastaFilePath <<- system.file("Proteomes", package = "ProteoMaker")
        } else {
          valid <- FALSE
          error_message <- paste("Error: File not available, please upload existing", param_name)
        }
      } else {
        proteomaker_config$fastaFilePath <<- system.file("Proteomes", package = "ProteoMaker")
      }
      output$PathToFasta_file <- renderText(basename(input_value))
    } else if (class == "boolean") {
      # Boolean check: Should be TRUE or FALSE
      if (!is.logical(input_value)) {
        valid <- FALSE
        error_message <- paste("Error: Value for", param_name, "should be a boolean (TRUE/FALSE).")
      }
    } else if (class == "list" || class == "string") {
      # Ensure the selected value is in the list of choices (if defined)
      if (!is.null(choices) && !input_value %in% choices) {
        valid <- FALSE
        error_message <- paste("Error: Value for", param_name, "is not a valid option. Available choices:", paste(choices, collapse = ", "))
      }
    } else if (class == "ptm") {
      return(list(valid = TRUE, message = paste(param_name, "is a PTM parameter.")))
    }
    # Overwrite params value when valid
    if (valid) {
      params <- parameters()
      params[[param_name]]$value <- input_value
      parameters(params)
    }

    # Return validation result
    if (valid) {
      return(list(valid = TRUE, message = paste(param_name, "is valid.")))
    } else {
      return(list(valid = FALSE, message = error_message))
    }
  }

  # Separate validation for PTM parameters
  validate_ptm_params <- function(params, up_params) {
    print("Check and update ptm parameters")
    valid <- TRUE
    message <- NULL

    # Set default values
    params$PTMTypes$value <- NULL
    params$PTMTypesDistr$value <- NULL
    params$PTMTypesMass$value <- NULL
    params$ModifiableResidues$value <- NULL
    params$ModifiableResiduesDistr$value <- NULL


    # Check PTM parameters
    ptms <- up_params$PTMTypes$value[[1]]
    ptms <- ptms[!is.na(ptms)]
    if (is.null(ptms) || length(ptms) == 0) {
      return(list(valid = TRUE, message = "No PTMs defined."))
    }
    if (!is.null(ptms) & length(ptms) > 0) {
      # check distr
      if (!all(sort(names(up_params$PTMTypesDistr$value[[1]])) == sort(ptms))) {
        valid <- FALSE
        error_message <- "Error: PTMTypesDistr needs to be list of named vector/list containing the relative
          frequency of each PTM."
      } else if (!all(sort(names(up_params$PTMTypesMass$value[[1]])) == sort(ptms))) {
        valid <- FALSE
        error_message <- "Error: PTMTypesMass needs to be list of named vector/list containing the mass of each PTM."
      } else if (!all(sort(names(up_params$ModifiableResidues$value[[1]])) == sort(ptms))) {
        valid <- FALSE
        error_message <- "Error: ModifiableResidues needs to be list of named vector/list containing the amino acids
          that can be modified by each PTM."
      } else if (!all(sort(names(up_params$ModifiableResiduesDistr$value[[1]])) == sort(ptms))) {
        valid <- FALSE
        error_message <- "Error: ModifiableResiduesDistr needs to be list of named vector/list containing the relative
          frequency of each amino acid that can be modified by each PTM."
      } else if (!all(is.numeric(unlist(up_params$PTMTypesDistr$value)))) {
        valid <- FALSE
        error_message <- "Error: PTMTypesDistr needs to be list of named vector/list containing the relative
          frequency of each PTM."
      } else if (!all(is.numeric(unlist(up_params$PTMTypesMass$value)))) {
        valid <- FALSE
        error_message <- "Error: PTMTypesMass needs to be list of named vector/list containing the mass of each PTM."
      } else if (!all(is.character(unlist(up_params$ModifiableResidues$value)))) {
        valid <- FALSE
        error_message <- "Error: ModifiableResidues needs to be list of named vector/list containing the amino acids
          that can be modified by each PTM."
      } else if (!all(is.numeric(unlist(up_params$ModifiableResiduesDistr$value)))) {
        valid <- FALSE
        error_message <- "Error: ModifiableResiduesDistr needs to be list of named vector/list containing the relative
          frequency of each amino acid that can be modified by each PTM."
      }
    }

    # Overwrite parameter values when valid
    if (valid) {
      print("Overwriting ptms")
      params$PTMTypes$value <- up_params$PTMTypes$value
      params$PTMTypesDistr$value <- up_params$PTMTypesDistr$value
      params$PTMTypesMass$value <- up_params$PTMTypesMass$value
      params$ModifiableResidues$value <- up_params$ModifiableResidues$value
      params$ModifiableResiduesDistr$value <- up_params$ModifiableResiduesDistr$value
    }

    parameters(params)

  # return validation result
  if (valid) {
    return(list(valid = TRUE, message = paste("PTM parameters ares valid.")))
  } else {
    return(list(valid = FALSE, message = error_message))
  }
}

# set global parameters
observeEvent(input$run_stat, proteomaker_config$runStatTests <<- input$run_stat)
observeEvent(input$run_benchmarks, proteomaker_config$calcAllBenchmarks <<- input$run_benchmarks)

# enable/disable button when results are avilable
observeEvent(sim_available(), {
  print("switching")
  if (sim_available()) {
    shinyjs::enable("results")
    shinyjs::enable("result_subtable")
  } else {
    print("disabling")
    shinyjs::disable("results")
    shinyjs::disable("result_subtable")
  }
})

# Render text about PTM table
output$ptm_table <- renderText(ptm_info())


# Placeholder for displaying selected values (can be modified for actual simulation)
output$preview_text <- renderText({
  paste("Selected values: \n", paste(sapply(names(parameters()), function(param_name) {
    paste(param_name, "=", input[[param_name]])
  }), collapse = "\n"))
})

# Allow changing theme
# bs_themer()


# Provide default parameters for each PTM
observeEvent( input$ptm_types, {
  params <- parameters()
  # Provide different updates of distributions for ph, ac, me, de
  updateNumericInput(session, "ptm_fraction", value = 1)# / (1 + length(params$PTMTypes$value)))
  for(aa in params$ModifiableResidues$choices) {
    updateNumericInput(session, paste0("ptm_aa_", aa), value = 0)
  }
  if (input$ptm_types == "ph") {
    updateNumericInput(session, "ptm_aa_S", value = 0.86)
    updateNumericInput(session, "ptm_aa_T", value = 0.13)
    updateNumericInput(session, "ptm_aa_Y", value = 0.01)
  } else if (input$ptm_types == "ox") {
    updateNumericInput(session, "ptm_aa_M", value = 1)
  } else if (input$ptm_types == "ac") {
    updateNumericInput(session, "ptm_aa_K", value = 1)
  } else if (input$ptm_types == "me") {
    updateNumericInput(session, "ptm_aa_K", value = 0.75)
    updateNumericInput(session, "ptm_aa_R", value = 0.25)
  } else if (input$ptm_types == "de") {
    updateNumericInput(session, "ptm_aa_N", value = 0.5)
    updateNumericInput(session, "ptm_aa_Q", value = 0.5)
  }
})

# Add PTM config to params and update PTM details in "ptm_table" (textOutput)
observeEvent(input$add_ptm, {
  params <- parameters()
  # Add PTM to params
  ptm_name <- input$ptm_types
  if (ptm_name %in% unlist(params$PTMTypes$value)) {
    shinyalert("Error: PTM already defined.")
    return()
  }
  ptm_fraction <- input$ptm_fraction
  ptm_aa <- c()
  ptm_aa_freq <- c()
  for (aa in params$ModifiableResidues$choices) {
    tfreq <- input[[paste0("ptm_aa_", aa)]]
    if (tfreq > 0) {
      ptm_aa <- append(ptm_aa, aa)
      ptm_aa_freq <- append(ptm_aa_freq, tfreq)
    }
  }
  if (length(ptm_aa) > 0) {
    # Normalize all PTM fractions to 1
    ptm_aa_freq <- ptm_aa_freq / sum(ptm_aa_freq)

    # Update params
    ptms <- params$PTMTypes$value[[1]]
    params$PTMTypes$value <- list(mods = append(ptms, ptm_name))
    params$PTMTypesDistr$value[[1]][[ptm_name]] <- ptm_fraction
    params$PTMTypesMass$value[[1]][[ptm_name]] <- params$PTMTypesMass$choices[[which(params$PTMTypes$choices == ptm_name)]]
    params$ModifiableResidues$value[[1]][[ptm_name]] <- ptm_aa
    params$ModifiableResiduesDistr$value[[1]][[ptm_name]] <- ptm_aa_freq
    parameters(params)
    # Update PTM details in "ptm_table"
    update_ptm_info <- paste0(ptm_info(), ifelse(is.null(ptm_info()), "", "<br/>"),
                              paste0("PTM: ", ptm_name,
                                     " | Frequency ", ptm_fraction, " | Amino acids: ",
                                     paste0(sapply(seq_along(ptm_aa), function(i) paste(" ", ptm_aa[[i]], ": ",
                                                                                        ptm_aa_freq[[i]])),
                                            collapse = ",")
                              )
    )
    ptm_info(update_ptm_info)

  } else {
    shinyalert("Error: PTM must be defined for at least one amino acid.")
  }
})

# Remove PTM from params and update PTM details in "ptm_table" (textOutput)
observeEvent(input$del_ptm, {
  params <- parameters()
  # remove from params
  ptm_name <- input$ptm_types
  print(paste("Removing PTM", ptm_name))
  if (length(ptm_name) == 1) {
    ptms <- unlist(params$PTMTypes$value)
    params$PTMTypes$value <- list(mods = ptms[ptms != ptm_name])
    params$PTMTypesDistr$value[[1]][[ptm_name]] <- NULL
    params$PTMTypesMass$value[[1]][[ptm_name]] <- NULL
    params$ModifiableResidues$value[[1]][[ptm_name]] <- NULL
    params$ModifiableResiduesDistr$value[[1]][[ptm_name]] <- NULL
    parameters(params)

    # remove line of PTM from ptm_table
    t_ptm_info <- strsplit(ptm_info(), "<br/>")[[1]]
    print(t_ptm_info)
    line_to_remove <- grep(paste0("^PTM: ", ptm_name), t_ptm_info)
    print(line_to_remove)
    t_ptm_info <- t_ptm_info[-line_to_remove]

    ptm_info(paste0(t_ptm_info, collapse = "<br/>"))
  }
})

# Read yaml file using upload button
observeEvent(input$yaml_file, {
  print("Loading yaml file")
  req(input$yaml_file)
  isolate({
    params <- parameters()
    # trycatch for errors running yaml.load_file
    tryCatch({
      up_params <- yaml::yaml.load_file(input$yaml_file$datapath)
      if (is.null(up_params) || length(up_params) == 0) {
        shinyalert("The params object is empty or invalid.")
      }
      if (is.null(up_params$params)) {
        shinyalert("Yaml file main item should be 'params
                     '")
      }
      up_params <- up_params$params
      # Fill all fields
      all_valid(TRUE)
      for (up_param_name in names(up_params)) {
        # validate parameter values
        value <- up_params[[up_param_name]]$value
        checked <- validate_param_values(up_param_name, params[[up_param_name]], value)

        if (checked$valid)  {
          if (!is.null(value)) {
            if (params[[up_param_name]]$class == "file") {
              # update fileInput with file name (no path)
              output[[paste0(up_param_name,"_file")]] <- renderText({
                basename(value)
              })
              if (!is.null(value$datapath))
                value <- value$datapath
            } else if (params[[up_param_name]]$class == "boolean") {
              # update checkboxInput
              updateCheckboxInput(session, up_param_name, value = value)
            } else if (params[[up_param_name]]$class == "list" | params[[up_param_name]]$class == "string") {
              # update selectInput
              value[is.na(value)] <- "NA"
              updateSelectInput(session, up_param_name, selected = value)
            } else if (params[[up_param_name]]$class == "numeric"){
              # update textInput
              updateNumericInput(session, up_param_name, value = value)
            }
          }
        } else {
          print("Error: Invalid parameter value")
          shinyalert(checked$message)
          all_valid(FALSE)
        }
      }
      checked <- validate_ptm_params(params, up_params)
      if (checked$valid) {
        update_ptm_info <- NULL
        ptms <- up_params$PTMTypes$value[[1]]
        ptms[ptms == "NA"] <- NA
        ptms <- ptms[!is.na(ptms)]
        if (length(ptms) == 0)
          ptms <- NULL
        for (ptm_name in ptms) {
          print(paste("adding ptm:", ptm_name))
          ptm_fraction <- up_params$PTMTypesDistr$value[[1]][[ptm_name]]
          ptm_aa <- up_params$ModifiableResidues$value[[1]][[ptm_name]]
          ptm_aa_freq <- up_params$ModifiableResiduesDistr$value[[1]][[ptm_name]]
          # Update PTM details in "ptm_table"
          update_ptm_info <- paste0(ptm_info(), ifelse(is.null(ptm_info()), "", "<br/>"),
                                    paste0("PTM: ", ptm_name,
                                           " | Frequency ", ptm_fraction, " | Amino acids: ",
                                           paste0(sapply(seq_along(ptm_aa), function(i) paste(" ", ptm_aa[[i]], ": ",
                                                                                              ptm_aa_freq[[i]])),
                                                  collapse = ",")
                                    )
          )
          ptm_info(update_ptm_info)
        }
        if (!is.null(update_ptm_info))
          ptm_info(update_ptm_info)
      } else {
        print("Error: Invalid PTM parameter value")
        shinyalert(checked$message)
        all_valid(FALSE)
      }
    }, error = function(e) {
      shinyalert("Error reading the YAML file: ", e$message)
    })
  })

})

# Check and update parameters when any of the fields a changed
observe({
  all_valid(TRUE)
  params <- parameters()
  print("Check and update parameters")
  for (param_name in names(params)) {
    # validate parameter values
    value <- input[[param_name]]
    checked <- validate_param_values(param_name, params[[param_name]], value)
    if (!checked$valid) {
      print("Error: Invalid parameter value")
      shinyalert(checked$message)
      all_valid(FALSE)
    }
  }
})



# Download the parameter file
output$download_yaml <- downloadHandler(
  filename = function() {
    paste("parameters", "yaml", sep = ".")
  },
  content = function(file) {
    # Add values from the input fields to the params list
    tparams <- parameters()
    for (param_name in names(tparams)) {
      if (is.null(tparams[[param_name]]$value)) {
        tparams[[param_name]]$value <- input[[param_name]]
      }
    }
    tparams <- list(params=tparams)
    # Write the yaml file
    yaml::write_yaml(tparams, file)

  }
)

observeEvent(input$start_sim, isolate({
  message(" ---- Starting simulation ----")
  sim_available(FALSE)
  # message while running via progressbar
  progress <- shiny::Progress$new()

  progress$set(message = "Running simulation...")

  # Initialize empty lists for each category
  Param <- list(
    paramGroundTruth = list(),
    paramProteoformAb = list(),
    paramDigest = list(),
    paramMSRun = list(),
    paramDataAnalysis = list()
  )

  message("read parameters")
  params <- parameters()
  print(params$PTMTypesDistr)
  # Normalize PTM fractions to total of one
  if(!is.null(params$PTMTypesDistr$value)) {
    if (!all(is.na(params$PTMTypesDistr$value))) {
      total_fractions <- sum(unlist(params$PTMTypesDistr$value))
      for(i in 1:length(unlist(params$PTMTypesDistr$value)))
        params$PTMTypesDistr$value[[1]][[i]] <- params$PTMTypesDistr$value[[1]][[i]] / total_fractions
    }
  }


  # Iterate over each parameter and place it into the correct category based on "type"
  for (param_name in names(params)) {
    #print(param_name)
    param_info <- params[[param_name]]
    category <- param_info$type
    value <- param_info$value
    if (is.null(value))
      value <- param_info$default
    value[value == "NA"] <- NA
    if (!is.null(category)) {
      Param[[category]][[param_name]] <- value
    }
  }


  withCallingHandlers({
    shinyjs::html("sim_log", "")
    cat <- message
    print(proteomaker_config)
    res <- run_sims(Param, proteomaker_config, input$overwrite)
    sim_results(res)
    message("---- Simulation completed ----")
  },
  message = function(m) {
    shinyjs::html(id = "sim_log", html = paste0(m$message), add = TRUE)
  })


  sim_available(TRUE)

  ttable <- get_simulation(sim_results()[[1]]$Param, Config = proteomaker_config, stage = input$result_type)
  ttable <- ttable[[which(names(ttable) != "Param")]]
  enable("result_subtable")
  updateSelectInput(session, choices = names(ttable[[1]]), inputId = "result_subtable")
  ttable <- ttable[[1]]
  current_stage("DataAnalysis")
  if (input$result_subtable == "globalBMs") {
    # Create a data frame with item names and values as rows
    ttable <- do.call(rbind, lapply(names(ttable), function(name) {
      data.frame(
        Item = name,
        Value = paste(ttable[[name]], collapse = ", "),  # Make one value
        stringsAsFactors = FALSE
      )
    }))
  }
  sim_table(as.data.frame(ttable))

  on.exit(progress$close())

}))

# update sim_table from result_type
observeEvent(c(input$result_type, input$result_subtable), {
  if(sim_available()) {
    print("updating table")

    ttable <- get_simulation(sim_results()[[1]]$Param, Config = proteomaker_config, stage = input$result_type)
    ttable <- ttable[names(ttable) != "Param"]
    if (input$result_type == "DataAnalysis" & length(ttable) > 1) {
      ttable$globalBMs <- ttable$Benchmarks$globalBMs
      ttable$Benchmarks  <- NULL
    } else {
    ttable <- ttable[[1]]
    }
    if (!is.data.frame(ttable)) {
      if (!(input$result_type == current_stage())) {
        enable("result_subtable")
        updateSelectInput(session, choices = names(ttable), inputId = "result_subtable")
      }
      ttable <- ttable[[input$result_subtable]]
    } else {
      updateSelectInput(session, selected = "NA", choices = c("NA" = ""), inputId = "result_subtable")
      disable("result_subtable")
    }
    if (input$result_subtable == "globalBMs") {
      # Create a data frame with item names and values as rows
      ttable <- do.call(rbind, lapply(names(ttable), function(name) {
        data.frame(
          Item = name,
          Value = paste(ttable[[name]], collapse = ", "),  # Wrap the list in I() to keep it as a list column
          stringsAsFactors = FALSE
        )
      }))
    }
    sim_table(as.data.frame(ttable))
    current_stage(input$result_type)
  }
})

# Button to retrieve results
output$download_results <- downloadHandler(
  print("Updated Download Handler"),
  # Download file
  filename = function() {
    paste0("ProteoMaker_results_", input$result_type, ".csv")
  },
  content = function(file) {
    if (!is.null(sim_table())) {
      message(" ---- Retrieving results ----")
      # Turn lists in to strings
      out_data <- sim_table()
      out_data <- out_data %>%
        mutate(across(where(is.list),
                      \(x) vapply(x, function(el)
                        paste(el, collapse = "|"),
                        FUN.VALUE = character(1))))
      write.csv(out_data, file, row.names = FALSE)
    }
  }
)


output$results_table <- renderDataTable({
  if (!is.null(sim_table())) {
    datatable(sim_table(), options = list(pageLength = 10))
  }
})

output$download_all_results <- downloadHandler(

  filename = function() {
    paste0("ProteoMaker_all_results_", Sys.Date(), ".zip")
  },
  content = function(file) {
    if(sim_available()) {
      message(" ---- Downloading all results ----")
    # Create a temporary directory to store the results
    temp_dir <- tempdir()
    # Create a subdirectory for ProteoMaker results
    result_dir <- file.path(temp_dir, "ProteoMaker_results")
    dir.create(result_dir, showWarnings = FALSE)

    stages_all <- c("GroundTruth", "ProteoformAb",
                    "Digest", "MSRun", "DataAnalysis")

    # Save each simulation result as a CSV file
    for (i in stages_all) {
      ttable <- get_simulation(sim_results()[[1]]$Param, Config = proteomaker_config, stage = i)

      ttable <- ttable[names(ttable) != "Param"]
      if (i == "DataAnalysis" & length(ttable) > 1) {
        ttable$globalBMs <- ttable$Benchmarks$globalBMs
        ttable$globalBMs <- do.call(rbind, lapply(names(ttable$globalBMs), function(name) {
          data.frame(
            Item = name,
            Value = paste(ttable$globalBMs[[name]], collapse = ", "),  # Wrap the list in I() to keep it as a list column
            stringsAsFactors = FALSE
          )
        }))
        ttable$Benchmarks  <- NULL
      } else {
        ttable <- ttable[[1]]
      }
      if (!is.data.frame(ttable)) {
        for (n in names(ttable)) {
          print(head(ttable[[n]]))
          if (!is.null(ttable[[n]])) {
          ttable[[n]] <- ttable[[n]] %>%
            mutate(across(where(is.list),
                          \(x) vapply(x, function(el)
                            paste(el, collapse = "|"),
                            FUN.VALUE = character(1))))
          write.csv(ttable[[n]], file = file.path(result_dir, paste0("ProteoMaker_results_", i, "_", n, ".csv")), row.names = FALSE)
          }
        }
      } else {
        if (!is.null(ttable)) {

        ttable <- ttable %>%
          mutate(across(where(is.list),
                        \(x) vapply(x, function(el)
                          paste(el, collapse = "|"),
                          FUN.VALUE = character(1))))

        write.csv(ttable, file = file.path(result_dir, paste0("ProteoMaker_results_", i, ".csv")), row.names = FALSE)
        }
      }
    }

    # Create a zip file containing all results
    zip::zipr(file, files = list.files(result_dir, full.names = TRUE))
  }
  }
)

# File upload size limit
options(shiny.maxRequestSize=1*1024^3)

}

