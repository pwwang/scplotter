#' SCPlotterChat Class
#'
#' An R6 class that provides chat functionality for SCPlotter
#'
#' @import R6
#' @export
#' @importFrom utils help getFromNamespace
#' @importFrom tools parse_Rd
#' @importFrom rlang %||%
#' @examples
#' \donttest{
#' if (FALSE) {
#' provider <- tidyprompt::llm_provider_openai(api_key = Sys.getenv("OPENAI_API_KEY"))
#' chat <- SCPlotterChat$new(provider)
#' chat$ask("What are the tools to use?")
#' # Tool identified:  ListTools
#' # Available tools:
#' # -  ClonalOverlapPlot :  ClonalOverlapPlot
#' #    Plot the overlap of the clones in different samples/groups.
#' #
#' # -  EnrichmentPlot :  Enrichment Plot
#' #    This function generates various types of plots for enrichment analysis.
#' # ...
#'
#' chat$ask("Plot the default cell-cell communication plot for the cellphonedb_res dataset")
#' # Tool identified:  CCCPlot
#' # Data object identified:  cellphonedb_res
#' # Running tool:  CCCPlot
#'
#' chat$ask("do a heatmap instead")
#' # Tool identified:  CCCPlot
#' # Data object identified:  cellphonedb_res
#' # Running tool:  CCCPlot
#' }
#' }
SCPlotterChat <- R6::R6Class(
    "SCPlotterChat",
    private = list(
        provider = NULL,
        this_package = "scplotter",
        tools = list(),
        datasets = list(),
        data_name = NULL,
        data = NULL,
        verbose = FALSE,
        history = c(),

        get_tools = function() {
            if (length(private$tools) > 0) {
                return(private$tools)
            }
            # Get all functions in the package
            all_funcs <- getNamespaceExports(private$this_package)
            all_funcs <- setdiff(all_funcs, class(self))

            get_help_rd <- getFromNamespace(".getHelpFile", "utils")
            # Get title and description of all functions from their help pages
            for (func_name in all_funcs) {
                rd <- get_help_rd(help(func_name, package = private$this_package))
                for (el in rd) {
                    if (attr(el, "Rd_tag") == "\\title") {
                        title <- el
                    } else if (attr(el, "Rd_tag") == "\\description") {
                        description <- el
                    }
                }
                private$tools[[func_name]] <- c(
                    title = paste(unlist(title), collapse = " "),
                    description = paste(sapply(unlist(description), function(x) paste0("  ", x)), collapse = " ")
                )
            }

            private$tools$ListTools <- c(
                title = "List all available tools\n",
                description = "List all available tools that can be used to handle the chat request.\n"
            )

            private$tools$ListData <- c(
                title = "List all available data objects\n",
                description = "List all available data objects that can be used to handle the chat request.\n"
            )
        },

        get_datasets = function() {
            if (length(private$datasets) > 0) {
                return(private$datasets)
            }
            # Data in global Environment
            data_objects <- ls(envir = .GlobalEnv)
            for (data_object in data_objects) {
                if (!inherits(get(data_object, envir = .GlobalEnv), c("data.frame", "Seurat"))) {
                    next
                }
                private$datasets[[data_object]] <- paste(class(get(data_object, envir = .GlobalEnv)), collapse = ", ")
            }

            for (pkg in c(private$this_package, "Seurat", "SeuratObject", "scRepertoire")) {
                pkg_data <- utils::data(package = pkg)
                for (i in 1:nrow(pkg_data$results)) {
                    private$datasets[[paste0(pkg_data$results[, "Package"][i], "::", pkg_data$results[, "Item"][i])]] <- pkg_data$results[, "Title"][i]
                }
            }
        },

        get_tool_info = function(tool_name, package = NULL, exclude_arguments = c()) {
            get_help_rd <- getFromNamespace(".getHelpFile", "utils")
            rd <- get_help_rd(help(tool_name, package = if(is.null(package)) private$this_package else package))

            out <- list()
            # datatype <- NULL
            for (el in rd) {
                if (attr(el, "Rd_tag") == "\\title") {
                    out$title <- paste(unlist(el), collapse = " ")
                } else if (attr(el, "Rd_tag") == "\\description") {
                    out$description <- paste(unlist(el), collapse = " ")
                } else if (attr(el, "Rd_tag") == "\\usage") {
                    out$usage <- paste(unlist(el), collapse = "")
                } else if (attr(el, "Rd_tag") == "\\examples") {
                    out$examples <- paste(unlist(el), collapse = "")
                } else if (attr(el, "Rd_tag") == "\\value") {
                    out$value <- paste(unlist(el), collapse = "")
                } else if (attr(el, "Rd_tag") == "\\arguments") {
                    out$arguments <- list()
                    arg_i <- 0
                    for (argument in el) {
                        if (attr(argument, "Rd_tag") == "\\item") {
                            arg_i <- arg_i + 1
                            arg_name <- trimws(argument[[1]])
                            if (arg_name %in% exclude_arguments) {
                                next
                            }
                            if (identical(arg_name, "...")) {
                                # Extract function references from the description
                                desc_text <- paste(unlist(argument[[-1]]), collapse = "")
                                # Find all references to package functions
                                pkg_funcs <- list()

                                # Extract function references using regex
                                # Handle pkg::func
                                func_matches <- gregexpr("[a-zA-Z0-9._]+::[a-zA-Z0-9._]+", desc_text)
                                if (length(func_matches[[1]]) > 1 || func_matches[[1]][1] > 0) {
                                    func_names <- regmatches(desc_text, func_matches)[[1]]

                                    for (func_name in func_names) {
                                        parts <- strsplit(func_name, "::", fixed = TRUE)[[1]]
                                        pkg <- parts[1]
                                        func <- parts[2] # Remove () if present

                                        pkg_funcs[[func_name]] <- private$get_tool_info(
                                            func, package = pkg, exclude_arguments = c(names(out$arguments), "data", "...")
                                        )
                                    }
                                }

                                # Add package function arguments as nested info
                                out$arguments[["..."]] <- desc_text

                                # Insert package function info after where it is mentioned
                                if (length(pkg_funcs) > 0) {
                                    expanded_desc <- desc_text
                                    for (func_ref in names(pkg_funcs)) {
                                        # Find the position where function is mentioned
                                        func_pos <- gregexpr(func_ref, expanded_desc)[[1]]

                                        if (func_pos[1] > 0) {
                                            func_end <- func_pos + attr(func_pos, "match.length")[1]
                                            # Find the \n position after the function reference
                                            newline <- regexpr("\n", substring(expanded_desc, func_end))
                                            if (newline[1] > 0) {
                                                func_end <- func_end + newline[1]
                                            } else {
                                                func_end <- nchar(expanded_desc)
                                            }
                                            # Format the function info with indentation
                                            func_info <- pkg_funcs[[func_ref]]
                                            info_text <- ""
                                            if (!is.null(func_info$arguments)) {
                                                info_text <- paste0(info_text, "\n    [...] can be:\n    ")
                                                info_text <- paste0(info_text, gsub("\n", "\n    ", func_info$arguments, fixed = TRUE), "\n")
                                            }

                                            expanded_desc <- paste0(
                                                substr(expanded_desc, 1, func_end),
                                                info_text,
                                                substr(expanded_desc, func_end + 1, nchar(expanded_desc))
                                            )
                                        }
                                    }
                                    out$arguments[["..."]] <- expanded_desc
                                }
                            } else {
                                out$arguments[[arg_name]] <- paste(unlist(argument[[-1]]), collapse = "")
                                # Extract the argument type, which is specified at the beginning of the description
                                # e.g. "data: (data.frame) A data frame or matrix."
                                # if (grepl("^\\([^)]+\\)", argument[[2]][1])) {
                                #     datatype <- gsub("^\\(([^)]+)\\).*", "\\1", argument[[2]][1])
                                # }
                            }
                        }
                    }
                    out$arguments <- paste(sapply(names(out$arguments), function(arg) {
                        paste0("- ", arg, ": ", out$arguments[[arg]])
                    }), collapse = "\n")
                }
            }

            # attr(out, "datatype") <- datatype
            out
        },

        locate_data = function(prompt, datatype, verbose) {
            if (!is.null(private$history)) {
                history <- paste0(
                    "\n\n", "--- Chat History ---\n",
                    paste(private$history, collapse = "\n")
                )
            } else {
                history <- ""
            }

            wrapt_prompt <- tidyprompt::prompt_wrap(
                prompt,
                modify_fn = function(base_prompt) {
                    paste0(
                        "Based on the following prompt, identify the name of the data object that is mentioned in the prompt.\n",
                        "The name should be one of the available ones listed.\n",
                        "If no data object is found based on the prompt, use the last mentioned data object in the chat history.\n",
                        "If no data object is found in the chat history, just answer \"None\".\n\n",
                        "--- Prompt ---\n",
                        base_prompt, "\n\n",
                        "-- Available Data Objects ---\n",
                        paste(
                            sapply(names(private$datasets), function(data) {
                                paste0("- ", data, ": ", private$datasets[[data]])
                            }),
                            collapse = "\n"
                        ),
                        history
                    )
                },
                validation_fn = function(response) {
                    # Check if the response is a valid R variable name
                    if (!grepl("^[a-zA-Z.][a-zA-Z0-9_.:]*$", response)) {
                        return(tidyprompt::llm_feedback("Provide a valid R variable name or namespace::name format."))
                    }
                    if (identical(response, "None")) {
                        if (
                            length(private$history) > 0 &&
                            grepl("^Assistant: tool - \\w+; data - [a-zA-Z0-9._:]+;", private$history[length(private$history)])
                        ) {
                            return(tidyprompt::llm_feedback("There are data object names in the chat history. Provide a valid data object name."))
                        }
                        return(TRUE)
                    }
                    if (!response %in% names(private$datasets)) {
                        return(tidyprompt::llm_feedback("The data object is not available. It is not in the list of available ones."))
                    }
                    if (is.null(datatype)) {
                        return(TRUE)
                    }
                    # Check if the object is of the expected type
                    # obj <- get(response, envir = .GlobalEnv)
                    # if (!inherits(obj, datatype)) {
                    #     return(tidyprompt::llm_feedback(
                    #         paste0(
                    #             "The data object is not of the expected type: ", datatype, ".\n",
                    #             "Please provide a valid data object name."
                    #         )
                    #     ))
                    # }
                    # return(TRUE)
                }
            )

            tidyprompt::send_prompt(wrapt_prompt, private$provider, max_interactions = 3, verbose = verbose)
        },

        locate_tool = function(prompt, verbose) {
            if (!is.null(private$history)) {
                history <- paste0(
                    "\n\n", "--- Chat History ---\n",
                    paste(private$history, collapse = "\n")
                )
            } else {
                history <- ""
            }
            wrapt_prompt <- tidyprompt::prompt_wrap(
                prompt,
                modify_fn = function(base_prompt) {
                    paste0(
                        "Based on the following prompt and identify the tool that can be used to handle the request.\n",
                        "Please only answer with the name of the tool from the listed available ones. \n",
                        "If no proper tool is identified from the prompt, use the last mentioned tool in the chat history.\n",
                        "If no tool is found in the chat history, just answer \"None\".\n\n",
                        "--- Prompt ---\n",
                        base_prompt, "\n\n",
                        "--- Available Tools ---\n",
                        paste(
                            sapply(names(private$tools), function(tool) {
                                paste0(
                                    "- ", tool, ": ", private$tools[[tool]]["title"],
                                    private$tools[[tool]]["description"]
                                )
                            }),
                            collapse = "\n"
                        ),
                        history
                    )
                },
                validation_fn = function(response) {
                    if (
                        identical(response, "None") &&
                        length(private$history) > 0 &&
                        grepl("^Assistant: tool - \\w+;", private$history[length(private$history)])
                    ) {
                        return(tidyprompt::llm_feedback("There are tool names in the chat history. Provide a valid tool name."))
                    }
                    if (identical(response, "None") || response %in% names(private$tools)) {
                        return(TRUE)
                    }
                    return(tidyprompt::llm_feedback("Provide a valid tool name from the available ones."))
                }
            )

            tidyprompt::send_prompt(wrapt_prompt, private$provider, max_interactions = 3, verbose = verbose)
        },

        run_tool = function(prompt, tool_name, data_name, tool_info, verbose, add_to_history) {
            if (!is.null(private$history)) {
                history <- paste0(
                    "\n\n", "--- Chat History ---\n",
                    paste(private$history, collapse = "\n")
                )
            } else {
                history <- ""
            }

            wrapt_prompt <- tidyprompt::prompt_wrap(
                prompt,
                modify_fn = function(base_prompt) {
                    paste0(
                        "Based on the following prompt and the given tool information, generate the code to run the tool.\n",
                        "The tool or function to be used is: ", tool_name, ". The data object to be used is: ",
                        # remove the prefix package:: from the data name, because the object passed to object_to_use will be used
                        sub("^[a-zA-Z0-9._]+::", "", data_name), ".\n",
                        "Don't quote the data name when using it. The code should be valid R code.\n",
                        "Only answer with the code that is wrapped between between ```r and ``` to run the tool.\n",
                        "If there is not enough information in the prompt to run the tool, also refer to the chat history.\n\n",
                        "--- Prompt ---\n",
                        base_prompt, "\n\n",
                        "--- Tool Information ---\n",
                        paste(
                            sapply(names(tool_info), function(sec) {
                                paste0("- ", sec, "\n  ", gsub("\n", "\n  ", tool_info[[sec]], fixed = TRUE), "\n")
                            }),
                            collapse = ""
                        ),
                        history
                    )
                },
                validation_fn = function(response) {
                    # Check if the response is valid R code
                    if (grepl("^```r", response) && grepl("```$", trimws(response))) {
                        if (add_to_history) {
                            private$history <<- c(private$history, paste0("User: ", prompt))
                            private$history <<- c(
                                private$history,
                                paste0(
                                    "Assistant: tool - ", tool_name,
                                    "; data - ", data_name,
                                    "; code - ", trimws(gsub("^```r|```$", "", trimws(response)))
                                )
                            )
                        }
                        return(TRUE)
                    }
                    return(tidyprompt::llm_feedback("Provide valid R code wrapped between ```r and ``` to run the tool."))
                }
            )

            env <- new.env()
            if (is.null(private$data)) {
                if (grepl("::", data_name)) {
                    parts <- strsplit(data_name, "::", fixed = TRUE)[[1]]
                    pkg <- parts[1]
                    dn <- parts[2]
                    data(list = dn, package = pkg, envir = env)
                } else {
                    data(list = data_name, envir = env)
                }
            } else {
                # Use the data set that is set in the chat
                env[[private$data_name]] <- private$data
            }

            # We only need to return the plot object, not display it
            # Capture current device and open a null device to suppress plotting
            old_dev <- dev.cur()
            pdf(NULL)
            on.exit({
                dev.off()
                if (old_dev > 1) dev.set(old_dev)
            }, add = TRUE)

            wrapt_prompt <- tidyprompt::answer_using_r(
                wrapt_prompt,
                pkgs_to_use = "scplotter",
                objects_to_use = as.list(env),
                evaluate_code = TRUE,
                return_mode = "full"
            )

            tidyprompt::send_prompt(wrapt_prompt, private$provider, verbose = verbose)
        }
    ),

    public = list(

        #' @description Create a new instance of the SCPlotterChat class
        #' @param provider An LLM provider object
        #' See https://tjarkvandemerwe.github.io/tidyprompt/articles/getting_started.html#setup-an-llm-provider
        #' for more details on LLM providers
        #' @param verbose A logical value indicating whether to print verbose messages
        #' Default is FALSE. And this will override the verbose setting of the provider
        #' @return A new instance of the SCPlotterChat class
        initialize = function(provider, verbose = FALSE) {
            if (!inherits(provider, "LlmProvider")) {
                stop("The provider must be an instance of llm_provider-class from tidyprompt.")
            }
            private$provider <- provider
            private$verbose <- verbose
            private$get_tools()
            private$get_datasets()

            private$provider$verbose <- verbose
        },

        #' @description Clear the chat history
        #' @return NULL
        clear_history = function() {
            private$history <- NULL
        },

        #' @description Get the chat history
        #' @return The chat history
        get_history = function() {
            if (is.null(private$history)) {
                return(NULL)
            }
            private$history
        },

        #' @description Print the list of available tools
        #' @return NULL
        list_tools = function() {
            if (length(private$tools) == 0) {
                private$get_tools()
            }
            cat("Available tools:\n")
            for (tool in names(private$tools)) {
                cat("- ", tool, ": ", private$tools[[tool]]["title"])
                cat("  ", private$tools[[tool]]["description"])
            }
        },

        #' @description Print the list of available data objects that can be used
        #' @return NULL
        list_data = function() {
            if (length(private$datasets) == 0) {
                private$get_datasets()
            }
            cat("Available data objects:\n")
            for (data in names(private$datasets)) {
                cat("- ", data, ": ", private$datasets[[data]], "\n")
            }
        },

        #' @description Set the data to be analyzed
        #' @param data The data object to be set
        #' @param name The name of the data object (optional)
        #' @return NULL
        set_data = function(data, name = NULL) {
            if (is.null(name)) {
                private$data_name <- deparse(substitute(data))
            } else {
                private$data_name <- name
            }
            private$data <- data
        },

        #' @description Get the data to be analyzed
        #' @return The data object
        get_data = function() {
            private$data
        },

        #' @description Send a prompt to the chat interface and receive a response
        #' @param prompt A character string containing the user's query or instruction
        #' @param verbose A logical value indicating whether to print verbose messages
        #' Default is NULL, which will use the verbose setting of the SCPlotterChat object
        #' @param add_to_history A logical value indicating whether to add the prompt and response to the chat history
        #' Default is TRUE
        #' @return A response from the chat system
        ask = function(prompt, verbose = NULL, add_to_history = TRUE) {
            verbose <- verbose %||% private$verbose
            tool_name <- private$locate_tool(prompt, verbose)

            if (is.null(tool_name) || identical(tool_name, "None")) {
                stop("There is no tool that can handle the request. Please modify your request and try again.")
            }
            cat("Tool identified: ", tool_name, "\n")

            if (identical(tool_name, "ListTools")) {
                self$list_tools()
                return(invisible())
            }

            if (identical(tool_name, "ListData")) {
                self$list_data()
                return(invisible())
            }

            tool_info <- private$get_tool_info(tool_name)
            if (is.null(private$data)) {
                data_name <- private$locate_data(prompt, attr(tool_info, "datatype"), verbose)
                if (is.null(data_name) || identical(data_name, "None")) {
                    stop("No data set can be located from your request. Try again with a revised request or set the data manually with `chat$set_data(...)`.")
                }
                cat("Data object identified: ", data_name, "\n")
            } else {
                data_name <- private$data_name
                cat("Using preset data: ", data_name, "\n")
            }

            result <- private$run_tool(prompt, tool_name, data_name, tool_info,
                verbose = verbose, add_to_history = add_to_history
            )

            cat("Code ran:\n")
            cat(as.character(result$code), "\n")
            return(result$output$result)
        }
    )
)
