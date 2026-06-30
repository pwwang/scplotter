#' LLM-Powered Chat Interface for Single-Cell Visualization
#'
#' @description
#' `SCPlotterChat` is an R6 class that provides a conversational interface for
#' generating single-cell data visualizations through natural language. It
#' leverages Large Language Models (LLMs) via the \pkg{tidyprompt} package to
#' interpret user requests, select appropriate visualization functions,
#' identify relevant data objects, and generate executable R code — all
#' without requiring the user to know function names, parameter details, or
#' data formats.
#'
#' @section Architecture:
#' `SCPlotterChat` orchestrates a three-stage LLM pipeline for each user
#' request:
#'
#' **Stage 1 — Tool identification** (`locate_tool`): The user's prompt is
#' analyzed against a registry of all exported \pkg{scplotter} functions
#' (built automatically from package documentation). The LLM selects the most
#' appropriate visualization function, considering the current conversation
#' history for context (e.g., "do a heatmap instead" reuses the previous
#' tool).
#'
#' **Stage 2 — Data identification** (`locate_data`): The LLM identifies the
#' target data object from available sources:
#' \itemize{
#'   \item Data frames and Seurat objects in the global environment
#'   \item Datasets exported by \pkg{scplotter}, \pkg{Seurat},
#'         \pkg{SeuratObject}, and \pkg{scRepertoire}
#'   \item Manually set data via `$set_data()`
#' }
#'
#' **Stage 3 — Code generation and execution** (`run_tool`): The LLM
#' generates R code to call the selected visualization function with the
#' identified data object. The code is evaluated in an isolated environment,
#' and the resulting plot is returned. Validation ensures the generated code
#' is syntactically valid R.
#'
#' @section Conversation history:
#' `SCPlotterChat` maintains a conversation history across calls to `$ask()`.
#' Each interaction records the user's prompt, the selected tool, the data
#' object, and the generated R code. This enables contextual follow-up
#' requests such as "change the palette to Spectral" or "facet by condition"
#' — the LLM understands that these are refinements of the previous plot
#' rather than new independent requests.
#'
#' @section Tool registry:
#' On initialization, `SCPlotterChat` scans the \pkg{scplotter} package
#' namespace for all exported functions and parses their Rd documentation to
#' build a tool registry. Each tool entry contains the function's title,
#' description, usage signature, arguments (with descriptions), and usage
#' examples. For `...` arguments that reference other package functions
#' (e.g., `[plotthis::Heatmap()]`), the documentation is recursively
#' expanded to include those functions' arguments as well.
#'
#' Two special pseudo-tools are also registered:
#' \itemize{
#'   \item `ListTools` — Lists all available visualization functions
#'   \item `ListData` — Lists all available data objects
#' }
#' When the LLM selects either of these, the chat responds with a listing
#' rather than generating a plot.
#'
#' @section LLM provider setup:
#' `SCPlotterChat` requires an LLM provider object from the \pkg{tidyprompt}
#' package. Any provider supported by \pkg{tidyprompt} can be used:
#' OpenAI, DeepSeek, Ollama (local), and others. See
#' \url{https://tjarkvandemerwe.github.io/tidyprompt/articles/getting_started.html#setup-an-llm-provider}
#' for setup instructions.
#'
#' @import R6
#' @export
#' @importFrom utils help getFromNamespace
#' @importFrom tools parse_Rd
#' @importFrom rlang %||%
#' @examples
#' \donttest{
#' # Setup LLM provider (requires an API key)
#' provider <- tidyprompt::llm_provider_openai(
#'     parameters = list(model = "deepseek-v4-flash", stream = TRUE),
#'     url = "https://api.deepseek.com/chat/completions",
#'     api_key = Sys.getenv("OPENAI_API_KEY")
#' )
#'
#' # Create chat instance
#' chat <- SCPlotterChat$new(provider)
#'
#' # List available tools
#' chat$ask("What are the tools to use?")
#' # Tool identified:  ListTools
#' # Available tools:
#' # -  CCCPlot :  Cell-Cell Communication Plot
#' #    Visualizes ligand-receptor interaction inference results ...
#' # ...
#'
#' # List available data
#' chat$ask("What data is available?")
#'
#' # Generate a plot from natural language
#' chat$ask("Plot the default cell-cell communication plot for the cellphonedb_res dataset")
#' # Tool identified:  CCCPlot
#' # Data object identified:  cellphonedb_res
#' # Code ran: CCCPlot(cellphonedb_res)
#'
#' # Refine with conversational context
#' chat$ask("do a heatmap instead")
#' # Tool identified:  CCCPlot
#' # Data object identified:  cellphonedb_res
#' # Code ran: CCCPlot(cellphonedb_res, plot_type = "heatmap")
#'
#' # Add a title
#' chat$ask("Add a proper title to the plot")
#'
#' # Manually set data to avoid auto-detection
#' chat$set_data(scplotter::cellphonedb_res)
#' chat$ask("Make a dot plot")
#'
#' # Inspect conversation history
#' chat$get_history()
#'
#' # Clear history for a fresh conversation
#' chat$clear_history()
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
                    "\n\n", "Chat History:\n",
                    paste(private$history, collapse = "\n")
                )
            } else {
                history <- ""
            }

            wrapt_prompt <- tidyprompt::prompt_wrap(
                prompt,
                modify_fn = function(base_prompt) {
                    paste0(
                        'Objective: Identify the data object to be used.\n\n',
                        'Decision Process:\n',
                        '- If the user explicitly names a data object, output that data object name. If the name does not match exactly the available ones, use the one that is mostly related.\n',
                        '- Else, if the request appears to refine the previous output (e.g., "do X instead", "make it a heatmap/dot/bar/etc", "change to ...", "add ... to", "same plot but ..."), select the last dataset name mentioned in the chat history.\n',
                        '- Else, analyze the current user request; if there is a clear, unambiguous match to a dataset, select that data object.\n',
                        '- Else, use the last mentioned data object from the chat history.\n',
                        '- If no data object is found, respond with "None".\n\n',
                        'Response Format: Provide only the name of the selected data object, or "None" if no tool applies. The response should be unquoted.\n\n',
                        'User Request:\n',
                        base_prompt, "\n\n",
                        "Available Data Objects:\n",
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
                    return(TRUE)
                }
            )

            tidyprompt::send_prompt(wrapt_prompt, private$provider, max_interactions = 3, verbose = verbose)
        },

        locate_tool = function(prompt, verbose) {
            if (!is.null(private$history)) {
                history <- paste0(
                    "\n\n", "Chat History:\n",
                    paste(private$history, collapse = "\n")
                )
            } else {
                history <- ""
            }
            wrapt_prompt <- tidyprompt::prompt_wrap(
                prompt,
                modify_fn = function(base_prompt) {
                    paste0(
                        'Objective: Select the most appropriate tool to handle the user\'s request while preserving conversational context. If the user refines or changes how the previous result should be visualized (e.g., asks for a different plot type), continue with the last plotting tool used unless they explicitly name a different tool.\n\n',
                        'Decision Process:\n',
                        '- If the user explicitly names a tool, output that tool.\n',
                        '- Else, if the request appears to refine the previous output (e.g., "do X instead", "make it a heatmap/dot/bar/etc", "change to ...", "add ... to", "same plot but ..."), select the last tool mentioned in the chat history.\n',
                        '- Else, analyze the current user request; if there is a clear, unambiguous match to a single tool, select that tool.\n',
                        '- Else, use the last mentioned tool from the chat history.\n',
                        '- If no tool is found, respond with "None".\n\n',
                        'Response Format: Provide only the name of the selected tool, or "None" if no tool applies.\n\n',
                        'User Request:\n',
                        base_prompt, '\n\n',
                        'Available Tools:\n',
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
                    "\n\n", "Chat History:\n",
                    paste(private$history, collapse = "\n")
                )
            } else {
                history <- ""
            }

            wrapt_prompt <- tidyprompt::prompt_wrap(
                prompt,
                modify_fn = function(base_prompt) {
                    paste0(
                        'Objective: Generate the R code to run the specified tool using the specified data object based on the user\'s request and the provided tool information.\n\n',
                        'Decision Process:\n',
                        '- Analyze the user\'s request to understand what needs to be done with the specified tool and data object.\n',
                        '- When use the data object name, do not quote it.\n',
                        '- Refer to the provided tool information to understand the tool\'s usage, arguments, and examples.\n',
                        '- Consider any specific parameters or options mentioned in the user\'s request that need to be included in the tool call.\n',
                        '- If the user\'s request is ambiguous or lacks necessary details, refer to the chat history for additional context that may help clarify the intended use of the tool and data object.\n',
                        '- Construct the R code that correctly calls the tool with the specified data object, ensuring that all necessary arguments are included and correctly formatted.\n\n',
                        'Response Format: Provide only the valid R code wrapped between "```r" and "```" to run the tool.\n\n',
                        'Tool to be used: ', tool_name, '\n',
                        'Data object to be used: ', sub("^[a-zA-Z0-9._]+::", "", data_name), '\n\n',
                        "User Request:\n",
                        base_prompt, "\n\n",
                        "Tool Information:\n",
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

        #' @description Create a new instance of the SCPlotterChat class.
        #'   On initialization, the chat scans the \pkg{scplotter} package for all
        #'   exported visualization functions and builds a tool registry from their
        #'   documentation. It also discovers available data objects from the global
        #'   environment and from packages (\pkg{scplotter}, \pkg{Seurat},
        #'   \pkg{SeuratObject}, \pkg{scRepertoire}).
        #' @param provider An LLM provider object from the \pkg{tidyprompt} package
        #'   (e.g., created via \code{tidyprompt::llm_provider_openai()}). See
        #'   \url{https://tjarkvandemerwe.github.io/tidyprompt/articles/getting_started.html#setup-an-llm-provider}
        #'   for a full list of supported providers and setup instructions.
        #' @param verbose Logical; if `TRUE`, prints the full prompt and LLM
        #'   response for each interaction. Useful for debugging prompt engineering
        #'   or understanding how the LLM interprets requests. Default is `FALSE`.
        #'   Note: this overrides the `verbose` setting of the provider object.
        #' @return A new `SCPlotterChat` instance (R6 object).
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

        #' @description Clear the conversation history. Useful for starting a
        #'   fresh conversation without creating a new `SCPlotterChat` instance.
        #'   After clearing, the LLM will have no memory of previous requests.
        #' @return `NULL` (invisibly).
        clear_history = function() {
            private$history <- NULL
        },

        #' @description Retrieve the conversation history. Each entry records
        #'   the user prompt and the assistant's response (tool used, data object,
        #'   and generated R code). The history provides conversational context
        #'   for follow-up requests.
        #' @return A character vector of history entries, or `NULL` if no
        #'   interactions have occurred or the history has been cleared.
        get_history = function() {
            if (is.null(private$history)) {
                return(NULL)
            }
            private$history
        },

        #' @description Print a list of all available visualization tools
        #'   (exported \pkg{scplotter} functions) that the LLM can use. Each
        #'   entry shows the function name, its title, and a brief description
        #'   extracted from the package documentation.
        #' @return `NULL` (invisibly). Output is printed to the console.
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

        #' @description Print a list of all available data objects that the LLM
        #'   can use for visualization. Data sources include: objects in the
        #'   global environment (data frames and Seurat objects), and datasets
        #'   exported by \pkg{scplotter}, \pkg{Seurat}, \pkg{SeuratObject}, and
        #'   \pkg{scRepertoire}.
        #' @return `NULL` (invisibly). Output is printed to the console.
        list_data = function() {
            if (length(private$datasets) == 0) {
                private$get_datasets()
            }
            cat("Available data objects:\n")
            for (data in names(private$datasets)) {
                cat("- ", data, ": ", private$datasets[[data]], "\n")
            }
        },

        #' @description Manually set the data object to be used for
        #'   visualization. This bypasses the automatic data detection in
        #'   `$ask()` — when data is set via this method, all subsequent
        #'   `$ask()` calls will use this data regardless of what the LLM
        #'   identifies from the prompt. Set to `NULL` to restore automatic
        #'   detection.
        #' @param data The data object (e.g., a Seurat object, data frame, or
        #'   any object compatible with \pkg{scplotter} visualization functions).
        #'   Pass `NULL` to clear the preset data and revert to automatic
        #'   detection.
        #' @param name Optional character string to use as the data object's
        #'   name in generated code. If `NULL` (default), the name is inferred
        #'   from the `data` argument via `deparse(substitute())`.
        #' @return `NULL` (invisibly).
        set_data = function(data, name = NULL) {
            if (is.null(name)) {
                private$data_name <- deparse(substitute(data))
            } else {
                private$data_name <- name
            }
            private$data <- data
        },

        #' @description Retrieve the currently set data object (if any).
        #'   Returns `NULL` if no data has been manually set via `$set_data()`.
        #' @return The data object set via `$set_data()`, or `NULL` if
        #'   automatic data detection is active.
        get_data = function() {
            private$data
        },

        #' @description Send a natural language prompt to the chat interface.
        #'   This is the primary method for interacting with `SCPlotterChat`.
        #'   It executes a three-stage pipeline:
        #'   \enumerate{
        #'     \item **Tool identification** — The LLM selects the most
        #'           appropriate \pkg{scplotter} visualization function based
        #'           on the prompt and conversation history.
        #'     \item **Data identification** — The LLM identifies the target
        #'           data object from available sources (unless data has been
        #'           manually set via `$set_data()`).
        #'     \item **Code generation and execution** — The LLM generates
        #'           R code to call the selected function with the identified
        #'           data, which is then evaluated. The resulting plot is
        #'           returned.
        #'   }
        #'   If the selected tool is `ListTools` or `ListData`, the
        #'   respective listing is printed instead of generating a plot.
        #' @param prompt A character string containing the user's query or
        #'   instruction in natural language (e.g., "Plot a UMAP of the
        #'   pancreas data colored by cell type", "Make it a heatmap instead").
        #' @param verbose Logical; if `TRUE`, prints the full LLM prompt and
        #'   response for this interaction. Default is `NULL`, which falls
        #'   back to the `verbose` setting of the `SCPlotterChat` instance.
        #'   Use this to debug a single interaction without enabling verbose
        #'   mode globally.
        #' @param add_to_history Logical; if `TRUE` (default), this
        #'   interaction is recorded in the conversation history, enabling
        #'   the LLM to understand follow-up requests in context. Set to
        #'   `FALSE` for one-off queries that should not influence subsequent
        #'   interactions.
        #' @return The generated plot (typically a `ggplot` object), or
        #'   `NULL` invisibly for `ListTools`/`ListData` requests.
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
