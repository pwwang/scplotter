#' @keywords internal
#' @importFrom utils getFromNamespace
check_columns <- getFromNamespace("check_columns", "plotthis")

#' @keywords internal
#' @importFrom utils getFromNamespace
combine_plots <- getFromNamespace("combine_plots", "plotthis")

# #' Monkey patch a function from a namespace
# #' @keywords internal
# #' @param namespace The namespace where the function is located
# #' @param function_name The name of the function to be patched
# #' @param new_function The new function to be patched
# monkey_patch <- function(namespace, function_name, new_function) {
#     ns <- getNamespace(namespace)
#     unlockBinding(function_name, ns)
#     assign(function_name, new_function, envir = ns)
#     lockBinding(function_name, ns)
# }
