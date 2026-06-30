# LLM-Powered Chat Interface for Single-Cell Visualization

`SCPlotterChat` is an R6 class that provides a conversational interface
for generating single-cell data visualizations through natural language.
It leverages Large Language Models (LLMs) via the tidyprompt package to
interpret user requests, select appropriate visualization functions,
identify relevant data objects, and generate executable R code — all
without requiring the user to know function names, parameter details, or
data formats.

## Architecture

`SCPlotterChat` orchestrates a three-stage LLM pipeline for each user
request:

**Stage 1 — Tool identification** (`locate_tool`): The user's prompt is
analyzed against a registry of all exported scplotter functions (built
automatically from package documentation). The LLM selects the most
appropriate visualization function, considering the current conversation
history for context (e.g., "do a heatmap instead" reuses the previous
tool).

**Stage 2 — Data identification** (`locate_data`): The LLM identifies
the target data object from available sources:

- Data frames and Seurat objects in the global environment

- Datasets exported by scplotter, Seurat, SeuratObject, and scRepertoire

- Manually set data via `$set_data()`

**Stage 3 — Code generation and execution** (`run_tool`): The LLM
generates R code to call the selected visualization function with the
identified data object. The code is evaluated in an isolated
environment, and the resulting plot is returned. Validation ensures the
generated code is syntactically valid R.

## Conversation history

`SCPlotterChat` maintains a conversation history across calls to
`$ask()`. Each interaction records the user's prompt, the selected tool,
the data object, and the generated R code. This enables contextual
follow-up requests such as "change the palette to Spectral" or "facet by
condition" — the LLM understands that these are refinements of the
previous plot rather than new independent requests.

## Tool registry

On initialization, `SCPlotterChat` scans the scplotter package namespace
for all exported functions and parses their Rd documentation to build a
tool registry. Each tool entry contains the function's title,
description, usage signature, arguments (with descriptions), and usage
examples. For `...` arguments that reference other package functions
(e.g., `[plotthis::Heatmap()]`), the documentation is recursively
expanded to include those functions' arguments as well.

Two special pseudo-tools are also registered:

- `ListTools` — Lists all available visualization functions

- `ListData` — Lists all available data objects

When the LLM selects either of these, the chat responds with a listing
rather than generating a plot.

## LLM provider setup

`SCPlotterChat` requires an LLM provider object from the tidyprompt
package. Any provider supported by tidyprompt can be used: OpenAI,
DeepSeek, Ollama (local), and others. See
<https://tjarkvandemerwe.github.io/tidyprompt/articles/getting_started.html#setup-an-llm-provider>
for setup instructions.

## Methods

### Public methods

- [`SCPlotterChat$new()`](#method-SCPlotterChat-initialize)

- [`SCPlotterChat$clear_history()`](#method-SCPlotterChat-clear_history)

- [`SCPlotterChat$get_history()`](#method-SCPlotterChat-get_history)

- [`SCPlotterChat$list_tools()`](#method-SCPlotterChat-list_tools)

- [`SCPlotterChat$list_data()`](#method-SCPlotterChat-list_data)

- [`SCPlotterChat$set_data()`](#method-SCPlotterChat-set_data)

- [`SCPlotterChat$get_data()`](#method-SCPlotterChat-get_data)

- [`SCPlotterChat$ask()`](#method-SCPlotterChat-ask)

- [`SCPlotterChat$clone()`](#method-SCPlotterChat-clone)

------------------------------------------------------------------------

### `SCPlotterChat$new()`

Create a new instance of the SCPlotterChat class. On initialization, the
chat scans the scplotter package for all exported visualization
functions and builds a tool registry from their documentation. It also
discovers available data objects from the global environment and from
packages (scplotter, Seurat, SeuratObject, scRepertoire).

#### Usage

    SCPlotterChat$new(provider, verbose = FALSE)

#### Arguments

- `provider`:

  An LLM provider object from the tidyprompt package (e.g., created via
  [`tidyprompt::llm_provider_openai()`](https://kennispunttwente.github.io/tidyprompt/reference/llm_provider_openai.html)).
  See
  <https://tjarkvandemerwe.github.io/tidyprompt/articles/getting_started.html#setup-an-llm-provider>
  for a full list of supported providers and setup instructions.

- `verbose`:

  Logical; if `TRUE`, prints the full prompt and LLM response for each
  interaction. Useful for debugging prompt engineering or understanding
  how the LLM interprets requests. Default is `FALSE`. Note: this
  overrides the `verbose` setting of the provider object.

#### Returns

A new `SCPlotterChat` instance (R6 object).

------------------------------------------------------------------------

### `SCPlotterChat$clear_history()`

Clear the conversation history. Useful for starting a fresh conversation
without creating a new `SCPlotterChat` instance. After clearing, the LLM
will have no memory of previous requests.

#### Usage

    SCPlotterChat$clear_history()

#### Returns

`NULL` (invisibly).

------------------------------------------------------------------------

### `SCPlotterChat$get_history()`

Retrieve the conversation history. Each entry records the user prompt
and the assistant's response (tool used, data object, and generated R
code). The history provides conversational context for follow-up
requests.

#### Usage

    SCPlotterChat$get_history()

#### Returns

A character vector of history entries, or `NULL` if no interactions have
occurred or the history has been cleared.

------------------------------------------------------------------------

### `SCPlotterChat$list_tools()`

Print a list of all available visualization tools (exported scplotter
functions) that the LLM can use. Each entry shows the function name, its
title, and a brief description extracted from the package documentation.

#### Usage

    SCPlotterChat$list_tools()

#### Returns

`NULL` (invisibly). Output is printed to the console.

------------------------------------------------------------------------

### `SCPlotterChat$list_data()`

Print a list of all available data objects that the LLM can use for
visualization. Data sources include: objects in the global environment
(data frames and Seurat objects), and datasets exported by scplotter,
Seurat, SeuratObject, and scRepertoire.

#### Usage

    SCPlotterChat$list_data()

#### Returns

`NULL` (invisibly). Output is printed to the console.

------------------------------------------------------------------------

### `SCPlotterChat$set_data()`

Manually set the data object to be used for visualization. This bypasses
the automatic data detection in `$ask()` — when data is set via this
method, all subsequent `$ask()` calls will use this data regardless of
what the LLM identifies from the prompt. Set to `NULL` to restore
automatic detection.

#### Usage

    SCPlotterChat$set_data(data, name = NULL)

#### Arguments

- `data`:

  The data object (e.g., a Seurat object, data frame, or any object
  compatible with scplotter visualization functions). Pass `NULL` to
  clear the preset data and revert to automatic detection.

- `name`:

  Optional character string to use as the data object's name in
  generated code. If `NULL` (default), the name is inferred from the
  `data` argument via `deparse(substitute())`.

#### Returns

`NULL` (invisibly).

------------------------------------------------------------------------

### `SCPlotterChat$get_data()`

Retrieve the currently set data object (if any). Returns `NULL` if no
data has been manually set via `$set_data()`.

#### Usage

    SCPlotterChat$get_data()

#### Returns

The data object set via `$set_data()`, or `NULL` if automatic data
detection is active.

------------------------------------------------------------------------

### `SCPlotterChat$ask()`

Send a natural language prompt to the chat interface. This is the
primary method for interacting with `SCPlotterChat`. It executes a
three-stage pipeline:

1.  **Tool identification** — The LLM selects the most appropriate
    scplotter visualization function based on the prompt and
    conversation history.

2.  **Data identification** — The LLM identifies the target data object
    from available sources (unless data has been manually set via
    `$set_data()`).

3.  **Code generation and execution** — The LLM generates R code to call
    the selected function with the identified data, which is then
    evaluated. The resulting plot is returned.

If the selected tool is `ListTools` or `ListData`, the respective
listing is printed instead of generating a plot.

#### Usage

    SCPlotterChat$ask(prompt, verbose = NULL, add_to_history = TRUE)

#### Arguments

- `prompt`:

  A character string containing the user's query or instruction in
  natural language (e.g., "Plot a UMAP of the pancreas data colored by
  cell type", "Make it a heatmap instead").

- `verbose`:

  Logical; if `TRUE`, prints the full LLM prompt and response for this
  interaction. Default is `NULL`, which falls back to the `verbose`
  setting of the `SCPlotterChat` instance. Use this to debug a single
  interaction without enabling verbose mode globally.

- `add_to_history`:

  Logical; if `TRUE` (default), this interaction is recorded in the
  conversation history, enabling the LLM to understand follow-up
  requests in context. Set to `FALSE` for one-off queries that should
  not influence subsequent interactions.

#### Returns

The generated plot (typically a `ggplot` object), or `NULL` invisibly
for `ListTools`/`ListData` requests.

------------------------------------------------------------------------

### `SCPlotterChat$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SCPlotterChat$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# \donttest{
if (FALSE) {
# Setup LLM provider (requires an API key)
provider <- tidyprompt::llm_provider_openai(
    parameters = list(model = "deepseek-v4-flash", stream = TRUE),
    url = "https://api.deepseek.com/chat/completions",
    api_key = Sys.getenv("OPENAI_API_KEY")
)

# Create chat instance
chat <- SCPlotterChat$new(provider)

# List available tools
chat$ask("What are the tools to use?")
# Tool identified:  ListTools
# Available tools:
# -  CCCPlot :  Cell-Cell Communication Plot
#    Visualizes ligand-receptor interaction inference results ...
# ...

# List available data
chat$ask("What data is available?")

# Generate a plot from natural language
chat$ask("Plot the default cell-cell communication plot for the cellphonedb_res dataset")
# Tool identified:  CCCPlot
# Data object identified:  cellphonedb_res
# Code ran: CCCPlot(cellphonedb_res)

# Refine with conversational context
chat$ask("do a heatmap instead")
# Tool identified:  CCCPlot
# Data object identified:  cellphonedb_res
# Code ran: CCCPlot(cellphonedb_res, plot_type = "heatmap")

# Add a title
chat$ask("Add a proper title to the plot")

# Manually set data to avoid auto-detection
chat$set_data(scplotter::cellphonedb_res)
chat$ask("Make a dot plot")

# Inspect conversation history
chat$get_history()

# Clear history for a fresh conversation
chat$clear_history()
}
# }
```
