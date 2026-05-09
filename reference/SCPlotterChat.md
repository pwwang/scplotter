# SCPlotterChat Class

An R6 class that provides chat functionality for SCPlotter

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

Create a new instance of the SCPlotterChat class

#### Usage

    SCPlotterChat$new(provider, verbose = FALSE)

#### Arguments

- `provider`:

  An LLM provider object See
  https://tjarkvandemerwe.github.io/tidyprompt/articles/getting_started.html#setup-an-llm-provider
  for more details on LLM providers

- `verbose`:

  A logical value indicating whether to print verbose messages Default
  is FALSE. And this will override the verbose setting of the provider

#### Returns

A new instance of the SCPlotterChat class

------------------------------------------------------------------------

### `SCPlotterChat$clear_history()`

Clear the chat history

#### Usage

    SCPlotterChat$clear_history()

#### Returns

NULL

------------------------------------------------------------------------

### `SCPlotterChat$get_history()`

Get the chat history

#### Usage

    SCPlotterChat$get_history()

#### Returns

The chat history

------------------------------------------------------------------------

### `SCPlotterChat$list_tools()`

Print the list of available tools

#### Usage

    SCPlotterChat$list_tools()

#### Returns

NULL

------------------------------------------------------------------------

### `SCPlotterChat$list_data()`

Print the list of available data objects that can be used

#### Usage

    SCPlotterChat$list_data()

#### Returns

NULL

------------------------------------------------------------------------

### `SCPlotterChat$set_data()`

Set the data to be analyzed

#### Usage

    SCPlotterChat$set_data(data, name = NULL)

#### Arguments

- `data`:

  The data object to be set

- `name`:

  The name of the data object (optional)

#### Returns

NULL

------------------------------------------------------------------------

### `SCPlotterChat$get_data()`

Get the data to be analyzed

#### Usage

    SCPlotterChat$get_data()

#### Returns

The data object

------------------------------------------------------------------------

### `SCPlotterChat$ask()`

Send a prompt to the chat interface and receive a response

#### Usage

    SCPlotterChat$ask(prompt, verbose = NULL, add_to_history = TRUE)

#### Arguments

- `prompt`:

  A character string containing the user's query or instruction

- `verbose`:

  A logical value indicating whether to print verbose messages Default
  is NULL, which will use the verbose setting of the SCPlotterChat
  object

- `add_to_history`:

  A logical value indicating whether to add the prompt and response to
  the chat history Default is TRUE

#### Returns

A response from the chat system

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
provider <- tidyprompt::llm_provider_openai(api_key = Sys.getenv("OPENAI_API_KEY"))
chat <- SCPlotterChat$new(provider)
chat$ask("What are the tools to use?")
# Tool identified:  ListTools
# Available tools:
# -  ClonalOverlapPlot :  ClonalOverlapPlot
#    Plot the overlap of the clones in different samples/groups.
#
# -  EnrichmentPlot :  Enrichment Plot
#    This function generates various types of plots for enrichment analysis.
# ...

chat$ask("Plot the default cell-cell communication plot for the cellphonedb_res dataset")
# Tool identified:  CCCPlot
# Data object identified:  cellphonedb_res
# Running tool:  CCCPlot

chat$ask("do a heatmap instead")
# Tool identified:  CCCPlot
# Data object identified:  cellphonedb_res
# Running tool:  CCCPlot
}
# }
```
