# Changelog

## Version 0.7.5

- fix(DESCRIPTION): specify exact version for drieslab/GiottoClass
  dependency
- fix(DESCRIPTION): update GiottoClass and GiottoUtils dependency
  versions
- fix(DESCRIPTION): update GiottoClass dependency to use the latest
  version without a specific commit
- ci: update GiottoClass dependency to version 0.4.7
- ci: update GiottoClass dependency to version R4.4.0
- chore(DESCRIPTION): specify exact version for drieslab/GiottoClass
  dependency
- ci: update GiottoClass and GiottoUtils dependencies to specific
  versions
- ci: update GiottoClass dependency to use the latest version
- docs: add documentation section with link to online resources \[skip
  ci\]
- fix: enhance gene selection logic in MarkersPlot function for heatmap
  and dot plots
- fix(CCCPlot): fix rows_split_by and columns_split_by for heatmp with
  method = ‘interaction’
- docs: update LLM provider model and API endpoint in
  Visualizing_data_with_LLMs vignette
- docs: update example in FeatureStatPlot documentation to include row
  annotation palette
- fix: enhance handling of heatmap and dot plot types in FeatureStatPlot
  function
- fix(FeatureStatPlot): hide column names for heatmap with cell_type
  bars
- refactor(do_call): use do_call from plotthis
- feat: add linkedheatmap plot type to CCCPlot function with new
  parameters for ligand and receptor expressions
- docs(CCCPlot): enhance CCCPlot documentation: Update title, arguments,
  and descriptions for clarity and detail. Improve data format and
  method sections, and refine examples for better usability.
- docs(CellDimPlot/CellVelocityPlot): enhance documentation for
  CellDimPlot and CellVelocityPlot functions
- docs(CellStatPlot): enhance CellStatPlot documentation and argument
  descriptions
- docs(SCPlotterChat): enhance documentation for LLM-powered chat
  interface, including detailed descriptions of architecture,
  conversation history, tool registry, and LLM provider setup
- docs: update documentation for Clonal Composition, Length, Overlap,
  Residency, and Volume plots
- docs: enhance ClonalDiversityPlot and ClonalRarefactionPlot
  documentation for clarity and detail
- docs(ClonalGeneUsagePlot): enhance documentation for clarity and
  detail, including usage examples and parameter descriptions
- docs: update documentation for ClonalKmerPlot and ClonalPositionalPlot
- docs: enhance ClonalStatPlot documentation with detailed arguments and
  usage examples
- docs: enhance documentation for clonal data processing functions,
  improving clarity and detail
- docs: enhance ClustreePlot documentation for clarity and detail,
  including usage examples and parameter descriptions
- docs: enhance EnrichmentPlot documentation for clarity and detail,
  including comprehensive parameter descriptions and usage examples
- docs: enhance FeatureStatPlot documentation and internal dispatcher
- docs: enhance MarkersPlot documentation for clarity and detail
- docs: enhance documentation for spatial plot functions in Seurat and
  Giotto
- docs: update documentation links and examples for clarity and
  consistency across multiple plot functions

## Version 0.7.4

- feat(ClonalStatPlot): enhance group selection (argument: groups) with
  named vector/list
- feat(ClonalResidencyPlot): enhance group handling for scatter plots
  with named vectors/lists
- feat(MarkersPlot): enhance selection functionality with custom
  filtering conditions

## Version 0.7.3

- fix(ClonalResidencyPlot): update data selection to include unique
  split_by parameter for upset plots
- chore: update LLM provider configuration in vignette to use new model
  and API key
- fix: update LLM provider model to use openai/gpt-5.4-nano in vignette
- fix(ClonalResidencyPlot): remove unused title parameter and update
  function usage in documentation
- fix(ClonalResidencyPlot): refactor segment data handling for improved
  readability
- feat: implement default dimension reduction functions for Seurat
  objects
- ci: update package dependencies in workflow for improved compatibility
- feat(MarkersPlot): add options to show labels and customize
  significant markers in heatmap plots
- fix: fix setting default dimreduction when “DefaultDimReduc\<-\` is
  available
- fix(CellStatPlot): ensure rows_name is set correctly for heatmap
  arguments
- refactor(CellStatPlot): remove rows_name parameter from function and
  documentation
- chore(CellStatPlot): improve error messages for ‘group_by’ and
  ‘rows_by’ parameters in pie heatmap
- refactor(FeatureStatPlot): comment out row and column name annotation
  arguments

## Version 0.7.2

- docs: update title and description to include spatial data analysis
- refactor: make implementation of group level retaining consistant
  across all clonal functions
- fix: set sample identity for Seurat objects in ClonalLengthPlot
  function
- fix: adjust x-axis breaks in ClonalLengthPlot to nearest tens and
  handle missing breaks
- chore(ClonalLengthPlot): add default value for position_dodge_preserve

## Version 0.7.1

- fix(ClonalStatPlot): streamline xlab assignment for clarity in plots
- feat(ClonalResidencyPlot): add group_by_sep parameter for combining
  multiple group_by columns
- feat(ClonalResidencyPlot): add with_class parameter to control clonal
  class inclusion in plots
- feat(utils): implement do_call function for improved performance in
  function calls
- feat(CellStatPlot): refactor Heatmap call to use do_call for improved
  argument handling
- feat(ClonalStatPlot): refactor Heatmap call to use do_call for
  improved argument handling
- feat(CellStatPlot): enhance’rows_by’ and ‘rows_split_by’ parameters in
  heatmap plots
- feat(ClonalStatPlot): support subgroup_by for heatmap
- ci: update ggrepel package version to 0.9.5 in dependencies (0.9.6
  requires R 4.5)
- fix(ClonalVolumePlot): ensure grouping_levels are simplified correctly
  (fixing order not working)
