## Version 0.7.1

- fix(ClonalStatPlot): streamline xlab assignment for clarity in plots
- feat(ClonalResidencyPlot): add group_by_sep parameter for combining multiple group_by columns
- feat(ClonalResidencyPlot): add with_class parameter to control clonal class inclusion in plots
- feat(utils): implement do_call function for improved performance in function calls
- feat(CellStatPlot): refactor Heatmap call to use do_call for improved argument handling
- feat(ClonalStatPlot): refactor Heatmap call to use do_call for improved argument handling
- feat(CellStatPlot): enhance'rows_by' and 'rows_split_by' parameters in heatmap plots
- feat(ClonalStatPlot): support subgroup_by for heatmap
- ci: update ggrepel package version to 0.9.5 in dependencies (0.9.6 requires R 4.5)
- fix(ClonalVolumePlot): ensure grouping_levels are simplified correctly (fixing order not working)
