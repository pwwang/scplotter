# Names used inside dplyr's data-masking verbs, which codetools cannot resolve and therefore
# reports as "no visible binding for global variable". They are real columns: `.sig_metric` is
# created a few lines above its use in EnrichmentPlot().
utils::globalVariables(".sig_metric")
