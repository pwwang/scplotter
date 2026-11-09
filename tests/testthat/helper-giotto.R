# GiottoClass resolves its python environment through conda every time a gobject is loaded:
# `loadGiotto()` -> `set_giotto_python_path()` -> `.envname_to_pypath()` -> `.check_conda()`,
# and the last step stops with
#   "Unable to find a conda binary. Use `installGiottoEnvironment()` or install a custom conda."
# on any machine without conda, which includes every CI runner.
#
# Nothing in this suite needs python: the Giotto tests load minis with GiottoData and only plot
# them (no Leiden, no spatial enrichment). `giotto.use_conda = FALSE` makes
# `set_giotto_python_path()` return before it touches reticulate, so objects load without a
# python environment. Set once here rather than in each test, so a new test cannot reintroduce
# the failure by forgetting it.
options(giotto.use_conda = FALSE)
