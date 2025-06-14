SHELL = /bin/bash

readme:
	Rscript <(echo "devtools::build_readme()")

docs:
	Rscript <(echo "devtools::document(); pkgdown::build_site()")

install:
	Rscript <(echo "devtools::build(); devtools::install()")

test:
	Rscript <(echo "devtools::test()")

notebooks:
	jupyter nbconvert -y --to html notebooks/spatial/*.ipynb --output-dir=pkgdown/assets

notebook: notebooks
nb: notebooks

.PHONY: readme docs install test notebooks notebook nb
