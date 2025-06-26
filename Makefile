SHELL = /bin/bash

readme:
	Rscript <(echo "devtools::build_readme()")

docs:
	Rscript <(echo "devtools::document(); pkgdown::build_site()")

install:
	Rscript <(echo "devtools::build(); devtools::install()")

test:
	Rscript <(echo "devtools::test()")

# make notebooks EXECUTE=true to run the notebooks
notebooks:
	jupyter nbconvert $(if $(DEBUG),--execute) -y \
		--to html notebooks/*.ipynb notebooks/spatial/*.ipynb \
		--output-dir=pkgdown/assets $(if $(EXECUTE),--execute) \
		--template lab

notebook: notebooks
nb: notebooks

.PHONY: readme docs install test notebooks notebook nb
