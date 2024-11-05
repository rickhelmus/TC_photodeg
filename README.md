
## Overview

This repository accompanies the manuscript *“Comprehensive mass
spectrometry workflows to systematically elucidate transformation
processes of organic micropollutants: a case study on photodegradation
of four pharmaceuticals”* and contains the
[R](https://www.r-project.org/) and
[patRoon](https://rickhelmus.github.io/patRoon/) based data processing
scripts and algorithms that were used to elucidate transformation
products from photodegradation experiments.

The repository layout is as follows:

- `workflow.R`: the main processing script
- `workflow/` directory: script files that implement the main workflow
- `utils/` directory: scripts with utility functions used throughout the
  workflow
- `report_pd` directory: all support files needed for automatic report
  generation
- `misc/` directory: miscellaneous workflow input files (suspect lists,
  sample info etc)
- `output/` directory: output files generated throughout the workflow
  - specifically: `suspects_TPs-compounds.csv` and
    `suspects_TPs-formulas.csv`: suspect lists generated during TP
    suspect screening workflows
- `renv.lock`: an [renv](https://rstudio.github.io/renv) lockfile to
  recreate a reproducible R environment

## Used software

### General

| Software | Version | URL |
|----|----|----|
| R | 4.2.2 | <https://www.r-project.org/> |
| DataAnalysis | 4.4 | <https://www.bruker.com/> |
| OpenMS | 2.7 | <https://openms.de/> |
| GenForm | embededded in patRoon | <https://sourceforge.net/projects/genform/> |
| MetFrag | 2.5 | <http://ipb-halle.github.io/MetFrag/> |
| BioTransformer | 3.0.0 | <http://biotransformer.ca/> |
| Chemical Transformation Simulator | Accessed online July 19th 2024 | <https://qed.epa.gov/cts/> |
| MolView | 2.4 | <https://molview.org/> |
| CDK Depict | – | <https://www.simolecule.com/cdkdepict/depict.html> |
| OpenBabel | 3.1.1 | <https://github.com/openbabel/openbabel> |
| Inkscape | 1.2.2 | <https://inkscape.org/> |

### Used R packages

| package | version | url |
|:---|:---|:---|
| ChemmineR | 3.48.0 | <https://github.com/girke-lab/ChemmineR> |
| RColorBrewer | 1.1.3 | <https://cran.r-project.org/web/packages/RColorBrewer> |
| callr | 3.7.6 | <https://callr.r-lib.org> |
| cliqueMS | 1.7.1 | <http://cliquems.seeslab.net> |
| data.table | 1.15.4 | <https://r-datatable.com> |
| depict | 0.4.0 | <https://github.com/CDK-R/depict> |
| fastcluster | 1.2.6 | <https://danifold.net/fastcluster.html> |
| fmcsR | 1.38.0 | <https://github.com/girke-lab/fmcsR> |
| future | 1.33.2 | <https://future.futureverse.org> |
| future.apply | 1.11.2 | <https://future.apply.futureverse.org> |
| kableExtra | 1.4.0.4 | <http://haozhu233.github.io/kableExtra/> |
| knitr | 1.48 | <https://yihui.org/knitr/> |
| magick | 2.7.3 | <https://docs.ropensci.org/magick/> |
| mzR | 2.30.0 | <https://github.com/sneumann/mzR/> |
| pagedown | 0.20.1 | <https://github.com/rstudio/pagedown> |
| patRoon | 2.3.3 | <https://github.com/rickhelmus/patRoon> |
| processx | 3.8.4 | <https://processx.r-lib.org> |
| rcdk | 3.8.1 | <https://cdk-r.github.io/cdkr/> |
| reactable | 0.4.4 | <https://glin.github.io/reactable/> |
| renv | 1.0.7 | <https://rstudio.github.io/renv/> |
| rmarkdown | 2.27 | <https://github.com/rstudio/rmarkdown> |
| svglite | 2.1.3 | <https://svglite.r-lib.org> |
| withr | 3.0.0 | <https://withr.r-lib.org> |
