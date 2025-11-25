# Four decades of species-specific remotely sensed phenology

This project integrates remotely sensed phenology from Landsat series, enabled by the [Baysian Land Survace Phenology (BLSP)](https://github.com/ncsuSEAL/Bayesian_LSP) model, and tree species information obtained from the [Forest Inventory and Analysis (FIA) program](https://research.fs.usda.gov/programs/fia) to compile a unique phenology dataset that spans four decades (1985-present) and includes 10 major tree species across the Northestern and Midwestern United States. Along with long-term in-situ phenology observations at two research forests ([Harvard Forest, MA](https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=hf003) and [Hubbard Brook Experimental Forest, NH](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-hbr.51.14)), we investigate species-specific phenological dynamics and responses to climate variation.

![site-map](/out/site-map.png)


## The structure of this repository

Directory structure:

- `data`: Raw and analysis-ready data. Note that raw data is not version controlled.
- `out`: The output figures & tables.
- `pipe`: Intermediate results/data/vis...
- `src`: Source code. Files are named under the following rules to decouple data processing, modeling, and visualization:
    - `dat_dl_`: Data downloading.
    - `dat_pr_`: Data processing.
    - `hlp_`: Helper functions that may be used in multiple scripts.
    - `mod_`: Model and analysis.
    - `vis_`: Visualization.
    - `base.R`: This file contains project configuration variables, especially the file paths over multiple storage locations, and some common functions.

Key R packages include:
- [`blsp`](https://github.com/ncsuSEAL/Bayesian_LSP) (v1.5)
- [`terra`](https://cran.r-project.org/web/packages/terra/index.html) (v1.8.7)
- [`data.table`](https://cran.r-project.org/web/packages/data.table/index.html) (v1.17.8) 
- [`rjags`](https://cran.r-project.org/web/packages/rjags/index.html) (v4.16). Note that the [Just Another Gibbs Sampler (JAG; v4.3.0)](https://mcmc-jags.sourceforge.io/) needs to be installed to use `rjags` for Bayesian model fitting. 
