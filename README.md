# genomesizeR: Genome size prediction

### About the package

This R package uses statistical modelling on data from NCBI databases and provides three statistical methods for genome size prediction of a given taxon, or group of taxa. 

A straightforward weighted mean method identifies the closest taxa with available genome size information in the taxonomic tree and averages their genome sizes using weights based on taxonomic distance. A frequentist random effect model uses nested genus and family information to output genome size estimates. Finally, a third option provides predictions from a distributional Bayesian multilevel model which uses taxonomic information from genus all the way to superkingdom, therefore providing estimates and uncertainty bounds even for under-represented taxa.

All three methods use:

  - A list of queries; a query being a taxon or a list of several taxa.
  - A reference database containing all the known genome sizes, built from the NCBI databases, with associated taxa.
  - A taxonomic tree structure as built by the NCBI.

`genomesizeR` retrieves the taxonomic classification of input queries, estimates the genome size of each query, and provides 95% confidence intervals for each estimate.

### How to install

Prerequisites: [`R`](https://www.r-project.org/) with the already installed packages up-to-date, and [`git`](https://git-scm.com/downloads)

Run one of the commands below in an R console to install the package. We include four different installation methods, as some setups (for example, corporate networks) may block specific download mechanisms:

```
install.packages("remotes")
remotes::install_github("https://github.com/ScionResearch/genomesizeR")

- OR -

install.packages("remotes")
remotes::install_git("https://github.com/ScionResearch/genomesizeR")

- OR -

install.packages("devtools")
devtools::install_github("ScionResearch/genomesizeR")

- OR -

install.packages("pak")
pak::pkg_install("git::https://github.com/ScionResearch/genomesizeR")
```

You also need to download the archive containing the reference databases and the bayesian models from `zenodo.org`, using the `inborutils` package. You can change the `path` option to where you want to download the archive (default is current directory '.'):

```
remotes::install_github("inbo/inborutils")
inborutils::download_zenodo("10.5281/zenodo.13733183", path=".")
```

### Testing

[See test instructions for genomesizeR](TESTING.md)

### Get started

[Here is a simple tutorial using the default bayesian method on an example file](https://scionresearch.github.io/genomesizeR/articles/genomesizeR.html)

### How to contribute 

[Contribution guidelines for genomesizeR](CONTRIBUTING.md)
