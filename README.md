# genomesizeR: Genome size prediction

### About the package

This R package uses statistical modelling on data from NCBI databases and provides three statistical methods for genome size prediction of a given taxon, or group of taxa. 

A straightforward weighted mean method (`weighted-mean`) identifies the closest taxa with available genome size information in the taxonomic tree and averages their genome sizes using weights based on taxonomic distance. A frequentist random effect model uses nested genus and family information to output genome size estimates. Finally, a third option provides predictions from a distributional Bayesian multilevel model which uses taxonomic information from genus all the way to superkingdom, therefore providing estimates and uncertainty bounds even for under-represented taxa.

All three methods use:

  - A list of queries; a query being a taxon or a list of several taxa.
  - A reference database containing all the known genome sizes, built from the NCBI databases, with associated taxa.
  - A taxonomic tree structure as built by the NCBI.

`genomesizeR` retrieves the taxonomic classification of input queries, estimates the genome size of each query, and provides 95% confidence intervals for each estimate.

### How to install

Install from GitHub:

```
install.packages("remotes")
remotes::install_github("ScionResearch/genomesizeR")
```

Download the archive containing the reference databases and the bayesian models from `zenodo.org`, using the `inborutils` package. You can change the `path` option to where you want to download the archive (default is current directory '.'):

```
remotes::install_github("inbo/inborutils")
inborutils::download_zenodo("10.5281/zenodo.13733183", path=".")
```

### Simple example

Store the path to the archive containing the reference databases and the bayesian models:

```
refdata_archive_path = "path/to/genomesizeRdata.tar.gz"
```

Read the example input file from the package:

```
example_input_file = system.file("extdata", "example_input.csv", package = "genomesizeR")
```

Load the package:

```
library(genomesizeR)
```

Run the main function to get the estimated genome sizes (with the default method which is the bayesian method):

```
results = estimate_genome_size(example_input_file, refdata_archive_path, sep='\t', match_column='TAXID', output_format='input')
```

Plot the genome size histogram per sample:

```
plotted_df = plot_genome_size_histogram(results)
```

Plot the genome size histogram for one sample:

```
plotted_df = plot_genome_size_histogram(results, only_sample='16S_1')
```

Plot the genome size boxplot per sample:

```
plotted_df = plot_genome_size_boxplot(results)
```

Plot the genome size boxplot for one sample:

```
plotted_df = plot_genome_size_boxplot(results, only_sample='ITS_1')
```

Plot the simplified taxonomic tree with colour-coded estimated genome sizes:

```
plotted_df = plot_genome_size_tree(results, refdata_archive_path)
```
