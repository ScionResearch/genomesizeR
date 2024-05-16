# Genome size prediction tool

### Install prerequisites

```
install.packages('devtools')
install.packages("BiocManager")
BiocManager::install(c('biomformat', 'ggtree', 'ncbitax'))
install.packages("remotes")
remotes::install_github("raim/ncbitax")
library(devtools)
```

### Install from zip

```
devtools::install_local('genomesizeR.zip')
```

### Load package

```
library(genomesizeR)
```

### Read example tax_table input file from the package

```
example_input_file = system.file("extdata", "example_Taxonomy_cleaned.csv", package = "genomesizeR")
```

### Usage

```
estimate_genome_size(
  queries,
  format = "csv",
  sep = ",",
  match_column = NA,
  match_sep = ";",
  size_db = NA,
  taxonomy = NA,
  output_format = "input",
  method = "weighted_mean",
  ci_threshold = 0.2,
  prediction_variables = c("family", "genus"),
  n_cores = "max"
)
```

### Arguments

`queries`
Queries: path to file or table object

`format`
Query format if in a file ('csv' (default), 'tax_table' or 'biom' (taxonomy table files))

`sep`	
If 'csv' format, column separator

`match_column`	
If 'csv' format, the column containing match information (with one or several matches)

`match_sep`	
If 'csv' format and several matches in match column, separator between matches

`size_db`	
Path to the genome size reference database to use if not the default one

`taxonomy`	
Path to taxonomy database (NCBI taxdump) to use if not the default one

`output_format`	
Format in which the output should be. Default: "input" a data frame with the same columns as the input, with the added columns: "estimated_genome_size" and "estimated_genome_size_confidence_interval". Other formats available: "data.frame", a data frame with the same number of rows as the input, and 3 columns: "TAXID", "estimated_genome_size" and "estimated_genome_size_confidence_interval".

`method`	
Method to use for genome size estimation, 'weighted_mean' or 'hierarchical'

`ci_threshold`	
Threshold for the confidence interval as a proportion of the guessed size (e.g. 0.2 means that estimations with a confidence interval that represents more than 20% of the guessed size will be discarded)

`n_cores`	
Number of CPU cores to use (default is 'max': all available minus 1)


### Run the main function to get estimated genome sizes

```
results = estimate_genome_size(example_input_file, format='tax_table')
```

### Read example tax_table sample file that will be read by the plotting functions

```
example_sample_file = system.file("extdata", "example_ASV_counts_cleaned.csv", package = "genomesizeR")
```

### Plot genome size histogram without using sample data

```
plotted_df = plot_genome_size_histogram(results)
```

### Plot genome size histogram per sample, using sample data

```
plotted_df = plot_genome_size_histogram(results, sample_data=example_sample_file)
```

### Plot genome size histogram for one sample, using sample data

```
plotted_df = plot_genome_size_histogram(results, sample_data=example_sample_file, only_sample='X16S_CPP47_K2PF5_CGAATACTGACA')
```

### Plot genome size boxplot without using sample data

```
plotted_df = plot_genome_size_boxplot(results)
```

### Plot genome size boxplot per sample, using sample data

```
plotted_df = plot_genome_size_boxplot(results, sample_data=example_sample_file)
```

### Plot genome size boxplot for one sample, using sample data

```
plotted_df = plot_genome_size_boxplot(results, sample_data=example_sample_file, only_sample='X16S_CPP47_K2PF5_CGAATACTGACA')
```

### Plot simplified taxonomic tree with colour-coded estimated genome sizes

```
plotted_df = plot_genome_size_tree(results)
```
