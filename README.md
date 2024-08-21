# Genome size prediction tool

### Install prerequisites (TODO)

```
install.packages('devtools')
install.packages("BiocManager")
BiocManager::install(c('biomformat', 'ggtree', 'ncbitax'))
install.packages('pbapply')
install.packages('seqinr')
install.packages("merTools")
install.packages("remotes")
remotes::install_github("raim/ncbitax")
library(devtools)
```

### Install from zip

```
devtools::install_local('genomesizeR.zip')
```

### Install from GitHub

```
remotes::install_github("ScionResearch/genome_size_prediction")
```

### Load package

```
library(genomesizeR)
```

### Read example input file from the package

```
example_input_file = system.file("extdata", "example_input.csv", package = "genomesizeR")
```

### Store the path to the archive containing all the reference databases and the bayesian models

```
refdata_archive_path = "path/to/genomesizeRdata.tar.gz"
```


### Usage

```
estimate_genome_size(
  queries,
  refdata_archive_path,
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
  n_cores = "half"
)
```

### Arguments

`queries`
Queries: path to file or table object

`refdata_archive_path`
Path to the archive containing all the reference data and the bayesian models

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
Method to use for genome size estimation, 'weighted_mean' or 'lmm'

`ci_threshold`	
Threshold for the confidence interval as a proportion of the guessed size (e.g. 0.2 means that estimations with a confidence interval that represents more than 20% of the guessed size will be discarded)

`n_cores`	
Number of CPU cores to use (default is 'half': all available divided by two)


### Example

Run the main function to get estimated genome sizes:

```
results = estimate_genome_size(example_input_file, refdata_archive_path, format='csv', sep='\t', match_column='TAXID', output_format='input')
```

Plot genome size histogram per sample:

```
plotted_df = plot_genome_size_histogram(results)
```

Plot genome size histogram for one sample:

```
plotted_df = plot_genome_size_histogram(results, only_sample='16S_1')
```

Plot genome size boxplot per sample:

```
plotted_df = plot_genome_size_boxplot(results)
```

Plot genome size boxplot for one sample:

```
plotted_df = plot_genome_size_boxplot(results, only_sample='ITS_1')
```

Plot simplified taxonomic tree with colour-coded estimated genome sizes:

```
plotted_df = plot_genome_size_tree(results)
```
