# Genome size prediction tool

### Install from zip

```
devtools::install_local('genomesizeR.zip')
```

### Load package

```
library(genomesizeR)
```

### Read example dada2 input file from the package

```
example_input_file = system.file("extdata", "example_Taxonomy_cleaned.csv", package = "genomesizeR")
```

### Run the main function to get estimated genome sizes

```
results = estimate_genome_size(example_input_file, format='dada2')
```

### Read example dada2 sample file that will be read by the plotting functions

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
