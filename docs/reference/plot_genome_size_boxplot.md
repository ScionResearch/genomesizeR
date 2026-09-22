# Plot genome size boxplot

This function loads a result table from estimate_genome_size and plots a
boxplot of estimated genome sizes.

## Usage

``` r
plot_genome_size_boxplot(
  output_table,
  sample_data = NA,
  only_sample = NA,
  PA = F
)
```

## Arguments

- output_table:

  Result table from estimate_genome_size

- sample_data:

  The file or dataframe containing sample information (sample name and
  count)

- only_sample:

  The sample to plot the genome size boxplot for (default: NA, all
  samples)

- PA:

  Use presence/absence data instead of abundance data (all counts are
  set to 1)
