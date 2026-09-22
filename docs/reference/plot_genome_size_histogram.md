# Plot genome size histogram

This function loads a result table from estimate_genome_size and plots a
histogram of estimated genome sizes.

## Usage

``` r
plot_genome_size_histogram(
  output_table,
  sample_data = NA,
  only_sample = NA,
  bins = 50,
  PA = F
)
```

## Arguments

- output_table:

  Result table from estimate_genome_size()

- sample_data:

  The file or dataframe containing sample information (sample name and
  count)

- only_sample:

  The sample to plot the genome size histogram for (default: NA, all
  samples)

- bins:

  Histogram bin parameter

- PA:

  Use presence/absence data instead of abundance data (all counts are
  set to 1)
