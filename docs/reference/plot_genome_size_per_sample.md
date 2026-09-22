# Plot genome size per sample

This function loads a result table from estimate_genome_size_per_sample
and plots the estimated genome size of each sample with its confidence
interval.

## Usage

``` r
plot_genome_size_per_sample(per_sample_table, order_by_size = TRUE)
```

## Arguments

- per_sample_table:

  Result table from estimate_genome_size_per_sample()

- order_by_size:

  Order samples by estimated genome size (default: TRUE) rather than by
  sample name
