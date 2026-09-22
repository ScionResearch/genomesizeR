# Estimate genome sizes per sample

This function aggregates the genome size estimations of individual
queries (e.g. ASVs) into one genome size estimation per sample,
propagating the estimation uncertainty.

## Usage

``` r
estimate_genome_size_per_sample(
  output_table,
  otu_table,
  taxa_are_rows = TRUE,
  weights = "log_count",
  ci_level = 0.95,
  return_posterior = FALSE
)
```

## Arguments

- output_table:

  Result table from
  [`estimate_genome_size`](https://scionresearch.github.io/genomesizeR/reference/estimate_genome_size.md)
  run with `method='bayesian'` and `return_posterior=TRUE`.

- otu_table:

  Abundance table with one row per query (in the same order as the rows
  of output_table, or with row names matching the "ASVs" column of
  output_table when there is one) and one column per sample (e.g. the
  OTU table of a phyloseq object, as a matrix or data frame).

- taxa_are_rows:

  Whether queries are rows and samples are columns in otu_table
  (default: TRUE). If FALSE, otu_table is transposed.

- weights:

  How query abundances are used to weight the mean genome size of a
  sample: "log_count" (default): log(1 + count); "relative_abundance":
  count (proportional to the relative abundance of the query in the
  sample); "presence": 1 for every query present in the sample.

- ci_level:

  Confidence level of the interval computed from the sample-level
  posterior distribution (default: 0.95).

- return_posterior:

  If TRUE, the sample-level posterior distributions are returned in
  additional columns "posterior_1", ..., "posterior_N" (one column per
  draw).

## Value

A data frame with one row per sample and the columns "sample",
"estimated_genome_size" (mean of the sample-level posterior
distribution, in base pairs), "confidence_interval_lower",
"confidence_interval_upper", "n_queries_used" (queries with an
estimation and present in the sample), "n_queries_total" (queries
present in the sample) and "abundance_covered" (proportion of the sample
abundance coming from queries with an estimation).

## Details

It requires the full posterior predictive distribution of each query,
obtained by running
[`estimate_genome_size`](https://scionresearch.github.io/genomesizeR/reference/estimate_genome_size.md)
with the bayesian method and `return_posterior=TRUE`. For each sample
and each posterior draw, the mean genome size across queries is
computed, weighted by query abundance. The resulting sample-level
posterior distribution is summarised into a mean estimate and a
confidence interval. Queries without an estimation (e.g. taxon not
found) are ignored, and the proportion of the sample abundance that they
represent is reported in the column "abundance_covered".
