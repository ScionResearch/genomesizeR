# Estimate genome sizes

This function loads a query file or table and an archive containing
reference databases and bayesian models, and predicts genome sizes.

## Usage

``` r
estimate_genome_size(
  queries,
  refdata_path,
  format = "csv",
  sep = ",",
  match_column = NA,
  match_sep = ";",
  output_format = "input",
  method = "bayesian",
  ci_threshold = 0.3,
  n_cores = "half",
  return_posterior = FALSE,
  n_draws = NULL
)
```

## Arguments

- queries:

  Queries: path to csv or BIOM file, or R object used for input.

- refdata_path:

  Path to the downloadable archive containing the reference databases
  and the bayesian models.

- format:

  Input format: "csv" for csv file (default), "tax_table" for taxonomy
  table file or object as used in e.g. phyloseq, "biom" for BIOM file,
  "dataframe" for a table-style object (e.g. data.frame or matrix
  object), "vector" for a vector object.

- sep:

  If table-style or csv format, column separator (default: ",").

- match_column:

  If table-style or csv format, the column containing match information
  (with one or several matches).

- match_sep:

  If table-style or csv format and several matches in match column,
  separator between matches (default: ";").

- output_format:

  Format in which the output should be. Default: "input" a data frame
  with the same columns as the input, with the added columns: "TAXID",
  "estimated_genome_size", "confidence_interval_lower",
  "confidence_interval_upper", "genome_size_estimation_status",
  "model_used", as well as taxids at all ranks. Other formats available:
  "data.frame", a data frame with only the previously described columns,
  without the "taxids at all ranks" columns.

- method:

  Method to use for genome size estimation, 'bayesian' (default),
  'weighted_mean' or 'lmm'.

- ci_threshold:

  Threshold for the confidence interval as a proportion of the predicted
  size (e.g. 0.3 means that estimations with a confidence interval that
  represents more than 30% of the predicted size will be tagged in the
  output table).

- n_cores:

  Number of CPU cores to use (default is 'half': half of all available
  cores).

- return_posterior:

  Bayesian method only. If TRUE, the full posterior predictive
  distribution of each query is returned in additional columns
  "posterior_1", ..., "posterior_N" (one column per posterior draw, in
  base pairs), in addition to the summarised estimate and confidence
  interval. The distributions are then predicted jointly for all
  queries, so that queries sharing a taxon unseen by the model share the
  same sampled effects in each draw, as needed to propagate uncertainty
  correctly when aggregating queries. Queries estimated from a reference
  mean ("reference_mean" in "model_used") get a constant posterior equal
  to their estimate. Required to compute sample-level estimates with
  [`estimate_genome_size_per_sample`](https://scionresearch.github.io/genomesizeR/reference/estimate_genome_size_per_sample.md).

- n_draws:

  Bayesian method only. Number of posterior draws to use (default: NULL,
  all draws of the model). When return_posterior is TRUE and the models
  have different numbers of draws, the smallest number is used so that
  all queries have the same number of posterior columns.

## Value

A data frame with one row per query. If the input has row names (e.g.
ASV identifiers of a taxonomy table), they are kept in a column "ASVs".

## Details

Identical queries (same taxonomy for the "tax_table" and "biom" formats,
same value in match_column when one is given, same row otherwise) are
only computed once and get identical results.
