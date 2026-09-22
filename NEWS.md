# genomesizeR 1.1.0

## New features

* `estimate_genome_size()` gains a `return_posterior` argument (bayesian method only). When `TRUE`, the full posterior predictive distribution of each query is returned in additional columns `posterior_1`, ..., `posterior_N`, one per posterior draw. The distributions are predicted jointly for all queries, so that the correlation between the predictions of related queries is preserved when they are aggregated. The number of draws can be limited with the new `n_draws` argument.

* New function `estimate_genome_size_per_sample()`, which aggregates the posterior distributions of individual queries (e.g. ASVs) into one mean genome size estimation per sample, given an abundance table. Queries are weighted by log-transformed abundance by default (relative abundance and presence are also available). The interval returned for each sample mean accounts for the estimation uncertainty in the genome size of every query, and the proportion of each sample's abundance that could be estimated is reported.

* New function `plot_genome_size_per_sample()`, which plots the estimated mean genome size of each sample with its interval.

* The row names of the input (e.g. ASV identifiers of a taxonomy table) are now kept in an `ASVs` column of the result table.

## Improvements

* Identical queries (same taxonomy for the `tax_table` and `biom` formats, same value in `match_column` when one is given, same row otherwise) are only computed once and get identical results. On ASV tables this reduces the number of model predictions by one to two orders of magnitude.

## Bug fixes

* The lmm method no longer estimates queries whose family has no genome size reference. Such queries were previously predicted from the model intercept alone, which returned the overall mean of the reference data with a narrow interval and an `OK` status. They now get the status `No reference for family in lmm model` and no estimate.

* With the lmm method, queries estimated from a reference mean keep the interval derived from the reference standard error instead of having it overwritten by a model interval. The lmm confidence interval columns are returned as numeric.

## Documentation

* New sections on sample-level estimation in the *Get started* and *Method description* vignettes.

* The *Input and output* vignette documents the `ASVs` and `posterior_*` columns, the new lmm status, and the meaning of the `Bayesian model not found` status (the query does not belong to Bacteria, Archaea or Eukaryota).
