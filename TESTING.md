# Testing

Instructions to run tests for the `genomesizeR` package.

Prerequisites
- Install `devtools` and `testthat`:

```r
install.packages("devtools")
install.packages("testthat")
```

Download reference data archive and place in `inst/extdata/`, e.g.:

```r
inborutils::download_zenodo("10.5281/zenodo.13733183", path="inst/extdata/")
```

Run all tests

```r
devtools::test()
```

Run the three E2E tests for each method

```r
devtools::test_file("tests/testthat/test_estimate_genome_size_bayesian_small_e2e.R")
devtools::test_file("tests/testthat/test_estimate_genome_size_lmm_small_e2e.R")
devtools::test_file("tests/testthat/test_estimate_genome_size_weighted_mean_small_e2e.R")
```

Notes
- All tests require the reference archive `genomesizeRdata_v*.tar.gz` in `inst/extdata/`.
