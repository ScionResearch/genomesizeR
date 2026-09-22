# E2E Test: identical queries are computed once and get identical results,
# while each row keeps its own input columns

test_that("duplicated queries get identical estimates and keep their own input columns", {
  skip_on_cran()

  refdata_path <- system.file("extdata", "genomesizeRdata_v1.0.3.tar.gz", package = "genomesizeR")
  if (!file.exists(refdata_path)) {
    skip("Reference data not available")
  }

  # Repeat the first 4 queries of the small example in another sample with other counts
  input <- read.csv(system.file("extdata", "example_input_small.csv", package = "genomesizeR"))
  repeated <- input[1:4, ]
  repeated$COUNT <- repeated$COUNT * 10
  repeated$SAMPLE <- "DUP"
  input_file <- tempfile(fileext = ".csv")
  write.csv(rbind(input, repeated), input_file, row.names = FALSE)

  results <- estimate_genome_size(
    input_file,
    refdata_path,
    match_column = 'TAXID',
    method = 'bayesian',
    n_cores = 1,
    return_posterior = TRUE,
    n_draws = 50
  )

  n <- nrow(input)
  expect_equal(nrow(results), n + 4)

  # Input columns of the repeated rows are their own
  expect_equal(as.numeric(results$COUNT[(n + 1):(n + 4)]), repeated$COUNT)
  expect_equal(results$SAMPLE[(n + 1):(n + 4)], rep("DUP", 4))

  # Estimates and posteriors of repeated queries are identical (computed once)
  expect_equal(results$estimated_genome_size[(n + 1):(n + 4)], results$estimated_genome_size[1:4])
  posterior_columns <- grep("^posterior_[0-9]+$", names(results), value = TRUE)
  expect_equal(unname(as.matrix(results[(n + 1):(n + 4), posterior_columns])),
               unname(as.matrix(results[1:4, posterior_columns])))
})
