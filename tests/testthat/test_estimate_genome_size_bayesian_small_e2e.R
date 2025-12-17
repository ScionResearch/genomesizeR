# E2E Test: estimate_genome_size with Bayesian method on small example dataset

test_that("estimate_genome_size produces expected results with Bayesian method (small dataset)", {
  skip_on_cran()

  # Load small example data
  example_input <- system.file("extdata", "example_input_small.csv", package = "genomesizeR")
  expected_results <- read.csv(
    system.file("extdata", "example_output_bayesian_small.csv", package = "genomesizeR"),
    stringsAsFactors = FALSE
  )

  # Convert numeric columns that may have been read as character
  numeric_cols <- c("confidence_interval_lower", "confidence_interval_upper")
  for (col in numeric_cols) {
    if (col %in% names(expected_results)) {
      expected_results[[col]] <- as.numeric(as.character(expected_results[[col]]))
    }
  }

  # Get path to reference data
  refdata_path <- system.file("extdata", "genomesizeRdata_v1.0.3.tar.gz", package = "genomesizeR")

  # Skip test if reference data is not available
  if (!file.exists(refdata_path)) {
    skip("Reference data not available")
  }

  # Run estimation with bayesian method
  results <- estimate_genome_size(
    example_input,
    refdata_path,
    sep = ',',
    match_column = 'TAXID',
    output_format = 'input',
    method = 'bayesian',
    ci_threshold = 0.5,
    n_cores = 1
  )

  # Test 1: Check that results is a data frame
  expect_s3_class(results, "data.frame")

  # Test 2: Compare estimated genome sizes with 5% tolerance (Bayesian is stochastic)
  expect_equal(
    results$estimated_genome_size,
    expected_results$estimated_genome_size,
    tolerance = 0.05,
    label = "estimated_genome_size"
  )

  # Test 3: Check that all estimated genome sizes are positive
  expect_true(all(results$estimated_genome_size > 0))

  # Test 4: Check that confidence intervals make sense (lower <= estimate <= upper)
  valid_ci <- !is.na(results$confidence_interval_lower) &
    !is.na(results$estimated_genome_size) &
    !is.na(results$confidence_interval_upper)
  if (any(valid_ci)) {
    res_lower <- as.numeric(results$confidence_interval_lower[valid_ci])
    res_est <- as.numeric(results$estimated_genome_size[valid_ci])
    res_upper <- as.numeric(results$confidence_interval_upper[valid_ci])
    expect_true(all(res_lower <= res_est), label = "lower CI <= estimate")
    expect_true(all(res_est <= res_upper), label = "estimate <= upper CI")
  }

  # Test 5: Compare confidence intervals with 10% tolerance
  expect_equal(
    results$confidence_interval_lower,
    expected_results$confidence_interval_lower,
    tolerance = 0.1,
    label = "confidence_interval_lower"
  )

  expect_equal(
    results$confidence_interval_upper,
    expected_results$confidence_interval_upper,
    tolerance = 0.1,
    label = "confidence_interval_upper"
  )
})
