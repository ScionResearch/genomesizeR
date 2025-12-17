# E2E Test: estimate_genome_size with weighted_mean method on small example dataset

test_that("estimate_genome_size produces expected results with weighted_mean method (small dataset)", {
  skip_on_cran()

  # Load small example data
  example_input <- system.file("extdata", "example_input_small.csv", package = "genomesizeR")
  expected_results <- read.csv(
    system.file("extdata", "example_output_weighted_mean_small.csv", package = "genomesizeR"),
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

  # Run estimation with weighted_mean method
  results <- estimate_genome_size(
    example_input,
    refdata_path,
    sep = ',',
    match_column = 'TAXID',
    output_format = 'input',
    method = 'weighted_mean',
    ci_threshold = 0.5,
    n_cores = 1
  )

  # Test 1: Check that results is a data frame
  expect_s3_class(results, "data.frame")

  # Test 2: Compare estimated genome sizes (weighted_mean is deterministic, use very tight tolerance)
  # Only compare non-NA values
  non_na_indices <- !is.na(results$estimated_genome_size) & !is.na(expected_results$estimated_genome_size)
  if (any(non_na_indices)) {
    expect_equal(
      results$estimated_genome_size[non_na_indices],
      expected_results$estimated_genome_size[non_na_indices],
      tolerance = 1e-5,
      label = "estimated_genome_size"
    )
  }

  # Test 3: Compare confidence intervals (deterministic method)
  expect_equal(
    results$confidence_interval_lower[non_na_indices],
    expected_results$confidence_interval_lower[non_na_indices],
    tolerance = 1e-5,
    label = "confidence_interval_lower"
  )

  expect_equal(
    results$confidence_interval_upper[non_na_indices],
    expected_results$confidence_interval_upper[non_na_indices],
    tolerance = 1e-5,
    label = "confidence_interval_upper"
  )

  # Test 4: Exact match for categorical columns (non-NA values)
  if (any(non_na_indices)) {
    expect_equal(
      results$genome_size_estimation_status[non_na_indices],
      expected_results$genome_size_estimation_status[non_na_indices]
    )
  }

  # Test 5: Check that all non-NA estimated genome sizes are positive
  valid_estimates <- !is.na(results$estimated_genome_size)
  if (any(valid_estimates)) {
    expect_true(all(results$estimated_genome_size[valid_estimates] > 0))
  }

  # Test 6: Check that confidence intervals make sense (lower < estimate < upper)
  valid_ci <- !is.na(results$confidence_interval_lower) &
    !is.na(results$estimated_genome_size) &
    !is.na(results$confidence_interval_upper)
  if (any(valid_ci)) {
    expect_true(all(results$confidence_interval_lower[valid_ci] < results$estimated_genome_size[valid_ci]))
    expect_true(all(results$estimated_genome_size[valid_ci] < results$confidence_interval_upper[valid_ci]))
  }

  # Test 7: Check that model_used contains expected method (for non-NA values)
  if (any(valid_estimates)) {
    expect_true(all(grepl("weighted_mean|reference", results$model_used[valid_estimates])))
  }
})
