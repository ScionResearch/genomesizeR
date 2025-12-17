# E2E Test: estimate_genome_size with LMM method on small example dataset

test_that("estimate_genome_size produces expected results with LMM method (small dataset)", {
  skip_on_cran()

  # Load small example data
  example_input <- system.file("extdata", "example_input_small.csv", package = "genomesizeR")
  expected_results <- read.csv(
    system.file("extdata", "example_output_lmm_small.csv", package = "genomesizeR"),
    stringsAsFactors = FALSE
  )

  # Get path to reference data
  refdata_path <- system.file("extdata", "genomesizeRdata_v1.0.3.tar.gz", package = "genomesizeR")

  # Skip test if reference data is not available
  if (!file.exists(refdata_path)) {
    skip("Reference data not available")
  }

  # Run estimation with lmm method
  results <- estimate_genome_size(
    example_input,
    refdata_path,
    sep = ',',
    match_column = 'TAXID',
    output_format = 'input',
    method = 'lmm',
    ci_threshold = 0.5,
    n_cores = 1
  )

  # Test 1: Check that results is a data frame
  expect_s3_class(results, "data.frame")

  # Test 2: Compare estimated genome sizes with 5% tolerance (LMM model may vary slightly)
  expect_equal(
    results$estimated_genome_size,
    expected_results$estimated_genome_size,
    tolerance = 0.05,
    label = "estimated_genome_size"
  )

  # Test 3: Check that all estimated genome sizes are positive
  expect_true(all(results$estimated_genome_size > 0))

  # Test 4: Ensure lower CI < estimate < upper CI 
  valid_idx <- !is.na(results$confidence_interval_lower) & 
               !is.na(results$estimated_genome_size) & 
               !is.na(results$confidence_interval_upper)
  if (any(valid_idx)) {
    # Ensure columns are numeric
    res_lower <- as.numeric(results$confidence_interval_lower[valid_idx])
    res_est <- as.numeric(results$estimated_genome_size[valid_idx])
    res_upper <- as.numeric(results$confidence_interval_upper[valid_idx])
    
    # Check: lower CI should be less than or equal to estimate
    lower_check <- res_lower <= res_est
    expect_true(all(lower_check), 
                label = "lower CI <= estimate")
    
    # Check: upper CI should be greater than or equal to estimate
    upper_check <- res_est <= res_upper
    expect_true(all(upper_check), 
                label = "estimate <= upper CI")
  }

  # Test 5: Compare confidence interval bounds with expected (within 15% tolerance)
  exp_lower_num <- as.numeric(expected_results$confidence_interval_lower)
  exp_upper_num <- as.numeric(expected_results$confidence_interval_upper)
  res_lower_num <- as.numeric(results$confidence_interval_lower)
  res_upper_num <- as.numeric(results$confidence_interval_upper)
  
  # Only compare non-NA values
  valid_ci <- !is.na(res_lower_num) & !is.na(exp_lower_num)
  if (any(valid_ci)) {
    expect_equal(
      res_lower_num[valid_ci],
      exp_lower_num[valid_ci],
      tolerance = 0.15,
      label = "lower CI within 15% of expected"
    )
  }
  
  valid_ci_upper <- !is.na(res_upper_num) & !is.na(exp_upper_num)
  if (any(valid_ci_upper)) {
    expect_equal(
      res_upper_num[valid_ci_upper],
      exp_upper_num[valid_ci_upper],
      tolerance = 0.15,
      label = "upper CI within 15% of expected"
    )
  }
})
