# E2E Test: estimate_genome_size with return_posterior=TRUE, then estimate_genome_size_per_sample,
# on the small example dataset

test_that("estimate_genome_size returns the posterior and estimate_genome_size_per_sample aggregates it (small dataset)", {
  skip_on_cran()

  example_input <- system.file("extdata", "example_input_small.csv", package = "genomesizeR")
  refdata_path <- system.file("extdata", "genomesizeRdata_v1.0.3.tar.gz", package = "genomesizeR")
  if (!file.exists(refdata_path)) {
    skip("Reference data not available")
  }

  n_draws <- 200
  results <- estimate_genome_size(
    example_input,
    refdata_path,
    sep = ',',
    match_column = 'TAXID',
    output_format = 'input',
    method = 'bayesian',
    ci_threshold = 0.5,
    n_cores = 1,
    return_posterior = TRUE,
    n_draws = n_draws
  )

  posterior_columns <- grep("^posterior_[0-9]+$", names(results), value = TRUE)
  posterior <- as.matrix(results[posterior_columns])

  # One posterior column per draw, all numeric
  expect_length(posterior_columns, n_draws)
  expect_true(is.numeric(posterior))

  # Estimate and confidence interval are the summaries of the returned posterior
  estimated <- !is.na(results$estimated_genome_size)
  reference_mean <- !is.na(results$model_used) & results$model_used == "reference_mean"
  bayesian <- estimated & !reference_mean
  expect_true(any(bayesian))
  expect_equal(unname(rowMeans(posterior[estimated, , drop = FALSE])), results$estimated_genome_size[estimated],
               tolerance = 1e-6)
  expect_equal(unname(apply(posterior[bayesian, , drop = FALSE], 1, quantile, probs = 0.025)),
               results$confidence_interval_lower[bayesian],
               tolerance = 1e-6)

  # Reference-mean queries have a constant posterior
  if (any(reference_mean)) {
    expect_equal(unname(apply(posterior[reference_mean, , drop = FALSE], 1, function(x) length(unique(x)))),
                 rep(1, sum(reference_mean)))
  }

  # Build the abundance table (queries x samples) from the COUNT and SAMPLE columns of the example
  samples <- unique(results$SAMPLE)
  otu_table <- sapply(samples, function(s) ifelse(results$SAMPLE == s, as.numeric(results$COUNT), 0))
  rownames(otu_table) <- rownames(results)

  per_sample <- estimate_genome_size_per_sample(results, otu_table)

  expect_s3_class(per_sample, "data.frame")
  expect_equal(per_sample$sample, samples)
  expect_true(all(per_sample$n_queries_used <= per_sample$n_queries_total))
  expect_true(all(per_sample$abundance_covered >= 0 & per_sample$abundance_covered <= 1))
  valid <- !is.na(per_sample$estimated_genome_size)
  expect_true(all(per_sample$confidence_interval_lower[valid] <= per_sample$estimated_genome_size[valid]))
  expect_true(all(per_sample$estimated_genome_size[valid] <= per_sample$confidence_interval_upper[valid]))

  # A sample's estimate lies within the range of its queries' estimates (up to floating point error)
  for (s in samples[valid]) {
    in_sample <- results$SAMPLE == s & estimated
    expect_gte(per_sample$estimated_genome_size[per_sample$sample == s], min(results$estimated_genome_size[in_sample]) * (1 - 1e-9))
    expect_lte(per_sample$estimated_genome_size[per_sample$sample == s], max(results$estimated_genome_size[in_sample]) * (1 + 1e-9))
  }
})
