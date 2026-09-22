# E2E Test: the lmm method does not estimate queries whose family is absent from its reference,
# and reference-mean queries keep their own confidence interval

test_that("lmm method refuses families absent from the reference and keeps reference-mean intervals", {
  skip_on_cran()

  refdata_path <- system.file("extdata", "genomesizeRdata_v1.0.3.tar.gz", package = "genomesizeR")
  if (!file.exists(refdata_path)) {
    skip("Reference data not available")
  }

  # Coronaviridae: a family with genomes in RefSeq but absent from the lmm reference (cellular organisms only)
  # SARS-CoV-2: species with a reference mean
  # Mycobacterium: genus estimated with the lmm model
  queries <- data.frame(TAXID = c(11118, 2697049, 1763))

  results <- estimate_genome_size(
    queries,
    refdata_path,
    format = 'dataframe',
    match_column = 'TAXID',
    method = 'lmm',
    n_cores = 1
  )

  # Unknown family: no estimate, explicit status
  expect_true(is.na(results$estimated_genome_size[1]))
  expect_equal(results$genome_size_estimation_status[1], "No reference for family in lmm model")

  # Reference mean: the interval is not replaced by a model interval
  # (this species has no standard error in the reference, so its interval is NA)
  expect_equal(results$model_used[2], "reference_mean")
  expect_true(is.na(results$confidence_interval_lower[2]) ||
                results$confidence_interval_lower[2] > results$estimated_genome_size[2] / 2)
  expect_true(is.na(results$confidence_interval_upper[2]) ||
                results$confidence_interval_upper[2] < 2 * results$estimated_genome_size[2])

  # lmm estimate with a model interval around it
  expect_equal(results$model_used[3], "lmm|family/genus")
  expect_true(is.numeric(results$confidence_interval_lower))
  expect_lt(results$confidence_interval_lower[3], results$estimated_genome_size[3])
  expect_gt(results$confidence_interval_upper[3], results$estimated_genome_size[3])
})
