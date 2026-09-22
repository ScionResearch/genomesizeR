# Unit tests: aggregation of per-query posterior distributions into per-sample distributions

# 3 queries x 4 draws
posterior <- matrix(c(1, 2, 3, 4,
                      10, 20, 30, 40,
                      100, 200, 300, 400),
                    nrow = 3, byrow = TRUE,
                    dimnames = list(NULL, paste0("posterior_", 1:4)))

# 3 queries x 2 samples
counts <- matrix(c(5, 0,
                   0, 3,
                   7, 3),
                 nrow = 3, byrow = TRUE,
                 dimnames = list(c("ASV_1", "ASV_2", "ASV_3"), c("S1", "S2")))

test_that("log_count weights give the expected weighted means per draw", {
  res <- aggregate_posterior_by_sample(posterior, counts, weights = "log_count")

  w1 <- log1p(c(5, 0, 7))
  expected_s1 <- colSums(posterior * w1) / sum(w1)
  w2 <- log1p(c(0, 3, 3))
  expected_s2 <- colSums(posterior * w2) / sum(w2)

  expect_equal(dim(res$posterior), c(2, 4))
  expect_equal(unname(res$posterior["S1", ]), unname(expected_s1))
  expect_equal(unname(res$posterior["S2", ]), unname(expected_s2))
  expect_equal(unname(res$n_queries_used), c(2, 2))
  expect_equal(unname(res$n_queries_total), c(2, 2))
  expect_equal(unname(res$abundance_covered), c(1, 1))
})

test_that("relative_abundance and presence weights are computed as expected", {
  res_ra <- aggregate_posterior_by_sample(posterior, counts, weights = "relative_abundance")
  expect_equal(unname(res_ra$posterior["S1", ]), unname((5 * posterior[1, ] + 7 * posterior[3, ]) / 12))

  res_pa <- aggregate_posterior_by_sample(posterior, counts, weights = "presence")
  expect_equal(unname(res_pa$posterior["S2", ]), unname((posterior[2, ] + posterior[3, ]) / 2))

  expect_error(aggregate_posterior_by_sample(posterior, counts, weights = "foo"))
})

test_that("queries without posterior are ignored and reported in abundance_covered", {
  posterior_na <- posterior
  posterior_na[3, ] <- NA

  res <- aggregate_posterior_by_sample(posterior_na, counts, weights = "presence")

  expect_equal(unname(res$posterior["S1", ]), unname(posterior[1, ]))
  expect_equal(unname(res$posterior["S2", ]), unname(posterior[2, ]))
  expect_equal(unname(res$n_queries_used), c(1, 1))
  expect_equal(unname(res$n_queries_total), c(2, 2))
  expect_equal(unname(res$abundance_covered), c(5 / 12, 3 / 6))
})

test_that("samples without usable query get NA", {
  posterior_na <- posterior
  posterior_na[c(1, 3), ] <- NA

  res <- aggregate_posterior_by_sample(posterior_na, counts)

  expect_true(all(is.na(res$posterior["S1", ])))
  expect_equal(unname(res$abundance_covered["S1"]), 0)
  expect_false(any(is.na(res$posterior["S2", ])))
})

test_that("estimate_genome_size_per_sample summarises the sample-level posterior", {
  output_table <- data.frame(ASVs = c("ASV_1", "ASV_2", "ASV_3"),
                             estimated_genome_size = rowMeans(posterior),
                             model_used = "bayesian Bacteria",
                             stringsAsFactors = FALSE)
  output_table <- cbind(output_table, as.data.frame(posterior))

  res <- estimate_genome_size_per_sample(output_table, counts, weights = "presence", ci_level = 0.5,
                                         return_posterior = TRUE)

  expected_s2 <- (posterior[2, ] + posterior[3, ]) / 2
  expect_equal(res$sample, c("S1", "S2"))
  expect_equal(res$estimated_genome_size[2], mean(expected_s2))
  expect_equal(res$confidence_interval_lower[2], unname(quantile(expected_s2, 0.25)))
  expect_equal(res$confidence_interval_upper[2], unname(quantile(expected_s2, 0.75)))
  expect_equal(unname(unlist(res[2, paste0("posterior_", 1:4)])), unname(expected_s2))

  # Rows of otu_table are matched by name to the ASVs column, whatever their order
  res_shuffled <- estimate_genome_size_per_sample(output_table, counts[c(3, 1, 2), ], weights = "presence")
  expect_equal(res_shuffled$estimated_genome_size, res$estimated_genome_size)

  # Transposed abundance table
  res_t <- estimate_genome_size_per_sample(output_table, t(counts), taxa_are_rows = FALSE, weights = "presence")
  expect_equal(res_t$estimated_genome_size, res$estimated_genome_size)

  expect_error(estimate_genome_size_per_sample(output_table[, 1:3], counts), "return_posterior")
})
