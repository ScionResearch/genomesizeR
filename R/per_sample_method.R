
# Aggregate per-query posterior distributions into per-sample posterior distributions.
#
# posterior: numeric matrix, one row per query (ASV) and one column per posterior draw.
# counts: numeric matrix, one row per query (same order as posterior) and one column per sample.
# weights: how counts are converted to weights.
#
# Returns a list with:
#   posterior: matrix, one row per sample and one column per draw, holding the weighted mean genome size
#              across queries for that draw (NA for samples without any usable query)
#   n_queries_used: number of queries with a posterior and a non-zero count in each sample
#   n_queries_total: number of queries with a non-zero count in each sample
#   abundance_covered: proportion of the counts of each sample coming from queries with a posterior
aggregate_posterior_by_sample <- function(posterior, counts, weights='log_count') {

  if (nrow(posterior) != nrow(counts)) {
    stop("posterior and counts must have the same number of rows")
  }

  counts[is.na(counts)] = 0
  has_posterior = stats::complete.cases(posterior)

  if (weights == 'log_count') {
    w = log1p(counts)
  }
  else if (weights == 'relative_abundance') {
    w = counts
  }
  else if (weights == 'presence') {
    w = (counts > 0) * 1
  }
  else {
    stop("weights must be one of 'log_count', 'relative_abundance' or 'presence'")
  }

  n_queries_total = colSums(counts > 0)
  abundance_covered = colSums(counts[has_posterior, , drop=FALSE]) / colSums(counts)
  abundance_covered[colSums(counts) == 0] = NA

  w = w[has_posterior, , drop=FALSE]
  n_queries_used = colSums(w > 0)
  w_sum = colSums(w)
  w = sweep(w, 2, ifelse(w_sum > 0, w_sum, NA), '/')

  # Weighted mean across queries for each sample and each draw: (samples x queries) x (queries x draws)
  sample_posterior = t(w) %*% posterior[has_posterior, , drop=FALSE]
  sample_posterior[n_queries_used == 0, ] = NA
  colnames(sample_posterior) = colnames(posterior)
  rownames(sample_posterior) = colnames(counts)

  return(list(posterior=sample_posterior,
              n_queries_used=n_queries_used,
              n_queries_total=n_queries_total,
              abundance_covered=abundance_covered))
}


#' Estimate genome sizes per sample
#'
#' This function aggregates the genome size estimations of individual queries (e.g. ASVs)
#' into one genome size estimation per sample, propagating the estimation uncertainty.
#'
#' It requires the full posterior predictive distribution of each query, obtained by running
#' \code{\link{estimate_genome_size}} with the bayesian method and \code{return_posterior=TRUE}.
#' For each sample and each posterior draw, the mean genome size across queries is computed, weighted
#' by query abundance. The resulting sample-level posterior distribution is summarised into a mean
#' estimate and a confidence interval. Queries without an estimation (e.g. taxon not found) are ignored,
#' and the proportion of the sample abundance that they represent is reported in the column "abundance_covered".
#'
#' @param output_table Result table from \code{\link{estimate_genome_size}} run with \code{method='bayesian'}
#'                     and \code{return_posterior=TRUE}.
#' @param otu_table Abundance table with one row per query (in the same order as the rows of output_table,
#'                  or with row names matching the "ASVs" column of output_table when there is one) and one
#'                  column per sample (e.g. the OTU table of a phyloseq object, as a matrix or data frame).
#' @param taxa_are_rows Whether queries are rows and samples are columns in otu_table (default: TRUE).
#'                      If FALSE, otu_table is transposed.
#' @param weights How query abundances are used to weight the mean genome size of a sample:
#'                "log_count" (default): log(1 + count); "relative_abundance": count (proportional to the
#'                relative abundance of the query in the sample); "presence": 1 for every query present in the sample.
#' @param ci_level Confidence level of the interval computed from the sample-level posterior distribution (default: 0.95).
#' @param return_posterior If TRUE, the sample-level posterior distributions are returned in additional columns
#'                         "posterior_1", ..., "posterior_N" (one column per draw).
#' @return A data frame with one row per sample and the columns "sample", "estimated_genome_size" (mean of the
#'         sample-level posterior distribution, in base pairs), "confidence_interval_lower",
#'         "confidence_interval_upper", "n_queries_used" (queries with an estimation and present in the sample),
#'         "n_queries_total" (queries present in the sample) and "abundance_covered" (proportion of the sample
#'         abundance coming from queries with an estimation).
#' @importFrom stats quantile complete.cases
#' @export
estimate_genome_size_per_sample <- function(output_table, otu_table, taxa_are_rows=TRUE,
                                            weights='log_count', ci_level=0.95,
                                            return_posterior=FALSE) {

  posterior_columns = grep('^posterior_[0-9]+$', names(output_table), value=TRUE)
  if (length(posterior_columns) == 0) {
    stop("No posterior distribution found in output_table: run estimate_genome_size with method='bayesian' and return_posterior=TRUE")
  }
  if (ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be between 0 and 1")
  }

  counts = as.matrix(otu_table)
  if (!taxa_are_rows) {
    counts = t(counts)
  }
  if (!is.numeric(counts)) {
    stop("otu_table must contain numeric abundances")
  }
  if (is.null(colnames(counts))) {
    colnames(counts) = paste0('sample_', seq_len(ncol(counts)))
  }

  # Match queries and abundance rows, by name if possible, else by position
  if ('ASVs' %in% names(output_table) && !is.null(rownames(counts))) {
    idx = match(output_table$ASVs, rownames(counts))
    if (all(is.na(idx))) {
      stop("None of the 'ASVs' of output_table were found in the row names of otu_table")
    }
    if (any(is.na(idx))) {
      warning(sum(is.na(idx)), " queries of output_table were not found in otu_table and are ignored")
    }
    not_in_output = setdiff(rownames(counts), output_table$ASVs)
    if (length(not_in_output) > 0) {
      warning(length(not_in_output), " rows of otu_table have no estimation in output_table and are ignored")
    }
    counts = counts[ifelse(is.na(idx), 1, idx), , drop=FALSE]
    counts[is.na(idx), ] = 0
  }
  else if (nrow(counts) != nrow(output_table)) {
    stop("otu_table must have one row per query of output_table (", nrow(output_table), " rows), in the same order")
  }

  posterior = as.matrix(output_table[posterior_columns])
  mode(posterior) = 'numeric'

  cat("Computing genome sizes per sample", fill=T)
  aggregated = aggregate_posterior_by_sample(posterior, counts, weights=weights)

  probabilities = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  ci = t(apply(aggregated$posterior, 1, quantile, probs=probabilities, na.rm=TRUE))

  per_sample = data.frame(sample=colnames(counts),
                          estimated_genome_size=rowMeans(aggregated$posterior),
                          confidence_interval_lower=ci[, 1],
                          confidence_interval_upper=ci[, 2],
                          n_queries_used=aggregated$n_queries_used,
                          n_queries_total=aggregated$n_queries_total,
                          abundance_covered=aggregated$abundance_covered,
                          stringsAsFactors=F)
  row.names(per_sample) = NULL

  no_estimation = sum(is.na(per_sample$estimated_genome_size))
  if (no_estimation > 0) {
    warning(no_estimation, " samples have no query with an estimation")
  }

  cat("\n#############################################################################", fill=T)
  cat("# Genome size estimation per sample summary:", fill=T)
  cat('#\n')
  cat('# ', ncol(counts), " samples, ", no_estimation, " without estimation", fill=T)
  cat('#\n')
  print(summary(per_sample$estimated_genome_size))
  cat('\n# Proportion of sample abundance with an estimation:\n')
  print(summary(per_sample$abundance_covered))

  if (return_posterior) {
    sample_posterior = as.data.frame(aggregated$posterior)
    row.names(sample_posterior) = NULL
    per_sample = cbind(per_sample, sample_posterior)
  }

  return(per_sample)
}
