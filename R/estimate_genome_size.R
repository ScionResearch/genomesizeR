

.onLoad = function(libname, pkgname) {
  temp_dir = tempdir()
  refdata_dir = paste(temp_dir, 'refdata', sep = '/')
  extdata_dir = paste(temp_dir, 'refdata', 'extdata', sep = '/')
  taxonomy_dir = paste(temp_dir, 'taxdump', sep = '/')
  genome_size_db = paste(temp_dir, 'genome_size_db.csv', sep = '/')
  genome_size_db_for_lmm = paste(temp_dir, 'genome_size_db_for_lmm.csv', sep = '/')
  taxonomy_archive = paste(extdata_dir, 'taxdump.tar.gz', sep = '/')
  bayesian_model_bact = paste(extdata_dir, 'fits', 'm_superkingdom2.rds', sep = '/')
  bayesian_model_euka = paste(extdata_dir, 'fits', 'm_superkingdom2759.rds', sep = '/')
  bayesian_model_arch = paste(extdata_dir, 'fits', 'm_superkingdom2157.rds', sep = '/')
  genome_size_db_archive = paste(extdata_dir, 'genome_size_db.tar.gz', sep = '/')
  genome_size_db_for_lmm_archive = paste(extdata_dir, 'genome_size_db_for_lmm.tar.gz', sep = '/')

  assign('temp_dir', temp_dir, envir = topenv())
  assign('refdata_dir', refdata_dir, envir = topenv())
  assign('bayesian_model_bact', bayesian_model_bact, envir = topenv())
  assign('bayesian_model_euka', bayesian_model_euka, envir = topenv())
  assign('bayesian_model_arch', bayesian_model_arch, envir = topenv())
  assign('taxonomy_archive', taxonomy_archive, envir = topenv())
  assign('genome_size_db_archive', genome_size_db_archive, envir = topenv())
  assign('genome_size_db_for_lmm_archive', genome_size_db_for_lmm_archive, envir = topenv())
  assign('taxonomy_dir', taxonomy_dir, envir = topenv())
  assign('genome_size_db', genome_size_db, envir = topenv())
  assign('genome_size_db_for_lmm', genome_size_db_for_lmm, envir = topenv())
}


#' Build a frequentist random effect model using nested effects
#'
#' @param size_db Genome size reference database
#' @param effects Vector of nested effects to use in the formula e.g. c("family", "genus")
build_lmm_model <- function(size_db, effects) {
  f <- as.formula(
    paste("log(genome.size) ~ (1|",
    paste(effects, collapse = " / "),
    ")"))

  model <- tryCatch(
    {
      lmer(f, data = size_db)
    },
    error=function(cond) {
      message("Error building lmm model:")
      message(cond)
      return(NA)
    },
    warning=function(cond) {
      return(NA)
    }
  )
  return(model)
}


#' Read bayesian model from rds file
#'
#' @param superkingdom Target superkingdom, taxid or name
get_bayes_model <- function(superkingdom) {
  if (superkingdom == "Bacteria" | superkingdom == 2) {
    bmodel = readRDS(bayesian_model_bact)
  }
  else if (superkingdom == "Archaea" | superkingdom == 2157) {
    bmodel = readRDS(bayesian_model_arch)
  }
  else if (superkingdom == "Eukaryota" | superkingdom == 2759) {
    bmodel = readRDS(bayesian_model_euka)
  }
  else {
    stop("A valid superkingdom must be specified")
  }
  return(bmodel)
}

#' Untar reference data archive containing reference databases and bayesion models
#'
#' @param refdata_path Path to tar.gz archive
#' @importFrom utils untar
get_refdata <- function(refdata_path) {
  if (dir.exists(refdata_dir)) {
    unlink(refdata_dir, recursive = TRUE)
  }
  cat("Untarring reference data", fill=T)
  untar(refdata_path, exdir=refdata_dir)
  cat("Using reference data in:", refdata_dir, fill=T)
  return(refdata_dir)
}


#' Read taxonomy database
#'
#' @param taxonomy_path Path to taxonomy database file or NA
#' @importFrom utils untar
get_taxonomy <- function(taxonomy_path=NA) {
  if (is.na(taxonomy_path)) {
    if (dir.exists(taxonomy_dir)) {
      unlink(taxonomy_dir, recursive = TRUE)
    }
    cat("Untarring taxonomy", fill=T)
    untar(taxonomy_archive, exdir=taxonomy_dir)
    taxonomy_path = taxonomy_dir
  }
  cat("Using taxonomy:", taxonomy_path, fill=T)
  return(taxonomy_path)
}


#' Read genome size database
#'
#' @param genome_size_db_path Path to genome size database file or NA
#' @importFrom utils untar
get_genome_size_db <- function(genome_size_db_path=NA) {
  if (is.na(genome_size_db_path)) {
    if (dir.exists(genome_size_db)) {
      unlink(genome_size_db, recursive = TRUE)
    }
    cat("Untarring genome size reference database", fill=T)
    untar(genome_size_db_archive, exdir=temp_dir)
    genome_size_db_path = genome_size_db
  }
  cat("Using genome size reference database:", genome_size_db_path, fill=T)
  return(genome_size_db_path)
}


#' Read genome size database for lmm model
#'
#' @param genome_size_db_path Path to genome size database file or NA
#' @importFrom utils untar
get_genome_size_db_for_lmm <- function(genome_size_db_path=NA) {
  if (is.na(genome_size_db_path)) {
    if (dir.exists(genome_size_db_for_lmm)) {
      unlink(genome_size_db_for_lmm, recursive = TRUE)
    }
    cat("Untarring genome size reference database", fill=T)
    cat(genome_size_db_for_lmm_archive, fill=T)
    untar(genome_size_db_for_lmm_archive, exdir=temp_dir)
    genome_size_db_path = genome_size_db_for_lmm
  }
  cat("Using genome size reference database:", genome_size_db_path, fill=T)
  return(genome_size_db_path)
}


#' Estimate genome sizes
#'
#' This function loads a query file or table and an archive containing
#' reference databases and bayesian models, and predicts genome sizes.
#'
#' @param queries Queries: path to csv or BIOM file, or R object used for input.
#' @param refdata_path Path to the downloadable archive containing the reference databases and the bayesian models.
#' @param format Input format: "csv" for csv file (default), "tax_table" for taxonomy table file or object as used in e.g. phyloseq,
#' "biom" for BIOM file, "dataframe" for a table-style object (e.g. data.frame or matrix object), "vector" for a vector object.
#' @param sep If table-style or csv format, column separator (default: ",").
#' @param match_column If table-style or csv format, the column containing match information (with one or several matches).
#' @param match_sep If table-style or csv format and several matches in match column, separator between matches (default: ";").
#' @param output_format Format in which the output should be.
#'                      Default: "input" a data frame with the same columns as the input,
#'                      with the added columns: "TAXID", "estimated_genome_size", "confidence_interval_lower",
#'                      "confidence_interval_upper", "genome_size_estimation_status", "model_used", as well as taxids at all ranks.
#'                      Other formats available: "data.frame", a data frame with only the previously described columns, without the
#'                      "taxids at all ranks" columns.
#' @param method Method to use for genome size estimation, 'bayesian' (default), 'weighted_mean' or 'lmm'.
#' @param ci_threshold Threshold for the confidence interval as a proportion of the predicted size
#'                     (e.g. 0.3 means that estimations with a confidence interval that represents more than 30% of
#'                     the predicted size will be tagged in the output table).
#' @param n_cores Number of CPU cores to use (default is 'half': half of all available cores).
#' @importFrom utils read.csv
#' @importFrom pbapply pbapply
#' @importFrom lme4 lmer
#' @importFrom dplyr bind_rows
#' @importFrom merTools predictInterval
#' @importFrom biomformat read_biom observation_metadata
#' @importFrom stats as.formula na.omit predict median
#' @importFrom brms posterior_predict
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @export
estimate_genome_size <- function(queries, refdata_path,
                                 format='csv', sep=',', match_column=NA, match_sep=';',
                                 output_format='input',
                                 method='bayesian',
                                 ci_threshold=0.3,
                                 n_cores='half') {

  options(warn=1)

  cat("Reading queries", fill=T)
  if (format == 'biom') {
    biom_table = read_biom(queries)
    biom_table = observation_metadata(biom_table)
    query_table = biom_table
    for(i in 1:ncol(query_table)) {
      s = query_table[,i]
      query_table[,i] = substr(s, 4, nchar(s))
    }
    queries = as.data.frame(query_table, stringsAsFactors = F)
  }
  else if (format == 'csv') {
    query_table = read.csv(queries, stringsAsFactors = FALSE, sep=sep)
    queries = as.data.frame(query_table, stringsAsFactors = F)
  }
  else if (format == 'tax_table') {
    if (class(queries) == "character") {
      query_table = read.csv(queries, stringsAsFactors = FALSE, sep=sep)
      queries = as.data.frame(query_table, stringsAsFactors = F)
    }
    else {
      query_table = queries
      queries = as.data.frame(query_table, stringsAsFactors = F)
    }
  }
  else if (format == 'dataframe' || format == 'vector') {
    query_table = queries
    queries = as.data.frame(query_table, stringsAsFactors = F)
  }
  else {
    cat("Input format not supported", fill=T)
    return()
  }

  cat("Reading genome size reference database", fill=T)

  na_models = NA
  genusfamily_model = NA
  bayes_model_bact = NA
  bayes_model_euka = NA
  bayes_model_arch = NA

  refdata_path = get_refdata(refdata_path)

  if (method == 'bayesian') {
    bayes_model_bact = get_bayes_model('Bacteria')
    bayes_model_euka = get_bayes_model('Eukaryota')
    bayes_model_arch = get_bayes_model('Archaea')
  }

  if (method == 'lmm') {
    size_db_lmm = get_genome_size_db_for_lmm()
    size_db_lmm = read.table(size_db_lmm, sep=",", header=TRUE, stringsAsFactors=TRUE, fill=FALSE, row.names = NULL, check.names = FALSE)
    size_db_lmm$order = as.factor(size_db_lmm$order)
    size_db_lmm$family = as.factor(size_db_lmm$family)
    size_db_lmm$genus = as.factor(size_db_lmm$genus)
    size_db_lmm$species = as.factor(size_db_lmm$species)
    size_db_lmm$genome.size = trunc(round(size_db_lmm$TOTAL_LENGTH))
    genusfamily_size_db = na.omit(size_db_lmm[, c("genome.size", "family", "genus")])
    genusfamily_model = build_lmm_model(genusfamily_size_db, c("family", "genus"))
    na_models = c(0)
    if (typeof(genusfamily_model) == 'logical' && is.na(genusfamily_model)) {
      na_models[1] = 1
    }
  }
  full_size_db = get_genome_size_db()
  full_size_db = read.csv(full_size_db, sep='\t', quote="", stringsAsFactors = FALSE)

  taxonomy = get_taxonomy()
  cat("Reading taxonomy names", fill=T)
  names = getnames(taxonomy)
  cat("Reading taxonomy nodes", fill=T)
  nodes = getnodes(taxonomy)
  alltax = parseNCBITaxonomy(taxonomy)

  cat("Computing genome sizes", fill=T)

  if (n_cores == 'half') {
    n_cores = parallel::detectCores() / 2
  }
  cat("Using ", n_cores, " cores", fill=T)

  output_table = try(pbapply(queries, 1, method,
                         models=list('genusfamily_model'=genusfamily_model,
                                     'bayes_model_bact'=bayes_model_bact, 'bayes_model_euka'=bayes_model_euka,
                                     'bayes_model_arch'=bayes_model_arch),
                         na_models=na_models, size_db=full_size_db, taxonomy=taxonomy,
                         names=names, nodes=nodes, alltax=alltax, format=format, output_format=output_format, match_column=match_column,
                         match_sep=match_sep, ci_threshold=ci_threshold, cl=n_cores))

  cat('Done', fill=T)

  # Compute confidence interval if needed and handle formatting
  if (typeof(output_table) == 'character') {
    output_table = data.frame(t(data.frame(output_table)))
  }
  else {
    output_table = as.data.frame(bind_rows(output_table), stringsAsFactors = F)
  }
  if (method == 'lmm') {
    confidence_interval = exp(compute_confidence_interval_lmm(output_table, genusfamily_model, n_cores))
    output_table$confidence_interval_lower = as.numeric(confidence_interval$lwr)
    output_table$confidence_interval_upper = as.numeric(confidence_interval$upr)
    new_queries = output_table
    output_table = pbapply(new_queries, 1, ci_post_treat, ci_threshold=ci_threshold)
    output_table = as.data.frame(bind_rows(output_table), stringsAsFactors = F)
  }
  else if (method == 'bayesian') {
    output_table = pbapply(output_table, 1, ci_post_treat, ci_threshold=ci_threshold)
    output_table = as.data.frame(bind_rows(output_table), stringsAsFactors = F)
    output_table$confidence_interval_lower = as.numeric(output_table$confidence_interval_lower)
    output_table$confidence_interval_upper = as.numeric(output_table$confidence_interval_upper)
  }
  else {
    output_table$confidence_interval_lower = as.numeric(output_table$confidence_interval_lower)
    output_table$confidence_interval_upper = as.numeric(output_table$confidence_interval_upper)
    output_table$genome_size_estimation_distance = as.numeric(output_table$genome_size_estimation_distance)
  }

  output_table$estimated_genome_size = as.numeric(output_table$estimated_genome_size)

  # Drop NA columns if there are any
  all_na_cols = sapply(output_table, \(x) all(is.na(x)))
  output_table = output_table[!all_na_cols]

  # Rename rows
  row.names(output_table) = paste0('query_', 1:nrow(output_table))

  summary(output_table$estimated_genome_size)

  successful_estimations = nrow(output_table[output_table$genome_size_estimation_status=='OK', ])*100 / nrow(output_table)

  cat("\n#############################################################################", fill=T)
  cat("# Genome size estimation summary:", fill=T)
  cat('#\n')
  cat('# ', successful_estimations, "% estimations achieving required precision", fill=T)
  cat('#\n')
  print(summary(output_table$estimated_genome_size))
  cat('\n# Estimation status:')
  print(table(output_table$genome_size_estimation_status))

  # Minimal dataframe if dataframe output requested
  if (output_format == "data.frame") {
    output_table = output_table[c('LCA', 'estimated_genome_size',
                                  'confidence_interval_lower', 'confidence_interval_upper',
                                  'genome_size_estimation_status', 'model_used')]
    names(output_table)[names(output_table) == 'LCA'] = 'TAXID'
  }

  return(output_table)
}
