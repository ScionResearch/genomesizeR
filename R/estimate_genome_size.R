

.onLoad = function(libname, pkgname) {
  temp_dir = tempdir()
  taxonomy_dir = paste(temp_dir, 'taxdump', sep = '/')
  genome_size_db = paste(temp_dir, 'genome_size_db.csv', sep = '/')
  genome_size_db_for_hierarchical = paste(temp_dir, 'genome_size_db_for_hierarchical.csv', sep = '/')
  taxonomy_archive = system.file("extdata", 'taxdump.tar.gz', package = "genomesizeR")
  genome_size_db_archive = system.file("extdata", 'genome_size_db.tar.gz', package = "genomesizeR")
  genome_size_db_for_hierarchical_archive = system.file("extdata", 'genome_size_db_for_hierarchical.tar.gz', package = "genomesizeR")
  bayesian_model_bact = system.file("extdata/fits", 'm_superkingdom2.rds', package = "genomesizeR")
  bayesian_model_euka = system.file("extdata/fits", 'm_superkingdom2759.rds', package = "genomesizeR")
  bayesian_model_arch = system.file("extdata/fits", 'm_superkingdom2157.rds', package = "genomesizeR")
  assign('bayesian_model_bact', bayesian_model_bact, envir = topenv())
  assign('bayesian_model_euka', bayesian_model_euka, envir = topenv())
  assign('bayesian_model_arch', bayesian_model_arch, envir = topenv())
  assign('taxonomy_archive', taxonomy_archive, envir = topenv())
  assign('genome_size_db_archive', genome_size_db_archive, envir = topenv())
  assign('genome_size_db_for_hierarchical_archive', genome_size_db_for_hierarchical_archive, envir = topenv())
  assign('taxonomy_dir', taxonomy_dir, envir = topenv())
  assign('genome_size_db', genome_size_db, envir = topenv())
  assign('genome_size_db_for_hierarchical', genome_size_db_for_hierarchical, envir = topenv())
  assign('temp_dir', temp_dir, envir = topenv())
}


build_model <- function(size_db, effects) {

  f <- as.formula(
    paste("log(genome.size) ~ (1|",
    paste(effects, collapse = " / "),
    ")"))

  model <- tryCatch(
    {
      lmer(f, data = size_db)
    },
    error=function(cond) {
      message("Error building hierarchical model:")
      message(cond)
      return(NA)
    },
    warning=function(cond) {
      return(NA)
    }
  )
  return(model)
}


#' Read bayesian model from rds
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


#' Read taxonomy database
#'
#' @param taxonomy_path Path to taxonomy database file or NA
#' @importFrom utils untar
get_taxonomy <- function(taxonomy_path) {
  if (is.na(taxonomy_path)) {
    if ( ! dir.exists(taxonomy_dir)) {
      cat("Untarring taxonomy", fill=T)
      untar(taxonomy_archive, exdir=temp_dir)
    }
    taxonomy_path = taxonomy_dir
  }
  cat("Using taxonomy:", taxonomy_path, fill=T)
  return(taxonomy_path)
}


#' Read genome size database
#'
#' @param genome_size_db_path Path to genome size database file or NA
#' @importFrom utils untar
get_genome_size_db <- function(genome_size_db_path) {
  if (is.na(genome_size_db_path)) {
    if ( ! dir.exists(genome_size_db)) {
      cat("Untarring genome size reference database", fill=T)
      untar(genome_size_db_archive, exdir=temp_dir)
    }
    genome_size_db_path = genome_size_db
  }
  cat("Using genome size reference database:", genome_size_db_path, fill=T)
  return(genome_size_db_path)
}


#' Read genome size database for hierarchical model
#'
#' @param genome_size_db_path Path to genome size database file or NA
#' @importFrom utils untar
get_genome_size_db_for_hierarchical <- function(genome_size_db_path) {
  if (is.na(genome_size_db_path)) {
    if ( ! dir.exists(genome_size_db_for_hierarchical)) {
      cat("Untarring genome size reference database", fill=T)
      cat(genome_size_db_for_hierarchical_archive, fill=T)
      untar(genome_size_db_for_hierarchical_archive, exdir=temp_dir)
    }
    genome_size_db_path = genome_size_db_for_hierarchical
  }
  cat("Using genome size reference database:", genome_size_db_path, fill=T)
  return(genome_size_db_path)
}


#' Infer genome sizes
#'
#' This function loads a query file and predicts genome sizes.
#'
#' @param queries Queries: path to file or table object
#' @param format Query format if in a file ('csv' (default), 'tax_table' or 'biom' (taxonomy table files))
#' @param sep If 'csv' format, column separator
#' @param match_column If 'csv' format, the column containing match information (with one or several matches)
#' @param match_sep If 'csv' format and several matches in match column, separator between matches
#' @param size_db Path to the genome size reference database to use if not the default one
#' @param taxonomy Path to taxonomy database (NCBI taxdump) to use if not the default one
#' @param output_format Format in which the output should be.
#'                      Default: "input" a data frame with the same columns as the input,
#'                      with the added columns: "estimated_genome_size" and "estimated_genome_size_confidence_interval".
#'                      Other formats available: "data.frame", a data frame with the same number of rows as the input, and 3 columns:
#'                      "TAXID", "estimated_genome_size" and "estimated_genome_size_confidence_interval".
#' @param method Method to use for genome size estimation, 'bayesian' (default), 'weighted_mean' or 'hierarchical'
#' @param ci_threshold Threshold for the confidence interval as a proportion of the guessed size
#'                     (e.g. 0.2 means that estimations with a confidence interval that represents more than 20% of
#'                     the guessed size will be discarded)
#' @param n_cores Number of CPU cores to use (default is 'half': half of all available cores)
#' @importFrom utils read.csv
#' @importFrom pbapply pbapply
#' @importFrom seqinr s2c
#' @importFrom lme4 lmer
#' @importFrom dplyr bind_rows
#' @importFrom merTools predictInterval
#' @importFrom parallel detectCores
#' @importFrom biomformat read_biom observation_metadata
#' @importFrom stats as.formula na.omit predict median
#' @importFrom brms posterior_predict
#' @importFrom data.table fread
#' @import doParallel
# @import parabar
#' @export
estimate_genome_size <- function(queries, format='csv', sep=',', match_column=NA, match_sep=';',
                                 size_db=NA, taxonomy=NA, output_format='input', method='bayesian',
                                 ci_threshold=0.2, prediction_variables=c('family', 'genus'), n_cores='half') {

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
  else if (typeof(queries) == "character") {
    query_table = read.csv(queries, stringsAsFactors = FALSE, sep=sep)
    queries = as.data.frame(query_table, stringsAsFactors = F)
  }
  else {
    query_table = queries
    queries = as.data.frame(query_table, stringsAsFactors = F)
  }

  cat("Reading genome size reference database", fill=T)

  na_models = NA
  genusfamily_model = NA
  genusorder_model = NA
  familyorder_model = NA
  bayes_model_bact = NA
  bayes_model_euka = NA
  bayes_model_arch = NA

  if (method == 'bayesian') {
    bayes_model_bact = get_bayes_model('Bacteria')
    bayes_model_euka = get_bayes_model('Eukaryota')
    bayes_model_arch = get_bayes_model('Archaea')
  }

  if (method == 'hierarchical') {
    size_db_h = get_genome_size_db_for_hierarchical(size_db)
    size_db_h = read.table(size_db_h, sep=",", header=TRUE, na.strings="None", stringsAsFactors=TRUE,  quote="", fill=FALSE)
    size_db_h$order = as.factor(size_db_h$order)
    size_db_h$family = as.factor(size_db_h$family)
    size_db_h$genus = as.factor(size_db_h$genus)
    size_db_h$species = as.factor(size_db_h$species)
    size_db_h$genome.size = as.integer(as.character(size_db_h$genome_size))
    genusfamily_size_db = na.omit(size_db_h[, c("genome.size", "family", "genus")])
    genusorder_size_db = na.omit(size_db_h[, c("genome.size", "genus", "order")])
    familyorder_size_db = na.omit(size_db_h[, c("genome.size", "family", "order")])
    genusfamily_model = build_model(genusfamily_size_db, c('genus', 'family'))
    genusorder_model = build_model(genusorder_size_db, c('genus', 'order'))
    familyorder_model = build_model(familyorder_size_db, c('family', 'order'))
    na_models = c(0,0,0)
    if (typeof(genusfamily_model) == 'logical' && is.na(genusfamily_model)) {
      na_models[1] = 1
    }
    if (typeof(genusorder_model) == 'logical' && is.na(genusorder_model)) {
      na_models[2] = 1
    }
    if (typeof(familyorder_model) == 'logical' && is.na(familyorder_model)) {
      na_models[3] = 1
    }
  }
  full_size_db = get_genome_size_db(size_db)
  full_size_db = read.csv(full_size_db, sep='\t', quote="", stringsAsFactors = FALSE)

  taxonomy = get_taxonomy(taxonomy)
  cat("Reading taxonomy names", fill=T)
  names = getnames(taxonomy)
  cat("Reading taxonomy nodes", fill=T)
  nodes = getnodes(taxonomy)
  alltax = parseNCBITaxonomy(taxonomy)

  cat("Computing genome sizes", fill=T)
  #options(warn=2)

  if (n_cores == 'half') {
    n_cores = parallel::detectCores() / 2
  }
  cat("Using ", n_cores, " cores", fill=T)
  #cluster = parallel::makeCluster(n_cores)
  #parallel::clusterExport(cluster, c("alltax", "nodes", "names"), envir=environment())

  output_table = try(pbapply(queries, 1, method,
                         models=list('genusfamily_model'=genusfamily_model, 'genusorder_model'=genusorder_model,
                                     'familyorder_model'=familyorder_model,
                                     'bayes_model_bact'=bayes_model_bact, 'bayes_model_euka'=bayes_model_euka,
                                     'bayes_model_arch'=bayes_model_arch),
                         na_models=na_models, size_db=full_size_db, taxonomy=taxonomy,
                         names=names, nodes=nodes, alltax=alltax, format=format, output_format=output_format, match_column=match_column,
                         match_sep=match_sep, ci_threshold=ci_threshold, cl=n_cores))


  # parabar::set_option("progress_track", TRUE)
  # backend <- parabar::start_backend(cores = 4, cluster_type = "psock", backend_type = "async")
  # parabar::export(backend, c("method"), environment())
  # parabar::export(backend, c("bayesian"), .GlobalEnv)
  # print(parabar::peek(backend))
  # output_table = try(parabar::par_apply(backend, queries, 1, method,
  #                            models=list('genusfamily_model'=genusfamily_model, 'genusorder_model'=genusorder_model,
  #                                        'familyorder_model'=familyorder_model,
  #                                        'bayes_model_bact'=bayes_model_bact, 'bayes_model_euka'=bayes_model_euka,
  #                                        'bayes_model_arch'=bayes_model_arch),
  #                            na_models=na_models, size_db=full_size_db, taxonomy=taxonomy,
  #                            names=names, nodes=nodes, alltax=alltax, format=format, output_format=output_format, match_column=match_column,
  #                            match_sep=match_sep, ci_threshold=ci_threshold))

  cat('Done', fill=T)

  #parabar::stop_backend(backend)
  #parallel::stopCluster(cl = cluster)

  # Handle R formatting issues
  if (method == 'hierarchical') {
    output_table = as.data.frame(bind_rows(output_table), stringsAsFactors = F)
    confidence_interval = exp(compute_confidence_interval_hierarchical(output_table, genusfamily_model, n_cores))
    output_table$confidence_interval_lower = as.numeric(confidence_interval$lwr)
    output_table$confidence_interval_upper = as.numeric(confidence_interval$upr)
    new_queries = output_table
    output_table = pbapply(new_queries, 1, ci_post_treat, ci_threshold=ci_threshold)
    output_table = as.data.frame(bind_rows(output_table), stringsAsFactors = F)
  }
  else if (method == 'bayesian') {
    output_table = as.data.frame(bind_rows(output_table), stringsAsFactors = F)
    new_queries = output_table
    output_table = pbapply(new_queries, 1, ci_post_treat, ci_threshold=ci_threshold)
    output_table = as.data.frame(bind_rows(output_table), stringsAsFactors = F)
  }
  else {
    output_table = as.data.frame(t(as.data.frame(output_table, stringsAsFactors = F)), stringsAsFactors = F)
    output_table$estimated_genome_size_confidence_interval = as.numeric(output_table$estimated_genome_size_confidence_interval)
    output_table$genome_size_estimation_distance = as.numeric(output_table$genome_size_estimation_distance)
  }

  output_table$estimated_genome_size = as.numeric(output_table$estimated_genome_size)
  #output_table$data_density = as.numeric(output_table$data_density)

  summary(output_table$estimated_genome_size)

  successful_estimations = nrow(output_table[output_table$genome_size_estimation_status=='OK', ])*100 / nrow(output_table)

  cat("\n#############################################################################", fill=T)
  cat("# Genome size estimation summary:", fill=T)
  cat('#\n')
  cat('# ', successful_estimations, "% successful estimations", fill=T)
  cat('#\n')
  print(summary(output_table$estimated_genome_size))
  cat('\n# Estimation status:')
  print(table(output_table$genome_size_estimation_status))

  if (method == 'weighted_mean') {
    cat('\n# Estimation at taxonomic rank:')
    print(table(output_table$genome_size_estimation_rank))
    cat('\n# Estimation distance:', fill=T)
    cat('#\n')
    print(summary(output_table$genome_size_estimation_distance))
    cat("#############################################################################", fill=T)
  }

  # Reformat if needed
  if (output_format == "data.frame") {
    output_table = output_table[c('LCA', 'estimated_genome_size', 'estimated_genome_size_confidence_interval', 'genome_size_estimation_status')]
    names(output_table)[names(output_table) == 'LCA'] = 'TAXID'
  }

  return(output_table)
}
