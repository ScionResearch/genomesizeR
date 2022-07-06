
.onLoad = function(libname, pkgname) {
  temp_dir = tempdir()
  taxonomy_dir = paste(temp_dir, 'taxdump', sep = '/')
  genome_size_db = paste(temp_dir, 'genome_size_db.csv', sep = '/')
  taxonomy_archive = system.file("extdata", 'taxdump.tar.gz', package = "genomesizeR")
  genome_size_db_archive = system.file("extdata", 'genome_size_db.tar.gz', package = "genomesizeR")
  assign('taxonomy_archive', taxonomy_archive, envir = topenv())
  assign('genome_size_db_archive', genome_size_db_archive, envir = topenv())
  assign('taxonomy_dir', taxonomy_dir, envir = topenv())
  assign('genome_size_db', genome_size_db, envir = topenv())
  assign('temp_dir', temp_dir, envir = topenv())
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


#' Infer genome sizes
#'
#' This function loads a query file and predicts genome sizes.
#'
#' @param queries Queries: path to file or table object
#' @param format Query format if in a file ('csv' (default) or 'dada2' (taxonomy table file))
#' @param sep If 'csv' format, column separator
#' @param match_column If 'csv' format, the column containing match information (with one or several matches)
#' @param match_sep If 'csv' format and several matches in match column, separator between matches
#' @param size_db Path to the genome size reference database to use if not the default one
#' @param taxonomy Path to taxonomy database (NCBI taxdump) to use if not the default one
#' @param output_format Format in which the output should be.
#'                      Default: "input" a data frame with the same columns as the input,
#'                      with the added columns: "estimated_genome_size" and "estimated_genome_size_confidence_interval".
#'                      Other formats available: "data.frame", a data frame with the same number of rows as the input, and 2 columns:
#'                      "estimated_genome_size" and "estimated_genome_size_confidence_interval".
#' @importFrom utils read.csv
#' @importFrom pbapply pbapply
#' @importFrom seqinr s2c
#' @export
estimate_genome_size <- function(queries, format='csv', sep=',', match_column=NA, match_sep=',',
                                 size_db=NA, taxonomy=NA, output_format='input') {

  cat("Reading queries", fill=T)
  if (typeof(queries) == "character") {
    query_table = read.csv(queries, stringsAsFactors = FALSE, sep=sep)
  }
  else if (typeof(queries) == "data.frame") {
    query_table = queries
  }

  queries = as.data.frame(queries, stringsAsFactors = F)

  cat("Reading genome size reference database", fill=T)
  # Check memory usage (map if too big)
  size_db = get_genome_size_db(size_db)
  size_db = read.csv(size_db, sep='\t', quote="", stringsAsFactors = FALSE)

  taxonomy = get_taxonomy(taxonomy)
  cat("Reading taxonomy names", fill=T)
  names = getnames(taxonomy)
  cat("Reading taxonomy nodes", fill=T)
  nodes = getnodes(taxonomy)

  cat("Computing genome sizes", fill=T)
  method = weighted_mean
  output_table = pbapply(query_table, 1, method, size_db=size_db, taxonomy=taxonomy, names=names, nodes=nodes, format=format, match_column=match_column, match_sep=match_sep)

  cat('Done', fill=T)

  # Handle R formatting issues
  output_table = as.data.frame(t(as.data.frame(output_table, stringsAsFactors = F)), stringsAsFactors = F)
  output_table$estimated_genome_size = as.numeric(output_table$estimated_genome_size)
  output_table$estimated_genome_size = as.numeric(output_table$estimated_genome_size)
  output_table$estimated_genome_size_confidence_interval = as.numeric(output_table$estimated_genome_size_confidence_interval)
  output_table$estimated_genome_size_confidence_interval = as.numeric(output_table$estimated_genome_size_confidence_interval)

  summary(output_table$estimated_genome_size)
  successful_estimations = nrow(output_table[! is.na(output_table$estimated_genome_size), ])*100 / nrow(output_table)

  cat("\n#############################################################################", fill=T)
  cat("# Genome size estimation summary:", fill=T)
  cat('#\n')
  cat('# ', successful_estimations, "% successful estimations", fill=T)
  cat('#\n')
  print(summary(output_table$estimated_genome_size))
  cat("#############################################################################", fill=T)

  # Reformat if needed
  if (output_format == "data.frame") {
    output_table = output_table[c('estimated_genome_size', 'estimated_genome_size_confidence_interval')]
  }

  return(output_table)
}
