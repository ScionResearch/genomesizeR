
ci_post_treat <- function(query, ci_threshold) {
  if (!is.na(query['model_used']) && query['model_used'] == 'reference_mean') {
    query['genome_size_estimation_status'] = 'OK'
    query = as.data.frame(t(as.data.frame(query)))
    return(query)
  }
  confidence_interval_lwr = as.numeric(query['confidence_interval_lower'])
  confidence_interval_upr = as.numeric(query['confidence_interval_upper'])
  estimated_genome_size = as.numeric(query['estimated_genome_size'])
  if (is.na(estimated_genome_size)) {
    query = as.data.frame(t(as.data.frame(query)))
    return(query)
  }
  if (is.na(confidence_interval_lwr) || is.na(confidence_interval_upr)) {
    query['genome_size_estimation_status'] = 'Could not compute confidence interval'
    query = as.data.frame(t(as.data.frame(query)))
    return(query)
  }
  if (confidence_interval_upr - estimated_genome_size > ci_threshold*estimated_genome_size || estimated_genome_size - confidence_interval_lwr > ci_threshold*estimated_genome_size) {
    query['genome_size_estimation_status'] = 'Confidence interval to estimated size ratio > ci_threshold'
    query = as.data.frame(t(as.data.frame(query)))
    return(query)
  }
  query['genome_size_estimation_status'] = 'OK'

  query = as.data.frame(t(as.data.frame(query)))
  return(query)
}


compute_LCA <- function(taxids, nodes) {
  if (length(taxids) > 1) {
    parents = lapply(taxids, allparents, taxdir=NA, nodes=nodes)
    LCA = Reduce(intersect, parents)[1]
  }
  else {
    LCA = taxids
  }
  return(LCA)
}


read_match <- function(query, format, match_column, match_sep) {
  match = NA
  if (format == 'biom') {
    for(i in length(query):1) {
      if (nchar(query[i]) > 0) {
        match = query[i]
        break
      }
    }
  }
  else if (format == 'tax_table') {
    if ( ! is.na(query['Species'])) {
      if (grepl("__", query['Species'])) {
        species = strsplit(query['Species'], '__', fixed = T)[[1]][2]
        genus = strsplit(query['Genus'], '__', fixed = T)[[1]][2]
        match = paste(genus, species)
        names(match) = 'Species'
      }
      else {
        match = query['Species']
      }
    }
    else if ( ! is.na(query['Genus'])) {
      match = query['Genus']
    }
    else if ( ! is.na(query['Family'])) {
      match = query['Family']
    }
    else if ( ! is.na(query['Order'])) {
      match = query['Order']
    }
  }
  else if (is.na(match_column)) {
    match = query
  }
  else {
    match = query[match_column]
  }

  if ((format != 'tax_table') && ((format != 'biom')) && (length(query[match_column]) != 1)) {
    stop("Can't read matches, check if column name is correct")
  }

  # Remove special characters found in some formats
  if (! is.na(match)) {
    n = names(match)
    if (grepl("__", match)) {
      match = strsplit(match, '__', fixed = T)[[1]][2]
      names(match) = n
    }
    if (names(match)[1] != "Species") {
      if (grepl("_", match)) {
        match = strsplit(match, '_', fixed = T)[[1]][1]
        names(match) = n
      }
    }
    if (grepl("/", match)) {
      match = strsplit(match, '/', fixed = T)[[1]][1]
      names(match) = n
    }
    match = gsub("\\[|\\]", '', match)
    match = strsplit(match, match_sep, fixed=T)
    match = unlist(match)
  }

  return(match)
}


get_match_taxid <- function(match, taxo_names) {
  match_taxid = c()
  for (m_idx in 1:length(match)) {
    m = match[m_idx]

    # Convert names to taxids if needed
    if (typeof(m) == "character") {
      # Check that they are not actually numeric
      if (! is.na(suppressWarnings(as.numeric(m)))) {
        match_taxid = c(match_taxid, as.numeric(m))
      }
      # Real taxon name to convert
      else {
        # Find taxid in list of names and add it to vector
        if (length(taxo_names[taxo_names$name==m,]$id) > 0) {
          taxid = taxo_names[taxo_names$name==m,]$id[1]
        }
        else {
          taxid = NA
        }
        match_taxid = c(match_taxid, taxid)
      }
    }
    else {
      match_taxid = c(match_taxid, m)
    }
  }
  return(match_taxid)
}
