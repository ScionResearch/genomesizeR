

compute_LCA <- function(taxids, nodes) {
  if (length(taxids) > 1) {
    parents = lapply(taxids, allparents, taxdir=NA, nodes=nodes)
    LCA = Reduce(intersect, parents)[1]
  }
  else {
    LCA = NA
  }
  return(LCA)
}


read_match <- function(query, format, match_column, match_sep) {
  match = NA
  if (format == 'dada2') {
    if ( ! is.na(query['Species'])) {
      if (grepl("__", query['Species'])) {
        species = strsplit(query['Species'], '__', fixed = T)[[1]][2]
        genus = strsplit(query['Genus'], '__', fixed = T)[[1]][2]
        match = paste(genus, species)
        names(match) = 'Species'
        #print(('24')
        #print((match)
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

  # Remove special characters found in some formats
  if (! is.na(match)) {
    n = names(match)
    #print(("47")
    #print((match)
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
    #print(("64")
    #print((match)
  }

  return(match)
}


get_match_taxid <- function(match, taxo_names) {
  match_taxid = c()
  for (m_idx in 1:length(match)) {
    m = match[m_idx]
    #print('m')
    #print(m)
    # Convert names to taxids if needed
    if (typeof(m) == "character") {
      # Check that they are not actually numeric
      if (! is.na(suppressWarnings(as.numeric(m)))) {
        match_taxid = c(match_taxid, as.numeric(m))
      }
      # Real taxon name to convert
      else {
        ##print(('65')
        ##print((m)
        #print('79')
        #print(m)
        # Find taxid in list of names and add it to vector
        if (length(taxo_names[taxo_names$name==m,]$id) > 0) {
          taxid = taxo_names[taxo_names$name==m,]$id[1]
          #print('OUI, 84')
          #print(taxid)
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
  ##print((match_taxid)
  return(match_taxid)
}
