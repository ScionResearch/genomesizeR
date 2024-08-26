
weighted_mean <- function(query, models, na_models, size_db, taxonomy, names, nodes, alltax, format, output_format, match_column, match_sep, ci_threshold) {

  out = query
  # if one match:
  #   if match has a size and a CI: it's that size
  #   if a size but no CI: weighted mean up to 1st ancestor with CI
  #   if no size, weighted mean up to 1st ancestor with a size and CI (should take distance and/or density into account for confidence)

  #if (output_format == 'input') {
  #  out = query
  #}
  #else {
  #  out = structure(character(0), names=character(0))
  #}

  out['estimated_genome_size'] = NA
  out['confidence_interval_lower'] = NA
  out['confidence_interval_upper'] = NA
  out['genome_size_estimation_status'] = NA
  out['genome_size_estimation_rank'] = NA
  out['genome_size_estimation_distance'] = NA
  out['LCA'] = NA
  if (format == 'tax_table' || format == 'biom') {
    out['TAXID'] = NA
    #out['SCIENTIFIC_NAME'] = NA
  }

  match = read_match(query, format, match_column, match_sep)

  if ((length(match) == 1) && (is.na(match))) {
    out['genome_size_estimation_status'] = 'Match is NA'
    return(out)
  }

  match_taxid = get_match_taxid(match, names)
  if (any(is.na(match_taxid))) {
    cat("\nNCBI taxid not found for:", fill=T)
    cat(query, fill=T)
    out['genome_size_estimation_status'] = 'NCBI taxid not found'
    return(out)
  }

  # Compute LCA
  LCA = compute_LCA(match_taxid, nodes)

  if (format == 'tax_table' || format == 'biom') {
    if (is.na(LCA)) {
      out['TAXID'] = as.character(match_taxid)
    }
    else {
      out['TAXID'] = as.character(LCA)
    }
    #out['SCIENTIFIC_NAME'] = sciname(match_taxid, names=names)
  }
  out['LCA'] = as.character(LCA)

  # Initialize parent_sizes table with the same columns as size db
  parent_sizes = size_db[0, ]
  distances = c()
  p_idx = 1

  for (match in match_taxid) {
    # Get query's taxonomy
    parents = allparents(match, taxdir=NA, nodes=nodes)
    if (is.null(parents)) {
      cat("\nParent taxids not found for:", fill=T)
      cat(match, fill=T)
      out['genome_size_estimation_status'] = 'Parent taxids not found'
      return(out)
    }

    ranks = getrank(parents, taxdir=NA, nodes=nodes)
    if (is.null(ranks)) {
      cat("\nParent taxid ranks not found for:", fill=T)
      cat(match, fill=T)
      cat(parents, fill=T)
      out['genome_size_estimation_status'] = 'Parent taxid ranks not found'
      return(out)
    }
    for (i in 1:length(parents)) {
      out[ranks[[i]]] = parents[[i]]
    }

    for (i in 1:length(parents)) {
      parent_taxid = parents[[i]]
      # Store parent if not stored already
      parent_data = size_db[parent_taxid + 1, ]
      # stop if getting higher than order
      if (! 'order' %in% ranks[i:length(ranks)]) {
        out['genome_size_estimation_status'] = 'Not enough genome size references for close taxa'
        return(out)
      }

      if (! is.na(parent_data$MEAN_GENOME_SIZE) && (! parent_data$TAXID %in% parent_sizes$TAXID)) {
        parent_sizes[p_idx,] = parent_data
        distances = c(distances, i)
        p_idx = p_idx + 1
        # Break out of loop if more than one match
        if (length(match_taxid) == 1 && ! is.na(parent_data$GENOME_COUNT) && parent_data$GENOME_COUNT > 1 ) { #& parent_data$GENOME_DATA_DENSITY > 0.15) {
          # But if that parent has poor data density, exit
          #if (! is.na(parent_data$GENOME_DATA_DENSITY) && length(parent_data$GENOME_DATA_DENSITY) > 0 && parent_data$GENOME_DATA_DENSITY < 0.1) { # TODO tune threshold, compute entropy?
          #  out['estimated_genome_size'] = NA
          #  out['estimated_genome_size_confidence_interval'] = NA
          #  out['data_density'] = parent_data$GENOME_DATA_DENSITY
          #  return(out)
          #}
          estimation_rank = parent_data$TAXONOMIC_RANK
          break
        }
      }
      # Break out of loop if several matches and the LCA has been reached
      if (length(match_taxid) > 1 && parent_taxid == LCA && nrow(parent_sizes) > 0) {
        estimation_rank = parent_data$TAXONOMIC_RANK
        break
      }
    }
  }

  # Compute the weighted mean from each ref to compute from
  estimated_size = 0
  st_err = 0
  sum_weights = 0
  for (i in 1:length(distances)) {
    estimated_size = estimated_size + ( (1.0/(distances[i]+1)) * parent_sizes$MEAN_GENOME_SIZE[i] )
    weight = 1.0/(distances[i]+1)
    sum_weights = sum_weights + weight
    if (! is.na(parent_sizes$GENOME_COUNT[i]) && parent_sizes$GENOME_COUNT[i] > 1) {
      st_err_i = parent_sizes$STANDARD_ERROR_GENOME_SIZE[i]
      st_err = st_err + (weight * st_err_i)**2
    }
  }

  estimated_size = estimated_size / sum_weights
  standard_error = sqrt(st_err)
  Z = 1.96     # 95% CI

  # Compute confidence interval
  margin_of_error = Z * standard_error

  if (margin_of_error > ci_threshold*estimated_size) {
    out['genome_size_estimation_status'] = 'Confidence interval to estimated size ratio > ci_threshold'
    return(out)
  }

  out['estimated_genome_size'] = estimated_size
  out['confidence_interval_lower'] = estimated_size - margin_of_error
  out['confidence_interval_upper'] = estimated_size + margin_of_error
  out['genome_size_estimation_status'] = 'OK'
  out['model_used'] = 'weighted_mean'

  out['genome_size_estimation_rank'] = estimation_rank
  out['genome_size_estimation_distance'] = max(distances)

  return(out)
}
