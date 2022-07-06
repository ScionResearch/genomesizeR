

weighted_mean <- function(query, size_db, taxonomy, names, nodes, format, match_column, match_sep) {

  # if one match:
  #   if match has a size and a CI: it's that size
  #   if a size but no CI: weighted mean up to 1st ancestor with CI
  #   if no size, weighted mean up to 1st ancestor with a size and CI (should take distance and/or density into account for confidence)
  # if several matches: ideally that's better (test that it is). Might need to keep pre-computed db distances
  #   weighted mean up to LCA (or to LCA with data...)
  # TODO? add a special confidence index: CI*density or something

  # TODO fix multi match behaviour

  match = read_match(query, format, match_column, match_sep)

  if ((length(match) == 1) && (is.na(match))) {
    query['estimated_genome_size'] = NA
    query['estimated_genome_size_confidence_interval'] = NA
    return(query)
  }

  #print("matches:")
  #print(match)

  match_taxid = get_match_taxid(match, names)
  if (any(is.na(match_taxid))) {
    cat("\nNCBI taxid not found for:", fill=T)
    cat(query, fill=T)
    query['estimated_genome_size'] = NA
    query['estimated_genome_size_confidence_interval'] = NA
    return(query)
  }

  #print("match taxids:")
  #print(match_taxid)

  # Compute LCA
  LCA = compute_LCA(match_taxid, nodes)

  #print("LCA:")
  #print(LCA)

  # Initialize parent_sizes table with the same columns as size db
  parent_sizes = size_db[0, ]
  distances = c()
  p_idx = 1

  for (match in match_taxid) {
    # Get query's taxonomy
    parents = allparents(match_taxid, taxdir=NA, nodes=nodes)

    #print("parents:")
    #print(parents)

    for (i in 1:length(parents)) {
      parent_taxid = parents[[i]]
      #print("entering iteration for parent:")
      #print(parent_taxid)
      # Break out of loop if several matches and the LCA has been reached
      if (! is.na(LCA) && parent_taxid == LCA) {
        #print("breaking out of loop, LCA reached")
        break
      }

      # Store parent if not stored already
      parent_data = size_db[parent_taxid + 1, ]
      if (! is.na(parent_data$MEAN_GENOME_SIZE) && (! parent_data$TAXID %in% parent_sizes$TAXID)) {
        parent_sizes[p_idx,] = parent_data
        distances = c(distances, i)
        p_idx = p_idx + 1
        # Break out of loop if one match and enough data
        if (is.na(LCA) && ! is.na(parent_data$GENOME_COUNT) && parent_data$GENOME_COUNT > 1 ) { #&& parent_data$GENOME_DATA_DENSITY > 0.15) {
          # But if that parent has poor data density, exit
          if (! is.na(parent_data$GENOME_DATA_DENSITY) && parent_data$GENOME_DATA_DENSITY < 0.2) { # TODO tune threshold, compute entropy?
            #print("low genome density")
            #print(parent_data$GENOME_DATA_DENSITY)
            #print(parent_data$TAXONOMIC_RANK)
            query['estimated_genome_size'] = NA
            query['estimated_genome_size_confidence_interval'] = NA
            return(query)
          }
          else {
            break
          }
        }
      }
    }
  }

  #print("parent_sizes")
  #print(parent_sizes)
  #print("distances")
  #print(distances)

  # Compute the weighted mean from each ref to compute from
  estimated_size = 0
  st_err = 0
  sum_weights = 0
  for (i in 1:length(distances)) {
    ##print("computing - i")
    ##print(i)
    ##print("computing - estimated_size")
    ##print(estimated_size)
    ##print("computing - distances[i]+1")
    ##print(distances[i]+1)
    ##print("computing - parent_sizes$MEAN_GENOME_SIZE[i]")
    ##print(parent_sizes$MEAN_GENOME_SIZE[i])
    estimated_size = estimated_size + ( (1.0/(distances[i]+1)) * parent_sizes$MEAN_GENOME_SIZE[i] )
    weight = 1.0/(distances[i]+1)
    sum_weights = sum_weights + weight
    #print(query)
    #print(parent_sizes$GENOME_COUNT[i])
    if (! is.na(parent_sizes$GENOME_COUNT[i]) && parent_sizes$GENOME_COUNT[i] > 1) {
      st_err_i = parent_sizes$STANDARD_ERROR_GENOME_SIZE[i]
      st_err = st_err + (weight * st_err_i)**2
    }
  }
  ##print(estimated_size)
  ##print(sum_weights)

  estimated_size = estimated_size / sum_weights
  standard_error = sqrt(st_err)

  #print("estimated_size")
  #print(estimated_size)

  Z = 1.96     # 95% CI

  # Compute confidence interval
  confidence_interval = Z * standard_error

  if (confidence_interval > 0.2*estimated_size) {
    query['estimated_genome_size'] = NA
    query['estimated_genome_size_confidence_interval'] = NA
    return(query)
  }

  query['estimated_genome_size'] = estimated_size
  query['estimated_genome_size_confidence_interval'] = confidence_interval
  return(query)
}
