
compute_confidence_interval_lmm <- function(query, model, n_cores) {

  if (n_cores > 1) {
    cluster <- parallel::makeCluster(
      n_cores,
      type = "PSOCK"
    )
    doParallel::registerDoParallel(cl = cluster)
    ci =  predictInterval(model, newdata = query,  which = "full",  n.sims = 1000, include.resid.var = FALSE, level=0.95,  stat="mean", .parallel=T)
    parallel::stopCluster(cl = cluster)
  }
  else {
    ci =  predictInterval(model, newdata = query,  which = "full",  n.sims = 1000, include.resid.var = FALSE, level=0.95,  stat="mean")
  }
  return(ci)
}


lmm <- function(query, models, na_models, size_db, taxonomy, names, nodes, alltax, format, output_format, match_column, match_sep, ci_threshold) {

  genusfamily_model = models$genusfamily_model

  out = query
  out['estimated_genome_size'] = NA
  out['confidence_interval_lower'] = NA
  out['confidence_interval_upper'] = NA
  out['genome_size_estimation_status'] = NA
  out['model_used'] = NA
  out['LCA'] = NA
  if (format == 'tax_table' || format == 'biom') {
    out['TAXID'] = NA
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
  }
  out['LCA'] = as.character(LCA)

  ref_data = size_db[LCA + 1, ]

  out['species'] = NA
  out['genus'] = NA
  out['family'] = NA
  out['order'] = NA
  parents = allparents(LCA, taxdir=NA, nodes=nodes)
  #parents = ncbitax::get.parents(LCA, alltax)
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

  out = as.data.frame(t(as.data.frame(out)))

  if ((ref_data['INFO_NODE'] == 'True')) { #&& (! is.na(ref_data['MEAN_GENOME_SIZE']))) {
    estimated_size = ref_data['MEAN_GENOME_SIZE']
    out['model_used'] = 'reference_mean'
    # Compute confidence interval
    standard_error = sqrt(ref_data['STANDARD_ERROR_GENOME_SIZE'])
    Z = 1.96     # 95% CI
    margin_of_error = Z * as.numeric(standard_error)
    out['confidence_interval_lower'] = estimated_size - margin_of_error
    out['confidence_interval_upper'] = estimated_size + margin_of_error
    if ((!is.na(margin_of_error)) && (margin_of_error > ci_threshold*estimated_size)) {
      out['genome_size_estimation_status'] = 'Confidence interval to estimated size ratio > ci_threshold'
    }
    else if (is.na(margin_of_error)) {
      out['genome_size_estimation_status'] = 'OK but no confidence interval'
    }
    else {
      out['genome_size_estimation_status'] = 'OK'
    }
  }
  else if (! is.na(out$family)) {
    estimated_size = exp(predict(genusfamily_model, out, type="response", allow.new.levels=TRUE))
    out['model_used'] = 'lmm|family/genus'
  }
  else {
    out['genome_size_estimation_status'] = 'No reference and query too high in taxonomic tree to fit in model'
    return(out)
  }

  out['estimated_genome_size'] = estimated_size

  out = unlist(as.vector(out[1,]))

  return(out)
}
