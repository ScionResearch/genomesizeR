
compute_confidence_interval_hierarchical <- function(query, model, n_cores) {

  my.cluster <- parallel::makeCluster(
    n_cores,
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)
  ci =  predictInterval(model, newdata = query,  which = "full",  n.sims = 100, include.resid.var = FALSE, level=0.95,  stat="mean", .parallel=T)
  parallel::stopCluster(cl = my.cluster)

  return(ci)
}


hierarchical <- function(query, models, na_models, size_db, taxonomy, names, nodes, alltax, format, output_format, match_column, match_sep, ci_threshold) {

  genusfamily_model = models$genusfamily_model
  genusorder_model = models$genusorder_model
  familyorder_model = models$familyorder_model

  out = query
  out['estimated_genome_size'] = NA
  out['estimated_genome_size_confidence_interval'] = NA
  out['genome_size_estimation_status'] = NA
  out['model_used'] = NA
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
  #ranks = ncbitax::getRank(parents, alltax)
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
    confidence_interval = Z * as.numeric(standard_error)
    out['estimated_genome_size_confidence_interval'] = confidence_interval
    if ((!is.na(confidence_interval)) && (confidence_interval > ci_threshold*estimated_size)) {
      out['genome_size_estimation_status'] = 'Confidence interval to estimated size ratio > ci_threshold'
    }
    else if (is.na(confidence_interval)) {
      out['genome_size_estimation_status'] = 'OK but no confidence interval'
    }
    else {
      out['genome_size_estimation_status'] = 'OK'
    }
  }
  else if (! na_models[1] && ! is.na(out$genus) && ! is.na(out$family)) {
    estimated_size = exp(predict(genusfamily_model, out, type="response", allow.new.levels=TRUE))
    out['model_used'] = 'hierarchical_1|genus/family'
  }
  else if (! na_models[2] && ! is.na(out$genus) && ! is.na(out$order)) {
    estimated_size = exp(predict(genusorder_model, out, type="response", allow.new.levels=TRUE))
    out['model_used'] = 'hierarchical_1|genus/order'
  }
  else if (! na_models[3] && ! is.na(out$family) && ! is.na(out$order)) {
    estimated_size = exp(predict(familyorder_model, out, type="response", allow.new.levels=TRUE))
    out['model_used'] = 'hierarchical_1|family/order'
  }
  else {
    out['genome_size_estimation_status'] = 'No reference and query too high in taxonomic tree to fit in model'
    return(out)
  }

  out['estimated_genome_size'] = estimated_size

  out = unlist(as.vector(out[1,]))

  return(out)
}
