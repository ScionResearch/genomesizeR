
bayesian <- function(query, models, na_models, size_db, taxonomy, names, nodes, alltax, format, output_format, match_column, match_sep, ci_threshold,
                     n_draws=NULL, predict=TRUE) {

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
  out['class'] = NA
  out['phylum'] = NA
  out['superkingdom'] = NA

  parents = allparents(LCA, taxdir=NA, nodes=nodes)
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

  if (ref_data['INFO_NODE'] == 'True') { #&& (! is.na(ref_data['MEAN_GENOME_SIZE']))) {
    estimated_size = ref_data['MEAN_GENOME_SIZE']
    out['estimated_genome_size'] = estimated_size
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
      out['genome_size_estimation_status'] = 'Could not compute confidence interval'
    }
    else {
      out['genome_size_estimation_status'] = 'OK'
    }
  }
  else {

    # Get bayes model for the query's superkingdom
    if (2 %in% parents) {
      model = models$bayes_model_bact
      out['model_used'] = 'bayesian Bacteria'
    }
    else if (2759 %in% parents) {
      model = models$bayes_model_euka
      out['model_used'] = 'bayesian Eukaryota'
    }
    else if (2157 %in% parents) {
      model = models$bayes_model_arch
      out['model_used'] = 'bayesian Archeae'
    }
    else {
      cat("\nBayesian model not found for:", fill=T)
      cat(match, fill=T)
      out['genome_size_estimation_status'] = 'Bayesian model not found'
      return(out)
    }

    # When the posterior is predicted jointly for all queries afterwards
    # (see predict_posterior_jointly), only the model to use is needed here
    if (!predict) {
      out = unlist(as.vector(out[1,]))
      return(out)
    }

    draw_ids = NULL
    if (!is.null(n_draws)) {
      draw_ids = seq_len(n_draws)
    }
    pred = brms::posterior_predict(model, newdata = out, allow_new_levels=TRUE, draw_ids=draw_ids)

    pred = as.data.frame(pred) %>%
      mutate_all(function(x){exp(x)*10^7})

    probabilities <- c(0.025, 0.975)
    pred_quant <- pred %>%
      reframe(across(everything(), ~quantile(., probabilities))) %>%
      t() %>%
      as.data.frame()
    names(pred_quant) <- paste0("Q", probabilities)

    pred_mean <- pred %>%
      summarise_all(~mean(.)) %>%
      t()

    out['estimated_genome_size'] = pred_mean
    out['confidence_interval_lower'] = pred_quant$Q0.025
    out['confidence_interval_upper'] = pred_quant$Q0.975

  }

  out = unlist(as.vector(out[1,]))

  return(out)
}


#' Predict the posterior predictive distribution of all queries jointly
#'
#' One call to posterior_predict per bayesian model, on the distinct taxonomies of the queries
#' estimated with that model. Predicting jointly (rather than query by query) ensures that queries
#' sharing a taxon unseen by the model (e.g. the same new family) share the same sampled effects
#' in each draw, which is needed to propagate uncertainty correctly when aggregating queries.
#'
#' @param output_table Result table with columns "model_used" and the taxid columns of the model
#' @param models List of bayesian models
#' @param n_draws Number of posterior draws
#' @return Numeric matrix with one row per query and one column per draw (genome size in base pairs),
#'         NA for queries not estimated with a bayesian model
#' @noRd
predict_posterior_jointly <- function(output_table, models, n_draws) {
  model_by_name = c('bayesian Bacteria'='bayes_model_bact',
                    'bayesian Eukaryota'='bayes_model_euka',
                    'bayesian Archeae'='bayes_model_arch')
  ranks = c('phylum', 'class', 'order', 'family', 'genus')

  posterior = matrix(NA_real_, nrow=nrow(output_table), ncol=n_draws,
                     dimnames=list(NULL, paste0('posterior_', seq_len(n_draws))))

  for (model_used in names(model_by_name)) {
    rows = which(!is.na(output_table$model_used) & output_table$model_used == model_used)
    if (length(rows) == 0) {
      next
    }
    newdata = output_table[rows, intersect(ranks, names(output_table)), drop=FALSE]
    for (rank in setdiff(ranks, names(newdata))) {
      newdata[rank] = NA
    }
    newdata = newdata[ranks]
    taxonomy_key = do.call(paste, c(newdata, sep='|'))
    distinct = !duplicated(taxonomy_key)
    cat(paste0("Predicting posterior distributions for ", sum(distinct), " distinct taxonomies (", model_used, ")"), fill=T)

    pred = brms::posterior_predict(models[[model_by_name[model_used]]], newdata=newdata[distinct, , drop=FALSE],
                                   allow_new_levels=TRUE, draw_ids=seq_len(n_draws))
    if (ncol(pred) != sum(distinct)) {
      stop("Joint posterior prediction returned ", ncol(pred), " columns for ", sum(distinct), " taxonomies")
    }
    posterior[rows, ] = t(pred)[match(taxonomy_key, taxonomy_key[distinct]), , drop=FALSE]
  }

  return(exp(posterior) * 10^7)
}
