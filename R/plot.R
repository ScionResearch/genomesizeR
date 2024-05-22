
transform_for_plot <- function(input_table, sample_data=NA, only_sample=NA, PA=F) {

  if (! is.na(sample_data)) {
    if (typeof(sample_data) == "character") {
      samples = read.csv(sample_data, stringsAsFactors = FALSE)
    }
    else if (typeof(sample_data) == "data.frame") {
      samples = sample_data
    }
    else {
      stop("Error: can not read sample data")
    }
    ASVs = c()
    sample = c()
    count = c()
    estimated_genome_size = c()
    for (i in 1:length(colnames(samples))) {
      if (colnames(samples)[i] != 'X' && colnames(samples)[i] != 'ASVs') {
        samp = colnames(samples)[i]
        for (j in 1:length(rownames(samples))) {
          if (PA) {
            c = 1
          }
          else {
            c = samples[j,i]
          }
          gs = input_table[input_table$ASVs == samples$ASVs[j], ]$estimated_genome_size
          if (c>0 && !is.na(gs) && (is.na(only_sample) || (!is.na(only_sample) && only_sample==samp))) {
            ASVs = c(ASVs, samples$ASVs[j])
            sample = c(sample, samp)
            count = c(count, c)
            estimated_genome_size= c(estimated_genome_size, gs)
          }
        }
      }
    }
    df <- data.frame(ASVs, sample, count, estimated_genome_size)
  }
  else {
    df = input_table
    if ('SAMPLE' %in% colnames(output_table)) {
      sample = 'SAMPLE'
    }
    else {
      sample = 'sample'
    }
    if (!is.na(only_sample)) {
      df = df[df[,sample] == only_sample,]
    }
  }

  return(df)
}


#' Plot genome size histograms
#'
#' This function loads a result table from estimate_genome_size and plots the histograms of estimated genome sizes.
#'
#' @param output_table Result table from estimate_genome_size
#' @param sample_data The file or dataframe containing sample information (sample name and count)
#' @param only_sample The sample to plot the genome size histogram for (default: all samples)
#' @param bins Histogram bin parameter
#' @param PA Use presence/absence data instead of abundance data (all counts are set to 1)
#' @import ggplot2
#' @import utils
#' @export
plot_genome_size_histogram <- function(output_table, sample_data=NA, only_sample=NA, bins=50, PA=F) {

  estimated_genome_size = 'estimated_genome_size'
  if ('COUNT' %in% colnames(output_table)) {
    count = 'COUNT'
  }
  else {
   count = 'count'
  }
  if ('SAMPLE' %in% colnames(output_table)) {
    sample = 'SAMPLE'
  }
  else {
    sample = 'sample'
  }
  density = 'after_stat(density)'

  to_plot = transform_for_plot(output_table, sample_data=sample_data, only_sample=only_sample, PA=PA)
  if ((! PA) && (sample %in% colnames(output_table)) && (count %in% colnames(output_table))) {
    plot(ggplot(to_plot, aes_string(x=estimated_genome_size, y=density, weight=count, fill=sample)) +
         geom_histogram( color="#e9ecef", alpha=.2, position = 'identity', bins=bins) +
         labs(fill=""))
  }
  else if (sample %in% colnames(output_table)) {
    plot(ggplot(to_plot, aes_string(x=estimated_genome_size, y=density, fill=sample)) +
           geom_histogram( color="#e9ecef", alpha=.2, position = 'identity', bins=bins) +
           labs(fill=""))
  }
  else {
    plot(ggplot(to_plot, aes_string(x=estimated_genome_size, y=density)) +
           geom_histogram( color="#e9ecef", position = 'identity', bins=bins) +
           labs(fill=""))
  }

  return(to_plot)
}


#' Plot genome size boxplots
#'
#' This function loads a result table from estimate_genome_size and plots the boxplots of estimated genome sizes.
#'
#' @param output_table Result table from estimate_genome_size
#' @param sample_data The file or dataframe containing sample information (sample name and count)
#' @param only_sample The sample to plot the genome size boxplot for (default: all samples)
#' @param PA Use presence/absence data instead of abundance data (all counts are set to 1)
#' @import ggplot2
#' @import utils
#' @export
plot_genome_size_boxplot <- function(output_table, sample_data=NA, only_sample=NA, PA=F) {

  estimated_genome_size = 'estimated_genome_size'
  if ('COUNT' %in% colnames(output_table)) {
    count = 'COUNT'
  }
  else {
    count = 'count'
  }
  if ('SAMPLE' %in% colnames(output_table)) {
    sample = 'SAMPLE'
  }
  else {
    sample = 'sample'
  }

  to_plot = transform_for_plot(output_table, sample_data=sample_data, only_sample=only_sample, PA=PA)
  if ((! PA) && (sample %in% colnames(output_table)) && (count %in% colnames(output_table))) {
    plot(ggplot(to_plot, aes_string(x=sample, y=estimated_genome_size, weight=count, fill=sample)) +
           geom_boxplot() +
           ggtitle("Boxplot of the distribution of estimated genome sizes per sample"))
  }
  else if (sample %in% colnames(output_table)) {
    plot(ggplot(to_plot, aes_string(x=sample, y=estimated_genome_size, fill=sample)) +
           geom_boxplot() +
           ggtitle("Boxplot of the distribution of estimated genome sizes per sample"))
  }
  else {
    plot(ggplot(to_plot, aes_string(y=estimated_genome_size)) + geom_boxplot() +
           ggtitle("Boxplot of the distribution of estimated genome sizes"))
  }
  return(to_plot)
}


#' Plot genome size tree
#'
#' This function loads a result table from estimate_genome_size and plots the taxonomic tree with colour-coded estimated genome sizes.
#'
#' @param output_table Result table from estimate_genome_size
#' @import ggplot2
#' @import ncbitax
#' @import dplyr
#' @import ggtree
#' @export
plot_genome_size_tree <- function(output_table, taxonomy=NA) {

  taxid_column="LCA"
  if (is.na(taxonomy)) {
    taxonomy = get_taxonomy(taxonomy)
  }
  tax = parseNCBITaxonomy(taxonomy)
  to_plot = output_table[!is.na(as.numeric(output_table[,taxid_column])) & !is.na(output_table$estimated_genome_size), ]
  to_plot$label = to_plot[,taxid_column]
  to_plot$label = as.character(to_plot$label)
  to_plot = to_plot[,c('label', 'estimated_genome_size')]
  to_plot = unique(to_plot)
  taxids = to_plot$label

  # Look for parent nodes stuck with their children:
  # List of lists of parents and recombine with the ones needed at each iteration
  parents = list()
  i = 1
  for (taxid in taxids) {
    taxid_parents = as.vector(ncbitax::get.parents(taxid, tax))
    parents[[i]] = taxid_parents
    #cat('i:', i, 'taxid:', taxid, '\n')
    i = i+1
  }

  keep_taxid_idx = 1:length(taxids)
  for (i in 1:length(taxids)) {
    taxid = taxids[i]
    taxids_to_use_idx = keep_taxid_idx[keep_taxid_idx != i]
    parents_to_use = parents[taxids_to_use_idx]
    parents_to_use = unlist(parents_to_use)
    parents_to_use = unique(parents_to_use)
    if (taxid %in% parents_to_use) {
      keep_taxid_idx = keep_taxid_idx[keep_taxid_idx != i]
    }
  }

  to_plot = to_plot[keep_taxid_idx,]
  to_plot$sci_name_for_tree = ncbitax::getName(c(to_plot$label), tax)
  tn = tax2newick(taxids, tax)
  tree = full_join(tn$tree, to_plot, by='label')

  plot(ggtree(tree, branch.length="none", layout='circular', aes(color=estimated_genome_size)) +
    scale_color_gradient2(low = "red", midpoint = median(to_plot$estimated_genome_size), mid = "blue", high = "green", space="Lab") +
    geom_tiplab(aes(label=sci_name_for_tree)) +
    geom_nodelab(aes(label=sci_name_for_tree), geom='label'))

  return(to_plot)
}
