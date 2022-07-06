


transform_for_plot <- function(input_table, sample_data=NA, only_sample=NA, PA=F) {

  if (! is.na(sample_data)) {
    if (typeof(sample_data) == "character") {
      samples = read.csv(sample_data, stringsAsFactors = FALSE)
    }
    else if (typeof(sample_data) == "data.frame") {
      samples = sample_data
    }
    else {
      print("Error: can not read sample data")
      return(1) # TODO error handling
    }
    ASVs = c()
    sample = c()
    count = c()
    estimated_genome_size = c()
    for (i in 1:length(colnames(samples))) {
      if (colnames(samples)[i] != 'X' && colnames(samples)[i] != 'ASVs') {
        samp = colnames(samples)[i]
        for (j in 1:length(rownames(samples))) {
          #if (PA) {
          #  c = 1
          #}
          #else {
          c = samples[j,i]
          #}
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
  count = 'count'
  sample = 'sample'
  density = '..density..'

  to_plot = transform_for_plot(output_table, sample_data=sample_data, only_sample=only_sample, PA=PA)

  if (is.na(sample_data)) {
    plot(ggplot(to_plot, aes_string(x=estimated_genome_size, y=density)) +
         geom_histogram( color="#e9ecef", position = 'identity', bins=bins) +
         labs(fill=""))
  }
  else {
    if (! PA) {
      plot(ggplot(to_plot, aes_string(x=estimated_genome_size, y=density, weight=count, fill=sample)) +
           geom_histogram( color="#e9ecef", alpha=.2, position = 'identity', bins=bins) +
           labs(fill=""))
    }
    else {
      plot(ggplot(to_plot, aes_string(x=estimated_genome_size, y=density, fill=sample)) +
             geom_histogram( color="#e9ecef", alpha=.2, position = 'identity', bins=bins) +
             labs(fill=""))
    }
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
  count = 'count'
  sample = 'sample'

  to_plot = transform_for_plot(output_table, sample_data=sample_data, only_sample=only_sample, PA=PA)

  if (is.na(sample_data)) {
    plot(ggplot(to_plot, aes_string(y=estimated_genome_size)) + geom_boxplot() +
      ggtitle("Boxplot of the distribution of estimated genome sizes"))
  }
  else {
    if (! PA) {
      plot(ggplot(to_plot, aes_string(x=sample, y=estimated_genome_size, weight=count, fill=sample)) +
             geom_boxplot() +
             ggtitle("Boxplot of the distribution of estimated genome sizes per sample"))
    }
    else {
      plot(ggplot(to_plot, aes_string(x=sample, y=estimated_genome_size, fill=sample)) +
             geom_boxplot() +
             ggtitle("Boxplot of the distribution of estimated genome sizes per sample"))
    }
  }
  return(to_plot)
}


# Plot genome sizes
#
# This function loads a result table from estimate_genome_size and plots the results.
#
# @param output_table Result table from estimate_genome_size
# @param sample_column The column containing sample information
# @param sample The sample to plot the genome size density for (default: all samples)
# @export
#plot_genome_size_results <- function(output_table, sample_column=NA, sample=NA) {
#  plot_genome_size_density(output_table, sample_column=sample_column)
#  plot_genome_size_boxplot(output_table, sample_column=sample_column)
#}

