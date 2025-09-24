---
title: 'genomesizeR: An R package for genome size prediction'
tags:
- R
- molecular ecology
authors:
  - name: Celine Mercier
    email: celine.mercier@scionresearch.com
    correspondance: "yes"
    orcid: "0000-0002-4782-1530"
    affiliation: '1'
  - name: Joane Elleouet
    orcid: "0000-0002-9597-3360"
    affiliation: '1'
  - name: Sean Husheer
    orcid: "0000-0002-0422-2386"
    affiliation: "1"
  - name: Loretta Garrett
    orcid: "0000-0001-7687-6795"
    affiliation: '1'
  - name: Steve A Wakelin
    orcid: "0000-0002-1167-8699"
    affiliation: '1'
affiliations:
  - name: Scion Group, Bioeconomy Science Institute, Titokorangi Drive, Private Bag 3020, Rotorua 3046, New Zealand
    index: 1
date: 04 August 2025
bibliography: paper.bib

---

# Summary

The genome size of organisms present in an environment can provide many insights into evolutionary and ecological processes at play in that environment. The genomic revolution has enabled a rapid expansion of our knowledge of genomes in many living organisms, and most of that knowledge is classified and readily available in the databases of the National Center for Biotechnology Information (NCBI). The `genomesizeR` tool leverages the wealth of taxonomic and genomic information present in NCBI databases to infer the genome size of Archeae, Bacteria, or Eukaryote organisms identified at any taxonomic level.

This R package provides three statistical methods for genome size prediction of a given taxon, or group of taxa. A straightforward 'weighted mean' method identifies the closest taxa with available genome size information in the taxonomic tree, and averages their genome sizes using weights based on taxonomic distance. A frequentist random effect model uses nested genus and family information to output genome size estimates. Finally, a third option provides predictions from a distributional Bayesian multilevel model which uses taxonomic information from genus all the way to superkingdom, therefore providing estimates and uncertainty bounds even for under-represented taxa.

`genomesizeR` retrieves the taxonomic classification of input queries, estimates the genome size of each query, and provides 95% confidence intervals for each estimate. Some plotting functions are also provided to visualise the results.

# Statement of need

The size of genomes and its evolution can provide important insights into evolutionary and ecological processes influencing both species and the environments they inhabit. The shedding of unnecessary genetic elements and their associated biosynthetic pathways, for example, is a common phenomenon observed in organisms with a high degree of host symbiosis [@moran2002microbial; @brader2014metabolic; @vandenkoornhuyse2007active]. Among many others, these findings demonstrate the opportunities associated with including genome size as a key trait in studies on communities to provide insights into ecological and evolutionary processes.

However, characterizing genome size remains challenging. The exponentially growing genome databases are an inexpensive resource unlocking a myriad of research opportunities, but genome size estimates for many taxa found in environmental samples are missing from public databases, or fully unknown. The evolutionary rule that phylogenetically related organisms share genetic similarities can be exploited, and genome size for taxa with unknown genome size can be statistically inferred from related taxa with known genome size, using taxonomy as a proxy for phylogeny. Another challenge is the precision of identification: some taxa can only be identified at high taxonomic levels. Statistical methods can also be used to infer their genome size range from databases. To our knowledge, there is no convenient and fast way to obtain genome size estimates with uncertainty bounds for any organism.

Using the increased prevalence of whole-genome information for all organisms, we have therefore developed `genomesizeR`, allowing the inference of genome size of many queries at once, based on taxonomic information and available genome data from the NCBI.

# Methods

## NCBI database filtering and processing

The reference database is built by querying all genome metadata information from the curated NCBI RefSeq database [@OLeary2016-kw]. This raw database is then filtered and prepared to include more pre-computed information to be used by the package.

## Bayesian method

The reference database of genome sizes was split by superkingdom (Bacteria, Archeae, Eukaryotes). A distributional Bayesian linear hierarchical model using the `brm` function from the `brms` package [@burkner2021brms] was fitted to each superkingdom dataset. The general model structure is outlined below and corresponds exactly to the most complex model, implemented for the Bacteria superkingdom. This general model was simplified by dropping the class group effect in the standard deviation model for the Eukaryote superkingdom, and dropping both the class and phylum group effect in the standard deviation model for the Archeae superkingdom. The latter is therefore not addressed using a distributional model, as the response variance has no predictor. The model is as follows:

\begin{gather*}
log(G_i) \sim \mathcal{N}(\mu_i, \sigma_{i}^2)
\end{gather*}

where $G_i$ is the genome size of species $i$ in the units of 10 Mbp. The model uses taxonomic levels as predictors, and is described in more detail in the package vignettes.

The estimation process uses Stan's Hamiltonian Monte Carlo algorithm with the U-turn sampler.

Posterior predictions are obtained using the predict function from the `brms` package, and 95% credible intervals are obtained using 2.5% and 97.5% quantiles from the posterior distribution.

## Frequentist method

A frequentist linear mixed-effects model (LMM) using the `lmer` function from the `lme4` package [@bates2015lme4] was fitted to the NCBI database of species with known genome sizes. The model is as follows:

\begin{gather*}
log(G_i) =  \alpha_0 + \alpha_{genus_{g[i]}} +  \alpha_{family_{f[i]}} + e_i \\
\end{gather*}
where $\alpha_0$ is the overall mean, $\alpha_{genus_{g[i]}}$ and $\alpha_{family_{f[i]}}$ are random effect of genus and family for genus $g[i]$ and family $f[i]$ and $e_i$ is the residual error of observation $i$.

## Weighted mean method

The weighted mean method computes the genome size of a query by averaging the known genome sizes of surrounding taxa in the taxonomic tree, with a weighted system where further neighbours have less weight in the computed mean.

# Method validation and comparison

The strengths and limitations of each method are outlined in \autoref{table:method_comp}. The weighted mean method is less reliable but can be used on queries with several potential taxonomic matches. The Bayesian method is the most reliable method especially for quantifying uncertainty around estimated means, and obtaining estimates for taxa that are not well represented at low ranks in the NCBI database.

| | CI estimation	| Model information	| Behaviour with well-studied organisms	| Query is a list of several taxa	| Minimum number of references needed for estimation |
| -- | -- | -- | -- | -- | -- |
| Bayesian | very reliable | any rank | + | + | 1 |
| LMM | mostly reliable | up to family level | + | + | 1 |
| Weighted mean | unreliable | up to order level | ++ | ++ | 2 |

: Comparison of method behaviour and applicability \label{table:method_comp}

# Availability

- Project name: genomesizeR
- Project home page: https://github.com/ScionResearch/genomesizeR 
- Operating system(s): Platform independent
- Programming language: R
- License: GNU General Public License

# Acknowledgements

The authors declare that they have no conflict of interest. Funding for this research came from the Tree-Root-Microbiome programme, which is funded by MBIE’s Endeavour Fund and in part by the New Zealand Forest Growers Levy Trust (C04X2002). We make no warranties regarding the accuracy or integrity of the Data. We accept no liability for any direct, indirect, special, consequential or other losses or damages of whatsoever kind arising out of access to, or the use of the Data. We are in no way to be held responsible for the use that you put the Data to. You rely on the Data entirely at your own risk.

# References
