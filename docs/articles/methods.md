# Method description

## Bayesian method

The NCBI database of species with known genome sizes was split by
superkingdom (Bacteria, Archeae, Eukaryotes). A distributional Bayesian
linear hierarchical model using the `brm` function from the `brms`
package was fitted to each superkingdom dataset. The general model
structure is outlined below and corresponds exactly to the most complex
model, implemented for the Bacteria superkingdom. This general model was
simplified by dropping the class group effect in the standard deviation
model for the Eukaryote superkingdom, and dropping both the class and
phylum group effect in the standard deviation model for the Archeae
superkingdom. The latter is therefore not addressed using a
distributional model, as the response variance has no predictor.

$$\begin{matrix}
{log\left( G_{i} \right) \sim \mathcal{N}\left( \mu_{i},\sigma_{i}^{2} \right)}
\end{matrix}$$

where $G_{i}$ is the genome size of species $i$ in the units of 10 Mbp.
The model uses predictors for both mean and standard deviation. The mean
is modelled as follows:

$$\begin{matrix}
{\mu_{i} = \alpha_{0} + \alpha_{genus_{g{\lbrack i\rbrack}}} + \alpha_{family_{f{\lbrack i\rbrack}}} + \alpha_{order_{o{\lbrack i\rbrack}}} + \alpha_{class_{c{\lbrack i\rbrack}}} + \alpha_{phylum_{p{\lbrack i\rbrack}}}} \\
{\alpha_{genus_{g{\lbrack i\rbrack}}} \sim \mathcal{N}\left( 0,\sigma_{genus}^{2} \right)} \\
{\alpha_{family_{f{\lbrack i\rbrack}}} \sim \mathcal{N}\left( 0,\sigma_{family}^{2} \right)} \\
{\alpha_{order_{o{\lbrack i\rbrack}}} \sim \mathcal{N}\left( 0,\sigma_{order}^{2} \right)} \\
{\alpha_{class_{c{\lbrack i\rbrack}}} \sim \mathcal{N}\left( 0,\sigma_{class}^{2} \right)} \\
{\alpha_{phylum_{p{\lbrack i\rbrack}}} \sim \mathcal{N}\left( 0,\sigma_{phylum}^{2} \right)} \\

\end{matrix}$$

The following prior distributions are used:

$$\begin{matrix}
{\alpha_{0} \sim \mathcal{N}(0,5)} \\
{\left( \sigma_{genus},\sigma_{family},\sigma_{order},\sigma_{class},\sigma_{phylum},s_{class},s_{phylum} \right) \sim \mathcal{N}^{+}(0,1)} \\

\end{matrix}$$

Differences in genome size variability was observed among taxa,
therefore the model also adds predictors to the standard deviation of
the response. The standard deviation is modelled as follows:

$$\begin{matrix}
{log\left( \sigma_{i} \right) = \lambda_{0} + \lambda_{class_{c{\lbrack i\rbrack}}} + \lambda_{phylum_{p{\lbrack i\rbrack}}}} \\
{\lambda_{class_{c{\lbrack i\rbrack}}} \sim \mathcal{N}\left( 0,s_{class}^{2} \right)} \\
{\lambda_{phylum_{p{\lbrack i\rbrack}}} \sim \mathcal{N}\left( 0,s_{phylum}^{2} \right)} \\

\end{matrix}$$

with priors

$$\begin{matrix}
{\lambda_{0} \sim \mathcal{N}(0,1)} \\
{\left( s_{class},s_{phylum} \right) \sim \mathcal{N}^{+}(0,1)} \\

\end{matrix}$$

$\mathcal{N}^{+}$ is the normal distribution truncated to positive
values.
$g\lbrack i\rbrack$,$f\lbrack i\rbrack$,$o\lbrack i\rbrack$,$c\lbrack i\rbrack$
and $p\lbrack i\rbrack$ are respectively the index for the genus,
family, order, class, and phylum of entry $i$ in the species-level
database. Note that taxonomic groups are naturally nested within each
other and the indices are designed to be unique to the particular
taxonomic group it represents.

The estimation process uses Stan’s Hamiltonian Monte Carlo algorithm
with the U-turn sampler.

Posterior predictions are obtained using the `predict` function from the
`brms` package, and 95% credible intervals are obtained using 2.5% and
97.5% quantiles from the posterior distribution.

Queries corresponding to identified species with an available genome
size estimate in the NCBI database get allocated the genome size value
of the database (averaged at the species level) and 95% confidence
intervals are calculated based on the standard error of the mean of all
genome sizes available for that species in the processed NCBI database.

#### Sample-level estimation

The full posterior predictive distribution of each query can also be
kept (`return_posterior` option) and aggregated per sample with
`estimate_genome_size_per_sample`. Posterior predictions are drawn
jointly for all queries, so that the correlation between the predictions
of related queries is preserved when they are aggregated. For each
sample and each posterior draw, the genome sizes of the queries present
in the sample are averaged, weighted by $log\left( 1 + n_{i} \right)$
where $n_{i}$ is the abundance of query $i$ in the sample (weights based
on relative abundance or presence are also available). The mean of the
resulting sample-level posterior distribution is the estimated genome
size of the sample, and its 2.5% and 97.5% quantiles give the 95%
credible interval, which therefore accounts for the estimation
uncertainty of every query. Queries with a species-level genome size
from the database keep this value in every draw; queries without an
estimate are ignored, and the proportion of the sample abundance they
represent is reported.

## Frequentist method

A frequentist linear mixed-effects model using the `lmer` function from
the `lme4` package was fitted to the NCBI database of species with known
genome sizes. The model is as follows:

$$\begin{matrix}
{log\left( G_{i} \right) = \alpha_{0} + \alpha_{genus_{g{\lbrack i\rbrack}}} + \alpha_{family_{f{\lbrack i\rbrack}}} + e_{i}} \\

\end{matrix}$$ where $\alpha_{0}$ is the overall mean,
$\alpha_{genus_{g{\lbrack i\rbrack}}}$ and
$\alpha_{family_{f{\lbrack i\rbrack}}}$ are random effect of genus and
family for genus $g\lbrack i\rbrack$ and family $f\lbrack i\rbrack$ and
$e_{i}$ is the residual error of observation $i$.

The estimation process using the restricted maximum likelihood method
(REML). A prediction interval is computed using the `predictInterval`
function from the `merTools` package. As higher nested levels (order,
class) are not taken into account in the model, predictions are not
produced for queries above the family level.

## Weighted mean method

The weighted mean method computes the genome size of a query by
averaging the known genome sizes of surrounding taxa in the taxonomic
tree, with a weighted system where further neighbours have less weight
in the computed mean. The identification of related taxa is limited to
levels below and including order.

For queries relating to well-characterised species where many genetic
studies have been performed, such as model organisms, this might lead to
more precise predictions than the two other methods. This method can
also perform better than the others if your queries consist of lists of
taxa (for example, an output of *blastn* where several matches can be
obtained for each query). Otherwise, we suggest using one of the other
methods, as the confidence intervals calculated are less reliable for
the weighted mean method.

#### Pseudocode describing the weighted mean method computation:

    # 1. Build table of parent information

    for each match (a query can be a list of matches e.g. blast result):
      get the taxid associated
      get the list of parents and their ranks
      for each parent :
        if higher than order rank: 
          STOP and return 'Not enough genome size references for close taxa'
        if genome size associated with that taxon in the reference database:
          store it in a 'parent information table'
          store the distance as well (number of nodes between query and this taxon)
          if we have more than one genome that contributed to the data we have in the 'parent information table':
            STOP the iteration through parents here
        if the query is a list of matches and we have reached their common ancestor:
          STOP the iteration through parents here

    # 2. Compute weighted mean from the table of parent information

    estimated_size = 0
    estimated_standard_error = 0
    sum_weights = 0
    for each line of parent information:
      weight = 1.0/(distance+1)
      sum_weights = sum_weights + weight
      estimated_size = estimated_size + weight * genome_size
      estimated_standard_error = estimated_standard_error + (weight * precomputed_standard_error)**2
    estimated_size = estimated_size / sum_of_all_weights
    standard_error = square_root(estimated_standard_error)
    Z = 1.96 # 95% CI
    margin_of_error = Z * standard_error
    confidence_interval_lower = estimated_size - margin_of_error
    confidence_interval_upper = estimated_size + margin_of_error
