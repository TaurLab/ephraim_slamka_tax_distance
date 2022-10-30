formulas
================

Various distance metrics have been used to assess beta-diversity across
sample populations in microbial ecology.

Some are shown below, where $D_{jk}$ denotes the calculated distance
between samples $j$ and $k$, where $x_{ij}$ and $x_{ik}$ refer to the
quantity of bacterial species/strain $i$ in pairwise samples $j$ and
$k$.

| Distance Metric | Definition                                                                                                                                                                                                  |
|:----------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Euclidean       | $\textrm{Deuc}_{jk} = \sqrt{\sum_i (x_{ij}-x_{ik})^2}$                                                                                                                                                      |
| Manhattan       | $\textrm{Dman}_{jk}=\sum_i \mid x_{ij}-x_{ik} \mid$                                                                                                                                                         |
| Morisita        | $\textrm{Dmor}_{jk} = 1 - \dfrac{2 \sum_i x_{ij} x_{ik}}{(\lambda_j +  \lambda_k) \sum_i x_{ij} \sum_i  x_{ik}}$, where $\lambda_j = \dfrac{\sum_i x_{ij} (x_{ij} - 1)}{\sum_i x_{ij} \sum_i (x_{ij} - 1)}$ |
| Morisita-Horn   | $\textrm{Dhorn}_{jk} = 1 - \dfrac{2 \sum_i x_{ij} x_{ik}}{(\lambda_j +  \lambda_k) \sum_i x_{ij} \sum_i  x_{ik}}$, where $\lambda_j = \sum_i {x_{ij}}^2/(\sum_i x_{ij})^2$                                  |

Other metrics calculate beta diversity by examining genetic related-ness
among microbial members. Unifrac calculates beta-diversity distance
between two populations using a phylogenetic tree, where $l_n$
represents the length between node $n$ and its parent, and $X_{nj}$ and
$X_{nk}$ are indicators (0 or 1) as descendants of node $n$ are absent
or present in samples $j$ and $k$, respectively.

| Distance Metric    | Definition                                                                                         |
|:-------------------|:---------------------------------------------------------------------------------------------------|
| Unweighted Unifrac | $\textrm{Dunifrac}_{jk} = \dfrac{\sum_n l_n \mid A_n - B_n \mid}{\sum_n l_n \max(X_{nj}, X_{nk})}$ |

There are limitations with all of the above metrics.

We used Horn index to calculate $\textrm{Dhorn}_{jk}$, representing the
Horn calculation across all taxonomic levels

Let $H_{ijk}$ represent the Horn index calculated at the taxonomic level
$i$ for samples $j$ and $k$.

$$
\textrm{Dtaxhorn}_{jk} = \dfrac{\sum_{i=2}^{N} (N-i-1) \textrm{Dhorn}_{ijk}}{\sum_{i=2}^{N} (N-i-1)}
$$

That is, the weighted average of Horn distances calculated for each
taxonomic level. Lower taxonomic levels are weighted more than higher
taxonomic levels. We omit the first level, Superkingdom, since its value
is ‘Bacteria’ for all microbial members in the sample data. Therefore
the weights across levels are as follows: Superkingdom: 0, Phylum: 7,
Class: 6, Order: 5, Family: 4, Genus: 3, Species: 2, ASV: 1.