# Broken_detailed_balance
Code for "Broken detailed balance and entropy production in the human brain" by Christopher W. Lynn, Eli J. Cornblath, Lia Papadopoulos, Maxwell A. Bertolero, and Danielle S. Bassett.

Data can be found here: https://www.dropbox.com/sh/p0tbnom0oum8f3d/AAC4eJsGGGAxLRdocq0_4gWaa?dl=0

The script "HCP_cluster_script.m" demonstrates how we performed the main analysis in the paper (as illustrated in Fig. 4), first using k-means clustering to identify coarse-grained neural states and then estimating entropy production. This script uses "kmeans_bisection.m" to perform the hierarchical clustering, "bootstrap_transitions.m" to perform the bootstrap sampling, and "entProd_transitions.m" to estimate the entropy production.

In order to creat the flux maps illustrated in Fig. 1, first run "fluxMap_bootstrap.m". Then use either "plot_fluxMap_rest.m" (for the rest data) or "plot_fluxMap_gambling.m" (for the gambling data). We also include "plot_fluxDivergence.m", which illustrates the flux divergences plotted in Fig. S3 in the Supporting Information.

For queries or issues, please contact Chris Lynn at cwlynn@princeton.edu.
