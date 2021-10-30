# RGMMBench

Quick links: [mixtools](https://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf) | [mclust](https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html)

## Showcase

This repository contains the code used to generate automatically summary data and figures of the manuscript "Gaussian mixtures in R". Its main purpose is to compare 
automatically computational and statistical performances of R packages estimating Gaussian mixture models (temporary in univariate dimension). Especially, we compare 
packages  **Rmixmod**, **mixtools**, **bgmm**, **mclust**, **EMCluster**, **GMKMcharlie**, **flexmix** and **DCEM**. Additionally, **otrimle** is provided to estimate parameters in case of outliers and **mixsmnsn** is dedicated for the estimation of skewed GMMs.

![Boxplot of the estimated parameters with four overlapping and unbalanced components](https://github.com/bastienchassagnol/RGMMBench/images/four_components_unbalanced_overlapping_boxplots.pdf)

<details>
    <summary>Click here to display the source code</summary>

```R
# load useful libraries and packages
library(ggplot2)
import::from(magrittr, "%>%", .into = "operators") 
import::from(rebmix, .except = c("AIC", "BIC", "split"))
library(mclust)
library(Rmixmod)


relevant_mixture_functions <- list ("otrimle"=list(name_fonction=em_otrimle, list_params=list()),
                                    "mixsmsn"=list(name_fonction=em_mixsmsn, list_params=list()),
                                    "em R" = list(name_fonction=emnmix, list_params=list()),
                                    "Rmixmod" = list(name_fonction=em_Rmixmod, list_params=list()),
                                    "mixtools" = list(name_fonction=em_mixtools, list_params=list()),                                     
                                    "bgmm"= list(name_fonction=em_bgmm, list_params=list()),
                                    "mclust" = list(name_fonction=em_mclust, list_params=list(prior = NULL)),
                                    "EMCluster" = list(name_fonction=em_EMCluster, list_params=list()),
                                    "GMKMcharlie"=list(name_fonction=em_GMKMcharlie, list_params=list()),
                                    "flexmix"= list(name_fonction=em_flexmix, list_params=list()),
                                    "DCEM"=list(name_fonction=em_DCEM, list_params=list()))

##################################################################
##      Compare computational performances of the packages      ##
##################################################################
four_components_statistical_performances <- benchmark_distribution_parameters(mixture_functions=relevant_mixture_functions,
                                                                             sigma_values=list("high OVL"= rep(2, 4)),
                                                                             mean_values=list(c(0, 4, 8, 12)),
                                                                             proportions = list("highly unbalanced"=c(0.1, 0.7, 0.1, 0.1)),
                                                                             skewness_values = list("null skewness"=rep(0, 4),
                                                                             Nbootstrap=200,  nobservations=c(2000)))
```

#################################################################
##  Save results (example with the four components simulation  ##
#################################################################

# save summary scores and distributions of the bootstrap simulations
openxlsx::write.xlsx(four_components_statistical_performances$local_scores,file = "tables/four_components_local_scores.xlsx", asTable = T)
openxlsx::write.xlsx(four_components_statistical_performances$global_scores,file = "tables/four_components_global_scores.xlsx", asTable = T)
openxlsx::write.xlsx(four_components_statistical_performances$distributions,file = "tables/four_components_distributions.xlsx", asTable = T)

# save boxplots associated to the distribution of the estimates
unbalanced_overlapping_boxplots <- four_components_computational_performances$plots$`2000_observations_UR_0.9_skewness_0_OVL_0.08_prop_outliers_0`
ggsave("images/four_components_unbalanced_overlapping_boxplots.pdf", unbalanced_overlapping_boxplots,
       width = 15, height = 14,dpi = 600)

</details>



## Install

To get the most recent version, open `R` and run:

```R
if(!require(devtools)) install.packages("devtools")
devtools::install_github("bastienchassagnol/RGMMBench")
```

The package is composed of four scripts: **main** contains the main script to load required libraries and executes auxiliary functions to reproduce figures and tables of the paper
Gaussian Mixtures in R, **mixture** enlists the functions used to simulate a Gaussian mixture and estimate its parameters using the EM algorithm, benchmark gathers the two functions
used to compare the computational and statistical performances of R packages in learning GMMs. Finally, visualisation displays three functions, two to compare graphically the computational
performances of the packages, and one representing the boxplots of the bootstrap simulations.

## Get inspired

ArXiv publication associated with the paper:

- Chassagnol et al. (2021). [Gaussian mixtures in R](https://doi.org/10.1007/s11192-020-03690-4), *The R journal*
- Scrucca et al. (2016). [mclust 5: Clustering, Classification and Density Estimation Using Gaussian Finite Mixture Models], *The R journal*

## Acknowledgments

Originally developed from an original course on mixture models and use of the EM algorithm for complex MLE estimation supplied by Gregory Nuel
