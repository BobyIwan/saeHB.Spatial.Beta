
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saeHB.Spatial.Beta

We designed this package to provide several functions for area-level
small area estimation under Spatial Simultaneous Autoregressive (SAR)
and Leroux Conditional Autoregressive (CAR) models, accommodating survey
design effect (DEFF) adjustments, using hierarchical Bayesian (HB)
method with Beta distribution for variables of interest. Some datasets
simulated by a data generation are also provided. The rjags package is
employed to obtain parameter estimates using Markov Chain Monte Carlo
(MCMC) algorithms. Model-based estimators involve the HB estimators
which include the mean estimation, the estimated model coefficients, the
random effect, and the random effect variance. For the reference, see
Rao and Molina (2015), Liu (2009), Liu et al. (2014), Kubacki and
Jedrzejczak (2016), Leroux et al. (2000), Chung and Datta (2020),
Anselin (1988), and Anselin and Morrison (2019).

## Author

Boby Iwan, Cucu Sumarni

## Maintainer

Boby Iwan <bobyiwanboby2122@gmail.com>

## Functions

- `betadeff_sar()` Estimates small area proportions using a Hierarchical
  Bayes (HB) method under a Spatial SAR Model with a Beta distribution,
  incorporating survey design effect (DEFF) adjustments.
- `beta_sar()` Estimates small area proportions using a Hierarchical
  Bayes (HB) method under a Spatial SAR Model with a Beta distribution
  without DEFF adjustments, by estimating the unknown precision
  parameter.
- `betadeff_lerouxcar()` Estimates small area proportions using a
  Hierarchical Bayes (HB) method under a Spatial Leroux CAR Model with a
  Beta distribution, incorporating survey design effect (DEFF)
  adjustments.
- `beta_lerouxcar()` Estimates small area proportions using a
  Hierarchical Bayes (HB) method under a Spatial Leroux CAR Model with a
  Beta distribution without DEFF adjustments, by estimating the unknown
  precision parameter.
- `betadeff_nonspatial()` Estimates small area proportions using a
  Hierarchical Bayes (HB) method under a Non-Spatial Model with a Beta
  distribution and Independent and Identically Distributed (IID) random
  effects, incorporating DEFF adjustments.
- `beta_nonspatial()` Estimates small area proportions using a
  Hierarchical Bayes (HB) method under a Non-Spatial Model with a Beta
  distribution and IID random effects without DEFF adjustments, by
  estimating the unknown precision parameter.
- `build_w()` A utility function to construct spatial weights matrices
  (contiguity, distance, or kernel) required for spatial modeling.
- `moran_test()` A diagnostic function to perform Moran’s I test for
  spatial autocorrelation.

## Installation

**System Requirement:** Since this package relies on `rjags` for MCMC
computations, you must first install the [JAGS (Just Another Gibbs
Sampler)](https://mcmc-jags.sourceforge.io/) software on your computer
before installing this package.

You can install the development version of saeHB.Spatial.Beta from
[GitHub](https://github.com/BobyIwan/saeHB.Spatial.Beta) with:

``` r
# install.packages("devtools")
devtools::install_github("BobyIwan/saeHB.Spatial.Beta")
```

Or, to include the vignette, use the following command:

``` r
devtools::install_github("BobyIwan/saeHB.Spatial.Beta", build_vignettes = TRUE)
```

## Example

This is a basic example of using the `betadeff_sar()` function to make
an estimate based on synthetic data in this package:

``` r
library(saeHB.Spatial.Beta)

# Load dataset and proximity matrix
data(databeta)
data(weight_mat)

# Fitting the Spatial SAR model
model_sar_deff <- betadeff_sar(
  formula = y ~ x1 + x2,
  deff = "deff",
  n_i = "n_i",
  proxmat = weight_mat,
  data = databeta
)
```

<img src="man/figures/README-example-1.png" alt="" width="100%" /><img src="man/figures/README-example-2.png" alt="" width="100%" /><img src="man/figures/README-example-3.png" alt="" width="100%" />

Extract the mean estimation for the areas:

``` r
head(model_sar_deff$est)
#>        Estimate   Est.Error  l-95% CI  u-95% CI
#> mu[1] 0.6904202 0.122276168 0.4491634 0.9066046
#> mu[2] 0.7358053 0.075425351 0.5852298 0.8758472
#> mu[3] 0.8315318 0.048186333 0.7266684 0.9207047
#> mu[4] 0.9098379 0.052079819 0.7889493 0.9734647
#> mu[5] 0.9147718 0.034302527 0.8348752 0.9695360
#> mu[6] 0.9933859 0.003904432 0.9845815 0.9986099
```

Extract the estimated model coefficients:

``` r
model_sar_deff$coefficient
#>          Estimate  Est.Error  l-95% CI u-95% CI     Rhat       ESS
#> beta[0] 2.1552971 0.18871070 1.7942579 2.469571 2.351008  63.83022
#> beta[1] 0.9846968 0.16170654 0.6613031 1.298958 1.822902  43.60133
#> beta[2] 0.9015832 0.09514377 0.7211005 1.099114 1.044164  51.33361
#> rho     0.7882047 0.10719963 0.5517969 0.959489 1.215377 536.84195
```

Extract the random effect for the areas:

``` r
model_sar_deff$randeff
#>          Estimate Est.Error    l-95% CI    u-95% CI
#> v[1]  -1.55295590 0.5415769 -2.46067324 -0.30183277
#> v[2]  -0.40668766 0.4868080 -1.27137753  0.68989860
#> v[3]  -0.56451854 0.3761887 -1.29802589  0.20921295
#> v[4]   1.25377103 0.6708254  0.02248222  2.62497249
#> v[5]   1.57656883 0.4135432  0.83389540  2.37519730
#> v[6]   1.07018762 0.6214890  0.08702998  2.30824241
#> v[7]  -2.25189587 0.3772757 -2.96683152 -1.43691527
#> v[8]  -2.22335204 0.4368014 -2.92910385 -1.40066119
#> v[9]  -0.34010596 0.4563517 -1.26326303  0.45044624
#> v[10]  0.64987774 0.6773952 -0.64682674  2.08807649
#> v[11]  2.04639123 0.6617645  0.97155382  3.78632657
#> v[12]  1.09010398 0.6111019 -0.11469939  2.24601657
#> v[13] -1.05056022 0.6390102 -2.11056767  0.54101461
#> v[14] -0.97210465 0.5683333 -2.50486357 -0.06298717
#> v[15] -0.29789316 0.6093787 -1.30150329  0.92165587
#> v[16]  1.03759858 0.5797622  0.17884612  2.48395933
#> v[17]  1.63925819 0.5189512  0.78898753  2.78170881
#> v[18]  0.69487104 0.9392959 -1.42223674  2.17441560
#> v[19] -0.84008891 0.5460784 -1.92505921  0.17527037
#> v[20] -0.26573043 0.5216119 -1.40333441  0.60022711
#> v[21]  0.51006130 0.5192477 -0.50860557  1.41488293
#> v[22]  0.05989928 0.2736546 -0.48650270  0.61845032
#> v[23]  1.45811532 0.6405493  0.54155041  2.95391925
#> v[24]  1.41442004 0.7454328  0.09578068  3.01440523
#> v[25] -0.77148674 0.7099261 -2.04790565  0.60050748
#> v[26]  0.47685012 0.7054332 -0.95129397  2.12708828
#> v[27]  0.26415565 0.7737525 -1.08873258  1.94640679
#> v[28]  1.04560152 0.4962706  0.14640639  2.03711501
#> v[29]  0.80995469 0.7839949 -0.28068586  2.44469111
#> v[30]  1.52865418 0.8541867  0.20885108  3.05165411
#> v[31]  0.40121995 0.6305752 -0.64229385  1.68319967
#> v[32]  0.15909805 0.5044961 -0.71972688  1.12836390
#> v[33]  0.14512840 0.5387659 -1.25674861  1.07377001
#> v[34]  0.87835185 0.6600124 -0.18512457  2.35154680
#> v[35]  2.64810355 0.8006772  1.33574917  4.21903674
#> v[36]  2.02274589 0.9177265  0.64723607  3.78524746
```

Extract the random effect variance for the areas:

``` r
model_sar_deff$refvar
#>           Estimate Est.Error  l-95% CI u-95% CI
#> a.var[1]  4.920737  46.68959 0.9311315 15.25944
#> a.var[2]  4.837588  46.79251 0.8967958 15.18369
#> a.var[3]  4.620072  46.82920 0.8333057 14.65618
#> a.var[4]  4.620072  46.82920 0.8333057 14.65618
#> a.var[5]  4.837588  46.79251 0.8967958 15.18369
#> a.var[6]  4.920737  46.68959 0.9311315 15.25944
#> a.var[7]  4.837588  46.79251 0.8967958 15.18369
#> a.var[8]  4.793385  46.96290 0.8680262 15.27067
#> a.var[9]  4.578069  47.00342 0.8075354 14.69694
#> a.var[10] 4.578069  47.00342 0.8075354 14.69694
#> a.var[11] 4.793385  46.96290 0.8680262 15.27067
#> a.var[12] 4.837588  46.79251 0.8967958 15.18369
#> a.var[13] 4.620072  46.82920 0.8333057 14.65618
#> a.var[14] 4.578069  47.00342 0.8075354 14.69694
#> a.var[15] 4.407630  47.05779 0.7676809 14.18318
#> a.var[16] 4.407630  47.05779 0.7676809 14.18318
#> a.var[17] 4.578069  47.00342 0.8075354 14.69694
#> a.var[18] 4.620072  46.82920 0.8333057 14.65618
#> a.var[19] 4.620072  46.82920 0.8333057 14.65618
#> a.var[20] 4.578069  47.00342 0.8075354 14.69694
#> a.var[21] 4.407630  47.05779 0.7676809 14.18318
#> a.var[22] 4.407630  47.05779 0.7676809 14.18318
#> a.var[23] 4.578069  47.00342 0.8075354 14.69694
#> a.var[24] 4.620072  46.82920 0.8333057 14.65618
#> a.var[25] 4.837588  46.79251 0.8967958 15.18369
#> a.var[26] 4.793385  46.96290 0.8680262 15.27067
#> a.var[27] 4.578069  47.00342 0.8075354 14.69694
#> a.var[28] 4.578069  47.00342 0.8075354 14.69694
#> a.var[29] 4.793385  46.96290 0.8680262 15.27067
#> a.var[30] 4.837588  46.79251 0.8967958 15.18369
#> a.var[31] 4.920737  46.68959 0.9311315 15.25944
#> a.var[32] 4.837588  46.79251 0.8967958 15.18369
#> a.var[33] 4.620072  46.82920 0.8333057 14.65618
#> a.var[34] 4.620072  46.82920 0.8333057 14.65618
#> a.var[35] 4.837588  46.79251 0.8967958 15.18369
#> a.var[36] 4.920737  46.68959 0.9311315 15.25944
```

## References

- Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd ed.).
  New Jersey: John Wiley & Sons, Inc. <doi:10.1002/9781118735855>.
- Liu, B. (2009). *Hierarchical Bayes estimation and empirical best
  prediction of small-area proportions*.
  <https://api.drum.lib.umd.edu/server/api/core/bitstreams/cb8e2cbf-441e-4f0f-b4b3-6182f3cf24de/content>.
- Liu, B., Lahiri, P., & Kalton, G. (2014). Hierarchical Bayes Modeling
  of Survey-Weighted Small Area Proportions. *Statistics Canada*.
  <https://www150.statcan.gc.ca/n1/pub/12-001-x/2014001/article/14030-eng.pdf>.
- Kubacki, J., & Jedrzejczak, A. (2016). Small Area Estimation of Income
  Under Spatial SAR Model. *Statistics in Transition New Series*,
  *17*(3), 365–390. <doi:10.59170/stattrans-2016-022>.
- Leroux, B. G., Lei, X., & Breslow, N. (2000). Estimation of Disease
  Rates in Small Areas: A New Mixed Model for Spatial Dependence.
  In M. E. Halloran & D. Berry (Eds.), *Statistical Models in
  Epidemiology, the Environment, and Clinical Trials* (Vol. 116,
  pp. 179–191). New York: Springer. <doi:10.1007/978-1-4612-1284-3_4>.
- Chung, H. C., & Datta, G. S. (2020). *Bayesian Hierarchical Spatial
  Models for Small Area Estimation* (Research Report Series).
  Washington, D.C.: U.S. Census Bureau.
  <https://www.census.gov/content/dam/Census/library/working-papers/2020/adrm/RRS2020-07.pdf>.
- Anselin, L. (1988). *Spatial Econometrics: Methods and Models*.
  Dordrecht: Springer Netherlands. <doi:10.1007/978-94-015-7799-1>.
- Anselin, L., & Morrison, S. (2019). *Spatial Weights as Distance
  Functions*.
  <https://spatialanalysis.github.io/lab_tutorials/Spatial_Weights_as_Distance_Functions.html>.
