
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saeHB.Spatial.Beta

Provides several functions and datasets for area-level Small Area
Estimation using the Hierarchical Bayesian (HB) method. Model-based
estimators are designed for variables of interest that follow a Beta
distribution (proportions bounded between 0 and 1). The package supports
spatial structures under the Simultaneous Autoregressive (SAR) model and
the Leroux Conditional Autoregressive (CAR) model. It also accommodates
survey design effect (DEFF) adjustments to handle complex survey data.
The `rjags` package is employed to obtain parameter estimates via Markov
Chain Monte Carlo (MCMC).For the reference, see Rao and Molina (2015)
<doi:10.1002/9781118735855>.

## Author

Boby Iwan, Cucu Sumarni

## Maintainer

Boby Iwan <bobyiwanboby2122@gmail.com>

## Functions

- `betaDeffSAR()` Estimates small area means using Spatial SAR Model
  with Beta distribution and Design Effect (DEFF) adjustments.
- `betaSAR()` Estimates small area means using Spatial SAR Model with
  Beta distribution without DEFF adjustments.
- `betaDeffLerouxCAR()` Estimates small area means using Spatial Leroux
  CAR Model with Beta distribution and Design Effect (DEFF) adjustments.
- `betaLerouxCAR()` Estimates small area means using Spatial Leroux CAR
  Model with Beta distribution without DEFF adjustments.
- `betaDeffNonSpatial()` Estimates small area means using a Non-Spatial
  Beta Model with Independent and Identically Distributed (IID) random
  effects and DEFF adjustments.
- `betaNonSpatial()` Estimates small area means using a Non-Spatial Beta
  Model without DEFF adjustments.
- `build_W()` A utility function to construct spatial weights matrices
  (contiguity, distance, or kernel) required for spatial modeling.
- `spatial_moran()` A diagnostic function to perform Moran’s I test for
  spatial autocorrelation.

## Installation

You can install the development version of saeHB.Spatial.Beta from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BobyIwan/saeHB.Spatial.Beta")
```

## Example

This is a basic example of using the `betaDeffSAR()` function to make an
estimate based on synthetic data in this package:

``` r
library(saeHB.Spatial.Beta)

# Load dataset and proximity matrix
data(dataBeta)
data(weight_mat)

# Fitting the Spatial SAR model
model_sarDeff <- betaDeffSAR(
  formula = y ~ x1 + x2,
  DEFF = "deff",
  n_i = "n_i",
  proxmat = weight_mat,
  data = dataBeta,
  n.adapt = 500
)
```

<img src="man/figures/README-example-1.png" alt="" width="100%" /><img src="man/figures/README-example-2.png" alt="" width="100%" /><img src="man/figures/README-example-3.png" alt="" width="100%" />

Extract the mean estimation for the areas:

``` r
head(model_sarDeff$Est)
#>        Estimate  Est.Error   l-95% CI  u-95% CI
#> mu[1] 0.8784156 0.04693805 0.75802237 0.9425771
#> mu[2] 0.6269277 0.10464048 0.42570608 0.8392672
#> mu[3] 0.5816371 0.07945872 0.42829699 0.7439266
#> mu[4] 0.2173992 0.06685184 0.10416101 0.3559373
#> mu[5] 0.2133118 0.07370978 0.07678379 0.3639901
#> mu[6] 0.9436485 0.03235324 0.87259401 0.9869702
```

Extract the estimated model coefficients:

``` r
model_sarDeff$coefficient
#>          Estimate  Est.Error  l-95% CI  u-95% CI     Rhat       ESS
#> beta[0] 1.9369675 0.26837727 1.3527652 2.3800517 1.531591  12.28722
#> beta[1] 0.7997211 0.15218007 0.5135726 1.0837170 1.026906  36.98647
#> beta[2] 0.7181547 0.09340642 0.5334150 0.8996733 1.264257 128.53339
#> rho     0.6856899 0.14729914 0.3389427 0.9094680 1.005473 186.02937
```

Extract the random effect for the areas:

``` r
model_sarDeff$randeff
#>          Estimate Est.Error    l-95% CI    u-95% CI
#> v[1]  -0.10783106 0.4694424 -1.22314383  0.65208038
#> v[2]  -0.83092581 0.4648043 -1.78006765  0.09412078
#> v[3]  -1.63847130 0.3707255 -2.41374318 -0.92099516
#> v[4]  -2.51754355 0.5581144 -3.43731502 -1.35424193
#> v[5]  -2.30098568 0.5887278 -3.74067841 -1.29130149
#> v[6]  -0.52233714 0.8005825 -1.72470315  1.05442508
#> v[7]  -1.66925454 0.4698649 -2.60536600 -0.82500371
#> v[8]  -0.57121330 0.5285002 -1.57143948  0.44252808
#> v[9]  -1.30153688 0.6162749 -2.51395763 -0.25232119
#> v[10] -1.67748450 0.6307478 -2.66634995 -0.06789242
#> v[11] -0.84583663 0.5816477 -1.80148401  0.47766972
#> v[12] -1.46822314 0.9819288 -3.29026773  0.74777123
#> v[13] -0.21620930 0.5482899 -1.24473605  0.84657183
#> v[14] -0.66745317 0.5388253 -1.71359029  0.44634915
#> v[15] -0.50636718 0.4323544 -1.40563726  0.34476162
#> v[16] -1.49723267 0.5563584 -2.61034980 -0.40560201
#> v[17]  0.21941788 0.6118077 -0.89912832  1.35873935
#> v[18]  0.17845158 0.6430400 -1.03012296  1.36311163
#> v[19]  1.19992376 0.4916931  0.23987611  2.21252226
#> v[20] -1.81351622 0.4833503 -2.71441080 -0.78701459
#> v[21] -1.55952188 0.4242196 -2.43860511 -0.65552459
#> v[22] -0.37897902 0.4390642 -1.16291765  0.46089258
#> v[23] -0.45795712 0.5899375 -1.39252273  0.82914422
#> v[24]  1.71089629 0.6958473  0.63023164  3.31862436
#> v[25]  1.04035629 0.6770139 -0.08068904  2.64834422
#> v[26]  0.20606932 0.5083907 -0.78327292  1.18082533
#> v[27] -0.79603822 0.6558185 -2.00496089  0.44080602
#> v[28]  0.17091684 0.7165652 -1.20628130  1.58143535
#> v[29]  0.76915414 0.5874883 -0.35246676  1.81390827
#> v[30]  1.18239101 0.6650039  0.03475275  2.48468866
#> v[31] -0.15064489 0.7903554 -1.67726815  1.28887993
#> v[32]  0.08943888 0.6771465 -1.34644607  1.09438958
#> v[33]  0.69779495 0.5941036 -0.21665746  2.08466342
#> v[34]  0.80000180 0.5564170 -0.07084046  2.07429692
#> v[35]  0.96171072 0.6598250 -0.25961727  2.54931781
#> v[36]  1.61171337 0.4533350  0.73045626  2.48363024
```

Extract the random effect variance for the areas:

``` r
model_sarDeff$refVar
#>           Estimate Est.Error  l-95% CI u-95% CI
#> a.var[1]  3.760232  43.28017 0.9646746 6.857557
#> a.var[2]  3.677417  43.38839 0.9380096 6.610790
#> a.var[3]  3.496129  43.42892 0.9079578 5.975525
#> a.var[4]  3.496129  43.42892 0.9079578 5.975525
#> a.var[5]  3.677417  43.38839 0.9380096 6.610790
#> a.var[6]  3.760232  43.28017 0.9646746 6.857557
#> a.var[7]  3.677417  43.38839 0.9380096 6.610790
#> a.var[8]  3.620138  43.56685 0.9215755 6.462746
#> a.var[9]  3.442653  43.61167 0.8892324 5.825479
#> a.var[10] 3.442653  43.61167 0.8892324 5.825479
#> a.var[11] 3.620138  43.56685 0.9215755 6.462746
#> a.var[12] 3.677417  43.38839 0.9380096 6.610790
#> a.var[13] 3.496129  43.42892 0.9079578 5.975525
#> a.var[14] 3.442653  43.61167 0.8892324 5.825479
#> a.var[15] 3.303048  43.67095 0.8639656 5.365353
#> a.var[16] 3.303048  43.67095 0.8639656 5.365353
#> a.var[17] 3.442653  43.61167 0.8892324 5.825479
#> a.var[18] 3.496129  43.42892 0.9079578 5.975525
#> a.var[19] 3.496129  43.42892 0.9079578 5.975525
#> a.var[20] 3.442653  43.61167 0.8892324 5.825479
#> a.var[21] 3.303048  43.67095 0.8639656 5.365353
#> a.var[22] 3.303048  43.67095 0.8639656 5.365353
#> a.var[23] 3.442653  43.61167 0.8892324 5.825479
#> a.var[24] 3.496129  43.42892 0.9079578 5.975525
#> a.var[25] 3.677417  43.38839 0.9380096 6.610790
#> a.var[26] 3.620138  43.56685 0.9215755 6.462746
#> a.var[27] 3.442653  43.61167 0.8892324 5.825479
#> a.var[28] 3.442653  43.61167 0.8892324 5.825479
#> a.var[29] 3.620138  43.56685 0.9215755 6.462746
#> a.var[30] 3.677417  43.38839 0.9380096 6.610790
#> a.var[31] 3.760232  43.28017 0.9646746 6.857557
#> a.var[32] 3.677417  43.38839 0.9380096 6.610790
#> a.var[33] 3.496129  43.42892 0.9079578 5.975525
#> a.var[34] 3.496129  43.42892 0.9079578 5.975525
#> a.var[35] 3.677417  43.38839 0.9380096 6.610790
#> a.var[36] 3.760232  43.28017 0.9646746 6.857557
```

## References

- Rao, J. N. K., & Molina, I. (2015). Small Area Estimation (2nd
  Edition). New Jersey: John Wiley and Sons,
  Inc. <doi:10.1002/9781118735855>.
- Kubacki, J., & Jedrzejczak, A. (2016). Small Area Estimation of Income
  Under Spatial SAR Model. Statistics in Transition New Series, Vol. 17,
  No. 3, pp. 365–390. <doi:10.59170/stattrans-2016-022>.
- Leroux, B. G., Lei, X., & Breslow, N. (2000). Estimation of Disease
  Rates in Small Areas: A New Mixed Model for Spatial Dependence.
  In M. E. Halloran & D. Berry (Eds.), Statistical Models in
  Epidemiology, the Environment, and Clinical Trials (Vol. 116,
  pp. 179–191). New York: Springer. <doi:10.1007/978-1-4612-1284-3_4>.
- Chung, H. C., & Datta, G. S. (2020). Bayesian Hierarchical Spatial
  Models for Small Area Estimation. Research Report Series. Washington,
  D.C.: U.S. Census Bureau.
