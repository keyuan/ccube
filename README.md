[![Build Status](https://travis-ci.org/keyuan/ccube.svg?branch=master)](https://travis-ci.org/keyuan/ccube)
 [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/keyuan/ccube?branch=master)](https://ci.appveyor.com/project/keyuan/ccube)

Ccube 
=====
`ccube` is an R package for clustering and estimating cancer cell fractions (CCF) of somatic variants (SNVs/SVs) from bulk whole genome/exome data. The method takes the reference and alternative allele read counts of called variants, corrects for copy number alterations and purity, then produces CCF estimates for all variants within the tumour sample. It identifies clusters of mutations, which can be used to determine the clonal architecture of the sample. 

Models
----
The package contains four Bayesian mixture models, all fitted with variational inference. 

- Ccube: Normal-Binomial mixture model, used for clustering and estimating CCFs of SNVs. Details can be found in this [manuscript](https://www.biorxiv.org/content/10.1101/484402v1.abstract). 
- CcubeSV: Normal-Binomial mixture model, similar with Ccube but modified for clustering and estimating CCFs of SVs. Details can be found in this [manuscript](https://www.biorxiv.org/content/10.1101/172486v1.abstract) and [repo](https://github.com/mcmero/SVclone).
- Student-t mixture model: Main model for purity estimation. Implements model described in this [paper](https://www.sciencedirect.com/science/article/pii/S0893608006001791?via%3Dihub).
- Normal mixture model: Alternative model for purity estimation, also used for code calibration and testing.

Installation  
----
```r
require(devtools)
devtools::install_github("keyuan/ccube")
```

Getting Started
----
```
vignettes/ccube.Rmd
```
