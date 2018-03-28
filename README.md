
STAMP Subset Testing and Analysis of Multiple Phenotypes
========================================================

Version 0.0

Authors: Andriy Derkach, Ruth Pfeiffer

Maintainer: Andriy Derkach (<andriy.derkach@nih.gov>)

Description: We develop a flexible procedure (STAMP) based on mixture models to perform region based meta-analysis of different phenotypes using data from different GWAS and identify subsets of associated phenotypes. To run variance-component (quadratic test) under mixture model use mixtureQuad.r and mixtureLinear.r to run for burden (linear) tests.

Depndens: R (&gt;= 3.2.1), psych, pscl, optimx, nleqslv, mvnfast, CompQuadForm

Installation
------------

Use [devtools](https://github.com/hadley/devtools) to install:

``` r
library(devtools)
install_github("matt-sutton/STAMP") #Update gitrepo after pull request
```

Example Usage
-------------

Import example

``` r
library(STAMP)
data(Example_SC11)
names(Example_SC11)
```

    ## [1] "ListBeta" "ListCor"  "ListSE"   "MafList"

``` r
length(Example_SC11$ListBeta)
```

    ## [1] 20

``` r
length(Example_SC11$ListBeta[[1]])
```

    ## [1] 210

Run variance-component (quadratic test) under mixture model:

``` r
mixtureQuad(1000,Example_SC11$ListBeta,Example_SC11$ListSE,Example_SC11$ListCor,Example_SC11$MafList)
```
