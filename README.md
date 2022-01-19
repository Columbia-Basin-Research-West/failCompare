
# failCompare  
<img src="man/figures/failCompare manual thumbnail.png" align="right" height="400" hspace="15" />  <br>
*failCompare* is a package for fitting and comparing failure-time/survival models from among the F distribution and [vitality](http://www.cbr.washington.edu/analysis/vitality) families. The package facilitates the visualization and ranks the performance of the various models using the [Skalski and Whitlock (2020)]((http://animalbiotelemetry.biomedcentral.com/articles/10.1186/s40317-020-00213-z)) goodness-of-fit measure. Dependencies include the [survival](https://cran.r-project.org/web/packages/survival/index.html), [flexsurv](https://cran.r-project.org/web/packages/flexsurv/index.html), and [Vitality](https://cran.r-project.org/web/packages/vitality/index.html) R packages. Additional features include: (1) handling of Type I and Type II right-censored data, (2) specialized plotting for examining model fit, and (3) generalized one-sample Kolmogorov-Smirnov tests of lack-of-fit ([Lilliefors 1967](https://www.tandfonline.com/doi/abs/10.1080/01621459.1967.10482916)).

**User manual:**
[failCompare (Version 1.0)](http://www.cbr.washington.edu/analysis/apps/failcompare)

The package was developed by Columbia Basin Research (CBR), an interdisciplinary research center in the [School of Aquatic and Fishery Sciences](http://www.fish.washington.edu/index.html) at the [University of Washington](http://www.washington.edu/) whose mission is to develop quantitative approaches and provide user-friendly tools to aid in the study of endangered salmonid stocks in the western North America.

Please visit the *failCompare* page on the [CBR website](http://www.cbr.washington.edu/analysis/apps/failcompare) to download the binary version of the package, view the current user manual, and access additional information on this and other software tools.

**Creators**: Steve Whitlock, Rebecca Buchanan, and Rich Townsend

#### Releases
-11/30/2021 on the CBR site <br>
-Coming soon to CRAN

**License:** [GPL-3](https://cran.r-project.org/web/licenses/GPL-3)


## Installation

#### Dependencies
*failCompare* depends upon a prior installation of the following survival, flexsurv, and vitality R packages.
Run the following code to install or update these packages
```
install.packages(c("survival", "flexsurv", "vitality"))
```

#### From binary 
Download (.zip) file here: 
[failCompare (Version 1.0)](http://www.cbr.washington.edu/sites/default/files/manuals/failCompare%20User%20Manual_0.pdf)

Select the option for "installing package(s) from local files..." or enter:
```
utils:::menuInstallLocal()
```
Navigate to and select the downloaded .zip folder.

#### From source code

If package `devtools` is installed:
```r
devtools::install_github("swhitCBR/failCompare")
```
Rtools must also be installed to build a package from the source code. Instructions for installing Rtools (on [Windows](https://cran.r-project.org/bin/windows/Rtools/) and [Mac OS X](https://cran.r-project.org/bin/macosx/tools/)) may be found on the CRAN main site.

## Usage

For detailed descriptions and examples, please refer to the user manual PDF: <br>
[failCompare Manual (Version 1.0)](http://www.cbr.washington.edu/sites/default/files/manuals/failCompare%20User%20Manual.pdf)


### Main Functions

#### fc_fit()
The model fitting routine used to fit one or a set of failure time models. If multiple models are specified an "fc_list" is created containing output from all model fits.

```r
# Where "failure_times" is a numeric vector of only positive values
# Fitting a 3-parameter Weibull model
mod_obj <- fc_fit(time=failure_times, model="weibull3")
```
The following code defines a list of four models to be ranked: (1) Gompertz; (2) 3-parameter Weibull; (3) Vitality 2009 ; and (4) Vitality 2013 models
```
mod_list <- fc_fit(time, model=c("gompertz","weibull3,"vitality.ku","vitality.4p"))

```
#### fc_rank()
Ranks the performance of the model using the Skalski and Whitlock (2020) goodness-of-fit measure.

```r
# Ranking the list of model
fc_rank(mod_list)
```

#### fc_test()
Provides lack-of-fit testing for a single model based on a Monte Carlo simulation with 10k iterations
```r
fc_test(time=failure_times,model=mod_obj,iters=10000)
```
## Association with *cbrATLAS* package

Once a failure model object (of class `fc_obj`) has been created it may be used to adjust for tag failure in a CJS model fit using the R pacakge [cbrATLAS](https://github.com/Columbia-Basin-Research-West/ATLAS) package. 

