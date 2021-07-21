# failCompare (Failure Time Model Fitting and Goodness of Fit Comparison)

(This package is still in a pre-release phase)

### Description 

A package for fitting and comparing among failure-time/survival models from the F distribution and vitality families, which facilitates visualizing andranking the performance of the model using the Skalski and Whitlock (2020) GOF measure.

### Installation

```r
devtools::install_github("swhitCBR/failCompare")
```

*Rtools required for building package from source

#### Main functions

```r
fc_fit(time, model,...)
```
Model fitting routine used to fit one or a set of failure time models. If a single model is specified a "fc_obj" is created, which can be used to adjust a CJS model in the forthcoming "ATLAS" package. If multiple models are specified a "failmod_list" is created containing output from all model fits. 

```r
# Defining list of three models to rank
mod_list = fc_fit(time, model=c("gompertz","weibull3,"vitality.ku"))

# Ranking the list of model
fc_rank(mod_list)
```
Ranks the performance of the model using the Skalski and Whitlock (2020) GOF measure.

### Additional resources

Vitality model information<br>
http://www.cbr.washington.edu/analysis/vitality

GOF measure definition<br>
[Skalski and Whitlock (2020)](http://animalbiotelemetry.biomedcentral.com/articles/10.1186/s40317-020-00213-z)

