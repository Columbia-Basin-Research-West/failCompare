# failmod package

### Description 

This repository contains the intial attempt to build and publish an R package for fitting and comparing the fit of failure time models.

### PENDING 

#### Workhorse functions

```r
fail_select(fmods,model)
```
input failmod_list, select one of the models from a failmod_list
returns a failmod_obj for correcting.
DONE!!!

```r
fail_compare(obj,group,test)
```
Helper function that use log-rank. F dist tests? to determine lots should be modeled separately.


#### Generic functions

```r
predict.failmod_obj(x,newdata,shift)
```
Input failmod_obj.Predict survival probability given new "times" data depending on the model. A temporal shift arguement also a good idea

```r
plot.failmod_obj()
```

Quick visualization of parmeteric model vs. the KM estimates
DONE!!!
```r
plot.failmod_list()
```
Quick visualization of up to three parmeteric model vs. the KM estimates. If rank table present, then default to displaying 3 top ranking models
DONE!!!
