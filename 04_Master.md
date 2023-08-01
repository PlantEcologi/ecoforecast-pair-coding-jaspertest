



This report was automatically generated with the R package **knitr**
(version 1.39).


```r
---
title: "Master script for postfire analysis"
output: html_document
---
```

```
## Error: <text>:6:0: unexpected end of input
## 4: ---
## 5: 
##   ^
```


### 1. Source functions, get data and plot

First we'll _source()_ (i.e. "run all code in") the scripts with the functions we made. Then we'll set the URL, read in the data with _download.NDVI()_, and plot it with _plot.NDVI()_.


```r
## Load required functions by running source() on the individual function files
if(file.exists("01_download.NDVI.R")) source("01_download.NDVI.R")
if(file.exists("02_plot.NDVI.R"))     source("02_plot.NDVI.R")
if(file.exists("03_negexp.R"))        source("03_negexp.R")

## Download NDVI data
URL = "https://raw.githubusercontent.com/jslingsby/BIO3019S_Ecoforecasting/master/data/modisdata.csv"
dat <- download.NDVI(URL)

# Convert "calendar_date" to postfire age in days since fire - assuming the first date in the times eries is the time of the fire 
dat$age <- (as.numeric(dat$calendar_date) - min(as.numeric(dat$calendar_date), na.rm = T))/365.25

## Plot overall NDVI time series
plot.NDVI(dat)
```

<img src="figure/04-Master-Rmdunnamed-chunk-1-1.png" title="plot of chunk unnamed-chunk-1" alt="plot of chunk unnamed-chunk-1" style="display: block; margin: auto;" />

<br>

Q1: This plot suggests that Fynbos greenness (NDVI) as observed from satellite saturates with time since fire. Why do you think it saturates rather than increasing linearly with time?

>*Answer 1:*

<br>

### 2. Fit models using Non-linear Least Squares (NLS)

Now we'll fit the simple and full negative exponential models using Non-linear Least Squares (NLS).

First the simpler model:


```r
## Simple model

# set parameters
par <- c(alpha = 0.2, gamma = 0.4, lambda = 0.5)

# fit model
fit_negexp <- nls(NDVI ~ alpha + gamma * (1 - exp(- age/lambda)),
                  data = dat, start = par, trace = F, 
                  control = nls.control(maxiter = 500))

# plot
plot.NDVI(dat = dat, fit = fit_negexp)
```

<img src="figure/04-Master-Rmdunnamed-chunk-2-1.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

<br>

And let's look at the model summary with parameter estimates


```r
# print model summary
summary(fit_negexp)
```

```
## 
## Formula: NDVI ~ alpha + gamma * (1 - exp(-age/lambda))
## 
## Parameters:
##        Estimate Std. Error t value Pr(>|t|)    
## alpha   0.25107    0.02887   8.695 1.04e-14 ***
## gamma   0.32371    0.02723  11.887  < 2e-16 ***
## lambda  1.17687    0.21396   5.500 1.84e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.07302 on 135 degrees of freedom
## 
## Number of iterations to convergence: 12 
## Achieved convergence tolerance: 3.924e-06
```

<br>

Now the full model:


```r
## Full model

# set parameters
par <- c(alpha = 0.2, gamma = 0.4, lambda = 0.5, A = 0.6, phi = 0)

# fit model
fit_negexpS <- nls(NDVI ~ alpha + gamma * (1 - exp(- age/lambda))
                   + A*sin(2*pi*age + (phi + pi/6*(3 - 1))), 
                   data = dat, start = par, trace = F, 
                   control = nls.control(maxiter = 500))

# plot
plot.NDVI(dat = dat, fit = fit_negexpS)
```

<img src="figure/04-Master-Rmdunnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />



```r
# print model summary
summary(fit_negexpS)
```

```
## 
## Formula: NDVI ~ alpha + gamma * (1 - exp(-age/lambda)) + A * sin(2 * pi * 
##     age + (phi + pi/6 * (3 - 1)))
## 
## Parameters:
##         Estimate Std. Error t value Pr(>|t|)    
## alpha   0.207522   0.024948   8.318 9.31e-14 ***
## gamma   0.364746   0.023926  15.245  < 2e-16 ***
## lambda  0.989154   0.126064   7.846 1.25e-12 ***
## A       0.063136   0.007114   8.875 4.12e-15 ***
## phi    -0.839167   0.111887  -7.500 8.10e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.05835 on 133 degrees of freedom
## 
## Number of iterations to convergence: 15 
## Achieved convergence tolerance: 6.939e-06
```

<br>

Lots more parameters...

Q2: How do the estimates for the common parameters compare?

>*Answer 2:*

<br>

### 3. Compare NLS models using ANOVA

Modelers often want to know which of a set of models are better. One way to do this when comparing nested* models using least squares is using analysis of variance (ANOVA). In this case the `anova()` function will take the model objects as arguments, and return an ANOVA testing whether the full model results in a significant reduction in the residual sum of squares (and thus is better at capturing the data), returning an F-statistic, Degrees of Freedom (the difference in the number of parameters between the models) and p-value.

*i.e. one model is a subset of the other, as in our case


```r
anova(fit_negexp, fit_negexpS)
```

```
## Analysis of Variance Table
## 
## Model 1: NDVI ~ alpha + gamma * (1 - exp(-age/lambda))
## Model 2: NDVI ~ alpha + gamma * (1 - exp(-age/lambda)) + A * sin(2 * pi * age + (phi + pi/6 * (3 - 1)))
##   Res.Df Res.Sum Sq Df  Sum Sq F value   Pr(>F)    
## 1    135    0.71976                                
## 2    133    0.45280  2 0.26696  39.207 4.12e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

<br>

Q3: Which model is better?

>*Answer 3:*

Q4: How many degrees of freedom are there in this ANOVA and why (i.e. what are they)?

>*Answer 4:*

<br>

### 4. Fit models using Maximum Likelihood Estimation (MLE)

First let's fit the simpler model:


```r
## Fit the simpler model using MLE

# set parameters
par <- c(alpha = 0.2, gamma = 0.4, lambda = 0.5)

# fit model
fit_negexpMLE <- fit.negexp.MLE(dat, par)

# plot
plot.NDVI(dat)
# add curve with MLE parameters
lines(dat$age, pred.negexp(fit_negexpMLE$par,dat$age), col = 'skyblue', lwd = 3)
```

<img src="figure/04-Master-Rmdunnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" style="display: block; margin: auto;" />



```r
fit_negexpMLE
```

```
## $par
##     alpha     gamma    lambda 
## 0.2510442 0.3237419 1.1767370 
## 
## $value
## [1] 359053.6
## 
## $counts
## function gradient 
##      118       NA 
## 
## $convergence
## [1] 0
## 
## $message
## NULL
```

<br>

Then the full model:


```r
## Fit the full model using MLE

# set parameters
par <- c(alpha = 0.2, gamma = 0.4, lambda = 0.5, A = 0.6, phi = 0)

# fit model
fit_negexpMLES <- fit.negexpS.MLE(dat, par)

# plot
plot.NDVI(dat)
# add curve with MLE parameters
lines(dat$age, pred.negexpS(fit_negexpMLES$par,dat$age), col = 'skyblue', lwd = 3)
```

<img src="figure/04-Master-Rmdunnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" style="display: block; margin: auto;" />


```r
fit_negexpMLES
```

```
## $par
##       alpha       gamma      lambda           A         phi 
##  0.20772317  0.36449293  0.98919689  0.06310554 -0.83741663 
## 
## $value
## [1] 225574.7
## 
## $counts
## function gradient 
##      914       NA 
## 
## $convergence
## [1] 0
## 
## $message
## NULL
```

<br>

### 5. Compare MLE models using Akaike's information criterion (AIC)

Note that we can't compare our MLE models using ANOVA because our custom functions do not return full model fits like the `nls()` function - only the parameter estimates, negative log-likelihoods and a few other diagnostics.

Another way to compare models (and probably the most common) is using the Akaike information criterion (AIC), which is an estimator of prediction error (i.e. relative quality) of statistical models for a given set of data. 

The formula for the Akaike information criterion is:

$AIC = 2K -2(ln(L))$

Where:

- $k$ = the number of estimated parameters in the model
- $L$ = maximum value of the likelihood function for the model

Since we have our negative log likelihoods (i.e. $-ln(L)$ in the formula above), we can calculate the AICs and compare them.


```r
AIC_simple = 6 + 2*fit_negexpMLE$value

AIC_simple
```

```
## [1] 718113.1
```

```r
AIC_full = 6 + 2*fit_negexpMLES$value

AIC_full
```

```
## [1] 451155.3
```

<br>

When comparing models, the lower the AIC the better, and in general a difference in AIC of 3 or more is analagous to the models being significantly different at an $\alpha$ of $p < 0.05$.


```r
AIC_simple - AIC_full
```

```
## [1] 266957.8
```

<br>

Q5: Is there a preferred model and if so, which one?

>*Answer 5:*

<br>

The nice thing about AIC is that the models you compare do not have to be nested like they do for ANOVA, as long as the data are the same. There are a few other constraints however... 

Here are the AIC scores for our pair of NLS models:


```r
AIC(fit_negexp, fit_negexpS)
```

```
##             df       AIC
## fit_negexp   4 -325.7135
## fit_negexpS  6 -385.6718
```

<br>

You'll notice that these are completely different to the AICs for the MLE models...

Q6: Why is it not okay to compare the AIC of these NLS models with the AIC of the MLE models? Hint: type `?AIC` into the R console and do some reading.

>*Answer 6:*

<br>
```

The R session information (including the OS info, R version and all
packages used):


```r
sessionInfo()
```

```
## R version 4.2.1 (2022-06-23)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_ZA.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_ZA.UTF-8        LC_COLLATE=en_ZA.UTF-8    
##  [5] LC_MONETARY=en_ZA.UTF-8    LC_MESSAGES=en_ZA.UTF-8   
##  [7] LC_PAPER=en_ZA.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_ZA.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
## [1] compiler_4.2.1 magrittr_2.0.3 tools_4.2.1    stringi_1.7.6 
## [5] highr_0.9      knitr_1.39     stringr_1.4.0  xfun_0.30     
## [9] evaluate_0.15
```

```r
Sys.time()
```

```
## [1] "2022-07-29 22:56:32 SAST"
```

