**corit** is a package which contains auxiliary functions for estimating time scale dependent correlations of irregularly sampled time series. The functions *DirectFiltering*, *IntegrandInterpolationMethod* and *InterpolationMethod* resample the irregular data to an equidistant spacing before or during filtering the time series. The Pearson correlation of the time series can be estimated using the function *CorIrregTimser* which allows the filtered time series as an optional output. To test the significance of the estimate, the function *CorQuantilesNullHyp* provides quantiles of correlations obtained from surrogate records, generated with the function *GenNullHypPair* as power-law time series with a spectral slope, estimated from the time series by using the function *estimateTimserSlopes*.   
The package belongs to the manuscript *Comparing methods for analysing time scale dependent correlations in irregularly sampled time series data* which was submitted to Computers & Geosciences. Please contact Maria Reschke (mreschke@awi.de) at the Alfred-Wegener-Institute, Helmholtz Centre for Polar and Marine Research, Germany, for more information.
 
**corit** can be installed directly from bitbucket:
```
if (!require("devtools")) {   
  install.packages("devtools")   
}   
devtools::install_bitbucket("ecus/corit")
``` 
After installation load the package:   
```
library("corit")
``` 
Using **corit** requires the package **zoo** (available on **CRAN**) which is automatically installed when installing **corit**, if it is not yet present. Note: There might be some conflicts when using functions of other packages. In this case, please specific the functions when using.


Code examples to illustrate the application of the R package ‘corit’

```
#(1) estimating the correlation of a pair of time series using linear interpolation and Gaussian filtering

#Assumes two time series with observations and observation times in a vector used to create a zoo-object.  
library(corit)
time.series1 <- zoo(observations1, order.by = observation.times1)	#create a zoo-object
time.series2 <- zoo(observations2, order.by = observation.times2)
Cor <- CorIrregTimser(
timser1 = time.series1,
timser2 = time.series2,
detr = FALSE,	#remove linear trend time series
method = "InterpolationMethod",
appliedFilter = "gauss",
fc = 1/200,	#cut-off frequency
dt = 10,	#time step for the interpolation
int.method = "linear",	#kind of interpolation
filt.output = FALSE)	#return filtered time series 
```

```
#(2) applying a significance test for the correlation estimate based on the correlation of independent noise

time.series1 <- zoo(observations1, order.by = observation.times1)
time.series2 <- zoo(observations2, order.by = observation.times2)
slopes <- estimateTimserSlopes(	#estimate spectral slopes of the time series
timeseries1 = time.series1,
timeseries2 = time.series2,
int.step = 1)	#time step of the interpolated time series
Quant <- CorQuantilesNullHyp(	#quantiles estimated based on surrogate correlations
timser1 = time.series1,
timser2 = time.series2,
beta.noise1 = slopes$s1,
beta.noise2 = slopes$s2,
detr = FALSE,
rep = 1000,	#repetition during Monte Carlo procedure
quant = c(0.05, 0.95),	#quantiles to be estimated
method = "InterpolationMethod",
appliedFilter = "gauss",
fc = 1/200,
dt = 10,
int.method = "linear")
```

```
#(3) code of the Monte Carlo procedure for an applied linear interpolation and Gaussian filtering

CorQuantilesNullHyp(timser1, timser2, beta.noise1, beta.noise2, detr, rep, 
quant, method, appliedFilter, fc, dt, int.method) 
{
n <- max(c(max(index(timser1)), max(index(timser2))))
time1 <- index(timser1)
time2 <- index(timser2)
corNullHypPair <- QuantCorPair <- list()
tmpCorNullHypPair <- numeric(length = rep)
for (p in 1:length(fc)) {
for (i in 1:rep) {
#generate surrogate time series and estimate the correlation
tmp <- GenNullHypPair(n, beta.noise1, beta.noise2, time1, time2, detr)
if (method == "InterpolationMethod") {
tmpCorNullHypPair[i] <- cor.test(InterpolationMethod(tmp$y1, 
fc[p], dt, n, int.method, appliedFilter, k), 
InterpolationMethod(tmp$y2, fc[p], dt, n, int.method, 
appliedFilter, k), method = "pearson", alternative = "two.sided", 
na.action = TRUE)$estimate
}
}
corNullHypPair[[p]] <- tmpCorNullHypPair
QuantCorPair[[p]] <- quantile(corNullHypPair[[p]], quant)
}
return(list(corPair = corNullHypPair, Quantile = QuantCorPair))
}
```