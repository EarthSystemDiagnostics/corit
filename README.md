**corit** is a package which contains auxiliary functions for estimating time scale dependent correlations of irregularly sampled time series. The functions *DirectFiltering*, *IntegrandInterpolationMethod* and *InterpolationMethod* resample the irregular data to an equidistant spacing before or during filtering the time series. The Pearson correlation of the time series can be estimated using the function *CorIrregTimser* which allows the filtered time series as an optional output. To test the significance of the estimate, the function *CorQuantilesNullHyp* provides quantiles of correlations obtained from surrogate records, generated with the function *GenNullHypPair* as power-law time series with a spectral slope, estimated from the time series by using the function *estimateTimserSlopes*.   
The package belongs to the manuscript ’Comparing methods for analysing time scale dependent correlations in irregularly sampled time series data’ which was submitted to Nonlinear Processes in Geophysics. Please contact Maria Reschke (mreschke@awi.de) at the Alfred-Wegener-Institute, Helmholtz Centre for Polar and Marine Research, Germany, for more information.
 
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
