**afetsdc** is a package which contains auxiliary functions for estimating time scale dependent correlations of irregularly sampled time series. The functions *DirectFiltering*, *IntegrandInterpolationMethod* and *InterpolationMethod* resample the irregular data to an equidistant spacing and filter the time series. The correlation itself can be estimated using standard routines working on regular data.   
The package belongs to the paper draft ’Comparison of methods for analyzing time scale dependent correlations of irregularly sampled time series’.
 
**afetsdc** can be installed directly from bitbucket:
```
if (!require("devtools")) {   
  install.packages("devtools")   
}   
devtools::install_bitbucket("ecus/afetsdc")
``` 
After installation load the package:   
```
library("afetsdc")
``` 
Using **afetsdc** requires the package **zoo** (available on **CRAN**) which is automatically installed when installing **afetsdc**, if it is not yet present.
