# stochasticcapnJAERE
This repository replicates the examples in Abbot et al. submitted to JAERE. All replications are directly derformable in R.

------------------------------------------------------------------------

0 Load functions and data for replication
====================================
To implement the V-approximation suggested in the manuscript, we provide the pre-defined R-functions. Before replication the following examples, these functions should be loaded at first.

```r
## load pre-defined functions
source('http://dl.dropbox.com/s/m1yo02onssypg6h/functions.R')

## R-packages used in the examples
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(rootSolve)
library(tidyverse)
```

1 Example 1: Optimized and non-optimized renewable resource management under convexity
====================================
This example is Example 1 in Pindyck (1984, pp. 296 - 297.). In the example,
