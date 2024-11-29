# Replicable Examples: Abbott, Fenichel, and Yun

This repository is to share the *R* data and code replicating the main results submitted to JAERE by Abbott, Fenichel, and Yun.

GitHub Repository maintained by: Seong Yun\
Department of Agricultural Economics\
Mississippi State University\
**<seong.yun@msstate.edu>**\
**<https://sites.google.com/site/yunsd2004/>**\
*Last updated: Nov. 25, 2024*

------------------------------------------------------------------------

Manuscript information
====================================

```r
  Abbott, K. Joshua, Eli P. Fenichel, and Seong, D. Yun, "Risky (Natural) Assets: Stochasticity and the Value of Natural Capital"
```

Example 1: Optimized and non-optimized renewable resource management under convexity
====================================

# Pindyck (1984)
This example is Example 1 in Pindyck (1984, pp. 296-297), which is Figure 2 and Figure 3 in the manuscript.

* Isoelastic demand function: $q(p) = b p^{-\eta}$
* Cost function: $c(s)=cs^{-\gamma}$

With $\eta=1/2$ and $\gamma=2$, the stock evolves according to a logistic drift with geometric Brownian motion as:

$$
ds = \left[ rs(1-s/K)-q \right] dt + \sigma s \ dZ
$$

where $s$ is a stock, $r$ is the intrinsic growth rate, $K$ is the carrying capacity, $q$ is a harvest level, $\sigma$ is stock volatility, and $dZ$ is an increment of a Wiener process.

This example has a closed form solution as:
$$
\text{Optimzed co-state: } V_s = \phi / s^2 \qquad \text{and}
$$
$$
\text{Economic program: } x(s) = q^*(s) = \frac{bs}{(\phi+c)^{(1/2)}} \text{,}
$$
where $\phi = \frac{2b^2+2b[b^2+c(r+\delta-\sigma^2)^2]^{1/2}}{(r+\delta-\sigma^2)^2}$.

To download the replication R-code,

``` r
    ## In R
    download.file("https://raw.github.com/ysd2004/stochasticcapnJAERE/main/replicationRcode/Figure2_and_Figure3.R", "Figure2_and_Figure3.R")
```

# Non-CBM stochasticity