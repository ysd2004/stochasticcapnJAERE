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

## Pindyck (1984) Example
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

## Non-CBM stochasticity
This example extends the GBM volatility to the demographic and environmental stochasticity in Sims, Horan, and Meadows (2018), which is Fgiure5 and Figure6 in the manuscript.

* GBM: $\sigma_{0}^{2} s^{2}$

* demographic and environmental: $\tilde{r}s \left( \frac{\tilde{M}+\tilde{C}}{\tilde{M}-\tilde{C}}  +(1-2 \mu) \frac{s}{\tilde{K}}\right) + \sigma_{M}^{2} \tilde{R}^2 s^{2} \left[ 1 - \frac{s}{K_{1}} \right]^2$

As shown in Table C.1, the economic parameters are the same in Figure 2 and Figure 3. The ecological parameters are:

* $R = 1$
* $\tilde{M} = 0.53$
* $\tilde{C} = 0.03$
* $\tilde{r} = 0.5$
* $\tilde{K} = 100$
* $K_{1} = 114.3$
* $\mu = 1$

To download the replication R-code,

``` r
    ## In R
    download.file("https://raw.github.com/ysd2004/stochasticcapnJAERE/main/replicationRcode/Figure5_and_Figure6.R", "Figure5_and_Figure6.R")
```

Example 2: A stock with non-convex drift and smooth stochasticity
====================================
This example is an extension of Fenichel and Abbott (2014) with the cubic growth function. In the deterministic GOM example, the followings are defined.

* Catch function: $h(s,x) = qsx^{\alpha}=qy^{\alpha}s^{\gamma \alpha +1}$ where $x(s)=ys^{\gamma}$
* Net benefit: $w(s)=pqy^{\alpha}s^{\gamma \alpha +1} - cys^{\gamma}$
* Logistic growth function: $f(s)=rs \left( 1-\frac{s}{K} \right)$

Supposing all other functions are the same, the stock evolves according to the cubic growth with geometric Brownian motion as:

$$
ds = \left(rs(t)\left(\frac{s(t)}{K_{1}} -1 \right)\left(1-\frac{s(t)}{K_2}\right)- h\left(x\left(s\left(t\right)\right),s\left(t\right)\right)     \right)dt+ \sigma s(t)dZ(t) \text{,}
$$

In the example, we used the following parameters in Fenichel and Abbott (2014).

* $r = 0.3847$
* $q = 0.0003172934$
* $p = 2.7$
* $c=153$
* $\alpha = 0.5436459$
* $\gamma = 0.7882$
* $y = 0.157455$
* $\delta = 0.02$
* $K = 359016000$
* $K_{1} = 0.13 K_{2} = 46672080$

To replicate Figure8, download the R-code with:

``` r
    ## In R
    download.file("https://raw.github.com/ysd2004/stochasticcapnJAERE/main/replicationRcode/Figure8.R", "Figure8.R")
```

In the code, there is the line to download the data of all brute force forward simulations (*simresult.rda*) detailed below. 

* *sto00.rda*: deterministic
* *sto01.rda*: stocastic with $\sigma=0.1$
* *sto02.rda*: stocastic with $\sigma=0.2$
* *sto03.rda*: stocastic with $\sigma=0.3$

