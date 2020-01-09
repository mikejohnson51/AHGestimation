
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AHGestimation

Using pre-processed observation data from the USGS, we can evaluate
learn some things about estimating power law fits from noisy data.

``` r
length(usgs_obs)
#> [1] 3544
```

Taking a random station we subset data occurring after 2010-01-01 and
that falls within the 2 year reoccurance interval (defined via NWM v20
reanalysis product).

``` r
index = 5
  
tmp = usgs_obs[[index]] %>% 
    filter(as.Date(date) > as.Date('2010-01-01')) %>% 
    filter(inChannel == TRUE) %>%
    filter(Ymean > 0) %>%
    filter(V > 0) %>%
    filter(TW > 0) %>%
    filter(Q > 0) %>%
    filter(is.finite(Ymean)) %>% dplyr::select(date,Q, TW, V, Y = Ymean)

dim(tmp)
#> [1] 63  5
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

## Single Relationship fits

Here we use the ahg estimate package to fit the Q-Y relationship using
OLS and NLS models:

``` r
ahg_estimate(Q = tmp$Q, Y = tmp$Y, allowance = .05)
#> NLS performs best for the Q-Y realtionship
#> $Y
#>         exp       coef nrmse   pb method
#> 1 0.5956429 0.06411443   8.8 -1.4    ols
#> 2 0.5068066 0.11060798   7.7  0.8    nls
```

Overall we see that the the NLS model provides a better fit (albeit
small) when measured both by nRMSE and pBais.

## Full Hydraulic fits

When we have data regarding all three hydraulic states (V,TW,Y) we can
ensure that the solutions found are physically valid (meets the
continuity constraint Q = Y*V*TW).

In this mode the OLS and NLS models are fit, and if continuity is not
met, then a Evolutionary approach is implemented. Doing so produces 3
unique fits for 3 variables (27 total combinations). These are crossed
to identify the best performing relationships that meet continuity at a
prescribed allowance:

``` r
ahg_estimate(Q = tmp$Q, V = tmp$V, TW = tmp$TW, Y = tmp$Y, allowance = .05)
#> NLS performs best for the Q-Y realtionship
#> NLS performs best for the Q-TW realtionship
#> NLS performs best for the Q-V realtionship
#> NLS does not meet continuity ...  👎
#> OLS does not meet continuity ...  👎
#> The best performing method (nls) is not physically valid ... 😢
#> Launching Evolutionary Algorithm... allownace (0.05) 💪
#> $Y
#>         exp       coef nrmse   pb method
#> 1 0.5956429 0.06411443   8.8 -1.4    ols
#> 2 0.5068066 0.11060798   7.7  0.8    nls
#> 3 0.2602169 0.52458233  15.1 19.6     GA
#> 
#> $TW
#>          exp     coef nrmse   pb method
#> 1 0.17789845 56.34749  21.7 -5.1    ols
#> 2 0.12695108 77.64203  20.8  0.1    nls
#> 3 0.09504582 93.88977  21.1  2.5     GA
#> 
#> $V
#>         exp       coef nrmse    pb method
#> 1 0.3816120 0.09803420   9.8  -4.1    ols
#> 2 0.4229717 0.08021054   9.2  -0.9    nls
#> 3 0.6326735 0.01930846  13.4 -18.6     GA
#> 
#> $g
#>    Var1 Var2 Var3    c1    c2 viable Verror TWerror Yerror tot_error
#> 1   nls   GA  ols 1.018 0.983   TRUE    9.8    21.1    7.7      38.6
#> 2    GA   GA   GA 0.951 0.988   TRUE   13.4    21.1   15.1      49.6
#> 3   nls  nls  nls 0.689 1.057  FALSE    9.2    20.8    7.7      37.7
#> 4   nls   GA  nls 0.833 1.025  FALSE    9.2    21.1    7.7      38.0
#> 5   nls  nls  ols 0.842 1.015  FALSE    9.8    20.8    7.7      38.3
#> 6   nls  ols  nls 0.500 1.108  FALSE    9.2    21.7    7.7      38.6
#> 7   ols  nls  nls 0.399 1.146  FALSE    9.2    20.8    8.8      38.8
#> 8   ols   GA  nls 0.483 1.114  FALSE    9.2    21.1    8.8      39.1
#> 9   nls  ols  ols 0.611 1.066  FALSE    9.8    21.7    7.7      39.2
#> 10  ols  nls  ols 0.488 1.104  FALSE    9.8    20.8    8.8      39.4
#> 11  ols   GA  ols 0.590 1.072  FALSE    9.8    21.1    8.8      39.7
#> 12  ols  ols  nls 0.290 1.197  FALSE    9.2    21.7    8.8      39.7
#> 13  ols  ols  ols 0.354 1.155  FALSE    9.8    21.7    8.8      40.3
#> 14  nls  nls   GA 0.166 1.266  FALSE   13.4    20.8    7.7      41.9
#> 15  nls   GA   GA 0.201 1.235  FALSE   13.4    21.1    7.7      42.2
#> 16  nls  ols   GA 0.120 1.317  FALSE   13.4    21.7    7.7      42.8
#> 17  ols  nls   GA 0.096 1.355  FALSE   13.4    20.8    8.8      43.0
#> 18  ols   GA   GA 0.116 1.323  FALSE   13.4    21.1    8.8      43.3
#> 19  ols  ols   GA 0.070 1.406  FALSE   13.4    21.7    8.8      43.9
#> 20   GA  nls  nls 3.267 0.810  FALSE    9.2    20.8   15.1      45.1
#> 21   GA   GA  nls 3.951 0.778  FALSE    9.2    21.1   15.1      45.4
#> 22   GA  nls  ols 3.993 0.769  FALSE    9.8    20.8   15.1      45.7
#> 23   GA   GA  ols 4.828 0.737  FALSE    9.8    21.1   15.1      46.0
#> 24   GA  ols  nls 2.371 0.861  FALSE    9.2    21.7   15.1      46.0
#> 25   GA  ols  ols 2.898 0.820  FALSE    9.8    21.7   15.1      46.6
#> 26   GA  nls   GA 0.786 1.020  FALSE   13.4    20.8   15.1      49.3
#> 27   GA  ols   GA 0.571 1.071  FALSE   13.4    21.7   15.1      50.2
#> 
#> $summary
#>   Var1 Var2 Var3    c1    c2 viable Verror TWerror Yerror tot_error
#> 1  nls   GA  ols 1.018 0.983   TRUE    9.8    21.1    7.7      38.6
#> 2  nls  nls  nls 0.689 1.057  FALSE    9.2    20.8    7.7      37.7
#> 3  nls  nls  nls 0.689 1.057  FALSE    9.2    20.8    7.7      37.7
#> 4  ols  ols  ols 0.354 1.155  FALSE    9.8    21.7    8.8      40.3
#>   condition
#> 1  physical
#> 2      best
#> 3       nls
#> 4       ols
```

In the above example we see that NLS was able to provide better fits the
OLS but neither NLS or OLS was able to provide physically valid
solutions (c1, c2, viable). While the GA approach was able to provide a
physically valid solution, its error was almost 10% higher then the
OLS/NLS methods.

However a combined approach of a NLS, OLS, and GA fit was able to
provide a physically valid result with only .9% more error the seen in
the best performing NLS method.

## Take home points of this paper:

> NLS is better then OLS for predicting single relationships in almost
> every case and should become the defacto approach.

> When fitting an entire system, OLS and NLS often provide results that
> are not physically valid.

> GA approaches can always find a physically valid solution but often
> introduce disproportionate error.

> Instead, using a mishmash of NLS, OLS, and GA fits (in cases of
> non-valid solutions), can find solutions with minimal error AND
> physical validity. We argue this approach is the proper way to fit AHG
> relationships.

> This process in formalized in the AHGestimation R packages.
