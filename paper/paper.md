---
title: 'AHGestimation: an R package for hybrid estimation of hydraulic geometry'
tags:
  - hydrology
  - AHG
  - optimization
  - hydraulics
authors:
  - name: J Michael Johnson
    orcid: 0000-0002-5288-8350
    affiliation: 1
  - name: Shahab Afshari
    orcid: 0000-0002-5166-6721
    affiliation: 2
  - name: Arash Modaresi Rad
    orcid: 0000-0002-6030-7923
    affiliation: 1
affiliations:
 - name: Lynker
   index: 1
 - name: UMass Amherst
   index: 2
date: 11 September 2023
bibliography: paper.bib
---

# Summary

### Background 

In the field of hydrology it is common to express the behavior of a river channel at a given cross section using power law equations relating top width (TW), mean depth (Y), and velocity (V) to a given discharge (Q). Collectively these equations define the "at a station hydraulic geometry" (AHG) [@leopold1953hydraulic].

$$
TW = a *Q^b
$$

$$
Y = c*Q^f
$$

$$
V = k*Q^m
$$

A set of AHG relations applies to within-bank flows at a specific cross section and assumes the channel characteristics do not significantly change with discharge. Under these assumptions, two continuity conditions determine if mass is conserved in the system. When either of these is violated, the AHG's equations will cause mass imbalances in computations.

$$
Q = TW·Y·V
$$

and therefore:

$$
\begin{eqnarray} b + f + m &=& \\ a·c·k &=&  \\ &=& 1  \end{eqnarray}
$$

### Statement of need

Large scale models simulating river flow are critical for forecasting water availability, drought, and flood inundation. These models must represent the shape of river channels in some generalized way. 

While hydraulic geometry relationships have been extensively studied, they remain unquantified for the majority of stream reaches across the country. Consequently, large-scale  models frequently use simple approximations that impact the accuracy of streamflow estimates [@hess-26-6121-2022; @johnson2023comprehensive] and flood forecasting [@zheng2018river; @maidment2014national; @johnson2019integrated; @fim]. At NOAA, these relationships have been based on trapezoidal geometries (e.g. [@wrfhydro]) derived from drainage area assumptions (e.g. [@blackburn2017development]). 

Other efforts have aimed to calculate, and synthesize river channel data at a large scale in the United States (e.g. [@enzminger_thomas_l_2023_7868764; @afshari_shahab_2019_2558565]) however each of these relied on traditional OLS fitting methods AND data preprocessing [@afshari2017statistical]. And while both efforts produced valuable data products, the software used is not shared.

This open source package is designed to assist work flows that are challenged by the following characteristics of hydraulic data:

1. Data is often distributed without consistent structure
2. Data is noisy and hard to fit with traditional methods
3. Data is tabular, making the development of efficient regional and continental datasets a challenge.


### Software

Estimating AHG's has predominately occurred on a location-by-location basis with site specific knowledge guiding the data used and the validation of the outputs. As efforts seek to extend to coverage of hydrologic models and flood inundation mapping (FIM), the ability to estimate these relations from disparate and often noisy datasets is increasingly important.

`AHGestimation` is an R package [@r-project] providing three capabilities:

1. Tools to estimate single and full system AHG relations using a robust estimation techniques that enforces flow continuity and minimizes total system error (`ahg_estimate`). This is accomplished by introducing a hybrid approach that supplements the traditional Ordinary Least Squares approach, with a  Nonlinear Least Square regression, Evolutionary algorithm (NSGA-2; [@mco]), and ensemble modeling approach.

2. Methods to filter outliers based on time (`date_filter`), mass conservation (`qva_filter`), and statistical detection (`mad_filter`, `nls_filter`).

3. Formalize, in code, many of the concepts derived in [@dingman2018field] that relate AHG coefficients and exponents to cross-section hydraulics and geometry. These include (`cross_section`, `compute_hydraulic_params`, `compute_n`). 

The package documentation includes several examples on the theory, design, and application of this tool set.

The first stable version of `AHGestimation` was made available in 2019 and applied to an aggregated dataset of USGS manual field measurements. Since then, it has been actively developed to better understand and quantify these fundamental relationships in the face of noisy, large, and disparate data sources. Applications of the software have been used to (1) better deliver actionable flood forecasts from the NOAA National Water Model [@johnson2022knowledge] (2) help the NOAA Office of Water Prediction's efforts to develop a 3D hydrography to improve national hydrologic and flood prediction and to (3) bolster the co-agency National Hydrologic Geospatial Fabric [[@referencefabric; @blodgett2021mainstems; @blodgett2023generating]].

# Example of use

`AHGestimation` is available on
[Gtihub](https://github.com/mikejohnson51/AHGestimation) and can be installed as follows:

``` r
#install.packages(remotes)
remotes::install.packages("mikejohnson51/AHGestimation")
```

This example shows how the package can be used to (1) remove data outliers based on time criteria and an NLS envelope (2) fit AHG parameters using our hybrid modeling approach and (3) estimate and plot the shape of the associated cross-section with an area-depth relation. The script to generate the plot can be found [here](image.R).

``` r 
(data = nwis  |>
  date_filter(10, keep_max = TRUE) |> 
  nls_filter(allowance = .5) )
  
#>  siteID       date           Q         Y        V       TW
#>  01096500 1987-04-07 317.1486864 3.7651765 1.624584 51.81600
#>  01096500 2009-06-10   4.9837651 0.3386667 0.445008 32.91840
#>  01096500 2010-03-25 122.6119475 2.2812075 1.106424 48.46320
#>  01096500 2010-04-06  67.1109274 1.4431347 1.039368 44.80560
#>  01096500 2010-09-16   2.1464170 0.3609009 0.283464 21.03120
#>  01096500 2011-03-10 124.0277899 2.2813252 1.094232 49.68240
#>  01096500 2011-03-10 118.9307574 2.2488293 1.063752 49.98720
#>  01096500 2011-06-01  15.0645626 0.6897725 0.850392 25.72512
#>  01096500 2011-12-01  49.2713138 1.4478000 0.932688 36.57600
  
ahg_fit = ahg_estimate(data)
t(ahg_fit[1,])

#> V_method  "nls"       
#> TW_method "nls"       
#> Y_method  "nls"       
#> viable    "TRUE"      
#> tot_error "0.352759"  
#> V_error   "0.1470922" 
#> TW_error  "0.1161101" 
#> Y_error   "0.08955664"
#> V_coef    "0.2822548" 
#> TW_coef   "18.17896"  
#> Y_coef    "0.1945348" 
#> V_exp     "0.3105107" 
#> TW_exp    "0.1850407" 
#> Y_exp     "0.5087359" 
#> condition "bestValid" 

shape = compute_hydraulic_params(ahg_fit[1,])

#>        r         p        d        R        bd        fd        md
#> 2.749318 0.6103574 5.427385 1.363727 0.1842508 0.5065641 0.3091851

cs = cross_section(r = shape$r,  TW = max(data$TW), Ymax = max(data$Y))

#>  ind         x            Y            A
#>    1  0.000000 3.5656665613 1.318953e+02
#>    2  1.040816 3.1871221498 1.130112e+02
#>    3  2.081633 2.8351885546 9.618829e+01
#>    4  3.122449 2.5090304653 8.127753e+01
#>    5  4.163265 2.2078034104 6.813466e+01
#>    6  5.204082 1.9306532333 5.662032e+01
#>    7  6.244898 1.6767155105 4.660008e+01
#>    8  7.285714 1.4451149046 3.794447e+01
#>    9  8.326531 1.2349644401 3.052910e+01
#>   10  9.367347 1.0453646866 2.423471e+01
```

![Faceted image with multiple views of the channel estimate.\label{fig:ahg-1}](paper-fig1.png)

As a proof of concept, this approach was applied to the synthetic rating curves in FIM4 data set ([@fim]). The package allowed the size of their rating curve database to be reduced by 99.68% and maintained average accuracy within 0.4% nRMSE or the source data. Further, the reduction to a consistent AHG formulation allows them to be more interoperable with efforts like [@enzminger_thomas_l_2023_7868764], or [@afshari_shahab_2019_2558565] while also beginning to provide key training data for more advanced prediction methods that seek to estimate the shape on non-measured river segments.

# Acknowledgements

The development of this package began in 2017 following the NOAA OWP Summer Institute and clear evidence channel shape may be a limiting factor in National Water Model Performance. 

The algorithm and implementation began as a graduate school project between friends at UC Santa Barbara and UMass Amherst and has since evolved to provide an open source utility for robust large scale data synthesis and evaluation. Funding from the National Science Foundation (Grants 1937099, 2033607) provided time to draft [@preprint] and apply an early version of this software to the Continental Flood Inundation Mapping synthetic rating curve dataset [@cfim]. Funding from the National Oceanic and Atmospheric Administration's Office of Water Prediction supported the addition of data filtering and hydraulic estimation, improved documentation, and code hardening. We are grateful to all involved.

# References
