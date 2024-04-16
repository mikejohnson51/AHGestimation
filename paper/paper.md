---
title: 'AHGestimation: An R package for computing robust, mass preserving hydraulic geometries and rating curves'
tags:
- hydrology
- AHG
- optimization
- hydraulics
date: "12 April 2024"
output:
  pdf_document: default
  html_document:
    df_print: paged
authors:
- name: J Michael Johnson
  orcid: "0000-0002-5288-8350"
  affiliation: 1
- name: Shahab Afshari
  orcid: "0000-0002-5166-6721"
  affiliation: 2
- name: Arash Modaresi Rad
  orcid: "0000-0002-6030-7923"
  affiliation: 1
bibliography: paper.bib
affiliations:
- name: Lynker, NOAA/NWS Office of Water Prediction
  index: 1
- name: University of Massachusetts Amherst
  index: 2
---

# Summary

### Background 

The behavior of a river channel at a given cross-section is often described using power law equations to relate top width (TW), mean depth (Y), and mean velocity (V) to discharge (Q). These equations collectively define the "at a station hydraulic geometry" (AHG), with coefficients (a, c, k) and exponents (b, f, m) for the Q-TW, Q-Y, and Q-V relationships [@leopold1953hydraulic]. While each relation can be studied independently, together they represent the full hydraulic system.

\begin{eqnarray} TW = a \cdot Q^b \\ Y = c \cdot Q^f \\ V = k \cdot Q^m \end{eqnarray}

Where: 

\begin{eqnarray} Q = TW \cdot Y \cdot V \end{eqnarray}

When using the traditional AHG fitting methods, power laws are estimated using Ordinary Least Squares Regression (OLS) on the logarithmic transformation of the variables, where the exponential of the intercept provides the power law coefficient, and the slope of the linear model provides the exponent. Several studies have highlighted the potential problem of using this approach given that when the log transformed values are back-transformed, the estimates are effectively medians, and not means of the estimates, resulting in a general low bias. AHG relations apply to within-bank flows at a specific location and assume the channel characteristics do not significantly change with discharge. Under these assumptions, two continuity conditions ensure mass conservation (equation 5-6). 

\begin{eqnarray} b + f + m = 1 \\ a \times c \times k = 1  \end{eqnarray}

Real-world data often fails to meet these conditions precisely, so an allowance is typically permitted. Violating these conditions leads to mass imbalances in computations, causing input streamflows to gain or lose mass in its translation to other states. Given the noisy nature of hydraulic data, tools are needed to prep input data and fit equations while balancing accuracy and mass preservation for hydrologic and hydraulic modeling.


### Statement of need

Hydrologic models that simulate streamflow are critical for forecasting water availability, drought, and flood inundation. A key aspect of these models is estimating the size and shape of channels, often achieved through hydraulic geometry relationships.

While extensively studied at local scales, these relationships remain unquantified for the majority of stream reaches globally, including in the United States. As a result, large-scale models often rely on incomplete approximations, leading to less accurate streamflow estimates [@hess-26-6121-2022; @johnson2023comprehensive; @fang2024] and flood forecasts [@zheng2018river; @maidment2014national; @johnson2019integrated; @fim]. For instance, the National Oceanic and Atmospheric Administration National Water Model [@cosgrove2023] uses trapezoidal geometries [@wrfhydro] that are in part derived from hydraulic geometry relationships and drainage area assumptions found in @bieger2015development, @bieger2016development, @blackburn2017development. Although these approximations are recognized as overly simplistic, advancing them requires integrating diverse observation systems and refining them into parameterized, mass-conserving relationships.

Several efforts have aimed to address this challenge in the United States, primarily relying on traditional OLS fitting methods and data preprocessing [@enzminger_thomas_l_2023_7868764; @afshari_shahab_2019_2558565; @afshari2017statistical]. However, these efforts are limited by the lack of shared software and source data, hindering the evolution, and interoperability, of their products.

### Software

AHG estimation has traditionally been conducted on a site-specific basis, with location-specific knowledge guiding data selection and output validation. However, as interest in large-scale model applications grows [@archfield2015], the importance of estimating these relationships from diverse and often noisy datasets becomes increasingly evident.

Towards this, `AHGestimation` is an R package [@r-project] providing three capabilities:

1. **AHG Estimation** (`ahg_estimate`): This tool provides robust estimation techniques to estimate single and full system AHG relations that enforce flow continuity and minimizes total system error. This is accomplished by introducing a hybrid approach that supplements the traditional OLS approach, with a  Nonlinear Least Square (NLS) regression, Evolutionary algorithm (NSGA-2) [@mco], and ensemble modeling approach.

2. **Outlier Filtering Methods** (`date_filter`, `qva_filter`, `mad_filter`, `nls_filter`): These methods allow users to filter outliers based on various criteria, including time, mass conservation, and statistical detection.

3. **Formalization of derived AHG Concepts** (`cross_section`, `compute_hydraulic_params`, `compute_n`): These functions formalize many of the concepts derived in @dingman2018field that relate AHG coefficients and exponents to cross-section hydraulics and geometry. 

The package documentation includes several examples on the theory, design, and application of these tools.

The first stable version of `AHGestimation` was made available in 2019 and was applied to an aggregated dataset of USGS manual field measurements. Since then, it has been actively developed to better understand and quantify these fundamental relationships in the face of noisy, large, and disparate data sources. Applications of the software have been used to (1) demonstrate how improved flood forecasts could be delivered from the NOAA/NWS National Water Model [@johnson2022knowledge], (2) help the NOAA/NWS Office of Water Prediction develop continental scale channel size and shape estimates to improve flood prediction and hydraulic routing, and (3) bolster the co-agency USGS/NOAA National Hydrologic Geospatial Fabric and Next Generation Water Resource Modeling Framework Hydrofabric efforts [@referencefabric; @blodgett2021mainstems; @blodgett2023generating; @nextgenhf].

# Example of use

`AHGestimation` is available on
[GitHub](https://github.com/mikejohnson51/AHGestimation) and can be installed as follows:

``` r
#install.packages(remotes)
remotes::install.packages("mikejohnson51/AHGestimation")
```

This example illustrates how the package can be utilized to:

1. **Remove Data Outliers**: Outliers are filtered based on time criteria and an NLS envelope.
2. **Fit AHG Parameters**: AHG parameters are estimated using the hybrid modeling approach.
3. **Estimate and Plot Cross-Section Shape**: The shape of the associated cross-section is estimated and plotted with an area-depth relation.

The script to generate the plot found in \autoref{fig:ahg-1} can be found [here](https://github.com/mikejohnson51/AHGestimation/blob/master/paper/image.R), and the `nwis` data object is exported with the package to provides field measurements taken at [USGS site 01096500 on the Nashua River at East Pepperell in Massachusetts](https://waterdata.usgs.gov/nwis/measurements/?site_no=01096500&agency_cd=USGS). 

``` r 
nwis

#> siteID       date           Q         Y        V       TW
#> 01096500 1984-11-14   9.7409954 0.5276645 0.652272 28.34640
#> 01096500 1985-01-04  11.8930757 0.6263473 0.682752 27.73680
#> 01096500 1985-05-06  10.8170356 0.5952226 0.563880 32.30880


nrow(nwis)
#>  245

# Keep only those observation made in the most recent 10 years, 
# and that fall withing the .5 nls envelope
(data = nwis  |>
  dplyr::rename(Q = Q_cms, Y  = Y_m, V= V_ms, TW = TW_m) |>
  date_filter(10, keep_max = TRUE) |> 
  nls_filter(allowance = 0.5) )

# data reduced to 80 observations based on filters
nrow(data)
#>  85
```

The reduced clean data can then be used to fit an AHG relation and compute a set of hydraulic parameters:

``` r
# Fit AHG relations
ahg_fit = ahg_estimate(data)
t(ahg_fit[1,])

#> V_method  "nls"      
#> TW_method "nls"      
#> Y_method  "nls"      
#> c1        "1.006"    
#> c2        "1.001"    
#> viable    "TRUE"     
#> tot_nrmse "0.3234122"
#> V_nrmse   "0.1337535"
#> TW_nrmse  "0.1009869"
#> Y_nrmse   "0.0886718"
#> V_coef    "0.2905399"
#> TW_coef   "18.23401" 
#> Y_coef    "0.1898499"
#> V_exp     "0.3101059"
#> TW_exp    "0.1756746"
#> Y_exp     "0.5155772"
#> condition "bestValid"

# Use the AHG relations to compute hydraulic parameters
shape = compute_hydraulic_params(ahg_fit[1,])

#> r    p    d    R    bd   fd   md
#> 2.93 0.60 5.70 1.34 0.18 0.51 0.31
```

Finally, the max width and depth, paired with the derived `r` coefficient can be used to generate a cross section:

``` r
# Use the max width, max depth, and derived `r` to generate a cross section
# x: is the relative distance from the left bank
# Y: is the associated depth
# A: is the area associated with depth Y

cs = cross_section(r = shape$r,  
                   TW = max(data$TW), 
                   Ymax = max(data$Y))

#>  ind         x            Y            A
#>    1  0.000000 3.5525892702 1.337069e+02
#>    2  1.040816 3.1514749668 1.136999e+02
#>    3  2.081633 2.7814289966 9.601412e+01
```

![Faceted image with multiple views of the channel estimate.\label{fig:ahg-1}](paper-fig1.png)

As a proof of concept, this approach was applied to the synthetic rating curves generated by NOAA/NWS OWP's inundation mapping software [@fim]. The package allowed the size of their rating curve database to be reduced by 99.68% while maintaining average accuracy within 0.4% NRMSE of the source data [@preprint]. This reduction to a consistent AHG formulation not only enhances interoperability with efforts such as those by @enzminger_thomas_l_2023_7868764 or @afshari_shahab_2019_2558565, but also lays the groundwork for providing essential training data for advanced prediction methods that seek to estimate the shape on non-measured river segments.

# Acknowledgements

The development of this package began following the 2017 NOAA/NWS Office of Water Prediction Summer Institute [@owpsi2017] where discussions highlighted the potential influence of channel shape representation on the performance of the National Water Model. 

The algorithm and implementation began as a graduate school project between friends at UC Santa Barbara and UMass Amherst and has since evolved to provide an open source utility for robust large scale data synthesis and evaluation. Funding from the National Science Foundation (Grants 1937099, 2033607) provided time to draft [@preprint] and apply an early version of this software to the Continental Flood Inundation Mapping synthetic rating curve dataset [@cfim]. Funding from the NOAA/NWS OWP supported the addition of data filtering and hydraulic geometry estimation, improved documentation, and code hardening. We are grateful to all involved.

# Disclaimer

The views expressed in this article do not necessarily represent the views of NOAA, the USGS, or the United States.

# References
