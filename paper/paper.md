---
title: 'AHGestimation: an R package for hybrid estimation of hydraulic geometry'
tags:
- hydrology
- AHG
- optimization
- hydraulics
date: "11 September 2023"
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
- name: Lynker
  index: 1
- name: UMass Amherst
  index: 2
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

AHG relations apply to within-bank flows at a specific cross section and assume the channel characteristics do not significantly change with discharge. Under these assumptions, two continuity conditions determine mass conservation. When either of these is violated, the AHG's equations will cause mass imbalances in computations.

$$
Q = TW·Y·V
$$

and therefore:

\begin{eqnarray} b + f + m &=& \\ a·c·k &=&  \\ &=& 1 \end{eqnarray}

### Statement of need

Large scale models simulating river flow are critical for forecasting water availability, drought, and flood inundation. These models must represent the size and shape of river channels in some generalized way. 

While hydraulic geometry relationships have been extensively studied, they remain unquantified for the majority of stream reaches across the country. Consequently, large-scale models frequently use generalized approximations that impact the accuracy of streamflow estimates [@hess-26-6121-2022; @johnson2023comprehensive] and flood forecasting [@zheng2018river; @maidment2014national; @johnson2019integrated; @fim]. At NOAA, National Water Model characteristics are based on trapezoidal geometries (e.g. [@wrfhydro]) that are in part derived from hydraulic geometry relationships and drainage area assumptions (e.g. @blackburn2017development). 

Other efforts have aimed to calculate, and synthesize river channel data at a large scale in the United States (e.g. [@enzminger_thomas_l_2023_7868764; @afshari_shahab_2019_2558565]) however each of these relied on traditional Ordinary Least Squares Regression (OLS) fitting methods _and_ data preprocessing [@afshari2017statistical]. And while both efforts produced valuable data products, the software used is either not shared difficult to use.

This open source package is designed to assist work flows that are challenged by the following characteristics of hydraulic data:

1. Data is often distributed without consistent structure
2. Data is noisy and hard to fit with traditional methods
3. Data is tabular, making the development of efficient regional and continental datasets a challenge.


### Software

AHG estimates have predominately occurred on a location-by-location basis with site specific knowledge guiding the data used and the validation of the outputs. As interest in large scale model applications increases, the need to estimate these relations from disparate and often noisy datasets is increasingly important.

`AHGestimation` is an R package [@r-project] providing three capabilities:

1. Tools to estimate single and full system AHG relations using a robust estimation techniques that enforces flow continuity and minimizes total system error (`ahg_estimate`). This is accomplished by introducing a hybrid approach that supplements the traditional OLS approach, with a  Nonlinear Least Square (NLS) regression, Evolutionary algorithm (NSGA-2; @mco), and ensemble modeling approach.

2. Methods to filter outliers based on time (`date_filter`), mass conservation (`qva_filter`), and statistical detection (`mad_filter`, `nls_filter`).

3. Formalize, in code, many of the concepts derived in [@dingman2018field] that relate AHG coefficients and exponents to cross-section hydraulics and geometry. These include (`cross_section`, `compute_hydraulic_params`, `compute_n`). 

The package documentation includes several examples on the theory, design, and application of this tool set.

The first stable version of `AHGestimation` was made available in 2019 and was applied to an aggregated dataset of USGS manual field measurements. Since then, it has been actively developed to better understand and quantify these fundamental relationships in the face of noisy, large, and disparate data sources. Applications of the software have been used to (1) demonstrate how improved flood forecasts could be delivered from the NWS/NOAA National Water Model [@johnson2022knowledge] (2)  help the NOAA/NWS Office of Water Prediction develop continental scale channel size and shape estimates to improve flood prediction and hydraulic routing and to (3) bolster the co-agency sponsored National Hydrologic Geospatial Fabric [@referencefabric; @blodgett2021mainstems; @blodgett2023generating].

# Example of use

`AHGestimation` is available on
[Gtihub](https://github.com/mikejohnson51/AHGestimation) and can be installed as follows:

``` r
#install.packages(remotes)
remotes::install.packages("mikejohnson51/AHGestimation")
```

This example shows how the package can be used to (1) remove data outliers based on time criteria and an NLS envelope (2) fit AHG parameters using a hybrid modeling approach and (3) estimate and plot the shape of the associated cross-section with an area-depth relation. The script to generate the plot can be found [here](image.R), and the `nwis` data object is exported with the package to provide sample data and contains the field measurements taken at [USGS site 01096500 on the Nashua River at East Pepperell in Massachusetts](https://waterdata.usgs.gov/nwis/measurements/?site_no=01096500&agency_cd=USGS). 

``` r 
nwis

#> siteID       date           Q         Y        V       TW
#> 01096500 1984-11-14   9.7409954 0.5276645 0.652272 28.34640
#> 01096500 1985-01-04  11.8930757 0.6263473 0.682752 27.73680
#> 01096500 1985-05-06  10.8170356 0.5952226 0.563880 32.30880
#> 01096500 1985-06-26   1.9453674 0.2400300 0.332232 24.38400
#> 01096500 1986-01-09  11.3550556 0.5920154 0.606552 31.69920
#> 01096500 1986-02-27  21.4358532 0.8574593 0.902208 27.73680
#> 01096500 1986-05-22   6.9093107 0.4806462 0.454152 31.69920
#> 01096500 1986-07-07   8.9198068 0.5457371 0.512064 32.00400

nrow(nwis)
#>  245

# Keep only those observation made in the last 10 year, 
# and that fall withing the .5 nls allowance
(data = nwis  |>
  date_filter(10, keep_max = TRUE) |> 
  nls_filter(allowance = .5) )

# data reduced to 80 observations based on filters
nrow(data)
#>  80
  
# Fit AHG relations
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

# Use the AHG relations to compute hydraulic parameters
shape = compute_hydraulic_params(ahg_fit[1,])

#>        r         p        d        R        bd        fd        md
#> 2.749318 0.6103574 5.427385 1.363727 0.1842508 0.5065641 0.3091851

# Use the max TW, max depth, and derived `r` to generate a cross section
# x is the relative distance from the left bank
# Y is the assocated depth
# A is the area assoicated with depth Y

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

As a proof of concept, this approach was applied to the synthetic rating curves generated by NWS/NOAA OWP's inundation mapping software ([@fim]). The package allowed the size of their rating curve database to be reduced by 99.68% and maintained average accuracy within 0.4% nRMSE or the source data [@preprint]. Further, the reduction to a consistent AHG formulation allows them to be more interoperable with efforts like [@enzminger_thomas_l_2023_7868764], or [@afshari_shahab_2019_2558565] while also beginning to provide key training data for more advanced prediction methods that seek to estimate the shape on non-measured river segments.

# Acknowledgements

The development of this package began in 2017 following the NWS/NOAA OWP Summer Institute and clear evidence channel shape may be a limiting factor in National Water Model performance. 

The algorithm and implementation began as a graduate school project between friends at UC Santa Barbara and UMass Amherst and has since evolved to provide an open source utility for robust large scale data synthesis and evaluation. Funding from the National Science Foundation (Grants 1937099, 2033607) provided time to draft [@preprint] and apply an early version of this software to the Continental Flood Inundation Mapping synthetic rating curve dataset [@cfim]. Funding from the NWS/NOAA OWP supported the addition of data filtering and hydraulic geometry estimation, improved documentation, and code hardening. We are grateful to all involved.

# References
