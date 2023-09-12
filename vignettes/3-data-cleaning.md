# Data Cleaning

This software also provides an optional capability based on mathematical concepts and literature to filter data outliers and inappropriate measurments. 

The observations of dicharge, width, depth and velocity from gauges although being extremly valuble require carfull consideration, and filttering due to reasons including but not limited to:

1- Changes in measurement instruments
2- Changes in location of measurements (slight deviations form cross sections can result in significant differences in geometry)
3- Changes in river due to flooding and etc. (river meandering)

As such the software comes with the following four filtering mechanisims.
Before, these filterings data are screened for nay negative, inf, and proportionaly large values
The filters are chosen based on logical bounds (in feet) found 
in historical records of bigest rivers in the U.S. (i.e., Hudson, Mississippi).
There are some gagues that are located in backwater zones and are affected by 
the coastal processes including  sea level changes, tides and storm surge 
that can extend hundreds of kilometers inland. These location often show negative
discharge vlaues.


## The **date_filter**

This filter is to select the time period of study. Selecting ranges greater than 10 years is not advised as river meandering, extreme floods, etc. can drastically change shape, location, and hydraulic geometry of a river. To use this filter one can call 

```r
df = usgs_obs[[index]] %>% 
  select(date = date, Y =Ymean, TW, V, Q)

hf_df = df %>% 
  date_filter(10, keep_max = TRUE) %>% 
```

## The **nls_filter**

This filter uses NLS fit to find data that are whithin 0.75 deviation from the fit and discards others as outlairs and then using the filtterd data the hybrid approuch is implemented to find the best fit. The choice for selection of the amount deviation allowed is given to the user.

```r
hf_df = df %>% 
  date_filter(10, keep_max = TRUE) %>% 
  nls_filter(.5)
```

## The **qva_filter**

This filter implements the approuch described in [HyG dataset](https://zenodo.org/record/7868764), where gages with measured Q are not equivilant of the product of measured velocity and measured area by more than 5% are discarded in process. Such gages are potenially near costal regions and are often tidally influenced and may exchibit filipped power law relation.

```r
hf_df = df %>% 
  date_filter(10, keep_max = TRUE) %>% 
  qva_filter()
```

## The **mad_filter**

This filter also implements the approuch described in [HyG dataset](https://zenodo.org/record/7868764), and uses and iterative fitting approuch to filter data in the following form

1- By appling a linear regression values of log-transformed TW, V, and D residuals inside the three median absolute deviation (MAD) range are kept. 
2- In a loop the above is reptead until there are no data outliers.
3- Based on P value of the regression statistical significance is tested and those with p-values >0.05 are excluded.

Then the filtered data can be passed to the solver to find FHG realtions with continuity constrain

```r
hf_df = df %>% 
  date_filter(10, keep_max = TRUE) %>% 
  qva_filter()
```

An example showing performance of each filtering method is shown below 

```r
df = usgs_obs[[index]] %>% 
  select(date = date, Y =Ymean, TW, V, Q)

hf_df = df %>% 
  date_filter(10, keep_max = TRUE) %>% 
  mad_filter() %>% 
  qva_filter() 

nls_df = df %>% 
  date_filter(10, keep_max = TRUE) %>% 
  nls_filter(.5)
  
x1 = fhg_estimate(df = df)
x2 = fhg_estimate(df = hf_df)
x3 = fhg_estimate(df = nls_df)
```
![fits](/man/figures/fits.png)