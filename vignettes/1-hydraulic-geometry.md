# **Channel Hydraulic Geometry and Rating Curves**

## The Need For Paramterizing Hydraulic Geoemtry Properties

Channel hydraulic geometry is a fundamental concept in fluvial geomorphology that describes the relationship between propreties of a river channel such as its size, shape, and dimensions and the flow of water within. The key charetristics that are extremly important in any hydraulic and hydrological modeling include channel width, depth, velocity, and cross-sectional shape. Several examples include applications in flood risk assessment, sediment transport and erosion, water resource management and etc. 

There are many measurments of width, depth, velocity recorded at USGS stations such as HYDRoacoustic dataset in support of the Surface Water Oceanographic Topography satellite mission ([HYDRoSWOT](https://data.usgs.gov/datacatalog/data/USGS:57435ae5e4b07e28b660af55)) and many NWIS stations. 

[Leopold & Maddock](https://books.google.com/books?hl=en&lr=&id=K496gE_YsJoC&oi=fnd&pg=PA1&dq=Leopold+and+Maddock+hydraulic+geometry&ots=FjsEDxo5Ci&sig=vysaOxNq0Y9RSDeZpASbHaQHs6E) first identified power law expressions can relate streamflow to top-width, depth, and velocity, hydrologists have been estimating ‘At-a-station Hydraulic Geometries’ (AHG) to describe average flow hydraulics. These simple but powerfull relations offer an alternative aprouch to reduce the size of the rating curve databases and reduce noiseness in data. 

This software uses the Feature based hydraulic geoemtry concept intorduced by [J.M. Johnson et al.]() that offers a computaionally efficent approuch to modeling hydraulic geometry and offers several methods to filtter data outliers, and inappropriate measurments.

## Hydraulic Geomtry Realtions

Under the steady-state continuity assumption we can describe the relationship between streamflow (Q), top-width (TW), depth (Y), and velocity (U) as powerlaw equations

TW = a*Q^b

Y = c*Q^f

U = k*Q^m

