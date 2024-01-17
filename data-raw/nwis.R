library(dplyr)
library(usethis)

nwis <- read.table(paste0('https://waterdata.usgs.gov/nwis/measurements?',
                          'site_no=01096500',
                          '&agency_cd=USGS',
                          '&format=rdb_expanded'),
           sep="\t", header=TRUE)[-1,] %>% 
  select(siteID = site_no, 
         date = measurement_dt, 
         Q = chan_discharge, 
         TW = chan_width, 
         A = chan_area, 
         V = chan_velocity) %>% 
  mutate(date = as.Date(date), 
         Q_cms = as.numeric(Q)* 0.028316847000000252,
         TW_m = as.numeric(TW)* 0.3048,
         A_m2 = as.numeric(A) * 0.092903,
         V_ms = as.numeric(V) * 0.3048,
         Y_m = A_m2 / TW_m) %>% 
  filter(complete.cases(.)) %>% 
  select(siteID, date, Q_cms,Y_m,V_ms,TW_m)

use_data(nwis, overwrite = TRUE)
