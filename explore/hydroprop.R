library(hydrofabric)
library(ggplot2) 
library(AHGestimation)

source = '/Users/mikejohnson/hydrofabric/v20.1'

compute_hydroprop = function(input, n = 10){
  
  range = seq(min(input$minZ), max(input$maxZ),  length.out = n)
  
  ll = list()
  
  for(i in 1:length(range)){
    
    d = range[i]
    
    tmp = filter(input, rZ <= d)
    
    ll[[i]] = data.frame(
      gid = input$gid[1],
      Y = d,
      wp = sum(sqrt(diff(tmp$relative_distance)^2 + diff(tmp$rZ)^2)),
      area = pmax(0, AUC(x = tmp$relative_distance, y = rep(d, nrow(tmp)), absolutearea = FALSE) - AUC(x = tmp$relative_distance, y = tmp$rZ, absolutearea = FALSE))
    ) %>% 
      mutate(hr = area/wp)
  }
  
  bind_rows(ll) 
}

xs  = glue('{source}/3D/cross-sections') %>% 
  open_dataset() %>% 
  filter(hy_id %in% c('wb-10000', "wb-10001")) %>% 
  select(hf_id = hy_id, cs_id, relative_distance, X, Y, Z) %>% 
  collect() %>% 
  group_by(hf_id, cs_id) %>% 
  mutate(gid = cur_group_id()) %>% 
  ungroup() %>% 
  group_by(hf_id) %>% 
  mutate(rZ = Z - min(Z),
         minZ = min(rZ),
         maxZ = max(rZ)) %>% 
  ungroup()

slope =  compute_channel_slope(extract_thalweg(xs))

hp = bind_rows(pblapply(split(x = xs, f = xs$gid), compute_hydroprop, n = 100)) %>% 
  left_join(select(xs, gid, hf_id, cs_id, minZ, maxZ), by = "gid", relationship = "many-to-many") %>% 
  left_join(slope, by = "hf_id", relationship = "many-to-many") %>% 
  distinct() %>% 
  mutate(n = .05)

hpt = open_dataset(glue('{source}/conus_net.parquet')) %>% 
  select(hf_id = id, lengthkm) %>% 
  inner_join(hp, by = 'hf_id', relationship = "many-to-many") %>% 
  collect() %>% 
  distinct() %>% 
  mutate(Q = (area * hr^(2/3) * sqrt(slope)) / n)

x = group_by(hpt, hf_id, Y) %>% 
  summarise(minQ  = median(Q, na.rm = TRUE), 
            meanQ = mean(Q,   na.rm = TRUE),
            maxQ  = max(Q,    na.rm = TRUE)) %>% 
  ungroup()

ggplot(data = x) + 
  geom_line(aes(y = minQ,  x = Y), color = "blue") + 
  geom_line(aes(y = meanQ, x = Y), color = "red") + 
  geom_line(aes(y = maxQ,  x = Y), color = "green") + 
  facet_wrap(~hf_id) + 
  labs(x = "Y", y = "Q")

(ahg = pblapply(split(x = select(x, Q = maxQ, Y = Y), f = x$hf_id), ahg_estimate) %>% 
  bind_rows(.id = 'hf_id') %>% 
  filter(method == "nls"))

