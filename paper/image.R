library(AHGestimation)
library(ggplot2)
library(patchwork)

# Processing
# 
nwis = AHGestimation::nwis |>
  dplyr::rename(Q = Q_cms, TW = TW_m, Y = Y_m, V = V_ms)

filter_data = nwis |>
  date_filter(10, keep_max = TRUE) |>
  nls_filter(allowance = .5) 

ahg_fit = ahg_estimate(filter_data)[1,]

shape = compute_hydraulic_params(ahg_fit)

cs = cross_section(r = shape$r,  TW = max(filter_data$TW), Ymax = max(filter_data$Y))

# Plotting
p1 = ggplot() +  
  geom_point(data = nwis, aes(x = Q, y = Y, col = "All data")) + 
  geom_point(data = filter_data, aes(x = Q, y = Y, col = "Filtered Data")) + 
  geom_line(data = filter_data, aes(x = Q, y = ahg_fit$Y_coef[1] * Q ^ (ahg_fit$Y_exp[1]), col = "Filtered Data")) +
  theme_light()  + 
  labs(title = 'Q-Y Relationship', color = "", y = "Y (m)", x = 'Q (cms)') + 
  scale_color_manual(values = c("Filtered Data" = "red", "All data" = "black")) + 
ggplot(data = nwis ) + 
  geom_point(data = nwis, aes(x = Q, y = TW, col = "All data")) + 
  geom_point(data = filter_data, aes(x = Q, y = TW, col = "Filtered Data")) + 
  geom_line(data = filter_data, aes(x = Q, y = ahg_fit$TW_coef[1] * Q ^ (ahg_fit$TW_exp[1]), col = "Filtered Data")) +
  theme_light() + 
  labs(title = 'Q-TW Relationship', color = "", y = "TW (m)", x = 'Q (cms)') + 
  scale_color_manual(values = c("Filtered Data" = "red", "All data" = "black")) + 
ggplot(data = nwis ) + 
  geom_point(data = nwis, aes(x = Q, y = V, col = "All data")) + 
  geom_point(data = filter_data, aes(x = Q, y = V, col = "Filtered Data")) + 
  geom_line(data = filter_data, aes(x = Q, y = ahg_fit$V_coef[1] * Q ^ (ahg_fit$V_exp[1]), col = "Filtered Data")) +
  theme_light()  + 
  labs(title = 'Q-V Relationship', color = "", y = "V (m/s)", x = 'Q (cms)') + 
  scale_color_manual(values = c("Filtered Data" = "red", "All data" = "black")) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom') 

p2 = ggplot(data = cs) + 
  geom_point(aes(x = x, y = Y)) + 
  geom_line(aes(x = x, y = Y)) + 
  theme_light() + 
  geom_hline(yintercept = max(cs$Y), color = "darkgreen") +
  geom_hline(yintercept = 2, color = "blue", linewidth = 2) +
  geom_vline(xintercept = mean(cs$x), color = "darkgreen") + 
  labs(x = "Relative Width (m)", y = "Y (m)", title = "Shape") +
ggplot(data = cs) + 
  geom_point(aes(x = A, y = Y)) + 
  geom_line(aes(x = A, y = Y)) + 
  geom_hline(yintercept = 2, color = "blue", linewidth = 2) +
  theme_light() +
  labs(Y = "Depth", x = "Cross Sectional Area (m²)", y = "", title = "Area to Depth") 


p3 = p1 + p2

ggsave(p3, filename = "paper/paper-fig1.png", height=8, width=8, units="in", dpi=600)


