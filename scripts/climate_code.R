#McLaughlin Climate Data
#Data gathered and cleaned by Jess Kowalski from McLaughlin Reserve/Knoxville Creek California via the Western Regional Climate Center; https://wrcc.dri.edu/cgi-bin/rawMAIN.pl?caucmc
#from 2012-2022
clim <- read_csv("data/clim_decade.csv")
#precip_in is the monthly average precipitation in inches (all other variables are clearly labelled)

library(weathermetrics)

ten_yr_mean_temp_F <- mean(clim$dailyav_temp_F); ten_yr_mean_temp_F

ten_yr_mean_precip_in <- mean(clim$precip_in); ten_yr_mean_precip_in

growing_season <- clim[94:102,]
growing_mean_precip_in <- mean(growing_season$precip_in); growing_mean_precip_in
growing_mean_precip_cm <- growing_mean_precip_in*2.54 ; growing_mean_precip_cm
growing_mean_temp_F <- mean(growing_season$dailyav_temp_F); growing_mean_temp_F
growing_mean_temp_C <- fahrenheit.to.celsius(growing_mean_temp_F, round = 2)

max(growing_season$dailyav_temp_F)
min(growing_season$dailyav_temp_F)
hist(growing_season$dailyav_temp_F)
hist(growing_season$precip_in)
