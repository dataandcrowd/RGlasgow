library(tidyverse)
library(sf)
library(automap)
#library(raster)
#library(gstat)


load("data/no2.RData")
stations <- read_sf("data/stations_10km.shp")
stations_df <- stations %>% filter (`F/R` == "Fixed")%>% st_set_geometry(NULL)
seoul <- read_sf("data/Seoul_City.shp") %>% as('Spatial') %>% fortify()


#-Merge Files for Summer 2013-##
no2_summer_df <- left_join(no2.sum.bk, stations_df, by = c("Station.ID" = "Station", "X" = "X", "Y" = "Y"))

# Some data exploration
no2_summer_df %>% 
  dplyr::select(starts_with("no2")) %>% 
  pivot_longer(everything(), names_to = "Date", values_to = "NO2") %>% 
  mutate(Date = str_replace(Date, "no2_", "")) -> no2_summer_df_tidy


summary(no2_summer_df_tidy)


# 3 day trend August 15-17th 2013
no2_summer_df_tidy %>% 
  slice(29:34) %>% #
  ggplot(aes(x= Date, y = NO2, fill= NO2)) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggsave("results/Hist_no2_hist_Aug.png", width = 5, height = 3.5)#, scale = 1.5)


# convert to shapefile points
no2_summer <- no2_summer_df %>% dplyr::select(X, Y, no2_8_15_day:no2_8_17_night)
#no2_summer <- st_as_sf(no2_summer, coords = c("X", "Y"), crs = 5181) %>% as('Spatial') fortify()

library(raster)
coordinates(no2_summer) <- ~X+Y
proj4string(no2_summer) <- CRS("+init=epsg:5181")


# Semivariogram
options(warn = -1) # don't print warnings
myVario <- list()

for(i in 1:6){
  myVario[[length(myVario)+1]] <- autofitVariogram(no2_summer[[i]] ~ 1, no2_summer, cutoff = 15000, width = 3000)
}

library(gridExtra)

p01 <- plot(myVario[[1]]$exp_var, myVario[[1]]$var_model, main = "Aug 15th\nDay hours")
p02 <- plot(myVario[[2]]$exp_var, myVario[[2]]$var_model, main = "Aug 15th\nNight hours")
p03 <- plot(myVario[[3]]$exp_var, myVario[[3]]$var_model, main = "Aug 16th\nDay hours")
p04 <- plot(myVario[[4]]$exp_var, myVario[[4]]$var_model, main = "Aug 16th\nNight hours")
p05 <- plot(myVario[[5]]$exp_var, myVario[[5]]$var_model, main = "Aug 17th\nDay hours")
p06 <- plot(myVario[[6]]$exp_var, myVario[[6]]$var_model, main = "Aug 17th\nNight hours")


varplot <- grid.arrange(p01, p02, p03, p04, p05, p06, nrow = 3)

ggsave("results/Semivariogram_no2_semvario_08_S2.png",varplot, width = 6, height = 8)#, scale = 1.5)



### Data Frame
seoul_grid <- data.frame(expand.grid(X = seq(min(no2_summer$X), max(no2_summer$X), length=200),
                                     Y = seq(min(no2_summer$Y), max(no2_summer$Y), length=200)))
coordinates(seoul_grid) <- ~X+Y
proj4string(seoul_grid) <- CRS("+init=epsg:5181")

#https://gis.stackexchange.com/questions/157279/saving-results-in-automap-r-package-for-time-series-data

##############
#--Kriging--##
##############
pred.model <- seoul_grid@coords
var.model <- seoul_grid@coords


kriging_result <- autoKrige(no2_summer@data[[6]]~X+Y,no2_summer,seoul_grid)
plot(kriging_result)


for(i in 1:6) {
  kriging_new <- autoKrige(no2_summer@data[[i]] ~X+Y,
                           no2_summer,
                           seoul_grid)
  xyz <- as.data.frame(kriging_new$krige_output$var1.pred)
  colnames(xyz) <- colnames(no2_summer@data)[i]
  pred.model <- cbind(pred.model, xyz)
}

## In ggplot!

##-- Add ColNames
colnames(pred.model) <- c("X", "Y", "aug15d", "aug15n", "aug16d", "aug16n", "aug17d", "aug17n")

##-- Find Mean and variance
stat <- pred.model %>% dplyr::select(-c(X,Y)) %>% 
  gather(factor_key = T) %>% 
  group_by(key) %>% summarise(mean= round(mean(value),1), median = round(median(value),1), 
                              sd= round(sd(value),1), max = max(value),min = min(value)) %>% rename(Hour = key)



##-- Plotting

krige_df <- pred.model %>% 
  pivot_longer(!c("X", "Y"), names_to = "Hour", values_to = "NO2") 

krige_df %>% 
  ggplot() +
  geom_tile(aes(x = X, y = Y, fill = NO2)) +
  scale_fill_distiller(palette = "Spectral", na.value = NA, limits = c(5,25), breaks = c(5,15,25)) +
  geom_contour(aes(x = X, y = Y, z = NO2),bins = 20, colour = "grey40", alpha = 0.7) +
  geom_path(data = seoul, aes(x = long, y = lat), color = 'black', size = 1) +
  geom_text(data = stat, aes(-Inf, -Inf, label = paste0("mean = " , mean)), hjust = -.1, vjust = -2, size = 3.5) +
  geom_text(data = stat, aes(-Inf, -Inf, label = paste0("sd = " , sd)), hjust = -.1, vjust = -1, size = 3.5) +
  facet_wrap(~ Hour, ncol = 8) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 20),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15)                                  
  ) 

# Export PNG
ggsave("results/NO2_kriged_pred_08_S2.png", width = 10, height = 2, dpi = 300)


# convert to Raster Bricks
krige <- rasterFromXYZ(pred.model, 
                       crs="+proj=tmerc +lat_0=38 +lon_0=127 +k=1 +x_0=200000 +y_0=500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
                       digits=5)

ras.road <- raster("data/road_10km_re.tif")  # Import raster
res.mgcv <- resample(krige, ras.road, method = "bilinear") # resample 
res.mgcv <- merge(ras.road, res.mgcv) # merge

# assign road
road_01 = road_02 = road_03 = road_04 = road_05 = road_06 = ras.road

# stack raster and remove individual raster files
road.stack <- stack(road_01, road_02, road_03, road_04, road_05, road_06)
rm(road_01, road_02, road_03, road_04, road_05, 
   road_06)

# add road ratio values to GAM raster
ratio.mid.aug <- no2.sum.ratio[29:34,]

for(i in 1:6){
  #road.stack[[i]] <- road.stack[[i]] * ratio.no2.sum$ratio[i]
  values(road.stack)[values(road.stack[[i]]) == 1] <- ratio.mid.aug$Back.Road.Ratio[i]
  values(road.stack)[values(road.stack[[i]]) == 2] <- ratio.mid.aug$Back.High.Ratio[i]
}

# add no2 and road values
r.poll.rd <- overlay(res.mgcv, road.stack, fun = function(x,y){ifelse(y != 0, x*y, x)})
names(r.poll.rd) <- c("feb15d", "feb15n", "feb16d", "feb16n", "feb17d", "feb17n")

ras <- xyFromCell(r.poll.rd, 1:ncell(r.poll.rd))
krige.df <- as.data.frame(r.poll.rd) 

##-- Find Mean and variance

ras.krige.stat <- data.frame(ras, krige.df)

stat1 <- ras.krige.stat %>% dplyr::select(-c(x,y)) %>% 
  gather(factor_key = T) %>% 
  group_by(key) %>% summarise(mean= round(mean(value),1), sd= round(sd(value),1), max = max(value),min = min(value)) %>% 
  rename(Hour = key)

#####
ras.krige.df <- data.frame(ras, krige.df) %>% 
  pivot_longer(!c("x", "y"), names_to = "Hour", values_to = "NO2") 

ras.krige.df %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = NO2)) +
  scale_fill_distiller(palette = "Spectral", na.value = NA, limits = c(10,50), breaks = c(10,20,30,40,50)) +
  geom_text(data = stat1, aes(-Inf, -Inf, label = paste0("mean = " , mean)), hjust = -.1, vjust = -2, size = 3.5) + 
  geom_text(data = stat1, aes(-Inf, -Inf, label = paste0("sd = " , sd)), hjust = -.1, vjust = -1, size = 3.5) + 
  geom_path(data = seoul, aes(x = long, y = lat), color = 'black', size = 1) +
  facet_wrap(~ Hour, ncol = 8) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 20),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15)                                  
  ) -> final


# Export PNG
ggsave("results/NO2_kriged_final_02_S2.png", final, width = 10, height = 2, dpi = 300)




