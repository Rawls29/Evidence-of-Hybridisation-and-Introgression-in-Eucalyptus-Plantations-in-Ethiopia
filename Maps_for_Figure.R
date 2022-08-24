rm(list=ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project/R Scripts')

####Libraries####
library("ggmap")
library(sf)
library(rgdal)
library(broom)

###Loading in metadata####
meta.path="eucalyptus_sample_data_2022.csv"
meta <- read.csv(meta.path)

####Plotting transects on map####
register_google(key="")#Need to input a key here to allow use of google maps
                  
#Mean of latitude and longitute for samples
llmeans <- sapply(meta[10:9], mean)

#Base maps
eth_map_2 <- get_map(location=llmeans,
                     maptype="satellite",
                     source="google",
                     zoom=8)

#Map of sample locations, uniform colour
ggmap(eth_map_2) +
  geom_point(data=meta, aes(x=Long., y=Lat.), color="red", size=2)+
  labs(y="Latitude", x="Longitude")

#Map of sample locations, coloured by transect
ggmap(eth_map_2) +
  geom_point(data=meta, aes(x=Long., y=Lat., colour=ï..Transect), size=2)+
  scale_colour_hue(name="Transect")+
  labs(y="Latitude", x="Longitude")

####Plot of Ethiopia####
shp <- readOGR(dsn="C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project/R Scripts/afr_g2014_2013_0", layer="afr_g2014_2013_0")
shp_tran <- spTransform(shp, CRS("+init=epsg:4326"))
tidy_shp <- tidy(shp_tran)

lat <- c(5, 5, 9, 9)
long <- c(36, 40, 40, 36)
lat <- as.data.frame(lat)
long <- as.data.frame(long)
study_area <- cbind(lat, long)

eth_map_3 <- get_map(location="Ethiopia",
                     maptype="satellite",
                     source="google",
                     zoom=5)

ggmap(eth_map_3, extent="normal", maprange=FALSE) +
  geom_polygon(data=tidy_shp,
               aes(long, lat, group=group),
               fill=NA, colour="black", alpha=0.2)+
  geom_polygon(data=study_area,
               aes(long, lat),
               fill=NA, colour="red")+
  theme_bw()+
  coord_map(projection = "mercator",
            xlim=c(30, 52),
            ylim=c(0, 20))+
  labs(y="Latitude", x="Longitude")
