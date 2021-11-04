library(tidyverse)
library(ggOceanMaps)
library(rgdal)
library(marmap)
library(raster)

# ABISS map

# bathymetry ####
lims <- c(-128, -123, 40, 50)
projection <- "+init=epsg:4326"
#gebcoPath <- "~/Dropbox/Harvard/ConsultingGig/map/GEBCO_2020_03_Mar_2021_b16ddbc41cd2/"
etopoPath <- "~/Dropbox/Harvard/ABISS/Analysis/map/etopo"
topo <- marmap::getNOAA.bathy(-130, -120, 40, 50,resolution=1,keep=TRUE,path=etopoPath)

rb <- raster_bathymetry(bathy = marmap::as.raster(topo),
                        depths = seq(0,4000,by=100), 
                        proj.out = projection, 
                        boundary = lims
)

bathym <- vector_bathymetry(rb)

#from https://maps.ngdc.noaa.gov/viewers/grid-extract/index.html
#coastal releif model, 3sec
rast <- raster("~/Dropbox/Harvard/ABISS/Analysis/map/exportImage.tiff")
rb <- raster_bathymetry(bathy = rast,
                        depths = seq(0,4000,by=100), 
                        proj.out = projection, 
                        boundary = lims
)
bathym <- vector_bathymetry(rb)

#locales ####
mclat = 45 + (50 + 54.2436/60)/60 # 45° 50' 54.2436" 
mclon = 124 + (53 + 43.9548/60)/60 #124° 53' 43.9548"
dt <- tibble(lat = c(mclat, 44.5691, 44+37.5/60, 45+33.3/60, 46+12.4/60),
             lon = c(mclon, 125.1481, 124+2.7/60, 123+55.1/60, 123+46.1/60)*-1,
             site = c("McArthur Promontory","OOI","NOAA","NOAA","NOAA"),
             id = c("McArthur\nPromontory","OOI:\nS. Hydrate Ridge","South Beach\ngauge","Garibaldi\ngauge","Astoria\ngauge"))
box1 <- tibble(lat = c(44.2,44.2,46.6,46.6),lon=c(-126,-123,-123,-126))
box2 <- tibble(lat = c(45.7,45.7,46,46),lon=c(-125.1,-124.7,-124.7,-125.1))

# map zoom 1 ####
basemap(limits = c(-135,-116,35,55)) +
  geom_spatial_polygon(data=box1,aes(x = lon, y = lat), fill = NA, col=1) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1, "cm"),width=unit(.7, "cm")) +
  theme(rect=element_rect(colour="black"),
        line=element_line(colour="black"),
        axis.text = element_text(color="black")) +
  labs(x=" ",y=" ") +
  coord_sf(expand = FALSE) 

# map zoom 2 ####
basemap(limits = c(-126.1,-122.9,44.2,46.6),#lims,
        shapefiles = list(land=ggOceanMapsData::dd_land,glacier=NULL,bathy=bathym),
        bathymetry = TRUE
        ) +
  #geom_spatial_polygon(data=box2,aes(x = lon, y = lat), fill = NA, col=1) +
  labs(x=" ",y=" ") +
  theme(legend.position = "left",
        #rect=element_rect(colour="black"),
        #line=element_line(colour="black"),
        axis.text = element_text(color="black")) +
  geom_spatial_point(data=dt, aes(x=lon, y=lat, shape=site), size=2.5, inherit.aes=FALSE) +
  geom_spatial_text(data=dt, aes(x=lon, y=lat, label=id), 
                    hjust=0,nudge_x=.1, nudge_y=-.1, 
                    size=3,
                    inherit.aes=FALSE) +
  coord_sf(expand = FALSE) 

# map zoom 3 ####
basemap(limits=c(-125.1,-124.7,45.7,46),
        shapefiles = list(land=ggOceanMapsData::dd_land,glacier=NULL,bathy=bathym),bathy.style = "poly_greys",
        bathymetry = TRUE
) +
  geom_spatial_point(data=dt[1,], aes(x=lon, y=lat, shape=site), col="red",inherit.aes=FALSE) +
  annotation_scale(location="br") +
  coord_sf(expand = FALSE) 

# scale bits

fake <- tibble(x=faithfuld$eruptions,y=faithfuld$waiting,
               depths = -faithfuld$density*110000+10)
bath_col_scale <- rev(colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08306B"))(31))

ggplot(fake,aes(x,y)) + 
  geom_raster(aes(fill=depths)) +
  scale_fill_gradientn(colours = bath_col_scale,name = "Depth (m)") 

"bathy_pb = colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08306B"))(x)"


