#stack bird range rasters (with gen time) and make plots

#rasterize range map

#see here: https://github.com/RS-eco/rasterSp

#could also use st_make_grid in sf or possibly turn raster into sf object to do this


devtools::install_github('RS-eco/rasterSp')
library(rasterSp)


rasterSp::rasterizeRange(dsn = list.files(paste0('~/Desktop/bird/'), 
                                          pattern = ".shp", full.names = TRUE), 
               id = "SCINAME", resolution = 0.5, save = TRUE, touches = TRUE,
               #seasonal = c(1,2), origin = 1, presence = c(1,2), 
               path = paste0('~/Desktop/processed/'))


lf <- list.files('~/Desktop/processed/', full.names = TRUE)
tt <- terra::rast(lf)
plot(tt)
