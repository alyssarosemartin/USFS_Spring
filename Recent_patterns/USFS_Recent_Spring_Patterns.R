library(raster)
library(sf)
library(sp)
library(exactextractr)
library(readr)
library(rnpn)
library(ggplot2)
library(tidyverse)
library(dplyr)
rm(list=ls())

#RasterFolderPRISM <-"~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/data_sources/rasters/PRISM_anomalies"
RasterFolderPRISM <-"~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/data_sources/rasters/PRISM_doy"
VectorFolder <- "~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/USFS/USFS_bound_buff"
OutputFolder <- "~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/USFS/Recent_Trends/Output"

setwd(RasterFolderPRISM)
for(year in 1981:2022){
  npn_download_geospatial('si-x:leaf_anomaly_prism', paste0(year,'-01-01'),output_path=paste0(year,"_si-x_fli_anom",".tif"))
}

for(year in 1981:2022){
  npn_download_geospatial('si-x:bloom_anomaly_prism', paste0(year,'-01-01'),output_path=paste0(year,"_si-x_fbi_anom",".tif"))
}


#Resample 2022 NCEP to match 2020 PRISM for Leaf if needed
setwd(RasterFolderPRISM)
FLI20 <- raster("2020_si-x_fli_anom.tif")
FLI21 <- raster("2021_si-x_leaf_anomaly_ncep.tif")
NAvalue(FLI20) <- -9999
NAvalue(FLI21) <- -9999
FLI21_4k <-resample(FLI21, FLI20, method='bilinear') 
plot(FLI21_4k)
writeRaster(FLI21_4k, "2021si-x_fli_anom_ncep_resample.tif", format="GTiff",overwrite=TRUE, NAflag=-9999)

#Resample 2022 NCEP to match 2020 PRISM for Bloom
FLI20 <- raster("2020_si-x_fbi_anom.tif")
FLI21 <- raster("2021_si-x_bloom_anomaly_ncep.tif")
NAvalue(FLI20) <- -9999
NAvalue(FLI21) <- -9999
FLI21_4k <-resample(FLI21, FLI20, method='bilinear') 
plot(FLI21_4k)
writeRaster(FLI21_4k, "2021si-x_fbi_anom_ncep_resample.tif", format="GTiff",overwrite=TRUE, NAflag=-9999)


#GET RASTER DATA TOGETHER
#Create Raster Stacks for leaf and bloom 1981-2021
setwd(RasterFolderPRISM)
FLI_PRISM <- stack(list.files(getwd(), pattern="fli")) #or pattern "leaf" for the doy values
NAvalue(FLI_PRISM) <- -9999
crs(FLI_PRISM) #check coordinate reference system
plot(FLI_PRISM[[39:42]], NAcol='blue')  #check out a few

FBI_PRISM <- stack(list.files(getwd(),pattern="fbi")) #or pattern "bloom" for the doy values
NAvalue(FBI_PRISM) <- -9999
crs(FBI_PRISM) #check coordinate reference system
plot(FBI_PRISM[[39:42]], NAcol='blue')  #check out a few

#GET VECTOR DATA TOGETHER
setwd(VectorFolder)

#Read in buffered boundaries as a special feature 
pl <- st_read("buffered_usfs_units.shp")
pli <- st_as_sf(pl)
crs(pli)
plit <- st_transform(pli, 4269) #code for NAD83
crs(plit)

pli_cent <- st_transform(pli, 32617) #code for UTM, per recs here https://gis.stackexchange.com/questions/43543/how-to-calculate-polygon-centroids-in-r-for-non-contiguous-shapes

#find centroid of polygons
centroids <- st_centroid(pli_cent)
plot(FLI_PRISM[[40]], NAcol='blue') 
plot(centroids$geometry, pch= 19, add=TRUE) #not working for me right now 

#this works, seems generally right - can get them to plot separately but not together
plot(plit$geometry)
plot(centroids$geometry, add=TRUE)

#ANALYZE VECTOR AND RASTER TOGETHER

#check that the raster and the vector line up for leaf
plot(FLI_PRISM[[40]], NAcol='blue') 
plot(plit, add=TRUE)

#check that the raster and the vector line up for bloom
plot(FBI_PRISM[[20]], NAcol='blue') 
plot(plit, add=TRUE)

#create a new column for each year with the DOY, weighted by how much of the polygon falls into each cell. 
setwd(OutputFolder)

#Extract weighted means for leaf
plit$FLI <- exact_extract(FLI_PRISM[[11:42]], plit, 'weighted_mean', weights=area(FLI_PRISM)) #for years 1991-2022
rownames(plit$FLI) <- pli$FORESTNAME
tab <- head(plit, n=5) #check this file
#write.csv(plit$FLI, file = 'USFS_FLI_values_geo_buff_PRISM.csv')
write.csv(plit$FBI, file = 'USFS_FLI_anom_geo_buff_PRISM.csv')

rm(plit)
#re run line 67 to remake the vector with no leaf info, then move to line 93

#Extract weighted means for bloom
plit$FBI <- exact_extract(FBI_PRISM[[1:41]], plit, 'weighted_mean', weights=area(FBI_PRISM))
rownames(plit$FBI) <- pli$FORESTNAME
tab <- head(plit, n=5)
#write.csv(plit$FLI, file = 'USFS_FBI_values_geo_buff_PRISM.csv')
write.csv(plit$FBI, file = 'USFS_FBI_anom_geo_buff_PRISM.csv')

#these two CSVs were created and combined and sent to Nathan Walker with USFS OSC in Jan 2022 for the prototype

#Call back up the CSVs
FLY <- readr::read_csv("USFS_FLI_values_geo_buff_PRISM.csv")
BLY <- readr::read_csv("USFS_FBI_values_geo_buff_PRISM.csv")



FLY2 <- as.data.frame(FLY)

FLY_Alleg <- filter(row.names(FLY2) %in% "Allegheny National Forest")
FLY_Alleg <- subset.matrix(FLY$X1 == "Allegheny National Forest")


filter(row.names(DF) %in% c("12a","13a"))

datawithoutVF = data[which(rownames(data) %nin% remove), ]



