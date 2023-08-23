library(raster)
library(sf)
library(sp)
library(exactextractr)
library(readr)
library(rnpn)
library(rgeos)
library(rgdal)
library(zoo)
library(dplyr)
library(leaflet)
rm(list=ls())

#This script replicates the Monahan et al 2016 methods (https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.1465) 
#to determine the location of the most recent 10, 20 and 30 year periods relative to the historic range of variability
#for each unit in the National Forest System
#Written by Alyssa Rosemartin, Aug 2023

RasterFolder <- "~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/data_sources/rasters/BEST"
VectorFolder <- "~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/USFS/USFS_bound"
OutputFolder <- "~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/USFS"

#GET RASTER DATA PREPPED
setwd(RasterFolder)

#Get the BEST SI-x First Leaf Dates raster data 1900-2021, using a loop
for(year in 1900:2021){
  npn_download_geospatial('si-x:average_leaf_best', paste0(year,'-01-01'),output_path=paste0(year,"_si-x_fli_best",".tif"))
}

#Get the BEST SI-x First Bloom Dates raster data 1900-2021, using a loop
for(year in 1900:2021){
  npn_download_geospatial('si-x:average_bloom_best', paste0(year,'-01-01'),output_path=paste0(year,"_si-x_fbi_best",".tif"))
}

# Read leaf and bloom raster data into a raster stack  (you may have to remove that aux.xml files that show up in the directory sometimes)
leaf_files <- list.files(getwd(),pattern="*fli_best.tif") 
FLI_stack <- stack(leaf_files) 
NAvalue(FLI_stack) <- -9999 #recode -9999 to NA
crs(FLI_stack) #check coordinate reference system
plot(FLI_stack[[118:119]])  #check out a few layers 

bloom_files <- list.files(getwd(),pattern="*fbi_best.tif")
FBI_stack <- stack(bloom_files) 
NAvalue(FBI_stack) <- -9999 #recode -9999 to NA
crs(FBI_stack) #check coordinate reference system
plot(FBI_stack[[119:122]])  #check out a few layers

#GET VECTOR DATA PREPPED
#This is the polygon shapefile for the USFS Units

setwd(VectorFolder)

#Simplify the complex polygons of the default USFS Unit polygons
#Create a 30km buffer around the polygons
#This is only done once, you can use the existing simplified and buffered version in the repo


poli <- st_read("AdminForest_Gen.shp")  #Read in vector data
ps <- st_as_sf(poli) #make it a special feature
ps1 <- st_transform(ps, 32198) #set to Quebec Lampert CRS since that is the native one, right?
ps1_simp <- st_simplify(ps1, dTolerance = 2000) #simplify the boundary
plot(ps1_simp) #check how it looks
st_area(ps1_simp) #determine area
ps_buff <-st_buffer(ps1_simp, dist=30000) #buffer at 30KM
plot(ps_buff) # check how it looks 
st_write(ps_buff, 'buffered_usfs_units.shp') #write it out as a shapefile

#Start here, shapefile already exists
#Read in buffered boundaries as a special feature 
pl <- st_read("buffered_usfs_units.shp") #Read in vector data
pli <- st_as_sf(pl) #make it a special feature
crs(pli) #check CRS
plit <- st_transform(pli, 4269) #set it to be code for NAD83 to match raster
crs(plit) #make sure that worked

#side effort, get polygon centroids for PRISM-based dashboard
centroids <- getSpPPolygonsLabptSlots(pl)

#ANALYZE VECTOR AND RASTER TOGETHER

#check that the raster and the vector line up for leaf
plot(FLI_stack[[95]]) 
plot(plit, add=TRUE)

#check that the raster and the vector line up for bloom
plot(FBI_stack[[95]]) 
plot(plit, add=TRUE)

#Create a new column for each year with the DOY, weighted by how much of the polygon falls into each cell. This takes ~20 mins each (leaf and bloom) to run for me. 
setwd(OutputFolder)

#Extract weighted means for leaf - eg, identify the day of year of leaf out for each year, for each polygon, weighted by how much
#of the SI-x raster data falls in the pixel
plit$FLI <- exact_extract(FLI_stack, plit, 'weighted_mean', weights=area(FLI_stack))
rownames(plit$FLI) <- pli$FORESTNAME
tab <-head(plit, n=5) #check this file- should have forest info + 122 columns w leaf dates
write.csv(plit$FLI, file = 'USFS_FLI_values_geo_buff_1900-2021_Aug22.csv') #write this file out

#remove and remake a clean plit and tab objects
#(because you don't want these objects with both leaf and bloom, you want 2 separate CSV files).
rm(plit, tab)
plit <- st_transform(pli, 4269) #code for NAD83

#Extract weighted means for bloom
plit$FBI <- exact_extract(FBI_stack, plit, 'weighted_mean', weights=area(FBI_stack))
rownames(plit$FBI) <- pli$FORESTNAME
tab <-head(plit, n=5) #check this file- should have forest info + 122 columns w bloom dates
write.csv(plit$FBI, file = 'USFS_FBI_values_geo_buff_1900-2021_Aug22.csv') #write this file out

#Call back up these CSVs
setwd(OutputFolder)
FLY <- readr::read_csv("USFS_FLI_values_geo_buff_1900-2021_Aug22.csv")
BLY <- readr::read_csv("USFS_FBI_values_geo_buff_1900-2021_Aug22.csv")

#Transpose the data frame, so each column is a forest or park UNIT (excluding the unit names/col 1 bc strings mess up transposing)
#You have to transpose because the roll apply that's coming up won't work with the years as columns, has to be years as rows.
FLI_forest <- t(as.matrix(FLY[c(2:123)]))  #make sure this period is correct, you are getting all your columns - 122 years in this case (obvi learned this the hard way)
FBI_forest <- t(as.matrix(BLY[c(2:123)]))


#TEN YEAR MOVING WINDOW MEANS
#Create 113 10y rolling window means for the 122 year period

#Apply moving window mean - this must be done where each FOREST is a column
FLI_10 <- rollapplyr(FLI_forest,10,mean,na.rm=TRUE, fill = NA)
FBI_10 <- rollapplyr(FBI_forest,10,mean,na.rm=TRUE, fill = NA)

# Drop NA rows - this is where we move from 122Y record to 113 moving windows (losing 9 rows)
FLI_10 <- na.omit(FLI_10)
FBI_10 <- na.omit(FBI_10)

#Transpose back, so that YEARS are columns
FLI_10_year <- t(FLI_10)
FBI_10_year <- t(FBI_10)

# Add forest names back to data 
rownames(FLI_10_year) <- plit$FORESTNAME
rownames(FBI_10_year) <- plit$FORESTNAME

#saved these out to troubleshoot the period - 122 years results in 113 moving windows
#not necessary otherwise to save these out.
write.csv(FLI_10_year, "FLI_10Y_Moving_Windows.csv")
write.csv(FBI_10_year, "FBI_10Y_Moving_Windows.csv")

#TWENTY YEAR MOVING WINDOW MEANS

#Apply moving window mean - this must be done where each FOREST is a column
FLI_20 <- rollapplyr(FLI_forest,20,mean,na.rm=TRUE, fill = NA)
FBI_20 <- rollapplyr(FBI_forest,20,mean,na.rm=TRUE, fill = NA)

# Drop NA rows
FLI_20 <- na.omit(FLI_20)
FBI_20 <- na.omit(FBI_20)

#Transpose back, so that YEARS are columns
FLI_20_year <- t(FLI_20)
FBI_20_year <- t(FBI_20)

# Add forest names to data 
rownames(FLI_20_year) <- plit$FORESTNAME
rownames(FBI_20_year) <- plit$FORESTNAME


#THIRTY YEAR MOVING WINDOW MEANS

#Apply moving window mean - this must be done where each FOREST is a column
FLI_30 <- rollapplyr(FLI_forest,30,mean,na.rm=TRUE, fill = NA)
FBI_30 <- rollapplyr(FBI_forest,30,mean,na.rm=TRUE, fill = NA)

# Drop NA rows
FLI_30 <- na.omit(FLI_30)
FBI_30 <- na.omit(FBI_30)

#Transpose back, so that YEARS are columns
FLI_30_year <- t(FLI_30)
FBI_30_year <- t(FBI_30)

# Add forest names to data 
rownames(FLI_30_year) <- plit$FORESTNAME
rownames(FBI_30_year) <- plit$FORESTNAME

#Fantastic, now you have all the rolling window means of leaf and bloom DOY for all the forest units
#This the dataset you are going to use to identify where the most recent 10, 20 or 30 year falls in the 
#distribution of all the rolling window means.

#Create a function to calculate a percentile for an observed value, using ecdf (empirical cumulative distribution function)
#for every row in the matrix
recent_percentile <- function(x) {
  percentile <- ecdf(x) # Define ecdf function for row x
  most_recent <- percentile(x[length(x)])  # Since x is a vector, not a matrix, we use length
  return(most_recent * 100) #gives us the position of the most recent rolling window mean as a percentile
}

# Apply this 'recent_percentile' function to the rolling window mean dataframes
fli_percentile_ten <- apply(FLI_10_year, MARGIN = 1, recent_percentile)
fbi_percentile_ten <- apply(FBI_10_year, MARGIN = 1, recent_percentile)

fli_percentile_twenty <- apply(FLI_20_year, MARGIN = 1, recent_percentile)
fbi_percentile_twenty <- apply(FBI_20_year, MARGIN = 1, recent_percentile)

fli_percentile_thirty <- apply(FLI_30_year, MARGIN = 1, recent_percentile)
fbi_percentile_thirty <- apply(FBI_30_year, MARGIN = 1, recent_percentile)


#Take these results as matrices
Leaf_Ten = as.matrix(fli_percentile_ten)
Leaf_Twenty = as.matrix(fli_percentile_twenty)
Leaf_Thirty = as.matrix(fli_percentile_thirty)

Bloom_Ten = as.matrix(fbi_percentile_ten)
Bloom_Twenty = as.matrix(fbi_percentile_twenty)
Bloom_Thirty = as.matrix(fbi_percentile_thirty)

#Combine them all into a single data frame (call up the rownames again)
results = data_frame(Leaf_Ten, Leaf_Twenty, Leaf_Thirty, Bloom_Ten, Bloom_Twenty, Bloom_Thirty)
rownames(results) <- plit$FORESTNAME

#Determine the average percentile ACROSS moving window sizes for leaf and bloom
#(eg, average of the 10, 20, 30 year window recent percentile result for leaf for each unit)
results$leaf_ave <- rowMeans(results[1:3], na.rm=TRUE)
results$bloom_ave <- rowMeans(results[4:6], na.rm=TRUE)
rownames(results) <- plit$FORESTNAME

#Write out the file with percentile values and leaf and bloom means
#This is a good spreadsheet to look at, sanity check and compare to manual calculations for some units to make sure it's right
write.csv(results, file = 'USFS_FLI_FBI_10_20_30_1900-2021_Aug22.csv')

#Restore the geometry column and forest names - because we want to plot the average percentile results
results <- cbind(rownames(results), data.frame(results, row.names=NULL))
colnames(results)[1] <- "FORESTNAME"

results2map = plit %>% 
  right_join(results, by = "FORESTNAME")

#Plot XY to check, the average across moving windows for leaf
plot(results2map$leaf_ave)

#Write the shapefile with the results and geography
setwd(VectorFolder)
st_write(results2map, 'usfs_results_thru2021_Aug22.shp')

#Pull back up to create summary table and then map
results2map <- st_read("usfs_results_thru2021_Aug22.shp")

#Create table to show summary results ordered by Region and then Forest Name
#There is some weirdness where the column names are messed up when you bring this file back in
#(eg bloom_ave becomes bloom_v)
results4table <- subset(results2map, select=c('REGION','FORESTNA','leaf_av','bloom_v')) %>%
  st_drop_geometry() %>%
  rename("Region" = "REGION",
         "Forest" = "FORESTNA",
         "Leaf %" = "leaf_av",
         "Bloom %" = "bloom_v") %>%
  arrange(Region, Forest) 

results4table$`Leaf %` <- round(results4table$`Leaf %`, digits = 2)
results4table$`Bloom %` <- round(results4table$`Bloom %`, digits = 2)

results4table$`Leaf %` <- results4table$`Leaf %` / 100
results4table$`Bloom %` <- results4table$`Bloom %` / 100

#write out this summary table for further formatting in Excel and display as image
setwd(OutputFolder)
write_csv(results4table, "average_percentiles_byForestbyRegion.csv")

#Create a map for leaf and a map for bloom

#Look at distribution of percentile values for leaf:
plot(results2map$leaf_av)

#Map prep - binning, colors, palette definitions for Leaf
bins <- c(0, 5, 25, 75, 95, 100)  
cb <- c('#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0')
pal <- colorBin(palette = cb, bins = bins, domain = results2map$leaf_av,  na.color = "transparent")

#Create map for Leaf
leaflet() %>% addTiles() %>%
  addPolygons(data = results2map, color = "black", fillColor = ~pal(results2map$leaf_av), fillOpacity = 0.8, weight = 1) %>%
  addLegend("bottomright", pal = pal, values = results2map$leaf_av, title = "Percentiles", opacity = 1)


#Look at distribution of percentile values for bloom:
plot(results2map$bloom_v)

#Map prep - binning, colors, palette definitions for Bloom
bins <- c(0, 5, 25, 75, 95, 100)  
cb <- c('#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0')
pal <- colorBin(palette = cb, bins = bins, domain = results2map$bloom_v,  na.color = "transparent")

#Create map for Bloom
leaflet() %>% addTiles() %>%
  addPolygons(data = results2map, color = "black", fillColor = ~pal(results2map$bloom_v), fillOpacity = 0.8, weight = 1)  %>%
  addLegend("bottomright", pal = pal, values = results2map$leaf_av, title = "Percentiles", opacity = 1)
