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

RasterFolder <- "your raster directory here"  #where the rasters will be downloaded
VectorFolder <- "your vector directory here"  #where you already have the vector/admin boundary file
OutputFolder <- "your output directory here"  #where the results will be output 

#GET RASTER DATA TOGETHER
setwd(RasterFolder)

#Get the BEST SI-x data 1901-2013 SI-x FLI dates, using a loop
for(year in 2010:2010){
  npn_download_geospatial('si-x:average_leaf_best', paste0(year,'-01-01'),output_path=paste0(year,"_si-x_fli_best",".tif"))
}

#Get the BEST SI-x data 1901-2013 SI-x FBI dates, using a loop
for(year in 1901:2013){
  npn_download_geospatial('si-x:average_bloom_best', paste0(year,'-01-01'),output_path=paste0(year,"_si-x_fbi_best",".tif"))
}

# Read in raster data with SI-x Day of Year (DOY) from 1901 to 2012 (you may have to remove that aux.xml files that show up in the directory sometimes)
leaf_files <- list.files(getwd(),pattern="*fli_best.tif") 
FLI_stack <- stack(leaf_files) 
NAvalue(FLI_stack) <- -9999
crs(FLI_stack) #check coordinate reference system
plot(FLI_stack[[95:97]], NAcol='blue')  #check out a few

bloom_files <- list.files(getwd(),pattern="*fbi_best.tif")
FBI_stack <- stack(bloom_files) 
NAvalue(FBI_stack) <- -9999
crs(FBI_stack) #check coordinate reference system
plot(FBI_stack[[100:103]], NAcol='blue')  #check out a few


#GET VECTOR DATA TOGETHER
setwd(VectorFolder)

# Create the buffered version of the USFS forest boundaries (30km)
poli <- st_read("AdminForest_Gen.shp")
ps <- st_as_sf(poli)
ps1 <- st_transform(ps, 32198)
ps1_simp <- st_simplify(ps1, dTolerance = 2000)
plot(ps1_simp)
st_area(ps1_simp)
ps_buff <-st_buffer(ps1_simp, dist=30000)
plot(ps_buff)
st_write(ps_buff, 'buffered_usfs_units.shp')

#Read in buffered boundaries as a special feature 
pl <- st_read("buffered_usfs_units.shp")
pli <- st_as_sf(pl)
crs(pli)
plit <- st_transform(pli, 4269) #code for NAD83
crs(plit)

#ANALYZE VECTOR AND RASTER TOGETHER

#check that the raster and the vector line up for leaf
plot(FLI_stack[[95]], NAcol='blue') 
plot(plit, add=TRUE)

#check that the raster and the vector line up for bloom
plot(FBI_stack[[95]], NAcol='blue') 
plot(plit, add=TRUE)

#create a new column for each year with the DOY, weighted by how much of the polygon falls into each cell. This takes ~20 mins to run for me. 
#Note - I think you may have to empty your environment and begin again to extract bloom values after leaf - you don't want a plit object with both leaf and bloom, you want the 2 CSV files.
setwd(OutputFolder)

#Extract weighted means for leaf
plit$FLI <- exact_extract(FLI_stack, plit, 'weighted_mean', weights=area(FLI_stack))
rownames(plit$FLI) <- pli$FORESTNAME
tab <-head(plit, n=5) #check this file
write.csv(plit$FLI, file = 'USFS_FLI_values_geo_buff.csv')

#Extract weighted means for bloom
plit$FBI <- exact_extract(FBI_stack, plit, 'weighted_mean', weights=area(FBI_stack))
rownames(plit$FBI) <- pli$FORESTNAME
write.csv(plit$FBI, file = 'USFS_FBI_values_geo_buff.csv')

#Call back up the CSVs
FLY <- readr::read_csv("USFS_FLI_values_geo_buff.csv")
BLY <- readr::read_csv("USFS_FBI_values_geo_buff.csv")

# Transpose the data frame, so each column is a forest or park UNIT (excluding the unit names/col 1 bc strings mess up transposing)
FLI_forest <- t(as.matrix(FLY[c(2:113)]))
FBI_forest <- t(as.matrix(BLY[c(2:113)]))


#TEN YEAR MOVING WINDOW MEANS

#Apply moving window mean - this must be done where each FOREST is a column
FLI_10 <- rollapplyr(FLI_forest,10,mean,na.rm=TRUE, fill = NA)
FBI_10 <- rollapplyr(FBI_forest,10,mean,na.rm=TRUE, fill = NA)

# Drop NA rows
FLI_10 <- na.omit(FLI_10)
FBI_10 <- na.omit(FBI_10)

#Transpose back, so that YEARS are columns
FLI_10_year <- t(FLI_10)
FBI_10_year <- t(FBI_10)

# Add forest names to data 
rownames(FLI_10_year) <- plit$FORESTNAME
rownames(FBI_10_year) <- plit$FORESTNAME


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


# Create a function using apply to calculate a percentile for an observed value, using ecdf, for every row in the matrix
recent_percentile <- function(x) {
  percentile <- ecdf(x) # Define ecdf function for row x
  most_recent <- percentile(x[length(x)])  # Since x is a vector, not a matrix, we use length
  return(most_recent * 100)
}

# Apply the 'recent_percentile' function to the rolling window mean dataframes
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

results$leaf_ave <- rowMeans(results[1:3], na.rm=TRUE)
results$bloom_ave <- rowMeans(results[4:6], na.rm=TRUE)
rownames(results) <- plit$FORESTNAME

#Write the file with percentile values and leaf and bloom means
write.csv(results, file = 'USFS_FLI_FBI_10_20_30.csv')

#Restore the geometry column and forest names
results <- cbind(rownames(results), data.frame(results, row.names=NULL))
colnames(results)[1] <- "FORESTNAME"

results2map = plit %>% 
  right_join(results, by = "FORESTNAME")

#Plot to check, the average across moving windows for leaf
plot(results2map$leaf_ave)

#Write the shapefile with the results and geography
setwd(VectorFolder)
st_write(results2map, 'usfs_results.shp')