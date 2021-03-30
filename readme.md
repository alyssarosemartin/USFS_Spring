# Spring Index - Recent Springs Relative to Historic Range of Variability

This repository contains the code for analyses of the Historic Range of Variability for Gridded Spring Index BEST First Leaf and First Bloom Layers.

## SIX HRV USFS

The file SIX_HRV.R produces a CSV file with the percentile value for recent springs relative to the historic range of variability (HRV), for each moving window (10, 20 and 30 years) and for each index (leaf and bloom) by National Forest, as well as means by index across moving windows. It also produces a shapefile of National Forests administrative boundaries with the values described above for each forest.


The code runs in R, with rnpn, raster, sf, sp, curl, rgdal, rgeos, readr, exactextractr, zoo, dplyr libraries loaded.


The input rasters are the Spring Index Leaf Leaf and Bloom (Average, BEST) layers for 1901-2012. The code downloads and names these files in your working directory. It also requires the shape file of the administrative boundaries of the USFS, available here:https://data.fs.usda.gov/geodata/edw/datasets.php?dsetCategory=boundaries
