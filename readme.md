# Spring Index - Patterns of Spring Onset at the USFS

This repository contains two scripts: 

## Recent Springs Relative to the 122 Year Historic Range of Variability

Folder: Long-term_changes

Script: USFS_Spring_LongTerm_Change.R

This script replicates the Monahan et al 2016 methods (https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.1465) to determine the location of the most recent 10, 20 and 30 year periods relative to the 122 year historic range of variability for each unit in the National Forest System.

The input rasters are the Spring Index Leaf Leaf and Bloom (Average, BEST) layers for 1900-2021. The code downloads and names these files in your working directory. 

The script produces a CSV file with the percentile value for recent springs relative to the historic range of variability (HRV), for each moving window (10, 20 and 30 years) and for each index (leaf and bloom) by National Forest, as well as means by index across moving windows. It also produces a shapefile of National Forests administrative boundaries with the values described above for each forest.

## Recent (1991-Present) Patterns in Spring Arrival

Folder: Recent_patterns

Script: USFS_Recent_Spring_Patterns.R
This script enables you to calculate the day of year of spring arrival 1991-present, based on PRISM data, National Forest System.

The input rasters are the Spring Index Leaf Leaf and Bloom (Average, PRISM) layers for 1991-present. The code downloads and names these files in your working directory. 






