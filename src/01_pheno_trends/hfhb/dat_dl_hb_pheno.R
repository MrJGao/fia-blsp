# ******************************************************************************
# Download phenology data from Hubbard Brook. Some code were adapted from their
# website 
# 
# Author: Xiaojie Gao
# Date: 2024-07-09
# ******************************************************************************
library(data.table)


# Package ID: knb-lter-hbr.51.14 Cataloging System:https://pasta.edirepository.org.
# Data set title: Hubbard Brook Experimental Forest: Routine Seasonal Phenology Measurements, 1989 - present.
# Data set creator:    - USDA Forest Service, Northern Research Station
# Contact:  Nina Lany -  USDA Forest Service, Northern Research Station  - nina.lany@usda.gov
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu

inUrl1 <- "https://pasta.lternet.edu/package/data/eml/knb-lter-hbr/51/14/9f623c83fb1da7595c6d2d498bde15df"
infile1 <- "data/raw/hb_pheno.csv"
try(download.file(inUrl1, infile1, method = "curl"))
if (is.na(file.size(infile1))) download.file(inUrl1, infile1, method = "auto")


