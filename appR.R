# create Shiny app
# Global VARIABLES CODES STARTS HERE

### load required libraries
pacman::p_load(raster, shiny, RColorBrewer, malariaAtlas, shinydashboard, rgdal, googledrive)
##
# load in shapefile
admin0shapefile <- shapefile('data/countries/admin2013_0.shp')

# download admin 0 raster
# https://drive.google.com/open?id=1xY9meFKuR4UqJvc_LwayI11zupvbj5ZQ
temp <- tempfile()
dl <- drive_download(as_id("1xY9meFKuR4UqJvc_LwayI11zupvbj5ZQ"), path = temp, overwrite = TRUE)
admin0raster <- raster(dl$local_path)

# grab a list of available rasters from the MAP API
raster_meta <- as.data.frame(listRaster())

# subset to remove rows where title = "Relative Abundance Africa"
raster_meta <- raster_meta[!raster_meta$title == "Relative Abundance Africa", ]

# turn returned variable as a list to feed into sourcing function
raster_name_list <- as.list(raster_meta$title)

# feed list into getRaster function
raster_list <- lapply(raster_name_list, getRaster)

# stack the surfaces
# error with surfaces due to differing extents
# loop through and get a vector of raster extents
# extents <- t(sapply(raster_list, function (x) as.vector(extent(x))))
# 
# get unique extents across all rasters
# u_extents <- unique(extents)

# loop through surfaces and generate info for each country; can't stack as there's 43 unique extents

for(i in 1:length(raster_list)){
  
  # grab the raster for that iteration
  r_i <- raster_list[[i]]
  
  # mask it by the extent of the admin 0 raster
  r_i <- crop(r_i, admin0raster)
  
  # check number of columns match (cell size is equal)
  logic_cell <- r_i@ncols == admin0raster@ncols
  
  # if the number of cells isn't equal, resample 
  if(logic_cell == FALSE){
    
    message(paste0("Resampling raster ", i, " to align with master grid"))
    
    r_i <- resample(r_i, admin0raster, method = "bilinear")
    
  }
  
  # generate the mean pixel value for all pixels within a country
  poly_mean <- zonal(r_i, admin0raster, fun = 'mean')
  
  # create a dataframe with the zonal output
  poly_mean <- as.data.frame(poly_mean)
  
  names(poly_mean) <- c("country_id_raster",
                        names(raster_list[[i]]))
  
  # convert to a binary surface indicating values for a country 
  poly_mean[[2]] <- ifelse(is.na(poly_mean[[2]]), 0, 1)
  
  # combine results across each raster
  if(i == 1){
    
    combined <- poly_mean
    
  } else {
    
    combined <- merge(combined, poly_mean, by = "country_id_raster")
    
  }
  
}

# add back ISO code 
iso3 <- as.data.frame(admin0shapefile[, c(1:2, 5)])

names(combined)[1] <- "GAUL_CODE"
combined <- merge(combined, iso3, by = "GAUL_CODE")
# combined$COUNTRY_ID.x <- NULL
# drop countries which are XXX
combined_trim <- combined[!combined$COUNTRY_ID.y == "XXX", ]

# write out combined dataset
write.csv(combined,
          'D:/Dropbox/PhD/collaboration projects/Wellcome Trust data re-use/SHINY-map-prize/data/combined_lookup.csv',
          row.names = FALSE)

## download-rasters.r
# andy south March 2019

# to download and save rasters for Africa from Malaria Atlas Project
# doesn't need to be run by the shiny app, run beforehand to save data

# needs packages declared in server.r

# seemingly had to convert sf back to sp to get bbox in format wanted by getRaster
extent_afr <- extent(c(-20,52), c(-38,38))
extent_afr <- bbox(extent_afr)

pfpr2_10_2015 <- getRaster(surface = "Plasmodium falciparum PR2-10", extent=extent_afr, year = 2015)
pf_incidence_2015 <- getRaster(surface = "Plasmodium falciparum Incidence", extent=extent_afr, year = 2015)
itn_2015 <- getRaster(surface = "Insecticide-treated bednet (ITN) coverage", extent=extent_afr, year = 2015)
time_to_city_2015 <- getRaster(surface = "A global map of travel time to cities to assess inequalities in accessibility in 2015", extent=extent_afr)

# could sample to reduce num pixels
# because this is all that are displayed by leaflet anyway
maxpixels <- 500000
time_to_city_2015 <- raster::sampleRegular(time_to_city_2015, maxpixels, asRaster = TRUE, useGDAL = TRUE)

#to make sure rasters are in memory
pfpr2_10_2015 <- readAll(pfpr2_10_2015)
pf_incidence_2015 <- readAll(pf_incidence_2015)
itn_2015 <- readAll(itn_2015)
time_to_city_2015 <- readAll(time_to_city_2015)

save(pfpr2_10_2015, file='data/rasters/pfpr2_10_2015.rda')
save(pf_incidence_2015, file='data/rasters/pf_incidence_2015.rda')
save(itn_2015, file='data/rasters/itn_2015.rda')
save(time_to_city_2015, file='data/rasters/time_to_city_2015.rda')

## PRE-PROCESS RASTER
# script to process MAP obtained raster files
pacman::p_load(raster, rgdal)

# load in lookup for raster paths
raster_lookup <- read.csv('data/raster_paths.csv',
                          stringsAsFactors = FALSE)

# read in a district raster
districts <- raster('data/getRaster/FAO_admin_1.tif')

raster_lookup$stats_path <- NA

# loop through each available surface, and generate district level statistics
for(i in 1:length(raster_lookup$path)){
  
  # get the path for the raster
  raster_path <- raster_lookup$path[[i]]
  
  # if the path isn't blank, read in the raster and generate zonal statistics
  if(raster_path != ""){
    
    # source raster
    raster_i <- raster(raster_path)
    raster_i[raster_i == -999] <- NA
    raster_i[raster_i == -9999] <- NA
    
    # if it's the same continent process:
    if(!is.null(intersect(extent(raster_i), extent(districts)))){
      
      # crop 'districts' by extent of raster_i
      districts_c <- crop(districts, raster_i)
      raster_i <- crop(raster_i, districts_c)
      
      # if the resolution of the two rasters isn't the same, resample the input raster to match the 5km district raster
      if(ncol(raster_i) != ncol(districts_c)){
        
        raster_i <- resample(raster_i, districts_c, method = "bilinear")
        
      } else {
        
        extent(districts_c) <- extent(raster_i)
        
      }
      
      # generate statistics
      raster_i_mean <- zonal(raster_i, districts_c, fun = 'mean', na.rm = TRUE)
      raster_i_min <- zonal(raster_i, districts_c, fun = 'min', na.rm = TRUE)
      raster_i_max <- zonal(raster_i, districts_c, fun = 'max', na.rm = TRUE)
      raster_i_sd <- zonal(raster_i, districts_c, fun = "sd", na.rm = TRUE)
      
      # merge into one dataframe
      raster_i_stats <- merge(raster_i_mean, raster_i_max, by = "zone")
      raster_i_stats <- merge(raster_i_stats, raster_i_min, by = "zone")
      raster_i_stats <- merge(raster_i_stats, raster_i_sd, by = "zone")
      
      # write out the merged dataframe
      # define a file name
      file_name_i <- substr(raster_path, 16, 1000)
      file_name_i <- gsub("\\.tiff", "_stats.csv", file_name_i)
      
      write.csv(raster_i_stats,
                file = paste0("data/processed/", file_name_i),
                row.names = FALSE)
      
      # add the filepath back into the raster_lookup dataframe
      raster_lookup$stats_path[i] <- paste0("data/processed/", file_name_i)
    }
  }
}

# write out the updated filepath document
write.csv(raster_lookup,
          'data/raster_stats_paths.csv',
          row.names = FALSE)

## admin-polys-simplify.r
# andy south 2/3/2019

# script to create faster admin polygons file for malaria atlas project shiny viewer

# create one polygons object that has both admin0 and admin1 boundaries in it.
# make district lookups simpler with consistent codes
# save as an R object rather than a shapefile for faster loading
# move to using package sf the modern replacement for sp that makes looking at polygon attributes easier.
# simplify polygons so object is smaller and display is faster, country maps look cleaner (polygon detail not needed)


# get iso country ids for africa from existing shapefile (then it may not be needed further)
# this could be replaced by a different method
library(raster)
admin_1 <- raster::shapefile('data/districts/admin_1.shp')
country_ids <- unique(admin_1$COUNTRY_ID)

# get polygons for admin0 and admin1 from as a spatialpolygonsdataframe
admin_0 <- raster::shapefile('data/countries/admin2013_0.shp')
admin_0 <- admin_0[admin_0$COUNTRY_ID %in% country_ids, ]

# convert to sf 
library(sf)
admin_1 <- admin_1[c(2:5, 1)]
names(admin_1) <- c("COUNTRY_ID",
                    "GAUL_CODE",
                    "ADMN_LEVEL",
                    "PARENT_ID",
                    "name")
combined <- rbind(admin_1, admin_0)
sf_africa <- sf::st_as_sf(combined)

# get rid of dodgy countrycodes
sf_africa <- sf_africa[sf_africa$country_id != 'XXX',]
# save as an R object that we can then read into the shiny app
# save(sf_africa, file='test.rda')
# then to plot the countries or districts just subset by admn_level
# ad0ago <- sf_africa[sf_africa$admn_level==0 & sf_africa$country_id=='AGO',]
# st_geometry needed to plot just the polygons
# plot(sf::st_geometry(ad0ago))

# sf_africa has spain in too !
# plot(sf::st_geometry(sf_africa))

#library(mapview)
#mapview(sf_africa)

#and has much more detail than we need e.g. fiddly islands.

#remove spain
sf_africa <- sf_africa[sf_africa$COUNTRY_ID!='ESP',]

# simplify
library(rmapshaper)
#default keeps 0.05 of points keep=0.05
#keep all polygons keep_shapes=TRUE
sf_afr_simp <- rmapshaper::ms_simplify(sf_africa, keep_shapes=TRUE)

save(sf_afr_simp, file='data/sf_afr_simp_fao.rda')

# to load
load('data/sf_afr_simp_fao.rda')

library(mapview)
mapview(sf_afr_simp, zcol='COUNTRY_ID', legend=FALSE)

ggplot(sf_afr_simp) +
  geom_sf(aes(fill=COUNTRY_ID))

# subset by country and then plot coloured named districts

#shorter name
sfafsi <- sf_afr_simp

#subset districts for one country
sfago <- sfafsi[sfafsi$admn_level==1 & sfafsi$country_id=='AGO',]

#datum=na needed to remove gridlines now, may not be in future

ggplot(sfago) +
  geom_sf(aes(fill=name)) + 
  coord_sf(datum = NA) + 
  theme_void()

# if I wanted to add district labels
ggplot(sfafsi) +
  geom_sf() +
  geom_text_repel(aes(x=X, y=Y, label=name))

###### LOOK UP GENERATION 
# load required libraries
pacman::p_load(raster, shiny, RColorBrewer, malariaAtlas, shinydashboard, rgdal, googledrive)
##
# load in shapefile
admin0shapefile <- shapefile('data/countries/admin2013_0.shp')

# download admin 0 raster
# https://drive.google.com/open?id=1xY9meFKuR4UqJvc_LwayI11zupvbj5ZQ
temp <- tempfile()
dl <- drive_download(as_id("1xY9meFKuR4UqJvc_LwayI11zupvbj5ZQ"), path = temp, overwrite = TRUE)
admin0raster <- raster(dl$local_path)

# grab a list of available rasters from the MAP API
raster_meta <- as.data.frame(listRaster())

# subset to remove rows where title = "Relative Abundance Africa"
raster_meta <- raster_meta[!raster_meta$title == "Relative Abundance Africa", ]

# turn returned variable as a list to feed into sourcing function
raster_name_list <- as.list(raster_meta$title)

# feed list into getRaster function
raster_list <- lapply(raster_name_list, getRaster)

# stack the surfaces
# error with surfaces due to differing extents
# loop through and get a vector of raster extents
# extents <- t(sapply(raster_list, function (x) as.vector(extent(x))))
# 
# get unique extents across all rasters
# u_extents <- unique(extents)

# loop through surfaces and generate info for each country; can't stack as there's 43 unique extents

for(i in 1:length(raster_list)){
  
  # grab the raster for that iteration
  r_i <- raster_list[[i]]
  
  # mask it by the extent of the admin 0 raster
  r_i <- crop(r_i, admin0raster)
  
  # check number of columns match (cell size is equal)
  logic_cell <- r_i@ncols == admin0raster@ncols
  
  # if the number of cells isn't equal, resample 
  if(logic_cell == FALSE){
    
    message(paste0("Resampling raster ", i, " to align with master grid"))
    
    r_i <- resample(r_i, admin0raster, method = "bilinear")
    
  }
  
  # generate the mean pixel value for all pixels within a country
  poly_mean <- zonal(r_i, admin0raster, fun = 'mean')
  
  # create a dataframe with the zonal output
  poly_mean <- as.data.frame(poly_mean)
  
  names(poly_mean) <- c("country_id_raster",
                        names(raster_list[[i]]))
  
  # convert to a binary surface indicating values for a country 
  poly_mean[[2]] <- ifelse(is.na(poly_mean[[2]]), 0, 1)
  
  # combine results across each raster
  if(i == 1){
    
    combined <- poly_mean
    
  } else {
    
    combined <- merge(combined, poly_mean, by = "country_id_raster")
    
  }
  
}

# add back ISO code 
iso3 <- as.data.frame(admin0shapefile[, c(1:2, 5)])

names(combined)[1] <- "GAUL_CODE"
combined <- merge(combined, iso3, by = "GAUL_CODE")
# combined$COUNTRY_ID.x <- NULL
# drop countries which are XXX
combined_trim <- combined[!combined$COUNTRY_ID.y == "XXX", ]

# write out combined dataset
write.csv(combined,
          'D:/Dropbox/PhD/collaboration projects/Wellcome Trust data re-use/SHINY-map-prize/data/combined_lookup.csv',
          row.names = FALSE)


#shinyApp('ui.R', 'server.R')
# UI
#### load required libraries ####
if(!require(raster)){
  install.packages("raster")
  library(raster)
}

if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}

if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if(!require(malariaAtlas)){
  install.packages("malariaAtlas")
  library(malariaAtlas)
}

if(!require(shinydashboard)){
  install.packages("shinydashboard")
  library(shinydashboard)
}

if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}

if(!require(shinyalert)){
  install.packages("shinyalert")
  library(shinyalert)
}

if(!require(shinyBS)){
  install.packages("shinyBS")
  library(shinyBS)
}

if(!require(shinythemes)){
  install.packages("shinythemes")
  library(shinythemes)
}

if(!require(shinycssloaders)){
  install.packages("shinycssloaders")
  library(shinycssloaders)
}

if(!require(shinyjs)){
  install.packages("shinyjs")
  library(shinyjs)
}

if(!require(mapview)){
  install.packages("mapview")
  library(mapview)
}

if(!require(leaflet)){
  install.packages("leaflet")
  library(leaflet)
}

# generate a list of countries for which MAP data exists
# fix some encoding and country issues
load('data/sf_afr_simp_fao.rda')
sf_afr_simp <- sf_afr_simp[sf_afr_simp$COUNTRY_ID != "XXX",]
sf_afr_simp <- sf_afr_simp[sf_afr_simp$COUNTRY_ID != "MYT",]	
sf_afr_simp$name[sf_afr_simp$GAUL_CODE == "16840"] <- "Goh-Djiboua"
sf_afr_simp$name[sf_afr_simp$GAUL_CODE == "818"] <- "Extreme-Nord"
country_names <- sf_afr_simp$name[sf_afr_simp$ADMN_LEVEL==0]	
country_names <- country_names[country_names != "Hala'ib triangle"]	
country_names <- country_names[country_names != "Ma'tan al-Sarra"]
country_names <- country_names[country_names != "Ilemi triangle"]
country_names <- country_names[country_names != "Abyei"]
country_names <- country_names[country_names != "Cape Verde"]
country_names <- country_names[country_names != "Djibouti"]
country_names <- country_names[country_names != "Seychelles"]
country_names <- sort(country_names)
country_names[7] <- "Cote d'Ivoire"

# define a UI use a fluid bootstrap layout
appCSS <- "
#loading-content {
position: absolute;
background: #344151;
opacity: 1;
z-index: 100;
left: 0;
right: 0;
height: 100%;
text-align: center;
color: #FFFFFF;
}
"
# ui PART OF THE SHINNNY STARTS HERE
ui <- navbarPage(
  "PMI-Malaria Initiative",
  tabPanel("Application",
           fluidPage(theme = shinytheme("flatly"),
                     useShinyalert(),
                     useShinyalert(),
                     useShinyjs(),
                     inlineCSS(appCSS),
                     
                     # Loading message
                     div(
                       id = "loading-content",
                       h2("Loading Application...")
                     ),
                     
                     hidden(
                       div(
                         id = "app-content",
                         p(" ")
                       )
                     ),
                     
                     # set a margin for the checkbox
                     tags$head(
                       tags$style(
                         HTML(".checkbox-inline {
                              margin-left: 0px;
                              margin-right: 10px;}
                              .checkbox-inline+.checkbox-inline {
                              margin-left: 0px;
                              margin-right: 10px;
                              }"))),
                    
                     # create a sidebar where the user can select a country, and districts (etc.)
                     # we may change this to a header once basic functionality is resolved
                     sidebarLayout(
                       
                       # sidebar panel for the inputs
                       sidebarPanel(
                         uiOutput("tab"),
                         br(),
                         "More information can be found in the help tab.",
                         br(),
                         br(),
                         
                         # country selection (only one country allowed at a time)
                         selectInput("country", "Country :",
                                     choices = country_names,
                                     selected = "Benin"),
                         
                         # hover-over tooltip
                         bsTooltip(id = "country",
                                   title = "Please select the country of interest, available districts will update based on this selection.",
                                   placement = "right", trigger = "hover", options = list(container = "body")),
                         
                         # choose raster layers
                         checkboxGroupInput("selected_raster", "Data to show and compare :",
                                            choices = list("Malaria in children (Falciparum)" = "Plasmodium falciparum PR2 10",
                                                           "Malaria incidence (Falciparum)" = "Plasmodium falciparum Incidence",
                                                           "Insecticide Treated Net distribution" = "Insecticide treated bednet  ITN  coverage",
                                                           "Travel time to nearest city" = "A global map of travel time to cities to assess inequalities in accessibility in 2015"), 
                                            selected = "Plasmodium falciparum PR2 10"),
                         
                         helpText("First layer is shown in map, other layers included in 'Output'"),
                         
                         # dynamic district selection
                         uiOutput("select_dist"),
                         
                         # hover-over tooltip
                         bsTooltip(id = "select_dist",
                                   title = "Please select the districts to feature within the comparison/ranking.",
                                   placement = "right", trigger = "hover", options = list(container = "body")),
                         
                         # button to generate stats
                         actionButton(inputId = "processStats", label = "Generate Report", class='butt'),
                         tags$head(tags$style(".butt{margin-bottom:5px;}")),
                         
                         # hover-over tooltip
                         bsTooltip(id = "processStats",
                                   title = "Run generation of statistics and ranking system. This will produce results which feature in the tabs to the right.",
                                   placement = "right", trigger = "hover", options = list(container = "body")),
                         
                         # download report button (defined in server)
                         uiOutput("downloadbutton")
                       ),
                       
                       mainPanel(
                         
                         tabsetPanel(id='main0', type = "tabs",
                                     #tabPanel(value ='tab1', title = "Selected country and districts", div(style = 'overflow-y:scroll;height:750px;',plotOutput("select_country", height = '750px', width = '750px'))),
                                     tabPanel(value ='tab1', title = "Map", div(style = 'overflow-y:scroll;height:750px;',leafletOutput("mapview_country_raster", height = '750px', width = '750px'))),
                                     #tabPanel(value ='tab1', title = "Map", leafletOutput("mapview_country_raster")),
                                     tabPanel(value ='tab3', title = "Table", DT::dataTableOutput("activetable")),
                                     tabPanel(value ='tab2', title = "Output Report", div(style = 'overflow-y:scroll;height:750px;',htmlOutput("report"))))
                       ) # end of main panel
                     ) # end of fluid page # end of sidebar layout
                         ) 
                       ), # end of tab panel
  
  # help main panel
  tabPanel("Help",
           # sub-panels within the 'help' tab
           tabsetPanel(type = 'tabs',
                       tabPanel(title='Help', includeMarkdown('help.md')),
                       tabPanel(title='About', includeMarkdown('about.md'))))
                     )










# THIS will Return a Shinny app Object
# RENDER_STATS_RMD
---
  title: "Summary statistics"
output: md_document
params:
  stats_list: NA
dist_select: NA
country_select: NA
always_allow_html: yes
geometry: "left = 3cm, right = 3cm, top = 2cm, bottom = 2cm"
---
  
  ```{r load_packages, echo = FALSE}

require(kableExtra)
require(knitr)
require(RColorBrewer)
require(plotfunctions)
require(sf)

country_select <- params$country_select
dist_select <- params$dist_select

dist_names <- unique(dist_select$NAME)

# create a lookup table for the variable names
var_name_orig <- c("A global friction surface enumerating land based travel speed for a nominal year 2015",
                   "A global map of travel time to cities to assess inequalities in accessibility in 2015",
                   "All cause fever",                                                                      
                   "Artemisinin based combination therapy  ACT  coverage",                                 
                   "Dominant Vectors",                                                                     
                   "Duffy Negativity Phenotype Frequency",                                                 
                   "G6PD Deficiency Allele Frequency",                                                     
                   "HbS  Sickle Haemoglobin  Allele Frequency",                                            
                   "Indoor residual spraying  IRS  coverage",                                              
                   "Insecticide treated bednet  ITN  coverage",                                            
                   "Malaria attributable fever as a proportion of all cause fever",                        
                   "Malaria attributable fever as a proportion of malaria positive fever",                 
                   "Malaria positive fever as a proportion of all cause fever",                            
                   "Non malarial fever",                                                                   
                   "Plasmodium falciparum Incidence",                                                      
                   "Plasmodium falciparum PR2 10",                                                         
                   "Plasmodium falciparum Spatial Limits",                                                 
                   "Plasmodium falciparum Support",                                                        
                   "Plasmodium falciparum Temperature Suitability",                                        
                   "Plasmodium knowlesi Risk",                                                             
                   "Plasmodium vivax Endemicity",                                                          
                   "Plasmodium vivax Relapse",                                                             
                   "Plasmodium vivax Spatial Limits",                                                      
                   "Plasmodium vivax Support",                                                             
                   "Plasmodium vivax Temperature Suitability",                                             
                   "Secondary Dominant Vectors")

var_name_new <- c("Global friction surface (2015)",
                  "Global accessibility (2015)",
                  "All cause fever",                                                                      
                  "ACT coverage",                                 
                  "Dominant Vectors",                                                                     
                  "Duffy Negativity",                                                 
                  "G6PD Deficiency",                                                     
                  "HbS Sickle",                                            
                  "IRS coverage",                                              
                  "ITN coverage",                                            
                  "Malaria attributable fever (prop ACF)",                        
                  "Malaria attributable fever (prop malaria +ve fever)",                 
                  "Malaria positive fever (proportion ACF)",                            
                  "Non malarial fever",                                                                   
                  "Plasmodium falciparum Incidence",                                                      
                  "Plasmodium falciparum PR2 10",                                                         
                  "Plasmodium falciparum Spatial Limits",                                                 
                  "Plasmodium falciparum Support",                                                        
                  "Plasmodium falciparum temp. suitability",                                        
                  "Plasmodium knowlesi Risk",                                                             
                  "Plasmodium vivax Endemicity",                                                          
                  "Plasmodium vivax Relapse",                                                             
                  "Plasmodium vivax Spatial Limits",                                                      
                  "Plasmodium vivax Support",                                                             
                  "Plasmodium vivax temp. suitability",                                             
                  "Secondary Dominant Vectors")

var_name_units <- c("Land-based travel speed (minutes)",
                    "Time to travel to the nearest city (minutes)",
                    "All-cause fever prevalence within the limits of stable Pf \ntransmission, in children under 5 years of age", 
                    "Percentage of cases of fever in under-5 year olds that were \ntreated with ACT",                                 
                    "Probability (0-1)",                                                                     
                    "Proportion of the population predicted to be Duffy negative",                                                 
                    "Allele frequency (median) for G6PD deficiency",                                                     
                    "Allele frequency (mean) of sickle haemoglobin alleles",                                            
                    "Percentage of the population protected by IRS",                                              
                    "Proportion of the population who were protected by ITNs",                                            
                    "Predicted malaria attributable fevers (MAF) as a proportion \nof all fevers, for children under 5",                        
                    "Predicted malaria attributable fevers as a proportion of \nmalaria positive fevers, for children under 5",                 
                    "Predicted malaria-positive fevers as a proportion of all \nfevers, for children under 5",                            
                    "Non-malarial fever prevalence is defined as the sum of the \nprevalence of febrile illness without a P. falciparum malaria \ninfection and the prevalence of febrile illness coincident \nwith, but not caused by, a P. falciparum malaria infection, \nfor children under 5", 
                    "All-age incidence rate (clinical cases per 1,000 population \nper annum) of Plasmodium falciparum malaria",
                    "Age-standardised parasite rate for Plasmodium falciparum \nmalaria for children two to ten years of age (PfPR2-10)",
                    "Plasmodium falciparum Spatial Limits",                                                 
                    "Number of days per year which could support infectious \nvectors",                                                        
                    "The index of suitability is a relative, not absolute, \nmeasure of the effect of temperature on vectorial capacity \nand reproductive number",
                    "Modelled relative risk of zoonotic Plasmodium knowlesi \ninfecting humans on a scale 0 (low) to 1.0 (high)",
                    "Predicted all-age parasite prevalence rate for Plasmodium \nvivax malaria",                                                          
                    "Predicted prevalence of Plasmodium vivax relapse",                                                             
                    "Plasmodium vivax Spatial Limits",                                                      
                    "Number of days per year which could support infectious \nvectors",                                                             
                    "The index of suitability is a relative, not absolute, measure \nof the effect of temperature on vectorial capacity \nand reproductive number",
                    "Secondary Dominant Vectors")

var_lookup <- data.frame(original = var_name_orig,
                         new = var_name_new,
                         units = var_name_units)

```

```{r caption_css, results = "asis", echo = FALSE}
cat("
    <style>
    caption {
    color: #0dc5c1;
    font-weight: bold;
    font-size: 1.20em;
    }
    </style>
    ")

```

## MaDD - Malaria Data by District Output Report
#### This document contains summary statistics for `r country_select$name`. 

#### Summary statistics for each selected variable appear below, and were generated using the "zonal" function within the "raster" R package.

</br>
  
  ```{r pop_stat_1, results = "axis", echo = FALSE}
knitr::opts_chunk$set(echo = TRUE)

stats_list <- params$stats_list

dist_select <- params$dist_select

# subset to get specified table
table_one <- stats_list[[1]]

# get variable names
variable_n <- names(table_one)[2]
variable_n <- substr(variable_n, 1, nchar(variable_n)-7)

# get units and shorter variable name
idx <- which(variable_n == var_lookup$original)
units_n <- as.character(var_lookup$units[idx])
variable_trim_n <- as.character(var_lookup$new[idx])

names(table_one) <- c("District",
                      "Mean",
                      "Max",
                      "Min",
                      "Standard Deviation")

# round so that values aren't ridiculously long
table_one$Mean <- round(table_one$Mean, 2)
table_one$Max <- round(table_one$Max, 2)
table_one$Min <- round(table_one$Min, 2)
table_one$`Standard Deviation` <- round(table_one$`Standard Deviation`, 2)

# replace negative values as NA
table_one$Mean[table_one$Mean < 0] <- NA
table_one$Max[table_one$Max < 0] <- NA
table_one$Min[table_one$Min < 0] <- NA
table_one$`Standard Deviation`[table_one$`Standard Deviation` < 0] <- NA

table_one$Mean[!complete.cases(table_one)] <- NA
table_one$Max[!complete.cases(table_one)] <- NA
table_one$Min[!complete.cases(table_one)] <- NA
table_one$`Standard Deviation`[!complete.cases(table_one)] <- NA

kable(table_one,
      row.names = FALSE,
      caption = variable_n,
      align = "c") %>% 
  add_footnote(label = units_n) %>%
  kable_styling(c("striped", "bordered")) %>% 
  column_spec(1, bold = TRUE, width = "7em") %>%
  column_spec(2:5, width = "5em") %>%
  kable_styling(full_width = FALSE,
                position = "left")

```

</br>
  
  ```{r pop_stats_2, results = "axis", echo = FALSE}
if(length(stats_list) > 1){
  
  # subset to get specified table
  table_two <- stats_list[[2]]
  
  # get variable names
  variable_n2 <- names(table_two)[2]
  variable_n2 <- substr(variable_n2, 1, nchar(variable_n2)-7)
  
  # get units and shorter variable name
  idx <- which(variable_n2 == var_lookup$original)
  units_n2 <- as.character(var_lookup$units[idx])
  variable_trim_n2 <- as.character(var_lookup$new[idx])
  
  names(table_two) <- c("District",
                        "Mean",
                        "Max",
                        "Min",
                        "Standard Deviation")
  
  # round so that values aren't ridiculously long
  table_two$Mean <- round(table_two$Mean, 2)
  table_two$Max <- round(table_two$Max, 2)
  table_two$Min <- round(table_two$Min, 2)
  table_two$`Standard Deviation` <- round(table_two$`Standard Deviation`, 2)
  
  # replace negative values as NA
  table_two$Mean[table_two$Mean < 0] <- NA
  table_two$Max[table_two$Max < 0] <- NA
  table_two$Min[table_two$Min < 0] <- NA
  table_two$`Standard Deviation`[table_two$`Standard Deviation` < 0] <- NA
  
  table_two$Mean[!complete.cases(table_two)] <- NA
  table_two$Max[!complete.cases(table_two)] <- NA
  table_two$Min[!complete.cases(table_two)] <- NA
  table_two$`Standard Deviation`[!complete.cases(table_two)] <- NA
  
  kable(table_two,
        row.names = FALSE,
        caption = variable_n2,
        align = "c") %>% 
    add_footnote(label = units_n2) %>%
    kable_styling(c("striped", "bordered")) %>% 
    column_spec(1, bold = TRUE, width = "7em") %>%
    column_spec(2:5, width = "5em") %>%
    kable_styling(full_width = FALSE,
                  position = "left")
  
}

```

</br>
  
  ```{r pop_stats_3, results = "axis", echo = FALSE}
if(length(stats_list) > 2){
  
  # subset to get specified table
  table_three <- stats_list[[3]]
  
  # get variable names
  variable_n3 <- names(table_three)[2]
  variable_n3 <- substr(variable_n3, 1, nchar(variable_n3)-7)
  
  # get units and shorter variable name
  idx <- which(variable_n3 == var_lookup$original)
  units_n3 <- as.character(var_lookup$units[idx])
  variable_trim_n3 <- as.character(var_lookup$new[idx])
  
  names(table_three) <- c("District",
                          "Mean",
                          "Max",
                          "Min",
                          "Standard Deviation")
  
  # round so that values aren't ridiculously long
  table_three$Mean <- round(table_three$Mean, 2)
  table_three$Max <- round(table_three$Max, 2)
  table_three$Min <- round(table_three$Min, 2)
  table_three$`Standard Deviation` <- round(table_three$`Standard Deviation`, 2)
  
  # replace negative values as NA
  table_three$Mean[table_three$Mean < 0] <- NA
  table_three$Max[table_three$Max < 0] <- NA
  table_three$Min[table_three$Min < 0] <- NA
  table_three$`Standard Deviation`[table_three$`Standard Deviation` < 0] <- NA
  
  table_three$Mean[!complete.cases(table_three)] <- NA
  table_three$Max[!complete.cases(table_three)] <- NA
  table_three$Min[!complete.cases(table_three)] <- NA
  table_three$`Standard Deviation`[!complete.cases(table_three)] <- NA
  
  kable(table_three,
        row.names = FALSE,
        caption = variable_n3,
        align = "c") %>% 
    add_footnote(label = units_n3) %>%
    kable_styling(c("striped", "bordered")) %>% 
    column_spec(1, bold = TRUE, width = "7em") %>%
    column_spec(2:5, width = "5em") %>%
    kable_styling(full_width = FALSE,
                  position = "left")
  
}
```

</br>
  
  ```{r pop_stats_4, results = "axis", echo = FALSE}
if(length(stats_list) > 3){
  
  # subset to get specified table
  table_four <- stats_list[[4]]
  
  # get variable names
  variable_n4 <- names(table_four)[2]
  variable_n4 <- substr(variable_n4, 1, nchar(variable_n4)-7)
  
  # get units and shorter variable name
  idx <- which(variable_n4 == var_lookup$original)
  units_n4 <- as.character(var_lookup$units[idx])
  variable_trim_n4 <- as.character(var_lookup$new[idx])
  
  names(table_four) <- c("District",
                         "Mean",
                         "Max",
                         "Min",
                         "Standard Deviation")
  
  # round so that values aren't ridiculously long
  table_four$Mean <- round(table_four$Mean, 2)
  table_four$Max <- round(table_four$Max, 2)
  table_four$Min <- round(table_four$Min, 2)
  table_four$`Standard Deviation` <- round(table_four$`Standard Deviation`, 2)
  
  # replace negative values as NA
  table_four$Mean[table_four$Mean < 0] <- NA
  table_four$Max[table_four$Max < 0] <- NA
  table_four$Min[table_four$Min < 0] <- NA
  table_four$`Standard Deviation`[table_four$`Standard Deviation` < 0] <- NA
  
  table_four$Mean[!complete.cases(table_four)] <- NA
  table_four$Max[!complete.cases(table_four)] <- NA
  table_four$Min[!complete.cases(table_four)] <- NA
  table_four$`Standard Deviation`[!complete.cases(table_four)] <- NA
  
  kable(table_four,
        row.names = FALSE,
        caption = variable_n4,
        align = "c") %>% 
    add_footnote(label = units_n4) %>%
    kable_styling(c("striped", "bordered")) %>% 
    column_spec(1, bold = TRUE, width = "7em") %>%
    column_spec(2:5, width = "5em") %>%
    kable_styling(full_width = FALSE,
                  position = "left")
  
}
```

```{r stat_plots, echo = FALSE, results = "hide", fig.height = 10, fig.width = 10}
# generate shapefile for the plot
if(length(stats_list) == 1){
  
  # first create an index
  dist_index <- match(dist_select$name, table_one$District)
  dist_select_v1 <- dist_select
  dist_select_v1$mean <- table_one$Mean[dist_index]
  
  # create a bespoke palette
  colour_p_v1 <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(length(dist_select$mean))
  dist_select_v1$colours <- colour_p_v1[as.numeric(cut(dist_select_v1$mean, breaks = (length(dist_select_v1$mean))))]
  
  country_select <- params$country_select
  
  plot(sf::st_geometry(country_select),
       col = "#d9d9d9")
  
  plot(sf::st_geometry(dist_select_v1),
       col = dist_select_v1$colours,
       lty = 3,
       add = TRUE)
  
  title(variable_trim_n, 
        cex.main = 1.2, col.main = "#0dc5c1",
        adj = 0.5)
  
  gradientLegend(valRange = c(min(dist_select_v1$mean, na.rm = TRUE), max(dist_select_v1$mean, na.rm = TRUE)), 
                 color = colour_p_v1, 
                 nCol = length(dist_select_v1$mean), 
                 pos = 0.5, 
                 side = 4,
                 length = 0.25, 
                 depth = 0.05, 
                 inside = TRUE, 
                 coords = FALSE,
                 pos.num = NULL, 
                 n.seg = 3, 
                 border.col = "black", 
                 dec = 2,
                 fit.margin = TRUE)
  
}

if(length(stats_list) == 2){
  
  dist_index <- match(dist_select$name, table_one$District)
  
  dist_select_v1 <- dist_select
  dist_select_v2 <- dist_select
  dist_select_v1$mean <- table_one$Mean[dist_index]
  dist_select_v2$mean <- table_two$Mean[dist_index]
  
  # create a bespoke palette
  colour_p_v1 <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(length(dist_select_v1$mean))
  colour_p_v2 <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(length(dist_select_v2$mean))
  dist_select_v1$colours <- colour_p_v1[as.numeric(cut(dist_select_v1$mean, breaks = (length(dist_select_v1$mean))))]
  dist_select_v2$colours <- colour_p_v2[as.numeric(cut(dist_select_v2$mean, breaks = (length(dist_select_v2$mean))))]
  
  # define plotting parameters
  par(mfrow = c(1,2),
      mar = c(5,4,4,2))
  
  country_select <- params$country_select
  
  # plot first stat
  plot(sf::st_geometry(country_select),
       col = "#d9d9d9")
  
  plot(sf::st_geometry(dist_select_v1),
       col = dist_select_v1$colours,
       lty = 3,
       add = TRUE)
  
  title(variable_trim_n, 
        cex.main = 1.2, col.main = "#0dc5c1",
        adj = 0)
  
  gradientLegend(valRange = c(min(dist_select_v1$mean, na.rm = TRUE), max(dist_select_v1$mean, na.rm = TRUE)), 
                 color = colour_p_v1, 
                 nCol = length(dist_select_v1$mean), 
                 pos = 0.5, 
                 side = 4,
                 length = 0.25, 
                 depth = 0.05, 
                 inside = TRUE, 
                 coords = FALSE,
                 pos.num = NULL, 
                 n.seg = 3, 
                 border.col = "black", 
                 dec = 2,
                 fit.margin = TRUE)
  
  # plot second stat
  plot(sf::st_geometry(country_select),
       col = "#d9d9d9")
  
  plot(sf::st_geometry(dist_select_v2),
       col = dist_select_v2$colours,
       lty = 3,
       add = TRUE)
  
  title(variable_trim_n2, 
        cex.main = 1.2, col.main = "#0dc5c1",
        adj = 0)
  
  gradientLegend(valRange = c(min(dist_select_v2$mean, na.rm = TRUE), max(dist_select_v2$mean, na.rm = TRUE)), 
                 color = colour_p_v2, 
                 nCol = length(dist_select_v2$mean), 
                 pos = 0.5, 
                 side = 4,
                 length = 0.25, 
                 depth = 0.05, 
                 inside = TRUE, 
                 coords = FALSE,
                 pos.num = NULL, 
                 n.seg = 3, 
                 border.col = "black", 
                 dec = 2,
                 fit.margin = TRUE)
  
}

if(length(stats_list) == 3){
  
  dist_index <- match(dist_select$name, table_one$District)
  
  dist_select_v1 <- dist_select
  dist_select_v2 <- dist_select
  dist_select_v3 <- dist_select
  dist_select_v1$mean <- table_one$Mean[dist_index]
  dist_select_v2$mean <- table_two$Mean[dist_index]
  dist_select_v3$mean <- table_three$Mean[dist_index]
  
  # create a bespoke palette
  colour_p_v1 <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(length(dist_select_v1$mean))
  colour_p_v2 <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(length(dist_select_v2$mean))
  colour_p_v3 <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(length(dist_select_v3$mean))
  dist_select_v1$colours <- colour_p_v1[as.numeric(cut(dist_select_v1$mean, breaks = (length(dist_select_v1$mean))))]
  dist_select_v2$colours <- colour_p_v2[as.numeric(cut(dist_select_v2$mean, breaks = (length(dist_select_v2$mean))))]
  dist_select_v3$colours <- colour_p_v3[as.numeric(cut(dist_select_v3$mean, breaks = (length(dist_select_v3$mean))))]
  
  # define plotting parameters
  par(mfrow = c(2,2),
      mar = c(5,4,4,2))
  
  country_select <- params$country_select
  
  # plot first stat
  plot(sf::st_geometry(country_select),
       col = "#d9d9d9")
  
  plot(sf::st_geometry(dist_select_v1),
       col = dist_select_v1$colours,
       lty = 3,
       add = TRUE)
  
  title(variable_trim_n, 
        cex.main = 1.2, col.main = "#0dc5c1",
        adj = 0)
  
  gradientLegend(valRange = c(min(dist_select_v1$mean, na.rm = TRUE), max(dist_select_v1$mean, na.rm = TRUE)), 
                 color = colour_p_v1, 
                 nCol = length(dist_select_v1$mean), 
                 pos = 0.5, 
                 side = 4,
                 length = 0.25, 
                 depth = 0.05, 
                 inside = TRUE, 
                 coords = FALSE,
                 pos.num = NULL, 
                 n.seg = 3, 
                 border.col = "black", 
                 dec = 2,
                 fit.margin = TRUE)
  
  # plot second stat
  plot(sf::st_geometry(country_select),
       col = "#d9d9d9")
  
  plot(sf::st_geometry(dist_select_v2),
       col = dist_select_v2$colours,
       lty = 3,
       add = TRUE)
  
  title(variable_trim_n2, 
        cex.main = 1.2, col.main = "#0dc5c1",
        adj = 0)
  
  gradientLegend(valRange = c(min(dist_select_v2$mean, na.rm = TRUE), max(dist_select_v2$mean, na.rm = TRUE)), 
                 color = colour_p_v2, 
                 nCol = length(dist_select_v2$mean), 
                 pos = 0.5, 
                 side = 4,
                 length = 0.25, 
                 depth = 0.05, 
                 inside = TRUE, 
                 coords = FALSE,
                 pos.num = NULL, 
                 n.seg = 3, 
                 border.col = "black", 
                 dec = 2,
                 fit.margin = TRUE)
  
  # plot third stat
  plot(sf::st_geometry(country_select),
       col = "#d9d9d9")
  
  plot(sf::st_geometry(dist_select_v3),
       col = dist_select_v3$colours,
       lty = 3,
       add = TRUE)
  
  title(variable_trim_n3, 
        cex.main = 1.2, col.main = "#0dc5c1",
        adj = 0)
  
  gradientLegend(valRange = c(min(dist_select_v3$mean, na.rm = TRUE), max(dist_select_v3$mean, na.rm = TRUE)), 
                 color = colour_p_v3, 
                 nCol = length(dist_select_v3$mean), 
                 pos = 0.5, 
                 side = 4,
                 length = 0.25, 
                 depth = 0.05, 
                 inside = TRUE, 
                 coords = FALSE,
                 pos.num = NULL, 
                 n.seg = 3, 
                 border.col = "black", 
                 dec = 2,
                 fit.margin = TRUE)
  
}

if(length(stats_list) == 4){
  
  dist_index <- match(dist_select$name, table_one$District)
  
  dist_select_v1 <- dist_select
  dist_select_v2 <- dist_select
  dist_select_v3 <- dist_select
  dist_select_v4 <- dist_select
  dist_select_v1$mean <- table_one$Mean[dist_index]
  dist_select_v2$mean <- table_two$Mean[dist_index]
  dist_select_v3$mean <- table_three$Mean[dist_index]
  dist_select_v4$mean <- table_four$Mean[dist_index]
  
  # create a bespoke palette
  colour_p_v1 <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(length(dist_select_v1$mean))
  colour_p_v2 <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(length(dist_select_v2$mean))
  colour_p_v3 <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(length(dist_select_v3$mean))
  colour_p_v4 <- colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))(length(dist_select_v4$mean))
  dist_select_v1$colours <- colour_p_v1[as.numeric(cut(dist_select_v1$mean, breaks = (length(dist_select_v1$mean))))]
  dist_select_v2$colours <- colour_p_v2[as.numeric(cut(dist_select_v2$mean, breaks = (length(dist_select_v2$mean))))]
  dist_select_v3$colours <- colour_p_v3[as.numeric(cut(dist_select_v3$mean, breaks = (length(dist_select_v3$mean))))]
  dist_select_v4$colours <- colour_p_v4[as.numeric(cut(dist_select_v4$mean, breaks = (length(dist_select_v4$mean))))]
  
  # define plotting parameters
  par(mfrow = c(2,2),
      mar = c(5,4,4,2))
  
  country_select <- params$country_select
  
  # plot first stat
  print(plot(sf::st_geometry(country_select),
             col = "#d9d9d9"))
  
  plot(sf::st_geometry(dist_select_v1),
       col = dist_select_v1$colours,
       lty = 3,
       add = TRUE)
  
  title(variable_trim_n, 
        cex.main = 1.2, col.main = "#0dc5c1",
        adj = 0.5)
  
  gradientLegend(valRange = c(min(dist_select_v1$mean, na.rm = TRUE), max(dist_select_v1$mean, na.rm = TRUE)), 
                 color = colour_p_v1, 
                 nCol = length(dist_select_v1$mean), 
                 pos = 0.5, 
                 side = 4,
                 length = 0.25, 
                 depth = 0.05, 
                 inside = TRUE, 
                 coords = FALSE,
                 pos.num = NULL, 
                 n.seg = 3, 
                 border.col = "black", 
                 dec = 2,
                 fit.margin = TRUE)
  
  # plot second stat
  plot(sf::st_geometry(country_select),
       col = "#d9d9d9")
  
  plot(sf::st_geometry(dist_select_v2),
       col = dist_select_v2$colours,
       lty = 3,
       add = TRUE)
  
  title(variable_trim_n2, 
        cex.main = 1.2, col.main = "#0dc5c1",
        adj = 0)
  
  gradientLegend(valRange = c(min(dist_select_v2$mean, na.rm = TRUE), max(dist_select_v2$mean, na.rm = TRUE)), 
                 color = colour_p_v2, 
                 nCol = length(dist_select_v2$mean), 
                 pos = 0.5, 
                 side = 4,
                 length = 0.25, 
                 depth = 0.05, 
                 inside = TRUE, 
                 coords = FALSE,
                 pos.num = NULL, 
                 n.seg = 3, 
                 border.col = "black", 
                 dec = 2,
                 fit.margin = TRUE)
  
  # plot third stat
  plot(sf::st_geometry(country_select),
       col = "#d9d9d9")
  
  plot(sf::st_geometry(dist_select_v3),
       col = dist_select_v3$colours,
       lty = 3,
       add = TRUE)
  
  title(variable_trim_n3, 
        cex.main = 1.2, col.main = "#0dc5c1",
        adj = 0)
  
  gradientLegend(valRange = c(min(dist_select_v3$mean, na.rm = TRUE), max(dist_select_v3$mean, na.rm = TRUE)), 
                 color = colour_p_v3, 
                 nCol = length(dist_select_v3$mean), 
                 pos = 0.5, 
                 side = 4,
                 length = 0.25, 
                 depth = 0.05, 
                 inside = TRUE, 
                 coords = FALSE,
                 pos.num = NULL, 
                 n.seg = 3, 
                 border.col = "black", 
                 dec = 2,
                 fit.margin = TRUE)
  
  # plot fourth stat
  plot(sf::st_geometry(country_select),
       col = "#d9d9d9")
  
  plot(sf::st_geometry(dist_select_v4),
       col = dist_select_v4$colours,
       lty = 3,
       add = TRUE)
  
  title(variable_trim_n4, 
        cex.main = 1.2, col.main = "#0dc5c1",
        adj = 0)
  
  gradientLegend(valRange = c(min(dist_select_v4$mean, na.rm = TRUE), max(dist_select_v4$mean, na.rm = TRUE)), 
                 color = colour_p_v4, 
                 nCol = length(dist_select_v4$mean), 
                 pos = 0.5, 
                 side = 4,
                 length = 0.25, 
                 depth = 0.05, 
                 inside = TRUE, 
                 coords = FALSE,
                 pos.num = NULL, 
                 n.seg = 3, 
                 border.col = "black", 
                 dec = 2,
                 fit.margin = TRUE)
  
}

```
# SERVER CODE STARTS HERE

#### load required libraries ####
if(!require(raster)){
  install.packages("raster")
  library(raster)
}

if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}

if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if(!require(malariaAtlas)){
  install.packages("malariaAtlas")
  library(malariaAtlas)
}

if(!require(shinydashboard)){
  install.packages("shinydashboard")
  library(shinydashboard)
}

if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}

if(!require(shinyalert)){
  install.packages("shinyalert")
  library(shinyalert)
}

if(!require(shinyBS)){
  install.packages("shinyBS")
  library(shinyBS)
}

if(!require(shinythemes)){
  install.packages("shinythemes")
  library(shinythemes)
}

if(!require(knitr)){
  install.packages("knitr")
  library(knitr)
}

if(!require(kableExtra)){
  install.packages("kableExtra")
  library(kableExtra)
}

if(!require(shinyjs)){
  install.packages("shinyjs")
  library(shinyjs)
}

if(!require(plotfunctions)){
  install.packages("plotfunctions")
  library(plotfunctions)
}

if(!require(sf)){
  install.packages("sf")
  library(sf)
}

if(!require(mapview)){
  #devtools::install_github("r-spatial/mapview@develop")
  install.packages("mapview")
  library(mapview)
}

if(!require(leaflet)){
  install.packages("leaflet")
  library(leaflet)
}
#for interactive table
if(!require(DT)){
  install.packages("DT")
  library(DT)
}


# read in MAP availability lookup table
lookup <- read.csv('data/combined_lookup.csv', sep = ',', check.names = FALSE)

# read in the processed data lookup table
lookup_processed <- read.csv('data/raster_stats_paths.csv', stringsAsFactors = FALSE)

# load simplified admin polygons
load('data/sf_afr_simp_fao.rda')
# correct some encoding issues
sf_afr_simp$name[sf_afr_simp$GAUL_CODE == "16840"] <- "Goh-Djiboua"
sf_afr_simp$name[sf_afr_simp$GAUL_CODE == "818"] <- "Extreme-Nord"
sf_afr_simp$name[sf_afr_simp$GAUL_CODE == "66"] <- "Cote d'Ivoire"

# raster layers for Africa downloaded, simplified and saved in download-rasters.r
load('data/rasters/pfpr2_10_2015.rda')
load('data/rasters/pf_incidence_2015.rda')
load('data/rasters/time_to_city_2015.rda')
load('data/rasters/itn_2015.rda')

# to allow checking for map changes
country_id_last <- FALSE
raster_id_last <- FALSE
lastmap <- FALSE

# get the country_id (e.g. CIV) for selected country name
get_country_id <- function(country_name) {
  
  country_name <- as.character(country_name)
  
  # fix for encoding issue
  if(country_name == "Cote d'Ivoire"){
    
    country_id <- "CIV"
    
  } else {
    
    country_id <- sf_afr_simp$COUNTRY_ID[sf_afr_simp$name==country_name]
    
  }
  
  country_id <- as.character(country_id)
  
}

# get district names for selected country id
get_dist_names <- function(country_id) {
  
  country_id <- as.character(country_id)
  dist_names <- sf_afr_simp$name[sf_afr_simp$COUNTRY_ID==country_id & sf_afr_simp$ADMN_LEVEL==1] 
  
}

# define the server logic
# SERVER CODES STARTS HERE
function(input, output, session) {
  
  # Hide the loading message when the rest of the server function has executed
  hide(id = "loading-content", anim = TRUE, animType = "fade")    
  show("app-content")
  
  url <- a("MAP Homepage", href="https://map.ox.ac.uk/")
  
  output$tab <- renderUI({
    
    tagList(" The PMI visuals is part of the Data Lake project which uses the ATLAS/MaDD project allows easy interaction for data provided by the Malaria Atlas Project:", url)
    
  })
  
  # create dynamic reactive list of districts, per input country
  output$select_dist <- renderUI({
    
    # get the country_id (e.g. CIV) for selected country name
    country_id <- get_country_id(input$country)
    
    # get names of districts in this country
    selected_dist <- get_dist_names(country_id)    
    
    checkboxGroupInput("selected_dist", "Select first-level administrative division (min 2):",
                       choices = selected_dist,
                       selected = selected_dist,
                       inline = TRUE)
  })
  
  # mapview interactive leaflet map plot
  output$mapview_country_raster <- renderLeaflet({
    
    # get the country_id (e.g. CIV) for selected country name
    country_id <- get_country_id(input$country)
    raster_id <- input$selected_raster[1]
    
    # exit function if country and layer haven't changed
    # to reduce waiting time for plot changes
    if(!is.null(raster_id) & !is.null(raster_id_last)){
      if(country_id==country_id_last & raster_id==raster_id_last)
        return(lastmap)      
    }
    
    # subset the country (includes districts)
    sf_cntry <- sf_afr_simp[sf_afr_simp$COUNTRY_ID==country_id & sf_afr_simp$ADMN_LEVEL==1,]
    
    # add country boundaries to the plot first
    m <- mapview(sf_cntry,
                 color = 'black',
                 lwd = 2,
                 legend = FALSE,
                 alpha.regions = 0, 
                 zcol = 'name')    
    
    # add the first selected raster to the plot
    if(!is.null(input$selected_raster)){
      
      switch(input$selected_raster[1],
             "Plasmodium falciparum Incidence" =  m <- m + mapView(pf_incidence_2015,
                                                                   col.regions = colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))),
             "Plasmodium falciparum PR2 10" = m <- m + mapView(pfpr2_10_2015,
                                                               col.regions = colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))),
             "Insecticide treated bednet  ITN  coverage" = m <- m + mapView(itn_2015,
                                                                            col.regions = colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu"))),
             # changed breaks to show more detail at the values in malaria countries
             "A global map of travel time to cities to assess inequalities in accessibility in 2015" = m <- m + mapview(time_to_city_2015, 
                                                                                                                        col.regions = colorRampPalette(brewer.pal(brewer.pal.info["YlGnBu",1], "YlGnBu")),
                                                                                                                        at = rev(c(0,60,120,180,240,400,600,800,1200,1600,2400,3200,6400,10000))))
    }
    
    # record current ids so can check if they change above
    # set global vars, possibly bad practice
    country_id_last <<- country_id
    raster_id_last <<- raster_id    
    
    # set extent of map to the selected country 
    bbox <- as.vector(sf::st_bbox(sf_cntry))      
    m <- leaflet::fitBounds(m@map, bbox[1], bbox[2], bbox[3], bbox[4])
    
    # save last map to a global var, possibly dodgy
    # to allow avoiding changing map if it hasn't changed
    lastmap <<- m
  })  
  
  ##########################################
  # interactive table to show district stats
  output$activetable = DT::renderDataTable({
    
    # get country and raster_ids
    # sets reactive dependence on these
    country_id <- get_country_id(input$country)
    raster_ids <- input$selected_raster
    
    #message("in activetable raster_ids :",raster_ids)
    
    sf_cntry <- sf_afr_simp[sf_afr_simp$COUNTRY_ID==country_id,]
    
    # subset selected districts
    sf_dist_select <- sf_cntry[sf_cntry$name %in% input$selected_dist,] 
    
    if (length(input$selected_raster) == 0){
      
      return(NULL)
      
    } else { 
      
      # code to get data copied from processStats
      
      # sub to get names aligned
      lookup_processed$surface_name <- str_replace_all(lookup_processed$surface_name, '\\.', ' ')
      
      # define empty vector to populate with stats tables
      stats_list <- NULL
      
      for(i in 1:length(input$selected_raster)){
        
        # grab csv with the stats
        # 1. get a path to the stats csv
        stats_i_idx <- which(lookup_processed$surface_name == input$selected_raster[[i]])
        stats_i_path <- lookup_processed$stats_path[stats_i_idx]
        
        # 2. read in the csv 
        stats_i <- read.csv(stats_i_path, stringsAsFactors = FALSE)
        
        # 3. subset csv to only retain selected_dist
        stats_i_sub <- stats_i[stats_i$zone %in% sf_dist_select$GAUL_CODE, ]
        
        # add a variable which is the district name
        index <- match(stats_i_sub$zone, sf_dist_select$GAUL_CODE)
        stats_i_sub$District <- sf_dist_select$name[index]
        sf_dist_select$mean <- stats_i$mean[index]
        
        # reorder stats_i_sub
        stats_i_sub <- stats_i_sub[c(6, 2:5)]
        names(stats_i_sub) <- c('District',
                                paste0(input$selected_raster[[i]], " (mean)"),
                                paste0(input$selected_raster[[i]], " (max)"),
                                paste0(input$selected_raster[[i]], " (min)"),
                                paste0(input$selected_raster[[i]], " (sd)"))
        
        #stats_list[[i]] <- stats_i_sub
        #andy new code to put just means from diff layers into one table
        #get district name from first layer
        if (i==1) stats_layer_means <- stats_i_sub[c(1)]
        #add means and ranks from other layers
        
        #BEWARE this relies on district names being the same and in same order
        # means
        stats_layer_means <- cbind(stats_layer_means, stats_i_sub[2])   
        # ranks
        # reverse priority for malaria & travel layers (high malaria or travel=high priority=low rank)
        #BEWARE we should create column in lookup_processed specifying whether a high value = high priority
        #if(length(grep("Plasmodium",input$selected_raster[[i]]))>0)
        #ITN is only layer where high=good and priority doesn't need reversing
        if(length(grep("ITN",input$selected_raster[[i]]))>0)
          priority <- rank(stats_i_sub[2])
        else
          priority <- rank(-stats_i_sub[2])
        
        
        stats_layer_means <- cbind(stats_layer_means, priority)  
        
        
      }      
    }
    
    #to go in the table  
    #stats_i_sub     
    #stats_layer_means
    
    #round data same as in reports, not first column
    stats_layer_means[-1] <- round(stats_layer_means[-1],2)
    
    # replace negative values as NA
    # try keeping in the negative values it exposes issues with the data that need to be addressed
    #stats_layer_means[stats_layer_means[] < 0] <- NA
    
    DT::datatable(stats_layer_means, 
                  rownames=FALSE, 
                  #fillContainer = FALSE,
                  options = list(paging=FALSE))
    
  })    
  
  # observeEvent for "processStats"
  observeEvent(input$processStats, {
    
    output$downloadbutton <- renderUI({
      downloadButton('download', 'Download File ...')
    })
    
    # check for district selection inputs   
    show("download")
    
    # check for max four inputs   
    if(length(input$selected_dist) < 2){
      
      shinyalert("Oops!", "Please select at least 2 districts to compare", type = "warning")
      
    }    
    
    # check for max four variable inputs   
    else if (length(input$selected_raster) == 0){
      
      shinyalert("Oops!", "Please select a raster", type = "warning")
      
    } else {
      
      updateTabsetPanel(session=session, inputId = 'main0', selected = 'tab2')
      
      # sub to get names aligned
      lookup_processed$surface_name <- str_replace_all(lookup_processed$surface_name, '\\.', ' ')
      
      # define empty vector to populate with stats tables
      stats_list <- NULL
      
      # get the country_id (e.g. CIV) for selected country name
      country_id <- get_country_id(input$country)
      
      sf_cntry <- sf_afr_simp[sf_afr_simp$COUNTRY_ID==country_id,]
      
      # subset selected districts
      sf_dist_select <- sf_cntry[sf_cntry$name %in% input$selected_dist,] 
      
      for(i in 1:length(input$selected_raster)){
        
        # grab csv with the stats
        # 1. get a path to the stats csv
        stats_i_idx <- which(lookup_processed$surface_name == input$selected_raster[[i]])
        stats_i_path <- lookup_processed$stats_path[stats_i_idx]
        
        # 2. read in the csv 
        stats_i <- read.csv(stats_i_path, stringsAsFactors = FALSE)
        
        # 3. subset csv to only retain selected_dist
        stats_i_sub <- stats_i[stats_i$zone %in% sf_dist_select$GAUL_CODE, ]
        
        # add a variable which is the district name
        index <- match(stats_i_sub$zone, sf_dist_select$GAUL_CODE)
        stats_i_sub$District <- sf_dist_select$name[index]
        sf_dist_select$mean <- stats_i$mean[index]
        
        # reorder stats_i_sub
        stats_i_sub <- stats_i_sub[c(6, 2:5)]
        names(stats_i_sub) <- c('District',
                                paste0(input$selected_raster[[i]], " (mean)"),
                                paste0(input$selected_raster[[i]], " (max)"),
                                paste0(input$selected_raster[[i]], " (min)"),
                                paste0(input$selected_raster[[i]], " (sd)"))
        
        stats_list[[i]] <- stats_i_sub
        
      }
      
      # generate tables as a markdown
      # copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "render_stats.rmd")
      file.copy("render_stats.rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(stats_list = stats_list,
                     dist_select = sf_dist_select,
                     country_select = sf_cntry)
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport,
                        output_file = file.path(tempdir(), "pop_stats.rmd"),
                        output_format = "md_document",
                        params = params,
                        envir = globalenv())
      
      getPage <- function() {
        
        return(includeMarkdown(file.path(tempdir(), "pop_stats.rmd")))
        
      }
      
      output$report <- renderUI({getPage()})
      
    }
  }
  )
  # download generated statistics as a html document
  output$download <- downloadHandler(
    
    filename = paste0("MAP_output_", Sys.Date(), ".html"),
    
    content = content <- function(file) {
      
      rmarkdown::render(input = file.path(tempdir(), "pop_stats.rmd"),
                        output_file = paste0(tempdir(), "/MAP_output_", Sys.Date(), ".html"),
                        output_format = "html_document")
      
      file.copy(paste0(tempdir(), "/MAP_output_", Sys.Date(), ".html"), file)
      
    },
    
    contentType = "text/html")
}
