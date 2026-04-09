library(terra)
library(sf)
library(tidyverse)
library(tibble)
library(rnaturalearth)
library(rnaturalearthdata)


################## DATA IMPORTATION ###############################

#Initial list of city - INDEX (C1, C2, C3...), NAME AND LON/LAT OF THE CENTER POINT
city_centers <- read.csv2('data_raw/city_centers.csv')
city_centers <- city_centers %>%  #remplace North America with 'NA'
  mutate(continent = replace_na(continent, "NA"))

city_centers_vect <- city_centers %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  vect()


##Building height from copernicus - 100m resolution
building_height <- rast('data_raw/building_height.tif')

#Land boundaries and mollweide transformation
land <- vect(ne_download(scale=10, type="land", category="physical", returnclass="sf"))
crs_moll <- "+proj=moll +datum=WGS84 +units=m +no_defs"
land_moll <- project(land, crs_moll)

#Human density mosaic, downloaded on copernicus as 2025 worldwide data with a 100m resolution
human_density <- rast('data_raw/human_density.tif')

#Human biomass data
metrics_data <- read.csv2('data_raw/metrics_data.csv') #compilation of several data
children_body_mass <- read.csv2('data_raw/children_body_mass.csv')

#Total biomass - urban wildlife
urban_biomass <- read.csv2('data_raw/urban_biomass.csv')
urban_animal_body_mass <- read.csv2('data_raw/urban_animal_body_mass.csv')


################### FUNCTIONS TO CALCULATE ALL METRICS ####################################

#Create a buffer zone around the geographic city center to find the highest density point
bbox_square <- function(city_point, buff){
  bbox <- st_bbox(city_point)
  bbox[c("xmin", "xmax")] <- c(bbox["xmin"] - buff, bbox["xmin"] + buff)
  bbox[c("ymin", "ymax")] <- c(bbox["ymin"] - buff, bbox["ymin"] + buff)
  return(st_as_sfc(st_bbox(bbox, crs = st_crs(point))))
}


#Function to extract complexity metrics from building heights with resolution in meters
GEE_complexity_city <- function(city) {

  mat <- as.matrix(city)/1000 #convert into matrix and building height into kilometers
  res_km <- res(city)/1000 #convert resolution into kilometers
  
  nb_values <- as.numeric(sum(!is.na(mat)))
  planar_area <- nb_values * (res_km[1]*res_km[2])
  na_surface <- as.numeric(sum(is.na(mat))) * (res_km[1]*res_km[2]) #surface with NA
  mat_na <- replace(mat, is.na(mat), 0) # replace na by 0 and remove the na _ surface
  surface_area <- sum(sp:::surfaceArea.matrix(mat_na, cellx = res_km[1], #sum of all cell surface area without NA 
                                              celly = res_km[2], 
                                              byCell = TRUE), na.rm = TRUE) - na_surface
  
  return(data.frame(
    planar_area = planar_area,
    surface_area = surface_area,    
    rugosity = surface_area / planar_area,    
    volume_per_unit_area = sum(mat, na.rm = TRUE) / planar_area,
    height_range = diff(range(mat, na.rm = TRUE))))}


#Function to extract human population metrics
human_population_metrics <- function(density){
  mat <- as.matrix(density)
  return(sum(mat, na.rm = TRUE))}


#Function that calculate human biomass, take 1 data line in enter
total_human_biomass <- function(city_metric){
  country <- city_metric$country
  metric_data <- metrics_data[metrics_data$code == country,]
  if (nrow(metric_data) == 0) {return(NA)}
  
  human_biomass <- 
    as.numeric(city_metric$pop_abundance) * as.numeric(metric_data$X0.4)/100 * as.numeric(children_body_mass$body_mass[1]) +
    as.numeric(city_metric$pop_abundance) * as.numeric(metric_data$X05.oct)/100 * as.numeric(children_body_mass$body_mass[2]) +
    as.numeric(city_metric$pop_abundance) * as.numeric(metric_data$nov.17)/100 * as.numeric(children_body_mass$body_mass[3]) +
    as.numeric(city_metric$pop_abundance) * as.numeric(metric_data$X18.)/100 * as.numeric(metric_data$sex_ratio)/(as.numeric(metric_data$sex_ratio)+100) * as.numeric(metric_data$BMI_men) * as.numeric(metric_data$mean_men_height_cm) +
    as.numeric(city_metric$pop_abundance) * as.numeric(metric_data$X18.)/100 * (1-(as.numeric(metric_data$sex_ratio)/(as.numeric(metric_data$sex_ratio)+100))) * as.numeric(metric_data$BMI_women) * as.numeric(metric_data$mean_women_height_cm)
  
  return(human_biomass)
}


#Function that calculate total biomass as human biomass + urban wildlife biomass
animal_biomass <- function(data){

  country <- data$country
  urban_biomass_subset <- urban_biomass[urban_biomass$code == country,] 
  if (nrow(urban_biomass_subset) == 0){return(NA)}
  
  total_biomass <- data$human_biomass + as.numeric(data$pop_abundance) * 
    (as.numeric(urban_biomass_subset$pigeon_per_capita) * as.numeric(urban_animal_body_mass$body_mass[1])+
       as.numeric(urban_biomass_subset$cat_per_capita) * as.numeric(urban_animal_body_mass$body_mass[3])+
       as.numeric(urban_biomass_subset$dog_per_capita) * as.numeric(urban_animal_body_mass$body_mass[4])) +
    as.numeric(data$pop_abundance) * as.numeric(urban_animal_body_mass$body_mass[2])
  
  return(total_biomass)
}




#############  Crop and calculate metrics ####################

#exception for which central point must be the geographic center point
except <- c('C23', 'C32','C35', 'C48', 'C69', 'C80', 'C84', 'C111', 'C120', 'C128', 'C155', 'C170', 'C201', 'C247', 'C311', 'C345', 'C365', 'C399')

cities_samples_complexity_population <- function(city){
  
  #City center as the highest density point
  city_point <- city_centers_vect[city_centers_vect$city_idx == city]
  city_point_m <- project(city_point, crs_moll)
  
  if (city %in% except){city_cent <- city_point_m}
  else{
    city_dens <- crop(human_density, bbox_square(city_point_m, buff=10000))
    city_center <- xyFromCell(city_dens, which.max(values(city_dens)))
    city_cent <- vect(city_center, type = "points", crs = crs_moll)
  }
  
  #define a bigger box around the highest density point
  box <- bbox_square(city_cent, buff=25000)
  city_land_full <- crop(building_height, box)
  city_land <- crop(city_land_full, land_moll, mask = TRUE)
  #plot(city_land)
  #plot(city_point_m, add=TRUE, col='red')
  #plot(city_cent, add=TRUE, col='yellow')
  
  #create the grid of 100 cells of 2kms side
  grid_box <- st_as_sfc(st_bbox(st_buffer(st_as_sf(city_cent), dist = 10000)))
  grid <- st_as_sf(st_make_grid(grid_box, cellsize = 2000, square = TRUE))
  #plot(grid, add=TRUE)
  
  #Same transformation for human density data
  city_dens_full <- crop(human_density, box)
  city_dens <- crop(city_dens_full, land_moll, mask = TRUE)
  #plot(city_dens)
  
  #Random grid cell sampling, crop same 25 cells in building height and human density simultaneously
  cropped_building <- list()
  cropped_density <- list()
  index <- sample(1:nrow(grid), 100) #random order of sampling
  
  i <- 1
  j <- 1
  
  while (length(cropped_building) < 15) {
    n <- index[i]
    vals <- crop(city_land, grid[n,])
    
    if (all(!is.na(values(vals))) & global(vals == 0, "sum") < length(values(vals))/2) { #keep the sample only is there is any NA and those with less than half values with 0
      cropped_building[[j]] <- vals
      cropped_density[[j]] <- crop(city_dens, grid[n,])
      j <- j+1}
    i <- i + 1}
  
  ########### Once 25 samples cropped - metrics calculation ###########
  
  #Initial data with city and sample information
  data1 <- data.frame(
    idx = values(city_point)[1],
    city = values(city_point)[2],
    country = values(city_point)[3],
    continent = values(city_point)[4],
    sample = seq(1:15))
  
  #Complexity metrics calculations
  data2 <- lapply(cropped_building, GEE_complexity_city)
  data2 <- do.call(rbind, data2)
  
  #Human population density calculation
  data3 <- lapply(cropped_density, human_population_metrics)
  pop_abundance <- do.call(rbind, data3)
  
  data <- cbind(data1, data2, pop_abundance)
  
  #Biomass metrics calculation
  data$human_biomass <- total_human_biomass(data)
  data$total_biomass <- animal_biomass(data)
  
  return(data)}





################### RUN #################

#Initial data
CITIES <- city_centers$city_idx

#run the function
final_data_list <- lapply(seq_along(CITIES), function(i) {
  message("Processing city index: ", i)
  cities_samples_complexity_population(CITIES[[i]])
})

#Final data
final_data <- do.call(rbind, final_data_list)
#write.csv2(final_data, 'results/cities_metrics.csv')
