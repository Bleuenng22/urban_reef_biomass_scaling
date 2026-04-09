
###########################################
# libraries 
###########################################

library(tidyverse)
library(brms)
library(tidybayes)
library(rnaturalearth)
library(rnaturalearthdata)




cities_circle <- read.csv2('data_proc/urban_biom_vol.csv')
cities_circle$continent <- cities_circle$continent %>% #problem with north america as NA
  replace_na("NA")





###########################################
# Modelisation
###########################################

## Model implementation ##
fmod <- brms::bf(total_biomass ~ log(volume_per_unit_area) + (1 | continent/country/city), family = Gamma(link='log'))
mod <- brms::brm(fmod, data = cities_circle, iter = 3000, chains = 6, seed = 123)
#saveRDS(mod, file = "mod.rds")
#mod <- readRDS("outputs/mod.rds")


## Outputs ##
summary(mod, digits = 2)
brms::pp_check(mod, ndraws = 50)
fixef(mod)
tidybayes::median_hdci(brms::bayes_R2(mod, summary = FALSE))



###########################################
# general plots - trends
###########################################

ggplot(cities_circle)+
  geom_point(aes(volume_per_unit_area, total_biomass,
                 col = continent)) + theme_bw() +
  scale_color_viridis_d()#+geom_text(aes(label = City_centers$City), hjust = 1.1, cex = 2.5)


ggplot(cities_circle)+
  geom_point(aes(log(volume_per_unit_area), log(total_biomass),
                 col=continent)) + theme_bw() +
  scale_color_viridis_d()#+geom_text(aes(label = City_centers$City), hjust = 1.1, cex = 2.5)


## Average metric per city
cities_grouped <- cities_circle %>% 
  filter(total_biomass != 0) %>% 
  group_by(city) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
            country = unique(country),
            continent = unique(continent))

ggplot(cities_grouped, aes(volume_per_unit_area, total_biomass, col = country))+
  geom_point() + theme_bw() + scale_color_viridis_d() +
  geom_text(aes(label = city), hjust = 1.1, cex = 2.5) +
  theme(legend.position="none")

#select the proportion of city names to display
donnees_echantillon <- cities_grouped[seq(1, nrow(cities_grouped), by = 3), ]

ggplot(cities_grouped, aes(log(volume_per_unit_area), log(total_biomass), col = continent))+
  geom_point() + theme_bw() + scale_color_viridis_d() +
  geom_text(
    data = donnees_echantillon,
    aes(label = city),
    hjust = 1.1,
    size = 2.5
  )#+
  #theme(legend.position="none")





###########################################
# Map of the 400 cities
###########################################

land <- ne_countries(scale = 110, returnclass = "sf")
land_vect <- vect(land)
land_sf <- st_as_sf(land_vect)

city_centers <- read.csv2('data/city_centers.csv')
city_centers_sf <- city_centers %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

ext(city_centers_sf)
ggplot() +
  geom_sf(data = land_sf, fill = "grey90", color = "grey40") +
  geom_sf(data = city_centers_sf, color = "blue", size = 1) +
  coord_sf(xlim = c( -140, 160), ylim = c(-50, 70), expand = FALSE) +
  theme_minimal()


