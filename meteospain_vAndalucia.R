## this script gets data from meteospain https://emf.creaf.cat/news/meteospain_available/ 

library(meteospain)
library(ggplot2)
library(sf)
library(dplyr)

# function to translate coordinates from arc min/sec to decimals
DMS2decimal = function(geometry) {
  g = st_coordinates(geometry)
  print(g)
  lon = paste0(as.character(g[1]), "00")
  lat = paste0(as.character(g[2]), "00")
  lon_deg = as.numeric(substring(lon, 1,2))
  lon_min = as.numeric(substring(lon, 4,5))
  lon_sec = as.numeric(substring(lon, 6,7))
  lat_deg = as.numeric(substring(lat, 1,2))
  lat_min = as.numeric(substring(lat, 4,5))
  lat_sec = as.numeric(substring(lat, 6,7))
  p_lon = lon_deg - lon_min/60 - lon_sec/3600
  p_lat = lat_deg + lat_min/60 + lat_sec/3600
  return(st_sfc(st_point(c(p_lon, p_lat)), crs = 4326))
}

# retrieve meteo data from april-may 2022
ria_options <- ria_options(resolution = 'daily', start_date = as.Date('2022-04-14'), end_date = as.Date('2022-05-14'))
meteo_ria <- get_meteo_from('ria', ria_options)

#check which crs is used
st_crs(meteo_ria)

# check unique number of stations
num_stations <- length(unique(meteo_ria$station_id)) # 101

# check unique number of geometries
length(unique(meteo_ria$geometry)) # matches 101

# get total no. of observations
num_obs <- length(meteo_ria$timestamp)

# calculate number of measurements per station
num_perstation <- num_obs / num_stations # value of num_perstation is 31, so there are 31 measurements per station, so 1 every day of the month

# check which columns contain NA values
colSums(is.na(meteo_ria))

# now we want to calculate margin of error for each station, first create new df 
stations <- data.frame(matrix(ncol = num_stations, nrow = round(num_perstation))) # use 31 rows because 31 days, and no. of columns is the amount of unique stations
names <- unique(meteo_ria$station_name) # get the unique names of the stations
colnames(stations) <- names # name the columns after the stations

# now per station, collect the mean temperature values
for (name in colnames(stations) ) { # iterate over each station
  vals <- c()
  for (i in 1:length(meteo_ria$timestamp)) {
    if (meteo_ria$station_name[i] == name) {
      vals <- c(vals, meteo_ria$mean_temperature[i])
      }
    }
  # if there are less than 31 (num of observations per station) observations, add NA values to the df to get 31 rows (needs 31 rows to add it to the df)
  if (length(vals) < round(num_perstation)) {
    print(paste0("The following station did not include all observations:", name))
    m <- round(num_perstation) - length(vals) # this is the amount of missing values
    while (m > 0) {
      vals <- c(vals, NA)
      m = m-1
      }
    }
  stations[name] <- vals # add collected values to the corresponding column in the df
}

stations$timestamp <- unique(meteo_ria$timestamp) # get the 31 timestamps used
stations <- relocate(stations, timestamp) # move stations to the start of the df

summary(stations)
colSums(is.na(stations))

# now calculate some statistics per station (assuming normal t-distrib)
n <- length(stations$timestamp) # this is the number of observations
crit_val <- qt(0.975, df=n-1) # get quantile of t distrib, with n-1 degrees of freedom

temp_spain <- data.frame(matrix(ncol = 3, nrow = 0)) # new df with columns for name, geom, temp, margin of error, and each row is a station

for (name in names) {
  # calculate mean temperature per station from the 31 observations
  mean_temperature <- mean((unlist(stations[name])))
  
  # calculate st dev for t distrib and correspondingly the margin of error using the quantile
  std_dev <- sd((unlist(stations[name])))/sqrt(n) 
  temp_moe <- std_dev * crit_val 
  
  # save calculations
  results <- c(name, mean_temperature, temp_moe)
  temp_spain <- rbind(temp_spain, results)
}

colnames(temp_spain) <- c("station_name", "mean_temp", "temp_moe")

# now adding the geometry info to the data
info <- get_stations_info_from('ria', ria_options) # get the geometry info
info <- select(info, c(station_name, geometry)) # keep only the geometry and name columns of the data
info <- info %>% distinct(station_name, .keep_all = TRUE) # get rid of duplicate stations

alldata <- left_join(temp_spain, info, by = "station_name") # add the geometry info as a new column to the temp data, matching station names with geometry

for(i in 1:length(alldata$geometry)){
  alldata$geometry[i] <- DMS2decimal(alldata$geometry[i])
}

### now plotting

library(Vizumap)
library(sf)
library(cowplot)

# use one of four pre-prepared colour palettes
cmBivPal <- build_palette(name = "GreenBlue")
view(cmBivPal)

# below, read.uv creates a df that is usable to build a map with the function after
temperature <- read.uv(data = alldata, estimate = "mean_temp", error = "temp_moe")

# convert geo info to sf object
geodata <- st_as_sf(get_stations_info_from('ria', ria_options))

# convert coordinates into decimals
for(i in 1:length(geodata$geometry)){
  geodata$geometry[i] <- DMS2decimal(geodata$geometry[i])
}

# create buffer around the points so they can be used as spatial polygons
geodata <- st_buffer(geodata, dist = 3000)

# convert into S4 object
geodata <- as_Spatial(geodata)

# create map using the biv color palette and the poverty df
tempBivMap <- build_bmap(data = temperature, geoData = geodata, id = "station_name")
map <- view(tempBivMap) + geom_sf(data = andalucia$geometry, color=alpha("black", 0.7), fill = NA) + ggtitle("Meteorological stations in Andalucia, Spain")

# create legend
tempBivKey <- build_bkey(data = temperature)
legend <- view(tempBivKey)

# plot map with legend, plot_grid is from cowplot library to plot things next to each other
plot_grid(map, legend, labels = NULL, scale = c(1, 0.5)) 

# read in spain-andalucia shapefile
andalucia <- st_read("/home/merel/Documents/I-CISK/uncertainty/Meteo_data_Spain_Andalucia/Andalucia_regions/13_23_DemarcacionCEPS.shp") 

# check the projection
st_crs(andalucia) # so its ETRS89, need to convert to WGS84 to match the weather stations crs
andalucia <- st_transform(andalucia, crs = st_crs(meteo_ria)) # take the crs from meteo ria to transform the data into


# now translate meteo ria geometries to decimals in stead of arc seconds
for(i in 1:length(meteo_ria$geometry)){+
  meteo_ria$geometry[i] <- DMS2decimal(meteo_ria$geometry[i])
}

# have a look at the stations in context of Andalucia
ggplot() + geom_sf(data = andalucia$geometry) + geom_sf(data = meteo_ria$geometry) + ggtitle("Meteorological stations in Andalucia, Spain")

View(st_contains(unique(st_as_sf(meteo_ria$geometry)), st_as_sf(andalucia$geometry)))

a <- aggregate(meteo_ria$mean_temperature, by = andalucia, FUN = mean)
