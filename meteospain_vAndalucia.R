# This script gets meteorological data from the Andalucia region with package meteospain: https://emf.creaf.cat/news/meteospain_available/ 
# The data is subsequently plotted along with uncertainty measures using package Vizumap. 

### functions section ----------------------------------

# function to translate coordinates from arc min/sec to decimals (meteo ria provides geometries of stations using arcseconds in stead of normal decimals)
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

### data retrieval and wrangling section ----------------------------------

library(meteospain)
library(sf)
library(dplyr)

# retrieve meteo data for Andalucia from april-may 2022
ria_options <- ria_options(resolution = 'daily', start_date = as.Date('2022-04-01'), end_date = as.Date('2022-05-01'))
meteo_ria <- get_meteo_from('ria', ria_options)

# now translate meteo ria geometries to decimals instead of arc seconds
for(i in 1:length(meteo_ria$geometry)){
    meteo_ria$geometry[i] <- DMS2decimal(meteo_ria$geometry[i])
}

#check which crs is used
st_crs(meteo_ria)

# read in spain-andalucia shapefile
andalucia <- st_read("/home/merel/Documents/I-CISK/uncertainty/Meteo_data_Spain_Andalucia/Andalucia_regions/13_23_DemarcacionCEPS.shp") 

# check the projection
st_crs(andalucia) # so its ETRS89, need to convert to WGS84 to match the weather stations crs
andalucia <- st_transform(andalucia, crs = st_crs(meteo_ria)) # take the crs from meteo ria to transform the data into

# check unique number of stations in data retrieved from meteospain
num_stations <- length(unique(meteo_ria$station_id)) # 101

# check unique number of geometries
length(unique(meteo_ria$geometry)) # matches 101

# get total no. of observations
num_obs <- length(meteo_ria$timestamp)

# calculate number of measurements per station
num_perstation <- num_obs / num_stations # 1 measurement for each station per day

# check which columns contain NA values
colSums(is.na(meteo_ria))

# now we want to calculate margin of error for each station -- create new df for this first
stations <- data.frame(matrix(ncol = num_stations, nrow = round(num_perstation))) # use 31 rows because 31 days of observations, and no. of columns is the amount of unique stations
names <- unique(meteo_ria$station_name) # get the unique names of the stations
colnames(stations) <- names # name the columns after the stations

# now per station, collect the mean temperature values and add them to the df created above
for (name in colnames(stations) ) { # iterate over each station
  vals <- c()
  for (i in 1:length(meteo_ria$timestamp)) {
    if (meteo_ria$station_name[i] == name) {
      vals <- c(vals, meteo_ria$mean_temperature[i]) # add mean temperature values to vector
      }
    }
  # if there are less than 31 (num of observations per station) observations collected, add NA values to the df to get 31 rows (so dims match)
  if (length(vals) < round(num_perstation)) {
    m <- round(num_perstation) - length(vals) # this is the amount of missing values
    print(paste0("The station", name, "missed", m, "observations"))
    while (m > 0) {
      vals <- c(vals, NA)
      m = m-1
      }
    }
  stations[name] <- vals # save collected temperature values to the corresponding column in the df
}

# add timestamp to the df
stations$timestamp <- unique(meteo_ria$timestamp)
stations <- relocate(stations, timestamp) # move timestamp to the 1st column of the df

# inspect values
summary(stations)
colSums(is.na(stations))

# now calculate some statistics per station (assuming normal t-distrib)
n <- length(stations$timestamp) # this is the number of observations
crit_val <- qt(0.975, df=n-1) # get quantile of t distrib, with n-1 degrees of freedom

temp_spain <- data.frame(matrix(ncol = 3, nrow = 0)) # new df with columns for name, geom, temp, margin of error -- each row represents one weather station

for (name in names) {
  # calculate mean temperature per station from the 31 observations
  mean_temperature <- mean((unlist(stations[name])))
  
  # calculate st dev for t distrib and correspondingly the margin of error using the quantile
  std_dev <- sd((unlist(stations[name])))/sqrt(n) 
  temp_moe <- std_dev * crit_val 
  
  # calculate no. of st devs difference from the mean
  
  # save calculations
  results <- c(name, mean_temperature, temp_moe)
  temp_spain <- rbind(temp_spain, results)
}

colnames(temp_spain) <- c("station_name", "mean_temp", "temp_moe")

# now adding the geometry info to the data
info <- get_stations_info_from('ria', ria_options) # get the geometry info
info <- select(info, c(station_id, station_name, geometry)) # keep only the geometry and name columns of the data

temp_spain <- left_join(temp_spain, info, by = "station_name") # add the geometry info as a new column to the temp data, matching station names with geometry

for(i in 1:length(temp_spain$geometry)){
  temp_spain$geometry[i] <- DMS2decimal(temp_spain$geometry[i])
}

# now get all data since 2001 and get all April values per station

ria_options2 <- ria_options(resolution = 'monthly', start_date = as.Date('2001-04-01'), end_date = as.Date('2022-05-01'))
meteo_ria_long <- get_meteo_from('ria', ria_options2) # this takes a few mins to run
meteoSubsetApril <- meteo_ria_long[grep(".....04-01", meteo_ria_long$timestamp), ] # select only the observations from April (the dots are so it only looks at the month and not the year)

# TODO 
a <- aggregate(mean_temperature~station_id+station_name, meteoSubsetApril, FUN = mean) # this is the mean for April per station

### plotting section ----------------------------------
library(ggplot2)
library(Vizumap)
library(cowplot)

# plot the stations with Andalucia
ggplot() + geom_sf(data = andalucia$geometry) + geom_sf(data = meteo_ria$geometry) + ggtitle("Meteorological stations in Andalucia, Spain")

# plot single station for checking which region it is in
ggplot() + geom_sf(data = andalucia$geometry) + geom_sf(data = meteo_ria$geometry[1]) # so this point is definitely in cordoba region, need to find function that returns this

# use one of four pre-prepared colour palettes for the uncertainty plotting
cmBivPal <- build_palette(name = "GreenBlue")
view(cmBivPal)

# below, read.uv creates a df that is usable to build a map with the function after
temperature <- read.uv(data = temp_spain, estimate = "mean_temp", error = "temp_moe")

# convert geo info to sf object
geodata <- st_as_sf(get_stations_info_from('ria', ria_options))

# convert coordinates into decimals
for(i in 1:length(geodata$geometry)){
  geodata$geometry[i] <- DMS2decimal(geodata$geometry[i])
}

# create buffer around the points so they can be used as spatial polygons, and then convert into S4 object as required for plotting with vizumap
buffered_geodata <- st_buffer(geodata, dist = 3000) |> as_Spatial()

# create map using the biv color palette and the poverty df
tempBivMap <- build_bmap(data = temperature, geoData = buffered_geodata, id = "station_name", palette = "CyanMagenta")
map <- view(tempBivMap) + geom_sf(data = andalucia$geometry, color=alpha("black", 0.7), fill = NA) + ggtitle("Meteorological stations in Andalucia, Spain")

# create legend
tempBivKey <- build_bkey(data = temperature, palette = "CyanMagenta")
legend <- view(tempBivKey)

# plot map with legend, plot_grid is from cowplot library to plot things next to each other
plot_grid(map, legend, labels = NULL, scale = c(1, 0.5)) 

# intersecting the polygons with the stations
pol = andalucia$geometry 
pts = geodata$geometry
ints = st_intersects(pol, pts) # this returns a list with the geometries and then the number of the points that lie in this geometry

meanpolytemps <- c()
for (i in 1:length(ints)) {
  polytemp <- c()
  for (point in ints[[i]]) {
    temp <- temp_spain$mean_temp[point] 
    polytemp <- c(polytemp, temp) # these are the mean temperatures of the stations that are in polygon i
  }
  result <- mean(as.numeric(polytemp), na.rm = TRUE) # this is the value that we want to assign to the polygon in the end, the mean temperature of all the stations inside the polygon
  meanpolytemps <- c(meanpolytemps, result)
}

andalucia$mean_temperature <- meanpolytemps

p1 <- ggplot(andalucia) + geom_sf(aes(fill = mean_temperature)) + ggtitle("Mean temperatures in Andalucia, April 2022")

# now using yearly April means to assign values to polygons
diff <- setdiff(geodata$station_id, a$station_id) # check which station ids do not match. df 'a' contains available data per station id
geodata2 <- geodata[ ! geodata$station_id %in% diff, ] # delete rows with station ids that are not in the available meteo data
geodata2 <- left_join(geodata2, a) # this adds the mean temperature from all Aprils since 2001 from the 'a' (aggregated) df which was composed in lines above

# now we can create a new df with the mean temps per year for each polygon
yearmeanpolytemps <- c()
for (i in 1:length(ints)) {
  polytemp <- c()
  for (point in ints[[i]]) {
    temp <- geodata2$mean_temperature[point] 
    polytemp <- c(polytemp, temp) # these are the mean temperatures of the stations that are in polygon i
  }
  result <- mean(as.numeric(polytemp), na.rm = TRUE) # this is the value that we want to assign to the polygon in the end, the mean temperature of all the stations inside the polygon
  yearmeanpolytemps <- c(yearmeanpolytemps, result)
}

andalucia$yearly_mean_temp <- yearmeanpolytemps

p2 <- ggplot(andalucia) + geom_sf(aes(fill = yearly_mean_temp)) + ggtitle("Historic mean temperatures in Andalucia, April 2001-2022")

plot_grid(p1, p2)
