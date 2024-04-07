
### MAPPING STATES
# ---------------------------------------------------------------------------- #
library(ggplot2)
library(maps)
library(mapdata)
library(ggmap)

# Create a dataframe with sample latitudes and longitudes
locations <- data.frame(
  state = c("Nebraska", "Wyoming", "Colorado"),
  latitude = c(41.4925, 42.7550, 39.5501),
  longitude = c(-99.9018, -107.3025, -105.7821),
  kept = c(TRUE, TRUE, FALSE)
)

# make a color column
locations$color <- ifelse(locations$kept, "darkgreen", "red")

# Get map data for the specified states
states_map <- map_data("state")
subset_map <- subset(states_map, region %in% c("nebraska", "wyoming", "colorado"))

# Plot the map
ggplot() +
  geom_polygon(data = subset_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "white") +
  geom_point(data = locations, aes(x = longitude, y = latitude), color = locations$color, size = 3) +
  coord_fixed(1.3) +
  labs(title = "Map of Nebraska, Wyoming, and Colorado with Sample Points",
       x = "Longitude", y = "Latitude") +
  theme_minimal()
# ---------------------------------------------------------------------------- #


### MAPPING RASTERS
# ---------------------------------------------------------------------------- #
library(tidyverse)
library(readxl)
library(dplyr)
library(rayshader)
library(magick)
library(raster)
library(osmdata)
library(elevatr)
library(raster)
library(RColorBrewer)

# read in the samplesheet
sheet_path <- "/Users/nickminor/Documents/uwyo_experiments/0002/resources/samplesheet.xlsx"
sample_locs <- read_excel(sheet_path, trim_ws = TRUE) %>%
  select(`Sample ID`, Latitude, Longitude)

# make a bounding box for the samples
bbox <- make_bbox(lon = Longitude, lat = Latitude,
                  data = sample_locs, f = c(0.3, 0.3))
bbox_df <- data.frame(
  x = c(bbox["left"], bbox["right"]),
  y = c(bbox["top"], bbox["bottom"])
)
row.names(bbox_df) <- c("max", "min")
bbox_df["max", "y"] <- bbox_df["max", "y"] + 2

# pull a raster based on that bounding box dataframe
transect_raster <- bbox_df %>%
  get_elev_raster(locations =., prj = "+proj=longlat +datum=WGS84", z = 5)

# prepare a dataframe for ggplot
raster_df <- as.data.frame(transect_raster, xy = TRUE)
names(raster_df)[3] <- "Elevation"

# get state boundaries
states_map <- map_data("state")
subset_map <- subset(states_map, region %in% c("nebraska", "south dakota", "wyoming", "colorado"))

# generate the map
transect_map <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = Elevation)) +
  scale_fill_gradientn(colors = c("darkgreen", "green3", "khaki", "tan", "beige", "white")) +
  geom_polygon(data = subset_map,
               aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  coord_fixed() +
  geom_point(data = sample_locs,
             aes(x = Longitude, y = Latitude), color = "red", size = 2) +
  xlab("Latitute") +
  ylab("Longitude") +
  theme_minimal()
transect_map
# ---------------------------------------------------------------------------- #
