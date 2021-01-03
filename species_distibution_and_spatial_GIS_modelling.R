##### SPECIES DISTRIBUTION MODELLING ####

setwd("C:/TFE")
install.packages('raster') # Core raster GIS data package
install.packages('sf') # Core vector GIS data package
install.packages('sp') # Another core vector GIS package
install.packages('dismo') # Species Distribution Modelling
install.packages('rgdal') # Interface to the Geospatial Data Abstraction Library

library(rgdal)
library(raster)
library(sf)
library(sp)
library(dismo)

vignette('sdm')

tapir_IUCN <- st_read('data/iucn_mountain_tapir/data_0.shp')
print(tapir_IUCN)

# Load the data frame
tapir_GBIF <- read.delim('data/gbif_mountain_tapir.csv', 
                         stringsAsFactors=FALSE)
# Drop rows with missing coordinates
tapir_GBIF <- subset(tapir_GBIF, ! is.na(decimalLatitude) | ! is.na(decimalLongitude))
# Convert to an sf object 
tapir_GBIF <- st_as_sf(tapir_GBIF, coords=c('decimalLongitude', 'decimalLatitude'))
st_crs(tapir_GBIF) <- 4326
print(tapir_GBIF)


# Load some (coarse) country background data
sne110 <- st_read('data/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp')

# Create a modelling extent for plotting and cropping the global data.
model_extent <- extent(c(-85,-70,-5,12))

# Plot the species data over a basemap
plot(st_geometry(sne110), xlim=model_extent[1:2], ylim=model_extent[3:4], 
     bg='lightblue', col='ivory')
plot(st_geometry(tapir_IUCN), add=TRUE, col='grey', border=NA)
plot(st_geometry(tapir_GBIF), add=TRUE, col='red', pch=4, cex=0.6)
box()

# Load the data
bioclim_hist <- getData('worldclim', var='bio', res=10, path='data')
bioclim_2050 <- getData('CMIP5', var='bio', res=10, rcp=60, model='HD', year=50, path='data')

# Relabel the future variables to match the historical ones
names(bioclim_2050) <- names(bioclim_hist)

# Look at the data structure
print(bioclim_hist)


par(mfrow=c(3,1), mar=c(1,1,1,1))
# Create a shared colour scheme
breaks <- seq(-300, 350, by=20)
cols <- hcl.colors(length(breaks) - 1)
# Plot the historical and projected data
plot(bioclim_hist[[1]], breaks=breaks, col=cols)
plot(bioclim_2050[[1]], breaks=breaks, col=cols)
# Plot the temperature difference
plot(bioclim_2050[[1]] - bioclim_hist[[1]], col=hcl.colors(20, palette='Inferno'))

# Reduce the global maps down to the species' range
bioclim_hist_local <- crop(bioclim_hist, model_extent)
bioclim_2050_local <- crop(bioclim_2050, model_extent)

# Create a simple land mask
land <- bioclim_hist_local[[1]] >= 0
# How many points to create? We'll use the same as number of observations
n_pseudo <- nrow(tapir_GBIF)
# Sample the points
pseudo_dismo <- randomPoints(mask=land, n=n_pseudo, p=st_coordinates(tapir_GBIF))
# Convert this data into an sf object, for consistency with the
# next example.
pseudo_dismo <- st_as_sf(data.frame(pseudo_dismo), coords=c('x','y'), crs=4326)

# Create buffers around the observed points
nearby <- st_buffer(tapir_GBIF, dist=1)
too_close <- st_buffer(tapir_GBIF, dist=0.2)
# merge those buffers
nearby <- st_union(nearby)
too_close <- st_union(too_close)
# Find the area that is nearby but _not_ too close
nearby <- st_difference(nearby, too_close)
# Get some points within that feature in an sf dataframe
pseudo_nearby <- st_as_sf(st_sample(nearby, n_pseudo))

par(mfrow=c(1,2), mar=c(1,1,1,1))
# Random points on land
plot(land, col='grey', legend=FALSE)
plot(st_geometry(tapir_GBIF), add=TRUE, col='red')
plot(pseudo_dismo, add=TRUE)
# Random points within ~ 100 km but not closer than ~ 20 km
plot(land, col='grey', legend=FALSE)
plot(st_geometry(tapir_GBIF), add=TRUE, col='red')
plot(pseudo_nearby, add=TRUE)

# Use kfold to add labels to the data, splitting it into 5 parts
tapir_GBIF$kfold <- kfold(tapir_GBIF, k=5)
# Do the same for the pseudo-random points
pseudo_dismo$kfold <- kfold(pseudo_dismo, k=5)
pseudo_nearby$kfold <- kfold(pseudo_nearby, k=5)


# Get the coordinates of 80% of the data for training 
train_locs <- st_coordinates(subset(tapir_GBIF, kfold != 1))
# Fit the model
bioclim_model <- bioclim(bioclim_hist_local, train_locs)

par(mfrow=c(1,2))
plot(bioclim_model, a=1, b=2, p=0.9)
plot(bioclim_model, a=1, b=5, p=0.9)

bioclim_pred <- predict(bioclim_hist_local, bioclim_model)
# Create a copy removing zero scores to focus on within envelope locations
bioclim_non_zero <- bioclim_pred
bioclim_non_zero[bioclim_non_zero == 0] <- NA
plot(land, col='grey', legend=FALSE)
plot(bioclim_non_zero, col=hcl.colors(20, palette='Blue-Red'), add=TRUE)

test_locs <- st_coordinates(subset(tapir_GBIF, kfold == 1))
test_pseudo <- st_coordinates(subset(pseudo_nearby, kfold == 1))
bioclim_eval <- evaluate(p=test_locs, a=test_pseudo, model=bioclim_model, 
                         x=bioclim_hist_local)
print(bioclim_eval)


par(mfrow=c(1,2))
# Plot the ROC curve
plot(bioclim_eval, 'ROC', type='l')
# Find the maximum kappa and show how kappa changes with the model threshold
max_kappa <- threshold(bioclim_eval, stat='kappa')
plot(bioclim_eval, 'kappa', type='l')
abline(v=max_kappa, lty=2, col='blue')

# Apply the threshold to the model predictions
tapir_range <- bioclim_pred >= max_kappa
plot(tapir_range, legend=FALSE, col=c('grey','red'))
plot(st_geometry(tapir_GBIF), add=TRUE, pch=4, col='#00000022')

# Create a single sf object holding presence and pseudo-absence data.
# - reduce the GBIF data and pseudo-absence to just kfold and a presence-absence value
present <- subset(tapir_GBIF, select='kfold')
present$pa <- 1
absent <- pseudo_nearby
absent$pa <- 0
# - rename the geometry column of absent to match so we can stack them together.
names(absent) <- c('geometry','kfold','pa')
st_geometry(absent) <- 'geometry'
# - stack the two dataframes
pa_data <- rbind(present, absent)
head(pa_data)
tail(pa_data)

envt_data <- extract(bioclim_hist_local, pa_data)
pa_data <- cbind(pa_data, envt_data)
head(pa_data)

glm_model <- glm(pa ~ bio2 + bio4 + bio3 + bio1 + bio12, data=pa_data, 
                 family=binomial(link = "logit"),
                 subset=kfold != 1)
# Look at the variable significances - which are important
summary(glm_model)

# Response plots
response(glm_model, fun=function(x, y, ...) predict(x, y, type='response', ...))


## evaluate model
glm_pred <- predict(bioclim_hist_local, glm_model, type='response')

# Extract the test presence and absence
test_present <- st_coordinates(subset(pa_data, pa == 1 & kfold == 1))
test_absent <- st_coordinates(subset(pa_data, pa == 0 & kfold == 1))
glm_eval <- evaluate(p=test_present, a=test_absent, model=glm_model, 
                     x=bioclim_hist_local)
print(glm_eval)


max_kappa <- plogis(threshold(glm_eval, stat='kappa'))
print(max_kappa)

## model perfomance plots
par(mfrow=c(1,2))
# ROC curve and kappa by model threshold
plot(glm_eval, 'ROC', type='l')
plot(glm_eval, 'kappa', type='l')
abline(v=max_kappa, lty=2, col='blue')


par(mfrow=c(2,2))
# Modelled probability
plot(glm_pred, col=hcl.colors(20, 'Blue-Red'))
# Threshold map
glm_map <- glm_pred >= max_kappa
plot(glm_map, legend=FALSE, col=c('grey','red'))
# Future predictions
glm_pred_2050 <- predict(bioclim_2050_local, glm_model, type='response')
plot(glm_pred_2050, col=hcl.colors(20, 'Blue-Red'))
# Threshold map
glm_map_2050 <- glm_pred_2050 >= max_kappa
plot(glm_map_2050, legend=FALSE, col=c('grey','red'))

table(values(glm_map), values(glm_map_2050), dnn=c('hist', '2050'))



#### my own go at modelling species ####

# Check the species (Pale throated Sloth) without downloading - this shows the number of records
gbif('Bradypus', 'tridactylus', download=FALSE)

# Download the data
locs <- gbif('Bradypus', 'tridactylus')
locs <- subset(locs, ! is.na(lat) | ! is.na(lon))
# Convert to an sf object 
locs <- st_as_sf(locs, coords=c('lon', 'lat'))
st_crs(locs) <- 4326

### check basis of records for locatioNs, remove preserved and fossil specimens ###
View(locs)
sloth_obvs <- subset(locs, locs$basisOfRecord=="HUMAN_OBSERVATION") 

#### see where sloths are in the world ###
plot(st_geometry(sne110))
plot(st_geometry(sloth_obvs), add=TRUE, col='red', pch=4, cex=0.6)

## change model extent for map ###
par(mfrow=c(1,1))
model_extent <- extent(c(-90,-30,-5,10))
# Plot the species data over a basemap
plot(st_geometry(sne110), xlim=model_extent[1:2], ylim=model_extent[3:4], 
     bg='lightblue', col='ivory')
plot(st_geometry(sloth_obvs), add=TRUE, col='purple', pch=3, cex=0.6)
box()


#### looking at climate data ####

# Load the data
bioclim_hist <- getData('worldclim', var='bio', res=10, path='data')
bioclim_2050 <- getData('CMIP5', var='bio', res=10, rcp=60, model='HD', year=50, path='data')

# Relabel the future variables to match the historical ones
names(bioclim_2050) <- names(bioclim_hist)

# Look at the data structure
print(bioclim_hist)

par(mfrow=c(3,1), mar=c(1,1,1,1))
# Create a shared colour scheme
breaks <- seq(-300, 350, by=20)
cols <- hcl.colors(length(breaks) - 1)
# Plot the historical and projected data
plot(bioclim_hist[[1]], breaks=breaks, col=cols)
plot(bioclim_2050[[1]], breaks=breaks, col=cols)
# Plot the temperature difference
plot(bioclim_2050[[1]] - bioclim_hist[[1]], col=hcl.colors(20, palette='Inferno'))

# Reduce the global maps down to the sloth range
bioclim_hist_local <- crop(bioclim_hist, model_extent)
bioclim_2050_local <- crop(bioclim_2050, model_extent)

# Create a simple land mask
land <- bioclim_hist_local[[1]] >= 0
# How many points to create? We'll use the same as number of observations
n_pseudo <- nrow(sloth_obvs)
# Sample the points
pseudo_dismo <- randomPoints(mask=land, n=n_pseudo, p=st_coordinates(sloth_obvs))
# Convert this data into an sf object, for consistency with the
# next example.
pseudo_dismo <- st_as_sf(data.frame(pseudo_dismo), coords=c('x','y'), crs=4326)


# Create buffers around the observed points
nearby <- st_buffer(sloth_obvs, dist=1)
too_close <- st_buffer(sloth_obvs, dist=0.2)
# merge those buffers
nearby <- st_union(nearby)
too_close <- st_union(too_close)
# Find the area that is nearby but _not_ too close
nearby <- st_difference(nearby, too_close)
# Get some points within that feature in an sf dataframe
pseudo_nearby <- st_as_sf(st_sample(nearby, n_pseudo))

par(mfrow=c(1,1), mar=c(1,1,1,1))
# Random points on land
plot(land, col='grey', legend=FALSE)
plot(st_geometry(sloth_obvs), add=TRUE, col='red')
plot(pseudo_dismo, add=TRUE)
# Random points within ~ 100 km but not closer than ~ 20 km
plot(land, col='grey', legend=FALSE)
plot(st_geometry(sloth_obvs), add=TRUE, col='red')
plot(pseudo_nearby, add=TRUE)


### testing and training dataset ####

# Use kfold to add labels to the data, splitting it into 5 parts
sloth_obvs$kfold <- kfold(sloth_obvs, k=5)
# Do the same for the pseudo-random points
pseudo_dismo$kfold <- kfold(pseudo_dismo, k=5)
pseudo_nearby$kfold <- kfold(pseudo_nearby, k=5)

#### SPATIAL MODELLING IN R ####

install.packages('ncf')
install.packages('SpatialPack') # For clifford test
install.packages('spdep') # Spatial dependence
install.packages('nlme') # GLS
install.packages('spgwr')

library(ncf)
library(spdep)
library(SpatialPack)
library(nlme)
library(spgwr)



#get TIFs
rich <- raster('data/avian_richness.tif')
aet <- raster('data/mean_aet.tif')
temp <- raster('data/mean_temp.tif')
elev <- raster('data/elev.tif')


#plot tifs
par(mfrow=c(2,2), mar=c(2,2,2,2))
plot(rich, main='Avian species richness')
plot(aet, main='Mean AET')
plot(temp, main='Mean annual temperature')
plot(elev, main='Elevation')

# split the figure area into a two by two layout
par(mfrow=c(2,2))
# plot a histogram of the values in each raster, setting nice 'main' titles
hist(rich, main='Avian species richness')
hist(aet, main='Mean AET')
hist(temp, main='Mean annual temperature')
hist(elev, main='Elevation')


# Stack the data
data_stack <- stack(rich, aet, elev, temp)
print(data_stack)

data_spdf <- as(data_stack, 'SpatialPixelsDataFrame')
summary(data_spdf)

data_sf <- st_as_sf(data_spdf)
print(data_sf)

cellsize <- res(rich)[[1]]
template <- st_polygon(list(matrix(c(-1,-1,1,1,-1,-1,1,1,-1,-1), ncol=2) * cellsize / 2))

polygon_data <- lapply(data_sf$geometry, function(pt) template + pt)
data_poly <- st_sf(avian_richness = data_sf$avian_richness, 
                   geometry=polygon_data)
plot(data_poly['avian_richness'])


layout(matrix(1:3, ncol=3), widths=c(5,5,1))
plot(data_spdf['avian_richness'], col=hcl.colors(20), what='image')
plot(data_sf['avian_richness'], key.pos=NULL, reset=FALSE, main='', 
     pal=hcl.colors, cex=0.7, pch=20)
plot(data_spdf['avian_richness'], col=hcl.colors(20), what='scale')


# Create three figures in a single panel
par(mfrow=c(1,3))
# Now plot richness as a function of each environmental variable
plot(avian_richness ~ mean_aet, data=data_sf)
plot(avian_richness ~ mean_temp, data=data_sf)
plot(avian_richness ~ elev, data=data_sf)


## corelations ##
temp_corr <- modified.ttest(x=data_sf$avian_richness, y=data_sf$mean_temp, 
                            coords=st_coordinates(data_sf))

print(temp_corr)

elev_corr <- modified.ttest(x=data_sf$avian_richness, y=data_sf$elev, 
                            coords=st_coordinates(data_sf))
print(elev_corr)

et_corr <- modified.ttest(x=data_sf$avian_richness, y=data_sf$mean_aet, 
                          coords=st_coordinates(data_sf))
print(et_corr)

tempvselev_corr <- modified.ttest(x=data_sf$mean_temp, y=data_sf$elev, 
                                  coords=st_coordinates(data_sf))
print(tempvselev_corr)

tempvsaet_corr <- modified.ttest(x=data_sf$mean_temp, y=data_sf$mean_aet, 
                                 coords=st_coordinates(data_sf))
print(tempvsaet_corr)

elevvsaet_corr <- modified.ttest(x=data_sf$elev, y=data_sf$mean_aet, 
                                 coords=st_coordinates(data_sf))
print(elevvsaet_corr)


# Give dnearneigh the coordinates of the points and the distances to use
rook <- dnearneigh(data_sf, d1=0, d2=cellsize)
queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize)

print(rook)
head(rook, n=3)
print(queen)
head(queen, n=3)

# Store the neighbourhood cardinalities in data_sf
data_sf$card_rook <- card(rook)
data_sf$card_queen <- card(queen)
# Look at the count of rook and queen neighbours for each point
plot(data_sf[c('card_rook', 'card_queen')], key.pos=4)



# Recreate the neighbour adding 1km to the distance
rook <- dnearneigh(data_sf, d1=0, d2=cellsize + 1)
queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize + 1)
data_sf$card_rook <- card(rook)
data_sf$card_queen <- card(queen)
plot(data_sf[c('card_rook', 'card_queen')], key.pos=4)


### closest locations, to account for islands with no neighbours when dsitance is fixed e.g. Principe, Comoros ###
knn <- knearneigh(data_sf, k=8)
# We have to look at the `nn` values in `knn` to see the matrix of neighbours
head(knn$nn, n=3)

## spatial weights ###
queen <- nb2listw(queen)
## didn't work!! cannot create weighted lists for locations with no neighbours
### remove locations with no neighbours
# Polygon covering Mauritius
mauritius <- st_polygon(list(matrix(c(5000, 6000, 6000, 5000, 5000,
                                      0, 0, -4000, -4000, 0), 
                                    ncol=2)))
mauritius <- st_sfc(mauritius, crs=crs(data_sf, asText=TRUE))
# Remove the island cells with zero neighbours
data_sf <- subset(data_sf, card_rook > 0)
# Remove Mauritius
data_sf <- st_difference(data_sf, mauritius)

rook <- dnearneigh(data_sf, d1=0, d2=cellsize + 1)
queen <- dnearneigh(data_sf, d1=0, d2=sqrt(2) * cellsize + 1)
data_sf$card_rook <- card(rook)
data_sf$card_queen <- card(queen)
knn <- knearneigh(data_sf, k=8)

## calculating weights
rook <- nb2listw(rook, style='W')
queen <- nb2listw(queen, style='W')
knn <- nb2listw(knn2nb(knn), style='W')


## global spatial autocorrelation ###
moran_avian_richness <- moran.test(data_sf$avian_richness, rook)
print(moran_avian_richness)
geary_avian_richness <- geary.test(data_sf$avian_richness, rook)
print(geary_avian_richness)

## try for different neighbourhood methods ###
qmoran_avian_richness <- moran.test(data_sf$avian_richness, queen)
print(qmoran_avian_richness)
qgeary_avian_richness <- geary.test(data_sf$avian_richness, queen)
print(qgeary_avian_richness)

knnmoran_avian_richness <- moran.test(data_sf$avian_richness, knn)
print(knnmoran_avian_richness)
knngeary_avian_richness <- geary.test(data_sf$avian_richness, knn)
print(knngeary_avian_richness)



## local spatial autocorrelation ###

local_moran_avian_richness <- localmoran(data_sf$avian_richness, rook)
data_sf$local_moran_avian_richness <- local_moran_avian_richness[, 'Ii'] # Observed Moran's I
plot(data_sf['local_moran_avian_richness'], cex=0.6, pch=20)

data_sf$local_g_avian_richness <- localG(data_sf$avian_richness, rook)
plot(data_sf['local_g_avian_richness'], cex=0.6, pch=20)




### autoregressive models ###

# fit a simple linear model
simple_model <- lm(avian_richness ~ mean_aet + elev + mean_temp, data=data_sf)
summary(simple_model)

# fit a spatial autoregressive model : much slower, can take minutes 
sar_model <- errorsarlm(avian_richness ~ mean_aet + elev + mean_temp, 
                        data=data_sf, listw=queen)
summary(sar_model)

# extract the predictions from the model into the spatial data frame
data_sf$simple_fit <- predict(simple_model)
data_sf$sar_fit <- predict(sar_model)
# Compare those two predictions with the data
plot(data_sf[c('avian_richness','simple_fit','sar_fit')], 
     pal=function(n) hcl.colors(n, pal='Blue-Red'))


# extract the residuals from the model into the spatial data frame
data_sf$simple_resid <- residuals(simple_model)
data_sf$sar_resid <- residuals(sar_model)
plot(data_sf[c('avian_richness','simple_resid', 'sar_resid')], 
     pal=function(n) hcl.colors(n, pal='Blue-Red'), key.pos=4)


### correlograms ####

# add the X and Y coordinates to the data frame
data_xy <- data.frame(st_coordinates(data_sf))
data_sf$x <- data_xy$X
data_sf$y <- data_xy$Y


# calculate a correlogram for avian richness: a slow process!
rich.correlog <- with(data_sf, correlog(x, y, avian_richness, increment=cellsize, resamp=0))
par(mfrow=c(1,1))
plot(rich.correlog)

par(mfrow=c(1,2))
# convert three key variables into a data frame
rich.correlog <- data.frame(rich.correlog[1:3])
# plot the size of the distance bins
plot(n ~ mean.of.class, data=rich.correlog, type='o')
# plot a correlogram for shorter distances
plot(correlation ~ mean.of.class, data=rich.correlog, type='o', subset=mean.of.class < 5000)
# add a horizontal  zero correlation line
abline(h=0)


# Calculate correlograms for the residuals in the two models, assessment of model's control of autocorrelation 
simple.correlog <- with(data_sf, correlog(x, y, simple_resid, increment=cellsize, resamp=0))
sar.correlog <- with(data_sf, correlog(x, y, sar_resid, increment=cellsize, resamp=0))

# Convert those to make them easier to plot
simple.correlog <- data.frame(simple.correlog[1:3])
sar.correlog <- data.frame(sar.correlog[1:3])

# plot a correlogram for shorter distances
par(mfrow=c(1,1))
plot(correlation ~ mean.of.class, data=simple.correlog, type='o', subset=mean.of.class < 5000)
# add the data for the SAR model to compare them
lines(correlation ~ mean.of.class, data=sar.correlog, type='o', subset=mean.of.class < 5000, col='red')
# add a horizontal  zero correlation line
abline(h=0)



### generalised least squares ####

# fit the simple model
gls_simple <- gls(avian_richness ~ mean_aet + mean_temp + elev, data=data_sf)
summary(gls_simple)

# Define the correlation structure
cor_struct_gauss <- corGaus(value=c(range=650, nugget=0.1), form=~ x + y, fixed=TRUE, nugget=TRUE)
# Add that to the simple model - this might take a few minutes to run!
gls_gauss <- update(gls_simple,  corr=cor_struct_gauss)
summary(gls_gauss)

# map model predictions
data_sf$gls_simple_pred <- predict(gls_simple)
data_sf$gls_gauss_pred <- predict(gls_gauss)
plot(data_sf[c('gls_simple_pred','gls_gauss_pred')], key.pos=4, cex=0.6, pch=20)

#Extract the residuals
data_sf$gls_simple_resid <- resid(gls_simple)
data_sf$gls_gauss_resid <- resid(gls_gauss)

# Calculate correlograms for the residuals in the two models
simple.correlog <- with(data_sf, correlog(x, y, gls_simple_resid, increment=cellsize, resamp=0))
gauss.correlog <- with(data_sf, correlog(x, y, gls_gauss_resid, increment=cellsize, resamp=0))

# Convert those to make them easier to plot
simple.correlog <- data.frame(simple.correlog[1:3])
gauss.correlog <- data.frame(gauss.correlog[1:3])

# plot a correlogram for shorter distances
plot(correlation ~ mean.of.class, data=simple.correlog, type='o', subset=mean.of.class < 5000)
# add the data for the SAR model to compare them
lines(correlation ~ mean.of.class, data=gauss.correlog, type='o', subset=mean.of.class < 5000, col='red')
# add a horizontal  zero correlation line
abline(h=0)

### is worse than sar model! ##


### relationship between observed and predicted richness ####

par(mfrow=c(1,3), mar=c(3,3,1,1), mgp=c(2,1,0))
# set the plot limits to show the scaling clearly
lims <- c(0, max(data_sf$avian_richness))
plot(avian_richness ~ gls_simple_pred, data=data_sf, xlim=lims, ylim=lims)
abline(a=0, b=1, col='red', lwd=2)
plot(avian_richness ~ gls_gauss_pred, data=data_sf, xlim=lims, ylim=lims)
abline(a=0, b=1, col='red', lwd=2)
plot(avian_richness ~ sar_fit, data=data_sf, xlim=lims, ylim=lims)
abline(a=0, b=1, col='red', lwd=2)


