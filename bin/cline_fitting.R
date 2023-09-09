#!/usr/bin/env Rscript

# libraries
library(ggplot2)
library(sp)
library(geosphere)
library(hzar)
library(viridis)
library(dplyr)
library(readxl)
library(forcats) 

setwd("~/Documents/Dissertation Research/Passerina_hybrid_zone_changes/admixture_plotting")

###################### make hzar clines ######################


# load metadata
metadata <- read_excel("collected_buntings,all_nick.xlsx", trim_ws = TRUE)
sample_ids_paul<-read.csv("sample_ids.csv")
sample_coords_paul<-read.csv("collection_coords.csv")

# need to remove PJD from sample IDs in sample coords
sample_coords_paul$SampleID <- sub("PJD", "", sample_coords_paul$SampleID)

# need to extract first number from all sample ids
sample_ids_paul$SampleID <-regmatches(sample_ids_paul$SampleID, regexpr("[[:digit:]]*", sample_ids_paul$SampleID))

#removing individuals from sampleIDs that were removed frm vcf for missing too much data
sample_ids_paul <- sample_ids_paul %>% 
  filter(SampleID != 86 & SampleID != 142)

qrep1 <- h5read("passerina6_k2_0.hdf5","q") #my samples
qrep2 <- h5read("passerina6_k2_1.hdf5","q")
qrep3 <- h5read("passerina6_k2_2.hdf5","q")
q1_mat <- qrep1[,1,]
q2_mat <- qrep2[,1,]
q3_mat <- qrep3[,1,]
q_long <- rbind(q1_mat,q2_mat,q3_mat)
q_means <- apply(q_long,2,mean)

# merging sample id and gbs data
sample_ids_paul$q <- q_means

# merging dataframes by sample id
sample_ids_paul <- merge(sample_ids_paul,sample_coords_paul, by = "SampleID", sort=FALSE)

##could we also look at the frequency of interspecific ancestry in certain areas over time?
site_ancestry_paul <- aggregate(q ~ Longitude, data=sample_ids_paul, mean)
plot(q ~ Longitude, data=site_ancestry_paul)





## reading in distance dataset
distance_paul<-read.csv("collection_coords.csv")

# aggregating mean longitude by population
distances <- aggregate(Longitude ~ Population, data=distance_paul, mean)

plot(Population~Longitude, data = distances)

### converting longitudinal points to distance along transect
# first need to know mean latitude of sampling locations
mean_sample_lat <- mean(sample_ids_paul$Latitude)

# conversion for 1 deg longitude to km: 1 deg = 111.320*cos(latitude) km
long_degree_km = 111.320*cos(mean_sample_lat)
distances$distance <- long_degree_km*(abs((min(distances$Longitude)-distances$Longitude)))
distances <- distances %>% 
  select(Population, distance)

#merging genomic hybrid index for every sample with distances
sample_ids_paul_distance <- merge(sample_ids_paul,distances)
plot(q ~ distance, data=sample_ids_paul_distance)


## need to get a dataset of mean q and bird count at each site for both collection datasets
site_counts_paul <- aggregate(SampleID ~ distance, data=sample_ids_paul_distance, length)
site_ancestry_paul <- aggregate(q ~ distance, data=sample_ids_paul_distance, mean)
paul_cline <- merge(site_ancestry_paul,site_counts_paul)

# need to flip q values for one dataset so the cline will be in the same direction
plot(q ~ distance, data=paul_cline)

paul_cline$q <- 1-paul_cline$q
plot(q ~ distance, data=paul_cline)


# need to comment out all the models except the one chosen by AICc for easier plug n chug

### add an nSamples column with 1 for each because we are using samples not pops
man_admix$nSamples <- rep(1, nrow(man_admix))





### make data object of allele frequency data
paul_Q <-
  hzar.doMolecularData1DPops(paul_cline$distance,
                             paul_cline$q,
                             paul_cline$SampleID)

### Plot the associated observed frequency versus distance
hzar.plot.obsData(paul_Q)


### Construct a clineMetaModel object for use with hzar.first.fitRequest.old.ML
paul_Q_model_free_both <- hzar.makeCline1DFreq(paul_Q, scaling="free", tails="both")
paul_Q_model_free_none <- hzar.makeCline1DFreq(paul_Q, scaling="free", tails="none")
paul_Q_model_free_right <- hzar.makeCline1DFreq(paul_Q, scaling="free", tails="right")
paul_Q_model_free_left <- hzar.makeCline1DFreq(paul_Q, scaling="free", tails="left")
paul_Q_model_free_mirror <- hzar.makeCline1DFreq(paul_Q, scaling="free", tails="mirror")
# # 
paul_Q_model_fixed_both <- hzar.makeCline1DFreq(paul_Q, scaling="fixed", tails="both")
paul_Q_model_fixed_none <- hzar.makeCline1DFreq(paul_Q, scaling="fixed", tails="none")
paul_Q_model_fixed_right <- hzar.makeCline1DFreq(paul_Q, scaling="fixed", tails="right")
paul_Q_model_fixed_left <- hzar.makeCline1DFreq(paul_Q, scaling="fixed", tails="left")
paul_Q_model_fixed_mirror <- hzar.makeCline1DFreq(paul_Q, scaling="fixed", tails="mirror")
# # 
paul_Q_model_none_both <- hzar.makeCline1DFreq(paul_Q, scaling="none", tails="both")
paul_Q_model_none_none <- hzar.makeCline1DFreq(paul_Q, scaling="none", tails="none")
paul_Q_model_none_right <- hzar.makeCline1DFreq(paul_Q, scaling="none", tails="right")
paul_Q_model_none_left <- hzar.makeCline1DFreq(paul_Q, scaling="none", tails="left")
paul_Q_model_none_mirror <- hzar.makeCline1DFreq(paul_Q, scaling="none", tails="mirror")


### The intent of these methods is to assist the optimizer in exploring the model parameter space by instructing it to ignore models that are not interesting. For example, if all of the sampled localities are in a region 100km wide, then a cline width of 110km is probably not interesting. A cline width of 500km in that scenario would definitely not be interesting at all.
paul_Q_model_free_both <- hzar.model.addBoxReq(paul_Q_model_free_both,-100,900)
paul_Q_model_free_left <- hzar.model.addBoxReq(paul_Q_model_free_left,-100,900)
paul_S1_22a_model_free_right <- hzar.model.addBoxReq(paul_Q_model_free_right,-100,900)
paul_S1_22a_model_free_mirror <- hzar.model.addBoxReq(paul_Q_model_free_mirror,-100,900)
paul_S1_22a_model_free_none <- hzar.model.addBoxReq(paul_Q_model_free_none,-100,900)

paul_Q_model_fixed_both <- hzar.model.addBoxReq(paul_Q_model_fixed_both,-100,900)
paul_Q_model_fixed_left <- hzar.model.addBoxReq(paul_Q_model_fixed_left,-100,900)
paul_Q_model_fixed_right <- hzar.model.addBoxReq(paul_Q_model_fixed_right,-100,900)
paul_Q_model_fixed_mirror <- hzar.model.addBoxReq(paul_Q_model_fixed_mirror,-100,900)
paul_Q_model_fixed_none <- hzar.model.addBoxReq(paul_Q_model_fixed_none,-100,900)

paul_Q_model_none_both <- hzar.model.addBoxReq(paul_Q_model_none_both,-100,900)
paul_Q_model_none_left <- hzar.model.addBoxReq(paul_Q_model_none_left,-100,900)
paul_Q_model_none_right <- hzar.model.addBoxReq(paul_Q_model_none_right,-100,900)
paul_Q_model_none_mirror <- hzar.model.addBoxReq(paul_Q_model_none_mirror,-100,900)
paul_Q_model_none_none <- hzar.model.addBoxReq(paul_Q_model_none_none,-100,900)


###cline model fitting
##generate an hzar.fitRequest object suitable for hzar.doFit
paul_Q_model_free_bothFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_free_both, paul_Q, verbose=FALSE)
paul_Q_model_free_leftFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_free_left, paul_Q, verbose=FALSE)
paul_Q_model_free_rightFitR <-hzar.first.fitRequest.old.ML(model=paul_Q_model_free_right, paul_Q, verbose=FALSE)
paul_Q_model_free_mirrorFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_free_mirror, paul_Q, verbose=FALSE)
paul_Q_model_free_noneFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_free_none, paul_Q, verbose=FALSE)

paul_Q_model_fixed_bothFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_fixed_both, paul_Q, verbose=FALSE)
paul_Q_model_fixed_leftFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_fixed_left, paul_Q, verbose=FALSE)
paul_Q_model_fixed_rightFitR <-hzar.first.fitRequest.old.ML(model=paul_Q_model_fixed_right, paul_Q, verbose=FALSE)
paul_Q_model_fixed_mirrorFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_fixed_mirror, paul_Q, verbose=FALSE)
paul_Q_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_fixed_none, paul_Q, verbose=FALSE)

paul_Q_model_none_bothFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_none_both, paul_Q, verbose=FALSE)
paul_Q_model_none_leftFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_none_left, paul_Q, verbose=FALSE)
paul_Q_model_none_rightFitR <-hzar.first.fitRequest.old.ML(model=paul_Q_model_none_right, paul_Q, verbose=FALSE)
paul_Q_model_none_mirrorFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_none_mirror, paul_Q, verbose=FALSE)
paul_Q_model_none_noneFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_model_none_none, paul_Q, verbose=FALSE)

##set mcmc chain length and burn in
paul_Q_model_free_bothFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_free_bothFitR$mcmcParam$burnin <- 5e5
paul_Q_model_free_leftFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_free_leftFitR$mcmcParam$burnin <- 5e5
paul_Q_model_free_rightFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_free_rightFitR$mcmcParam$burnin <- 5e5
paul_Q_model_free_mirrorFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_free_mirrorFitR$mcmcParam$burnin <- 5e5
paul_Q_model_free_noneFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_free_noneFitR$mcmcParam$burnin <- 5e5

paul_Q_model_fixed_bothFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_fixed_bothFitR$mcmcParam$burnin <- 5e5
paul_Q_model_fixed_leftFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_fixed_leftFitR$mcmcParam$burnin <- 5e5
paul_Q_model_fixed_rightFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_fixed_rightFitR$mcmcParam$burnin <- 5e5
paul_Q_model_fixed_mirrorFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_fixed_mirrorFitR$mcmcParam$burnin <- 5e5
paul_Q_model_fixed_noneFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_fixed_noneFitR$mcmcParam$burnin <- 5e5

paul_Q_model_none_bothFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_none_bothFitR$mcmcParam$burnin <- 5e5
paul_Q_model_none_leftFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_none_leftFitR$mcmcParam$burnin <- 5e5
paul_Q_model_none_rightFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_none_rightFitR$mcmcParam$burnin <- 5e5
paul_Q_model_none_mirrorFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_none_mirrorFitR$mcmcParam$burnin <- 5e5
paul_Q_model_none_noneFitR$mcmcParam$chainLength <- 1e5
paul_Q_model_none_noneFitR$mcmcParam$burnin <- 5e5

##Run the optimizer using the parameters listed in the hzar.fitRequest given.
paul_Q_model_free_bothFit <- hzar.doFit(paul_Q_model_free_bothFitR)
paul_Q_model_free_leftFit <- hzar.doFit(paul_Q_model_free_leftFitR)
paul_Q_model_free_rightFit <- hzar.doFit(paul_Q_model_free_rightFitR)
paul_Q_model_free_mirrorFit <- hzar.doFit(paul_Q_model_free_mirrorFitR)
paul_Q_model_free_noneFit <- hzar.doFit(paul_Q_model_free_noneFitR)

paul_Q_model_fixed_bothFit <- hzar.doFit(paul_Q_model_fixed_bothFitR)
paul_Q_model_fixed_leftFit <- hzar.doFit(paul_Q_model_fixed_leftFitR)
paul_Q_model_fixed_rightFit <- hzar.doFit(paul_Q_model_fixed_rightFitR)
paul_Q_model_fixed_mirrorFit <- hzar.doFit(paul_Q_model_fixed_mirrorFitR)
paul_Q_model_fixed_noneFit <- hzar.doFit(paul_Q_model_fixed_noneFitR)

paul_Q_model_none_bothFit <- hzar.doFit(paul_Q_model_none_bothFitR)
paul_Q_model_none_leftFit <- hzar.doFit(paul_Q_model_none_leftFitR)
paul_Q_model_none_rightFit <- hzar.doFit(paul_Q_model_none_rightFitR)
paul_Q_model_none_mirrorFit <- hzar.doFit(paul_Q_model_none_mirrorFitR)
paul_Q_model_none_noneFit <- hzar.doFit(paul_Q_model_none_noneFitR)

# ### plot model to look for run stability and convergence and return the mcmc data with an added a log likelihood column.
par(mar=c(1,1,1,1))
plot(hzar.mcmc.bindLL(paul_Q_model_free_bothFit))
plot(hzar.mcmc.bindLL(paul_Q_model_free_leftFit))
plot(hzar.mcmc.bindLL(paul_Q_model_free_rightFit))
plot(hzar.mcmc.bindLL(paul_Q_model_free_mirrorFit))
plot(hzar.mcmc.bindLL(paul_Q_model_free_noneFit))

### group multiple fits of the same model and the same observation data into a single object
paul_Q_model_free_bothData <- hzar.dataGroup.add(paul_Q_model_free_bothFit)
paul_Q_model_free_leftData <- hzar.dataGroup.add(paul_Q_model_free_leftFit)
paul_Q_model_free_rightData <- hzar.dataGroup.add(paul_Q_model_free_rightFit)
paul_Q_model_free_mirrorData <- hzar.dataGroup.add(paul_Q_model_free_mirrorFit)
paul_Q_model_free_noneData <- hzar.dataGroup.add(paul_Q_model_free_noneFit)

paul_Q_model_fixed_bothData <- hzar.dataGroup.add(paul_Q_model_fixed_bothFit)
paul_Q_model_fixed_leftData <- hzar.dataGroup.add(paul_Q_model_fixed_leftFit)
paul_Q_model_fixed_rightData <- hzar.dataGroup.add(paul_Q_model_fixed_rightFit)
paul_Q_model_fixed_mirrorData <- hzar.dataGroup.add(paul_Q_model_fixed_mirrorFit)
paul_Q_model_fixed_noneData <- hzar.dataGroup.add(paul_Q_model_fixed_noneFit)

paul_Q_model_none_bothData <- hzar.dataGroup.add(paul_Q_model_none_bothFit)
paul_Q_model_none_leftData <- hzar.dataGroup.add(paul_Q_model_none_leftFit)
paul_Q_model_none_rightData <- hzar.dataGroup.add(paul_Q_model_none_rightFit)
paul_Q_model_none_mirrorData <- hzar.dataGroup.add(paul_Q_model_none_mirrorFit)
paul_Q_model_none_noneData <- hzar.dataGroup.add(paul_Q_model_none_noneFit)


### Generate a hzar.dataGroup object representing a fit of the null model to a hzar.obsData object
paul_Q_modelNull <- hzar.dataGroup.null(paul_Q)

##make list of cline models and null models
paul_Q_dGs <- list(cline_free_bothModel = paul_Q_model_free_bothData,
                  cline_free_leftModel = paul_Q_model_free_leftData,
                  cline_free_rightModel = paul_Q_model_free_rightData,
                  cline_free_mirrorModel = paul_Q_model_free_mirrorData,
                  cline_free_noneModel = paul_Q_model_free_noneData,
                  cline_fixed_bothModel = paul_Q_model_fixed_bothData,
                  cline_fixed_leftModel = paul_Q_model_fixed_leftData,
                  cline_fixed_rightModel = paul_Q_model_fixed_rightData,
                  cline_fixed_mirrorModel = paul_Q_model_fixed_mirrorData,
                  cline_fixed_noneModel = paul_Q_model_fixed_noneData,
                  cline_none_bothModel = paul_Q_model_none_bothData,
                  cline_none_leftModel = paul_Q_model_none_leftData,
                  cline_none_rightModel = paul_Q_model_none_rightData,
                  cline_none_mirrorModel = paul_Q_model_none_mirrorData,
                  cline_none_noneModel = paul_Q_model_none_noneData,
                           nullModel = paul_Q_modelNull)

### Collect optimizer output based on the same hzar.obsData object
paul_Q_oDG <- hzar.make.obsDataGroup(paul_Q_dGs)
 
### Set the names of the list of hzar.dataGroup objects contained in a hzar.obsDataGroup object using the names from either a named list of hzar.dataGroup objects or another hzar.obsDataGroup object.
paul_Q_oDG <- hzar.copyModelLabels(paul_Q_dGs,paul_Q_oDG)
 
### Plots a line representing the expected frequency versus distance for the given object. For hzar.dataGroup and hzar.obsDataGroup objects, plots the observed data backing the model. For hzar.obsDataGroup objects, plots the maximum likelihood cline for each model.
hzar.plot.cline(paul_Q_oDG)
 
### Calculate the AIC or corrected AIC score table for the given hzar.obsDataGroup object. There will be one score generated for each model associated with this object.
print(hzar.AICc.hzar.obsDataGroup(paul_Q_oDG))
#                             AICc
# cline_free_bothModel    38.74405
# cline_free_leftModel    32.61257
# cline_free_rightModel   32.76373
# cline_free_mirrorModel  32.74017
# cline_free_noneModel    27.38734
# cline_fixed_bothModel   32.19472
# cline_fixed_leftModel   26.99131
# cline_fixed_rightModel  26.47339
# cline_fixed_mirrorModel 26.98859
# cline_fixed_noneModel   22.23681
# cline_none_bothModel    32.27241
# cline_none_leftModel    26.98895
# cline_none_rightModel   26.88876
# cline_none_mirrorModel  26.99258
# cline_none_noneModel    22.23561
# nullModel               56.57411


#### plot clines
pdf(file="paul_q_cline.pdf")
hzar.plot.fzCline(paul_Q_model_fixed_noneData, fzCol = rgb(68/255, 1/255, 84/255, 0.3), pch = 20, col=rgb(68/255, 1/255, 84/255, 0.7), xlab = "Distance (km)", ylab = "Q value")
dev.off()
print(paul_Q_model_fixed_noneData$ML.cline$param.free$center)
#[1] 0.8665297
print(paul_Q_model_fixed_noneData$ML.cline$param.free$width)
#[1] 20.0455
print(hzar.getLLCutParam(paul_Q_model_fixed_noneData,c("center","width")))
# center2LLLow center2LLHigh width2LLLow width2LLHigh
#   -8.183703      5.566193    10.53658     47.68713


############## Plot all the clines together ##############
### remake all the models with the points subtracted around the centre
paul_centered <- paul_cline
paul_centered$distance <- paul_cline$dist-(paul_Q_model_fixed_noneData$ML.cline$param.free$center)

paul_Q_centered <-
  hzar.doMolecularData1DPops(paul_centered$distance,
                             paul_centered$q,
                             paul_centered$SampleID)


hzar.plot.obsData(paul_Q_centered)


### make ClineMetaModel object
paul_Q_centered_model_fixed_none<- hzar.makeCline1DFreq(paul_Q_centered, scaling="none", tails="none")

### set BoxReq
paul_Q_centered_model_fixed_none <- hzar.model.addBoxReq(paul_Q_centered_model_fixed_none,-100,60)

### generate hzar fitRequests
paul_Q_centered_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_centered_model_fixed_none, paul_Q_centered, verbose=FALSE)

### set mcmc chainLength and burnin
paul_Q_centered_model_fixed_noneFitR$mcmcParam$chainLength <- 1e5
paul_Q_centered_model_fixed_noneFitR$mcmcParam$burnin <- 5e5

### Run the optimizer using the parameters listed in the hzar.fitRequest given.
paul_Q_centered_model_fixed_noneFit <- hzar.doFit(paul_Q_centered_model_fixed_noneFitR)

### group multiple fits of the same model and the same observation data into a single object
paul_Q_centered_model_fixed_noneData <- hzar.dataGroup.add(paul_Q_centered_model_fixed_noneFit)

### get the range of parameter values that are within two log likelihood units of the maximum likelihood for center and width
print(hzar.getLLCutParam(paul_Q_centered_model_fixed_noneData,c("center","width")))


### MANNY BOEHM CODE
# hzar.overPlot.fzCline doesn't contain any code for creating plots :(
# so, let's extract info from 'hzar.dataGroup' objects in the same way
# that hzar.plot.fzCline does

# --------------------
# 95% CI

# hzar.plot.fzCline() fits a line (and the 95% CI) to the points by first 
# generating a series of 107 x coordinates bounded by the current plot window
# later, the y coordinates are estimated by accessing 
# the function describing the cline in the dataGroup object
## sometimes this dies and makes the polygons real small but if I rerun this code a couple times they go backto normal
# the numbers are from hzar, no relation to my observations being 108 in caor
xSeries <- seq(from = par("usr")[1], to = par("usr")[2], 
               length.out = 109)
if (par("xaxs") == "r") 
  xSeries <- xSeries[2:108]

# hzar.plot.fzCline() plots the 95% CI using graphics::polygon() 
# you can customize the colours here
paul_fzCline = hzar.getCredParamRed(paul_Q_centered_model_fixed_noneData)
paul_fzCoor <- paul_fzCline$fzCline(xSeries)

### plot the polygons and lines
# the stats::line() function uses the x and y coordinates from above to draw the line

pdf(file = "centered_clines.pdf")

# first create an empty plot bounded by the biggest dataset
hzar.plot.obsData(paul_Q_centered_model_fixed_noneData, col = "transparent")

#add polygons and lines
polygon(x = c(paul_fzCoor$x, rev(paul_fzCoor$x)), y = c(paul_fzCoor$yMin, rev(paul_fzCoor$yMax)), border = rgb(68/255, 1/255, 84/255, 0.75), col=rgb(68/255, 1/255, 84/255, 0.3))
lines(x = xSeries, y = paul_Q_centered_model_fixed_noneData$ML.cline$clineFunc(xSeries),  col=rgb(68/255, 1/255, 84/255), lwd=2)
dev.off()   

# do a zoomed in version on just the centers
pdf(file = "centered_clines_zoom.pdf")
