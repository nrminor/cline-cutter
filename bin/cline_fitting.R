#!/usr/bin/env Rscript

###* In this script I model genomic clines for Passerina buntings collected by me (Paul Dougherty)
###* from 2018 to 2021 and collected by Matt Carling to visualize how the hybrid zone has changed
###* in the time between our sampling efforts. This script was heavily based on a script by Libby
###* Natola for Natola et al. 2023

# set working directory
setwd("~/Documents/UBC/Bioinformatics/rnrb/clines/")

# libraries
library(ggplot2)
library(geosphere)
library(hzar)
library(viridis)
library(dplyr)
library(forcats) 
library(sp)

setwd("~/Documents/Dissertation Research/Passerina_hybrid_zone_changes/admixture_plotting")

###################### make hzar clines ######################

#### reading in sampling IDs
## Matt's data
sample_ids_matt<-read.csv("historic_sampleIDs.csv")
sample_coords_matt<-read.csv("historic_passerina_collection_coords.csv")
sample_coords_matt<-read.csv("historic_passerina_collection_coords,matts_birds.csv")
#* need to remove samples that are in sampleIDs but not coord list and samples 
#* in coord list that aren't in sample IDs, I think the merge function does this automatically

# loading and adding Entropy estimates of hybrid index to data frame
qrep1 <- h5read("passerina_historic_k2_0.hdf5","q") #Matt's samples
qrep2 <- h5read("passerina_historic_k2_1.hdf5","q")
qrep3 <- h5read("passerina_historic_k2_2.hdf5","q")
q1_mat <- qrep1[,1,]
q2_mat <- qrep2[,1,]
q3_mat <- qrep3[,1,]
q_long <- rbind(q1_mat,q2_mat,q3_mat)
q_means <- apply(q_long,2,mean)

# merging sample id and gbs data
sample_ids_matt$q <- q_means

# merging dataframes by sample id
sample_ids_matt <- merge(sample_ids_matt,sample_coords_matt, by = "SampleID", sort=FALSE)


## my data
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

# plotting mean hybrid index score at each collection site
site_ancestry_matt <- aggregate(q ~ Longitude, data=sample_ids_matt, mean)
plot(q ~ Longitude, data=site_ancestry_matt)

site_ancestry_paul <- aggregate(q ~ Longitude, data=sample_ids_paul, mean)
points(q ~ Longitude, data=site_ancestry_paul, col="red")



## reading in distance datasets

distance_matt<-read.csv("historic_passerina_collection_coords,matts_birds.csv")
distance_paul<-read.csv("collection_coords.csv")

# for both datasets, aggregating mean longitude by population
distance_matt <- aggregate(Longitude ~ Population, data=distance_matt, mean)
distance_paul <- aggregate(Longitude ~ Population, data=distance_paul, mean)

distances <- rbind(distance_paul, distance_matt)

plot(Population~Longitude, data = distances)

### converting longitudinal points to distance along transect
# first need to know mean latitude of sampling locations
mean_sample_lat <- mean(c(sample_ids_paul$Latitude, sample_ids_matt$Latitude))

# conversion for 1 deg longitude to km: 1 deg = 111.320*cos(latitude) km
long_degree_km = 111.320*cos(mean_sample_lat)
distances$distance <- long_degree_km*(abs((min(distances$Longitude)-distances$Longitude)))
distances <- distances %>% 
  select(Population, distance)


#merging genomic hybrid index for every sample with distances
sample_ids_paul_distance <- merge(sample_ids_paul,distances)

matt_distances_samples <- merge(distance_matt,distances)
matt_distances_samples <- matt_distances_samples %>% 
  select(Longitude,distance)
matt_distances_samples <- aggregate(Longitude ~ distance, data=matt_distances_samples, mean)

sample_ids_matt_distance <- merge(sample_ids_matt, matt_distances_samples)

plot(q ~ distance, data=sample_ids_matt_distance)
points(q ~ distance, data=sample_ids_paul_distance, col="red")


## need to get a dataset of mean q and bird count at each site for both collection datasets
site_counts_matt <- aggregate(SampleID ~ distance, data=sample_ids_matt_distance, length)
site_ancestry_matt <- aggregate(q ~ distance, data=sample_ids_matt_distance, mean)
matt_cline <- merge(site_ancestry_matt,site_counts_matt)

site_counts_paul <- aggregate(SampleID ~ distance, data=sample_ids_paul_distance, length)
site_ancestry_paul <- aggregate(q ~ distance, data=sample_ids_paul_distance, mean)
paul_cline <- merge(site_ancestry_paul,site_counts_paul)

# need to flip q values for one dataset so the cline will be in the same direction
plot(q ~ distance, data=paul_cline)

paul_cline$q <- 1-paul_cline$q
plot(q ~ distance, data=paul_cline)

plot(q ~ distance, data=matt_cline)
points(q ~ distance, data=paul_cline, col="red")







# need to comment out all the models except the one chosen by AICc for easier plug n chug

### add an nSamples column with 1 for each because we are using samples not pops
man_admix$nSamples <- rep(1, nrow(man_admix))

### make data object of allele frequency data
paul_Q <-
  hzar.doMolecularData1DPops(paul_cline$distance,
                             paul_cline$q,
                             paul_cline$SampleID)
matt_Q <-
  hzar.doMolecularData1DPops(matt_cline$distance,
                             matt_cline$q,
                             matt_cline$SampleID)



### Plot the associated observed frequency versus distance
hzar.plot.obsData(paul_Q)
hzar.plot.obsData(matt_Q)


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



### add an nSamples column with 1 for each because we are using samples not pops
matt_admix$nSamples <- rep(1, nrow(matt_admix))

### make data object of allele frequency data
matt_Q <-
  hzar.doMolecularData1DPops(matt_cline$dist,
                             matt_cline$q,
                             matt_cline$SampleID)

### Plot the associated observed frequency versus distance
hzar.plot.obsData(matt_Q)


### Construct a clineMetaModel object for use with hzar.first.fitRequest.old.ML
matt_Q_model_free_both <- hzar.makeCline1DFreq(matt_Q, scaling="free", tails="both")
matt_Q_model_free_none <- hzar.makeCline1DFreq(matt_Q, scaling="free", tails="none")
matt_Q_model_free_right <- hzar.makeCline1DFreq(matt_Q, scaling="free", tails="right")
matt_Q_model_free_left <- hzar.makeCline1DFreq(matt_Q, scaling="free", tails="left")
matt_Q_model_free_mirror <- hzar.makeCline1DFreq(matt_Q, scaling="free", tails="mirror")

matt_Q_model_fixed_both <- hzar.makeCline1DFreq(matt_Q, scaling="fixed", tails="both")
matt_Q_model_fixed_none <- hzar.makeCline1DFreq(matt_Q, scaling="fixed", tails="none")
matt_Q_model_fixed_right <- hzar.makeCline1DFreq(matt_Q, scaling="fixed", tails="right")
matt_Q_model_fixed_left <- hzar.makeCline1DFreq(matt_Q, scaling="fixed", tails="left")
matt_Q_model_fixed_mirror <- hzar.makeCline1DFreq(matt_Q, scaling="fixed", tails="mirror")

matt_Q_model_none_both <- hzar.makeCline1DFreq(matt_Q, scaling="none", tails="both")
matt_Q_model_none_none <- hzar.makeCline1DFreq(matt_Q, scaling="none", tails="none")
matt_Q_model_none_right <- hzar.makeCline1DFreq(matt_Q, scaling="none", tails="right")
matt_Q_model_none_left <- hzar.makeCline1DFreq(matt_Q, scaling="none", tails="left")
matt_Q_model_none_mirror <- hzar.makeCline1DFreq(matt_Q, scaling="none", tails="mirror")


##The intent of these methods is to assist the optimizer in exploring the model parameter space by instructing it to ignore models that are not interesting. For example, if all of the sampled localities are in a region 100km wide, then a cline width of 110km is probably not interesting. A cline width of 500km in that scenario would definitely not be interesting at all.
matt_Q_model_free_both <- hzar.model.addBoxReq(matt_Q_model_free_both,-100,900)
matt_Q_model_free_left <- hzar.model.addBoxReq(matt_Q_model_free_left,-100,900)
matt_Q_model_free_right <- hzar.model.addBoxReq(matt_Q_model_free_right,-100,900)
matt_Q_model_free_mirror <- hzar.model.addBoxReq(matt_Q_model_free_mirror,-100,900)
matt_Q_model_free_none <- hzar.model.addBoxReq(matt_Q_model_free_none,-100,900)

matt_Q_model_fixed_both <- hzar.model.addBoxReq(matt_Q_model_fixed_both,-100,900)
matt_Q_model_fixed_left <- hzar.model.addBoxReq(matt_Q_model_fixed_left,-100,900)
matt_Q_model_fixed_right <- hzar.model.addBoxReq(matt_Q_model_fixed_right,-100,900)
matt_Q_model_fixed_mirror <- hzar.model.addBoxReq(matt_Q_model_fixed_mirror,-100,900)
matt_Q_model_fixed_none <- hzar.model.addBoxReq(matt_Q_model_fixed_none,-100,900)

matt_Q_model_none_both <- hzar.model.addBoxReq(matt_Q_model_none_both,-100,900)
matt_Q_model_none_left <- hzar.model.addBoxReq(matt_Q_model_none_left,-100,900)
matt_Q_model_none_right <- hzar.model.addBoxReq(matt_Q_model_none_right,-100,900)
matt_Q_model_none_mirror <- hzar.model.addBoxReq(matt_Q_model_none_mirror,-100,900)
matt_Q_model_none_none <- hzar.model.addBoxReq(matt_Q_model_none_none,-100,900)


#### cline model fitting
### generate an hzar.fitRequest object suitable for hzar.doFit
matt_Q_model_free_bothFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_free_both, matt_Q, verbose=FALSE)
matt_Q_model_free_leftFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_free_left, matt_Q, verbose=FALSE)
matt_Q_model_free_rightFitR <-hzar.first.fitRequest.old.ML(model=matt_Q_model_free_right, matt_Q, verbose=FALSE)
matt_Q_model_free_mirrorFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_free_mirror, matt_Q, verbose=FALSE)
matt_Q_model_free_noneFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_free_none, matt_Q, verbose=FALSE)

matt_Q_model_fixed_bothFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_fixed_both, matt_Q, verbose=FALSE)
matt_Q_model_fixed_leftFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_fixed_left, matt_Q, verbose=FALSE)
matt_Q_model_fixed_rightFitR <-hzar.first.fitRequest.old.ML(model=matt_Q_model_fixed_right, matt_Q, verbose=FALSE)
matt_Q_model_fixed_mirrorFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_fixed_mirror, matt_Q, verbose=FALSE)
matt_Q_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_fixed_none, matt_Q, verbose=FALSE)

matt_Q_model_none_bothFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_none_both, matt_Q, verbose=FALSE)
matt_Q_model_none_leftFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_none_left, matt_Q, verbose=FALSE)
matt_Q_model_none_rightFitR <-hzar.first.fitRequest.old.ML(model=matt_Q_model_none_right, matt_Q, verbose=FALSE)
matt_Q_model_none_mirrorFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_none_mirror, matt_Q, verbose=FALSE)
matt_Q_model_none_noneFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_model_none_none, matt_Q, verbose=FALSE)

##set mcmc chain length and burn in
matt_Q_model_free_bothFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_free_bothFitR$mcmcParam$burnin <- 5e5
matt_Q_model_free_leftFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_free_leftFitR$mcmcParam$burnin <- 5e5
matt_Q_model_free_rightFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_free_rightFitR$mcmcParam$burnin <- 5e5
matt_Q_model_free_mirrorFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_free_mirrorFitR$mcmcParam$burnin <- 5e5
matt_Q_model_free_noneFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_free_noneFitR$mcmcParam$burnin <- 5e5

matt_Q_model_fixed_bothFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_fixed_bothFitR$mcmcParam$burnin <- 5e5
matt_Q_model_fixed_leftFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_fixed_leftFitR$mcmcParam$burnin <- 5e5
matt_Q_model_fixed_rightFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_fixed_rightFitR$mcmcParam$burnin <- 5e5
matt_Q_model_fixed_mirrorFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_fixed_mirrorFitR$mcmcParam$burnin <- 5e5
matt_Q_model_fixed_noneFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_fixed_noneFitR$mcmcParam$burnin <- 5e5

matt_Q_model_none_bothFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_none_bothFitR$mcmcParam$burnin <- 5e5
matt_Q_model_none_leftFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_none_leftFitR$mcmcParam$burnin <- 5e5
matt_Q_model_none_rightFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_none_rightFitR$mcmcParam$burnin <- 5e5
matt_Q_model_none_mirrorFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_none_mirrorFitR$mcmcParam$burnin <- 5e5
matt_Q_model_none_noneFitR$mcmcParam$chainLength <- 1e5
matt_Q_model_none_noneFitR$mcmcParam$burnin <- 5e5

##Run the optimizer using the parameters listed in the hzar.fitRequest given.
matt_Q_model_free_bothFit <- hzar.doFit(matt_Q_model_free_bothFitR)
matt_Q_model_free_leftFit <- hzar.doFit(matt_Q_model_free_leftFitR)
matt_Q_model_free_rightFit <- hzar.doFit(matt_Q_model_free_rightFitR)
matt_Q_model_free_mirrorFit <- hzar.doFit(matt_Q_model_free_mirrorFitR)
matt_Q_model_free_noneFit <- hzar.doFit(matt_Q_model_free_noneFitR)

matt_Q_model_fixed_bothFit <- hzar.doFit(matt_Q_model_fixed_bothFitR)
matt_Q_model_fixed_leftFit <- hzar.doFit(matt_Q_model_fixed_leftFitR)
matt_Q_model_fixed_rightFit <- hzar.doFit(matt_Q_model_fixed_rightFitR)
matt_Q_model_fixed_mirrorFit <- hzar.doFit(matt_Q_model_fixed_mirrorFitR)
matt_Q_model_fixed_noneFit <- hzar.doFit(matt_Q_model_fixed_noneFitR)

matt_Q_model_none_bothFit <- hzar.doFit(matt_Q_model_none_bothFitR)
matt_Q_model_none_leftFit <- hzar.doFit(matt_Q_model_none_leftFitR)
matt_Q_model_none_rightFit <- hzar.doFit(matt_Q_model_none_rightFitR)
matt_Q_model_none_mirrorFit <- hzar.doFit(matt_Q_model_none_mirrorFitR)
matt_Q_model_none_noneFit <- hzar.doFit(matt_Q_model_none_noneFitR)

# ### plot model to look for run stability and convergence and return the mcmc data with an added a log likelihood column.
par(mar=c(1,1,1,1))
plot(hzar.mcmc.bindLL(matt_Q_model_free_bothFit))
plot(hzar.mcmc.bindLL(matt_Q_model_free_leftFit))
plot(hzar.mcmc.bindLL(matt_Q_model_free_rightFit))
plot(hzar.mcmc.bindLL(matt_Q_model_free_mirrorFit))
plot(hzar.mcmc.bindLL(matt_Q_model_free_noneFit))
# 
##group multiple fits of the same model and the same observation data into a single object
matt_Q_model_free_bothData <- hzar.dataGroup.add(matt_Q_model_free_bothFit)
matt_Q_model_free_leftData <- hzar.dataGroup.add(matt_Q_model_free_leftFit)
matt_Q_model_free_rightData <- hzar.dataGroup.add(matt_Q_model_free_rightFit)
matt_Q_model_free_mirrorData <- hzar.dataGroup.add(matt_Q_model_free_mirrorFit)
matt_Q_model_free_noneData <- hzar.dataGroup.add(matt_Q_model_free_noneFit)

matt_Q_model_fixed_bothData <- hzar.dataGroup.add(matt_Q_model_fixed_bothFit)
matt_Q_model_fixed_leftData <- hzar.dataGroup.add(matt_Q_model_fixed_leftFit)
matt_Q_model_fixed_rightData <- hzar.dataGroup.add(matt_Q_model_fixed_rightFit)
matt_Q_model_fixed_mirrorData <- hzar.dataGroup.add(matt_Q_model_fixed_mirrorFit)
matt_Q_model_fixed_noneData <- hzar.dataGroup.add(matt_Q_model_fixed_noneFit)

matt_Q_model_none_bothData <- hzar.dataGroup.add(matt_Q_model_none_bothFit)
matt_Q_model_none_leftData <- hzar.dataGroup.add(matt_Q_model_none_leftFit)
matt_Q_model_none_rightData <- hzar.dataGroup.add(matt_Q_model_none_rightFit)
matt_Q_model_none_mirrorData <- hzar.dataGroup.add(matt_Q_model_none_mirrorFit)
matt_Q_model_none_noneData <- hzar.dataGroup.add(matt_Q_model_none_noneFit)


### Generate a hzar.dataGroup object representing a fit of the null model to a hzar.obsData object
matt_Q_modelNull <- hzar.dataGroup.null(matt_Q)

##make list of cline models and null models
matt_Q_dGs <- list(cline_free_bothModel = matt_Q_model_free_bothData,
                  cline_free_leftModel = matt_Q_model_free_leftData,
                  cline_free_rightModel = matt_Q_model_free_rightData,
                  cline_free_mirrorModel = matt_Q_model_free_mirrorData,
                  cline_free_noneModel = matt_Q_model_free_noneData,
                  cline_fixed_bothModel = matt_Q_model_fixed_bothData,
                  cline_fixed_leftModel = matt_Q_model_fixed_leftData,
                  cline_fixed_rightModel = matt_Q_model_fixed_rightData,
                  cline_fixed_mirrorModel = matt_Q_model_fixed_mirrorData,
                  cline_fixed_noneModel = matt_Q_model_fixed_noneData,
                  cline_none_bothModel = matt_Q_model_none_bothData,
                  cline_none_leftModel = matt_Q_model_none_leftData,
                  cline_none_rightModel = matt_Q_model_none_rightData,
                  cline_none_mirrorModel = matt_Q_model_none_mirrorData,
                  cline_none_noneModel = matt_Q_model_none_noneData,
                  nullModel = matt_Q_modelNull)

##Collect optimizer output based on the same hzar.obsData object
matt_Q_oDG <- hzar.make.obsDataGroup(matt_Q_dGs)

##Set the names of the list of hzar.dataGroup objects contained in a hzar.obsDataGroup object using the names from either a named list of hzar.dataGroup objects or another hzar.obsDataGroup object.
matt_Q_oDG <- hzar.copyModelLabels(matt_Q_dGs,matt_Q_oDG)

##Plots a line representing the expected frequency versus distance for the given object. For hzar.dataGroup and hzar.obsDataGroup objects, plots the observed data backing the model. For hzar.obsDataGroup objects, plots the maximum likelihood cline for each model.
hzar.plot.cline(matt_Q_oDG)

##Calculate the AIC or corrected AIC score table for the given hzar.obsDataGroup object. There will be one score generated for each model associated with this object.
print(hzar.AICc.hzar.obsDataGroup(matt_Q_oDG))
# AICc
# cline_free_bothModel    37.84759
# cline_free_leftModel    32.93323
# cline_free_rightModel   32.77541
# cline_free_mirrorModel  33.01367
# cline_free_noneModel    28.41212
# cline_fixed_bothModel   32.92816
# cline_fixed_leftModel   37.27875
# cline_fixed_rightModel  28.50004
# cline_fixed_mirrorModel 29.98919
# cline_fixed_noneModel   32.87428
# cline_none_bothModel    32.53676
# cline_none_leftModel    43.87432
# cline_none_rightModel   28.25540
# cline_none_mirrorModel  29.36315
# cline_none_noneModel    39.47224
# nullModel               78.41916



### make em purty, get the range of parameter values that are within two log likelihood units of the maximum likelihood for center and width

min(hzar.AICc.hzar.obsDataGroup(paul_Q_oDG))
min(hzar.AICc.hzar.obsDataGroup(matt_Q_oDG))


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

pdf(file="matt_q_cline.pdf")
hzar.plot.fzCline(matt_Q_model_fixed_noneData,  fzCol = rgb(33/255, 145/255, 140/255, 0.3), pch = 20, col=rgb(33/255, 145/255, 140/255, 0.7), xlab = "Distance (km)", ylab = "Q value")
dev.off()
print(matt_Q_model_fixed_noneData$ML.cline$param.free$center)
#[1] -41.10436
print(matt_Q_model_fixed_noneData$ML.cline$param.free$width)
#[1] 138.0059
print(hzar.getLLCutParam(matt_Q_model_fixed_noneData,c("center","width")))
# center2LLLow center2LLHigh width2LLLow width2LLHigh
# 1    -69.69375      2.401899    79.31715     228.5547


## plotting the two clines together:

# first create an empty plot bounded by the biggest dataset
hzar.plot.obsData(matt_Q_model_fixed_noneData, col = "transparent", xlab = "Distance (km)", ylab = "Hybrid Index")

xSeries <- seq(from = par("usr")[1], to = par("usr")[2], 
               length.out = 109)
if (par("xaxs") == "r") 
  xSeries <- xSeries[2:108]

#add polygons and lines
paul_fzCline = hzar.getCredParamRed(paul_Q_model_fixed_noneData)
paul_fzCoor <- paul_fzCline$fzCline(xSeries)
polygon(x = c(paul_fzCoor$x, rev(paul_fzCoor$x)), y = c(paul_fzCoor$yMin, rev(paul_fzCoor$yMax)), border = rgb(68/255, 1/255, 84/255, 0.5), col=rgb(68/255, 1/255, 84/255, 0.5))
lines(x = xSeries, y = paul_Q_model_fixed_noneData$ML.cline$clineFunc(xSeries), col=rgb(68/255, 1/255, 84/255), lwd=2)
hzar.plot.obsData(paul_Q_model_fixed_noneData, fzCol = rgb(68/255, 1/255, 84/255, 0.3), pch = 20, col=rgb(68/255, 1/255, 84/255, 0.7), add=TRUE)

matt_fzCline = hzar.getCredParamRed(matt_Q_model_fixed_noneData)
matt_fzCoor <- matt_fzCline$fzCline(xSeries)
polygon(x = c(matt_fzCoor$x, rev(matt_fzCoor$x)), y = c(matt_fzCoor$yMin, rev(matt_fzCoor$yMax)), border = rgb(253/255, 231/255, 37/255, 0.5), col=rgb(253/255, 231/255, 37/255, 0.5))
lines(x = xSeries, y = matt_Q_model_fixed_noneData$ML.cline$clineFunc(xSeries), col=rgb(253/255, 231/255, 37/255), lwd=2)
hzar.plot.obsData(matt_Q_model_fixed_noneData, fzCol = rgb(68/255, 1/255, 84/255, 0.3), pch = 20, col=rgb(253/255, 231/255, 37/255, 0.7), add=TRUE)

# adding a legend
matt_col <- rgb(68/255, 1/255, 84/255)
paul_col <- rgb(253/255, 231/255, 37/255)

legend(0,1, legend=c("2004-2007", "2018-2021"), 
       col=c(paul_col,matt_col), pch = 16, cex=.8)





############## Plot all the clines together ##############
### remake all the models with the points subtracted around the centre
paul_centered <- paul_cline
paul_centered$distance <- paul_cline$dist-(paul_Q_model_fixed_noneData$ML.cline$param.free$center)
matt_centered <- matt_cline
matt_centered$distance <- matt_cline$dist-(matt_Q_model_fixed_noneData$ML.cline$param.free$center)

paul_Q_centered <-
  hzar.doMolecularData1DPops(paul_centered$distance,
                             paul_centered$q,
                             paul_centered$SampleID)
matt_Q_centered <-
  hzar.doMolecularData1DPops(matt_centered$distance,
                             matt_centered$q,
                             matt_centered$SampleID)


hzar.plot.obsData(paul_Q_centered)
hzar.plot.obsData(matt_Q_centered)


### make ClineMetaModel object
paul_Q_centered_model_fixed_none<- hzar.makeCline1DFreq(paul_Q_centered, scaling="none", tails="none")
matt_Q_centered_model_fixed_none <- hzar.makeCline1DFreq(matt_Q_centered, scaling="none", tails="right")

### set BoxReq
paul_Q_centered_model_fixed_none <- hzar.model.addBoxReq(paul_Q_centered_model_fixed_none,-100,60)
matt_Q_centered_model_fixed_none <- hzar.model.addBoxReq(matt_Q_centered_model_fixed_none,-100,550)

### generate hzar fitRequests
paul_Q_centered_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(model=paul_Q_centered_model_fixed_none, paul_Q_centered, verbose=FALSE)
matt_Q_centered_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(model=matt_Q_centered_model_fixed_none, matt_Q_centered, verbose=FALSE)

### set mcmc chainLength and burnin
paul_Q_centered_model_fixed_noneFitR$mcmcParam$chainLength <- 1e5
paul_Q_centered_model_fixed_noneFitR$mcmcParam$burnin <- 5e5
matt_Q_centered_model_fixed_noneFitR$mcmcParam$chainLength <- 1e5
matt_Q_centered_model_fixed_noneFitR$mcmcParam$burnin <- 5e5

### Run the optimizer using the parameters listed in the hzar.fitRequest given.
paul_Q_centered_model_fixed_noneFit <- hzar.doFit(paul_Q_centered_model_fixed_noneFitR)
matt_Q_centered_model_fixed_noneFit <- hzar.doFit(matt_Q_centered_model_fixed_noneFitR)

### group multiple fits of the same model and the same observation data into a single object
paul_Q_centered_model_fixed_noneData <- hzar.dataGroup.add(paul_Q_centered_model_fixed_noneFit)
matt_Q_centered_model_fixed_noneData <- hzar.dataGroup.add(matt_Q_centered_model_fixed_noneFit)

### get the range of parameter values that are within two log likelihood units of the maximum likelihood for center and width
print(hzar.getLLCutParam(paul_Q_centered_model_fixed_noneData,c("center","width")))
print(hzar.getLLCutParam(matt_Q_centered_model_fixed_noneData,c("center","width")))


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
matt_fzCline = hzar.getCredParamRed(matt_Q_centered_model_fixed_noneData)
matt_fzCoor <- matt_fzCline$fzCline(xSeries)

### plot the polygons and lines
# the stats::line() function uses the x and y coordinates from above to draw the line

pdf(file = "centered_clines.pdf")

# first create an empty plot bounded by the biggest dataset
hzar.plot.obsData(matt_Q_centered_model_fixed_noneData, col = "transparent")

#add polygons and lines
polygon(x = c(matt_fzCoor$x, rev(matt_fzCoor$x)), y = c(matt_fzCoor$yMin, rev(matt_fzCoor$yMax)), border = rgb(253/255, 231/255, 37/255, 0.7), col=rgb(253/255, 231/255, 37/255, 0.7))
lines(x = xSeries, y = matt_Q_centered_model_fixed_noneData$ML.cline$clineFunc(xSeries), col=rgb(253/255, 231/255, 37/255), lwd=2)
polygon(x = c(paul_fzCoor$x, rev(paul_fzCoor$x)), y = c(paul_fzCoor$yMin, rev(paul_fzCoor$yMax)), border = rgb(68/255, 1/255, 84/255, 0.75), col=rgb(68/255, 1/255, 84/255, 0.3))
lines(x = xSeries, y = paul_Q_centered_model_fixed_noneData$ML.cline$clineFunc(xSeries),  col=rgb(68/255, 1/255, 84/255), lwd=2)
dev.off()   

# do a zoomed in version on just the centers
pdf(file = "centered_clines_zoom.pdf")
