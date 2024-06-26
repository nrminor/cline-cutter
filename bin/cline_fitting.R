#!/usr/bin/env Rscript

###################### make hzar clines ######################

# retrieve positional command line arguments
args <- commandArgs(TRUE)

options(warn = 1)

suppressPackageStartupMessages({

  # libraries
  require(ggplot2)
  require(sp)
  require(geosphere)
  require(hzar)
  require(viridis)
  require(readr)
  require(stringr)
  require(dplyr)
  require(readxl)
  require(forcats)
  require(rhdf5)
  require(SimDesign)

  # print package versions
  sessioninfo::session_info("attached")

})

#### INPUT PATHS ####
if (is.function(args)) {
  meta_path <- "samplesheet.xlsx"
  interactive <- TRUE
  print("Running interactively.")
} else if (length(args) == 0) {
  meta_path <- "samplesheet.xlsx"
  interactive <- FALSE
  print("Running non-interactively.")
} else {
  meta_path <- args[1]
  interactive <- FALSE
  print("Running non-interactively.")
}
if (!(meta_path %in% list.files("."))) {
  stop("A metadata file is required.")
}
hdf_5_files <- sort(list.files(".", ".hdf5"))
if (length(hdf_5_files) < 3) {
  stop("3 hdf5 replicates must be provided in the current working directory.")
}
hdf5_1 <- hdf_5_files[1]
hdf5_2 <- hdf_5_files[2]
hdf5_3 <- hdf_5_files[3]
sample_rep <- str_remove(hdf5_1, "_1.hdf5")


#### LOAD METADATA ####
#### ---------------------------------------------------------------------------
#### #

# identify downsample regime
if (grepl("random", hdf5_1)) {
  downsample_regime <- "random"
} else if (grepl("uneven", hdf5_1)) {
  downsample_regime <- "uneven"
} else {
  downsample_regime <- "even"
}
samples_path <- paste(sample_rep, ".txt", sep = "")

# filter to the current subset and remove space from first colname
metadata <- read_tsv(samples_path, col_names = "Sample ID",
                     show_col_types = FALSE, trim_ws = TRUE) %>%
  left_join(read_excel(meta_path, trim_ws = TRUE), by = "Sample ID",
            keep = FALSE) %>%
  rename(SampleID = `Sample ID`)

stopifnot(all(metadata$Latitude > 0) || all(metadata$Latitude < 0))
stopifnot(all(metadata$Longitude > 0) || all(metadata$Longitude < 0))

# rename to conform with Paul's conventions
metadata$SampleID <- regmatches(metadata$SampleID, regexpr("[[:digit:]]*",
                                                           metadata$SampleID))
sample_ids <- metadata %>%
  select(SampleID, Population)
sample_coords <- metadata %>%
  select(SampleID, Latitude, Longitude)

# --------------------------------------------------------------------------- #


#### LOAD HDF5 FILES WITH ANCESTRY/HYBRID SCORES ####
#### ---------------------------------------------------------------------------
#### #
qrep1 <- h5read(hdf5_1, "q")
qrep2 <- h5read(hdf5_2, "q")
qrep3 <- h5read(hdf5_3, "q")
q1_mat <- qrep1[, 1, ]
q2_mat <- qrep2[, 1, ]
q3_mat <- qrep3[, 1, ]
q_long <- rbind(q1_mat, q2_mat, q3_mat)
q_means <- apply(q_long, 2, mean)

# make sure we have the number of mean Q-scores that we expect (1 per sample)
stopifnot(length(q_means) == nrow(metadata))
stopifnot(length(q_means) == nrow(sample_ids))
stopifnot(length(q_means) == nrow(sample_coords))

# merging sample id and gbs data
sample_ids$q <- q_means

# merging dataframes by sample id
sample_ids <- merge(sample_ids, sample_coords, by = "SampleID", sort = FALSE)

## could we also look at the frequency of interspecific ancestry in certain
## areas over time?
site_ancestry <- aggregate(q ~ Longitude, data = sample_ids, mean)
pdf(paste(downsample_regime, "regime_q_by_longitude.pdf", sep = "_"),
    height = 6, width = 9)
plot(q ~ Longitude, data = site_ancestry, pch = 16)
grid()
dev.off() |> invisible()
# --------------------------------------------------------------------------- #


#### COLLATE DISTANCES ####
#### ---------------------------------------------------------------------------
#### #

# aggregating mean longitude by population
distances <- aggregate(Longitude ~ Population, data = metadata, mean)

if (interactive == TRUE) {
  plot(Population ~ Longitude, data = distances, pch = 16)
  grid()
}

### converting longitudinal points to distance along transect first need to
### know mean latitude of sampling locations
mean_sample_lat <- mean(sample_ids$Latitude)

# conversion for 1 deg longitude to km: 1 deg = 111.320*cos(latitude) km
long_degree_km <- 111.32 * cos(mean_sample_lat)
distances$distance <- long_degree_km * (abs((min(distances$Longitude) -
                                               distances$Longitude)))
distances <- distances %>%
  select(Population, distance)

# merging genomic hybrid index for every sample with distances
sample_ids_distance <- merge(sample_ids, distances, by = "Population")
if (interactive == TRUE) {
  plot(q ~ distance, data = sample_ids_distance)
  grid()
}

## need to get a dataset of mean q and bird count at each site for both
## collection datasets
site_counts <- aggregate(SampleID ~ distance, data = sample_ids_distance,
                         length)
site_ancestry <- aggregate(q ~ distance, data = sample_ids_distance, mean)
cline <- merge(site_ancestry, site_counts)

# need to flip q values for one dataset so the cline will be in the same
# direction
if (interactive == TRUE) {
  plot(q ~ distance, data = cline, xlab = "inter-population distance",
       ylab = "Q score", pch = 16)
  grid()
}

# cline$q <- 1 - cline$q plot(q ~ distance, data = cline)


#### CLINE MODELING STEPS
#### ---------------------------------------------------------------------------
#### # make data object of allele frequency data
Q <- hzar.doMolecularData1DPops(cline$distance, cline$q, cline$SampleID)

### Plot the associated observed frequency versus distance
if (interactive == TRUE) {
  hzar.plot.obsData(Q)
}

### Construct a clineMetaModel object for use with
### hzar.first.fitRequest.old.ML
{
  Q_model_free_both <- hzar.makeCline1DFreq(Q, scaling = "free", tails = "both")
  Q_model_free_none <- hzar.makeCline1DFreq(Q, scaling = "free", tails = "none")
  Q_model_free_right <- hzar.makeCline1DFreq(Q, scaling = "free", tails = "right")
  Q_model_free_left <- hzar.makeCline1DFreq(Q, scaling = "free", tails = "left")
  Q_model_free_mirror <- hzar.makeCline1DFreq(Q, scaling = "free", tails = "mirror")

  # #
  Q_model_fixed_both <- hzar.makeCline1DFreq(Q, scaling = "fixed", tails = "both")
  Q_model_fixed_none <- hzar.makeCline1DFreq(Q, scaling = "fixed", tails = "none")
  Q_model_fixed_right <- hzar.makeCline1DFreq(Q, scaling = "fixed", tails = "right")
  Q_model_fixed_left <- hzar.makeCline1DFreq(Q, scaling = "fixed", tails = "left")
  Q_model_fixed_mirror <- hzar.makeCline1DFreq(Q, scaling = "fixed", tails = "mirror")

  # #
  Q_model_none_both <- hzar.makeCline1DFreq(Q, scaling = "none", tails = "both")
  Q_model_none_none <- hzar.makeCline1DFreq(Q, scaling = "none", tails = "none")
  Q_model_none_right <- hzar.makeCline1DFreq(Q, scaling = "none", tails = "right")
  Q_model_none_left <- hzar.makeCline1DFreq(Q, scaling = "none", tails = "left")
  Q_model_none_mirror <- hzar.makeCline1DFreq(Q, scaling = "none", tails = "mirror")
}


### The intent of these methods is to assist the optimizer in exploring the
### model parameter space by instructing it to ignore models that are not
### interesting. For example, if all of the sampled localities are in a
### region 100km wide, then a cline width of 110km is probably not
### interesting. A cline width of 500km in that scenario would definitely
### not be interesting at all.
{
  Q_model_free_both <- hzar.model.addBoxReq(Q_model_free_both, -100, 900)
  Q_model_free_left <- hzar.model.addBoxReq(Q_model_free_left, -100, 900)
  S1_22a_model_free_right <- hzar.model.addBoxReq(Q_model_free_right, -100, 900)
  S1_22a_model_free_mirror <- hzar.model.addBoxReq(Q_model_free_mirror, -100, 900)
  S1_22a_model_free_none <- hzar.model.addBoxReq(Q_model_free_none, -100, 900)

  Q_model_fixed_both <- hzar.model.addBoxReq(Q_model_fixed_both, -100, 900)
  Q_model_fixed_left <- hzar.model.addBoxReq(Q_model_fixed_left, -100, 900)
  Q_model_fixed_right <- hzar.model.addBoxReq(Q_model_fixed_right, -100, 900)
  Q_model_fixed_mirror <- hzar.model.addBoxReq(Q_model_fixed_mirror, -100, 900)
  Q_model_fixed_none <- hzar.model.addBoxReq(Q_model_fixed_none, -100, 900)

  Q_model_none_both <- hzar.model.addBoxReq(Q_model_none_both, -100, 900)
  Q_model_none_left <- hzar.model.addBoxReq(Q_model_none_left, -100, 900)
  Q_model_none_right <- hzar.model.addBoxReq(Q_model_none_right, -100, 900)
  Q_model_none_mirror <- hzar.model.addBoxReq(Q_model_none_mirror, -100, 900)
  Q_model_none_none <- hzar.model.addBoxReq(Q_model_none_none, -100, 900)
}


### cline model fitting generate an hzar.fitRequest object suitable for
### hzar.doFit
{
  Q_model_free_bothFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_free_both, Q, verbose = FALSE)
  Q_model_free_leftFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_free_left, Q, verbose = FALSE)
  Q_model_free_rightFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_free_right, Q, verbose = FALSE)
  Q_model_free_mirrorFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_free_mirror, Q, verbose = FALSE)
  Q_model_free_noneFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_free_none, Q, verbose = FALSE)

  Q_model_fixed_bothFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_fixed_both, Q, verbose = FALSE)
  Q_model_fixed_leftFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_fixed_left, Q, verbose = FALSE)
  Q_model_fixed_rightFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_fixed_right, Q, verbose = FALSE)
  Q_model_fixed_mirrorFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_fixed_mirror, Q, verbose = FALSE)
  Q_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_fixed_none, Q, verbose = FALSE)

  Q_model_none_bothFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_none_both, Q, verbose = FALSE)
  Q_model_none_leftFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_none_left, Q, verbose = FALSE)
  Q_model_none_rightFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_none_right, Q, verbose = FALSE)
  Q_model_none_mirrorFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_none_mirror, Q, verbose = FALSE)
  Q_model_none_noneFitR <- hzar.first.fitRequest.old.ML(
    model = Q_model_none_none, Q, verbose = FALSE)
} |> quiet() |> invisible()

## set mcmc chain length and burn in
{
  Q_model_free_bothFitR$mcmcParam$chainLength <- 1e+05
  Q_model_free_bothFitR$mcmcParam$burnin <- 5e+05
  Q_model_free_leftFitR$mcmcParam$chainLength <- 1e+05
  Q_model_free_leftFitR$mcmcParam$burnin <- 5e+05
  Q_model_free_rightFitR$mcmcParam$chainLength <- 1e+05
  Q_model_free_rightFitR$mcmcParam$burnin <- 5e+05
  Q_model_free_mirrorFitR$mcmcParam$chainLength <- 1e+05
  Q_model_free_mirrorFitR$mcmcParam$burnin <- 5e+05
  Q_model_free_noneFitR$mcmcParam$chainLength <- 1e+05
  Q_model_free_noneFitR$mcmcParam$burnin <- 5e+05

  Q_model_fixed_bothFitR$mcmcParam$chainLength <- 1e+05
  Q_model_fixed_bothFitR$mcmcParam$burnin <- 5e+05
  Q_model_fixed_leftFitR$mcmcParam$chainLength <- 1e+05
  Q_model_fixed_leftFitR$mcmcParam$burnin <- 5e+05
  Q_model_fixed_rightFitR$mcmcParam$chainLength <- 1e+05
  Q_model_fixed_rightFitR$mcmcParam$burnin <- 5e+05
  Q_model_fixed_mirrorFitR$mcmcParam$chainLength <- 1e+05
  Q_model_fixed_mirrorFitR$mcmcParam$burnin <- 5e+05
  Q_model_fixed_noneFitR$mcmcParam$chainLength <- 1e+05
  Q_model_fixed_noneFitR$mcmcParam$burnin <- 5e+05

  Q_model_none_bothFitR$mcmcParam$chainLength <- 1e+05
  Q_model_none_bothFitR$mcmcParam$burnin <- 5e+05
  Q_model_none_leftFitR$mcmcParam$chainLength <- 1e+05
  Q_model_none_leftFitR$mcmcParam$burnin <- 5e+05
  Q_model_none_rightFitR$mcmcParam$chainLength <- 1e+05
  Q_model_none_rightFitR$mcmcParam$burnin <- 5e+05
  Q_model_none_mirrorFitR$mcmcParam$chainLength <- 1e+05
  Q_model_none_mirrorFitR$mcmcParam$burnin <- 5e+05
  Q_model_none_noneFitR$mcmcParam$chainLength <- 1e+05
  Q_model_none_noneFitR$mcmcParam$burnin <- 5e+05
}

## Run the optimizer using the parameters listed in the hzar.fitRequest
## given.
cat("\nFitting model labeled 'cline_free_bothModel':")
Q_model_free_bothFit <- hzar.doFit(Q_model_free_bothFitR)
cat("\nFitting model labeled 'cline_free_leftModel:")
Q_model_free_leftFit <- hzar.doFit(Q_model_free_leftFitR)
cat("\nFitting model labeled 'cline_free_rightModel':")
Q_model_free_rightFit <- hzar.doFit(Q_model_free_rightFitR)
cat("\nFitting model labeled 'cline_free_mirrorModel':")
Q_model_free_mirrorFit <- hzar.doFit(Q_model_free_mirrorFitR)
cat("\nFitting model labeled 'cline_free_noneModel':")
Q_model_free_noneFit <- hzar.doFit(Q_model_free_noneFitR)

cat("\nFitting model labeled 'cline_fixed_bothModel':")
Q_model_fixed_bothFit <- hzar.doFit(Q_model_fixed_bothFitR)
cat("\nFitting model labeled 'cline_fixed_leftModel':")
Q_model_fixed_leftFit <- hzar.doFit(Q_model_fixed_leftFitR)
cat("\nFitting model labeled 'cline_fixed_rightModel':")
Q_model_fixed_rightFit <- hzar.doFit(Q_model_fixed_rightFitR)
cat("\nFitting model labeled 'cline_fixed_mirrorModel':")
Q_model_fixed_mirrorFit <- hzar.doFit(Q_model_fixed_mirrorFitR)
cat("\nFitting model labeled 'cline_fixed_noneModel':")
Q_model_fixed_noneFit <- hzar.doFit(Q_model_fixed_noneFitR)

cat("\nFitting model labeled 'cline_none_bothModel':")
Q_model_none_bothFit <- hzar.doFit(Q_model_none_bothFitR)
cat("\nFitting model labeled 'cline_none_leftModel':")
Q_model_none_leftFit <- hzar.doFit(Q_model_none_leftFitR)
cat("\nFitting model labeled 'cline_none_rightModel':")
Q_model_none_rightFit <- hzar.doFit(Q_model_none_rightFitR)
cat("\nFitting model labeled 'cline_none_mirrorModel':")
Q_model_none_mirrorFit <- hzar.doFit(Q_model_none_mirrorFitR)
cat("\nFitting model labeled 'cline_none_noneModel':")
Q_model_none_noneFit <- hzar.doFit(Q_model_none_noneFitR)

# ### plot model to look for run stability and convergence and return the
# mcmc data with an added a log likelihood column.
pdf(
  paste(downsample_regime, "free_both_model_trace.pdf", sep = "_"),
  height = 11, width = 8.5
)
plot(hzar.mcmc.bindLL(Q_model_free_bothFit))
dev.off() |> invisible()
pdf(
  paste(downsample_regime, "free_left_model_trace.pdf", sep = "_"),
  height = 11, width = 8.5
)
plot(hzar.mcmc.bindLL(Q_model_free_leftFit))
dev.off() |> invisible()
pdf(
  paste(downsample_regime, "free_right_model_trace.pdf", sep = "_"),
  height = 11, width = 8.5
)
plot(hzar.mcmc.bindLL(Q_model_free_rightFit))
dev.off() |> invisible()
pdf(
  paste(downsample_regime, "free_mirror_model_trace.pdf", sep = "_"),
  height = 11, width = 8.5
)
plot(hzar.mcmc.bindLL(Q_model_free_mirrorFit))
dev.off() |> invisible()
pdf(
  paste(downsample_regime, "free_none_model_trace.pdf", sep = "_"),
  height = 11, width = 8.5
)
plot(hzar.mcmc.bindLL(Q_model_free_noneFit))
dev.off() |> invisible()

### group multiple fits of the same model and the same observation data
### into a single object
Q_model_free_bothData <- hzar.dataGroup.add(Q_model_free_bothFit)
Q_model_free_leftData <- hzar.dataGroup.add(Q_model_free_leftFit)
Q_model_free_rightData <- hzar.dataGroup.add(Q_model_free_rightFit)
Q_model_free_mirrorData <- hzar.dataGroup.add(Q_model_free_mirrorFit)
Q_model_free_noneData <- hzar.dataGroup.add(Q_model_free_noneFit)

Q_model_fixed_bothData <- hzar.dataGroup.add(Q_model_fixed_bothFit)
Q_model_fixed_leftData <- hzar.dataGroup.add(Q_model_fixed_leftFit)
Q_model_fixed_rightData <- hzar.dataGroup.add(Q_model_fixed_rightFit)
Q_model_fixed_mirrorData <- hzar.dataGroup.add(Q_model_fixed_mirrorFit)
Q_model_fixed_noneData <- hzar.dataGroup.add(Q_model_fixed_noneFit)

Q_model_none_bothData <- hzar.dataGroup.add(Q_model_none_bothFit)
Q_model_none_leftData <- hzar.dataGroup.add(Q_model_none_leftFit)
Q_model_none_rightData <- hzar.dataGroup.add(Q_model_none_rightFit)
Q_model_none_mirrorData <- hzar.dataGroup.add(Q_model_none_mirrorFit)
Q_model_none_noneData <- hzar.dataGroup.add(Q_model_none_noneFit)


### Generate a hzar.dataGroup object representing a fit of the null model
### to a hzar.obsData object
Q_modelNull <- hzar.dataGroup.null(Q)

## make list of cline models and null models
Q_dGs <- list(cline_free_bothModel = Q_model_free_bothData,
              cline_free_leftModel = Q_model_free_leftData,
              cline_free_rightModel = Q_model_free_rightData,
              cline_free_mirrorModel = Q_model_free_mirrorData,
              cline_free_noneModel = Q_model_free_noneData,
              cline_fixed_bothModel = Q_model_fixed_bothData,
              cline_fixed_leftModel = Q_model_fixed_leftData,
              cline_fixed_rightModel = Q_model_fixed_rightData,
              cline_fixed_mirrorModel = Q_model_fixed_mirrorData,
              cline_fixed_noneModel = Q_model_fixed_noneData,
              cline_none_bothModel = Q_model_none_bothData,
              cline_none_leftModel = Q_model_none_leftData,
              cline_none_rightModel = Q_model_none_rightData,
              cline_none_mirrorModel = Q_model_none_mirrorData,
              cline_none_noneModel = Q_model_none_noneData,
              nullModel = Q_modelNull)

### Collect optimizer output based on the same hzar.obsData object
Q_oDG <- hzar.make.obsDataGroup(Q_dGs) |> quiet()

### Set the names of the list of hzar.dataGroup objects contained in a
### hzar.obsDataGroup object using the names from either a named list of
### hzar.dataGroup objects or another hzar.obsDataGroup object.
Q_oDG <- hzar.copyModelLabels(Q_dGs, Q_oDG)


### Plots a line representing the expected frequency versus distance for the
### given object. For hzar.dataGroup and hzar.obsDataGroup objects, plots the
### observed data backing the model. For hzar.obsDataGroup objects, plots the
### maximum likelihood cline for each model.
pdf(paste(downsample_regime, "model_cline_comparison.pdf", sep = "_"),
    height = 8, width = 8)
all_clines_plot <- hzar.plot.cline(Q_oDG, pch = 16) ; remove(all_clines_plot)
dev.off() |> invisible()

### Calculate the AIC or corrected AIC score table for the given
### hzar.obsDataGroup object. There will be one score generated for each model
### associated with this object.
model_aic_df <- data.frame(hzar.AICc.hzar.obsDataGroup(Q_oDG))
model_aic_df$model <- row.names(model_aic_df)
model_aic_df <- model_aic_df %>%
  select(model, AICc)
write_tsv(model_aic_df, paste(downsample_regime, "regime_model_aic.tsv",
                              sep = "_"))

#### plot clines
pdf(file = paste(downsample_regime, "q_cline.pdf", sep = "_"))
par(mar = c(1, 1, 1, 1))
hzar.plot.fzCline(Q_model_fixed_noneData, fzCol = rgb(68/255, 1/255,
                                                      84/255, 0.3),
                  pch = 20, col = rgb(68/255, 1/255, 84/255, 0.7),
                  xlab = "Distance (km)", ylab = "Q value")
dev.off() |> invisible()
print(
  paste("Estimated cline center for the",
        downsample_regime,
        "downsampling regime:",
        Q_model_fixed_noneData$ML.cline$param.free$center,
        sep = " ")
)
print(
  paste("Estimated cline width for the",
        downsample_regime,
        "downsampling regime:",
        Q_model_fixed_noneData$ML.cline$param.free$width,
        sep = " ")
)
# print("Two log likelihood unit range for center and width estimates:")
# print(hzar.getLLCutParam(Q_model_fixed_noneData, c("center", "width")))


############## Plot all the clines together ############## remake all the
############## models with the points subtracted around the centre
centered <- cline
centered$distance <- cline$dist -
  (Q_model_fixed_noneData$ML.cline$param.free$center)

Q_centered <- hzar.doMolecularData1DPops(centered$distance, centered$q,
                                         centered$SampleID)
if (interactive == TRUE) {
  hzar.plot.obsData(Q_centered)
}

### make ClineMetaModel object
Q_centered_model_fixed_none <- hzar.makeCline1DFreq(Q_centered,
                                                    scaling = "none", tails = "none")

### set BoxReq
Q_centered_model_fixed_none <- hzar.model.addBoxReq(Q_centered_model_fixed_none,
                                                    -100, 60)

### generate hzar fitRequests
Q_centered_model_fixed_noneFitR <- hzar.first.fitRequest.old.ML(
  model = Q_centered_model_fixed_none, Q_centered, verbose = FALSE) |> quiet()

### set mcmc chainLength and burnin
Q_centered_model_fixed_noneFitR$mcmcParam$chainLength <- 1e+05
Q_centered_model_fixed_noneFitR$mcmcParam$burnin <- 5e+05

### Run the optimizer using the parameters listed in the hzar.fitRequest given.
Q_centered_model_fixed_noneFit <- hzar.doFit(Q_centered_model_fixed_noneFitR) |> quiet()

### group multiple fits of the same model and the same observation data into a
### single object
Q_centered_model_fixed_noneData <- hzar.dataGroup.add(Q_centered_model_fixed_noneFit)

### get the range of parameter values that are within two log likelihood units
### of the maximum likelihood for center and width
print(paste("Two log-likelihood units around cline centers and widths for the",
            downsample_regime, "downsampling regime:", sep = " "))
print(hzar.getLLCutParam(Q_centered_model_fixed_noneData, c("center", "width")))


### MANNY BOEHM CODE hzar.overPlot.fzCline doesn't contain any code for
### creating plots :( so, let's extract info from 'hzar.dataGroup' objects in
### the same way that hzar.plot.fzCline does

# -------------------- 95% CI

# hzar.plot.fzCline() fits a line (and the 95% CI) to the points by first
# generating a series of 107 x coordinates bounded by the current plot window
# later, the y coordinates are estimated by accessing the function describing
# the cline in the dataGroup object sometimes this dies and makes the polygons
# real small but if I rerun this code a couple times they go backto normal the
# numbers are from hzar, no relation to my observations being 108 in caor
xSeries <- seq(from = par("usr")[1], to = par("usr")[2], length.out = 109)
if (par("xaxs") == "r") {
  xSeries <- xSeries[2:108]
}

# hzar.plot.fzCline() plots the 95% CI using graphics::polygon() you can
# customize the colours here
fzCline <- hzar.getCredParamRed(Q_centered_model_fixed_noneData)
fzCoor <- fzCline$fzCline(xSeries)

### plot the polygons and lines the stats::line() function uses the x and y
### coordinates from above to draw the line

pdf(file = paste(downsample_regime, "centered_clines.pdf", sep = "_"))
hzar.plot.obsData(Q_centered_model_fixed_noneData, col = "transparent")
polygon(x = c(fzCoor$x, rev(fzCoor$x)), y = c(fzCoor$yMin, rev(fzCoor$yMax)),
        border = rgb(68/255, 1/255, 84/255, 0.75),
        col = rgb(68/255, 1/255, 84/255, 0.3))
lines(x = xSeries, y = Q_centered_model_fixed_noneData$ML.cline$clineFunc(xSeries),
      col = rgb(68/255, 1/255, 84/255), lwd = 2)
dev.off() |> invisible()
