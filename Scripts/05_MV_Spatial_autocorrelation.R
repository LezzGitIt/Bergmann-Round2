## Bergmann's rule across full-annual-cycle in Caprimulgids##
## Model validation (mv) Spatial autocorrelation SHORTENED 05 -- 
## This is a shortened script, showing just the things that were important / actually worked in terms of justifying that we didn't see spatial autocorrelation in our data

## Contents ##
## This script 1)	Generates empirical variograms & bubble plots to visualize spatial autocorrelation, 2)	Calculates Moran’s I & correlograms, 3) implements gls() models to control for SAC, varying input parameters, 4) determine how modeled parameter estimates change 

## Variograms 
# 1. Use the top models from each Spp * DV combo
# 2. Make create_variog() function, which cycles through maximum distance parameter (argument "max.dist")
# 3. Apply the function to empirical data, varying the 'option' argument
# 4. Generate empirical variogram plots

## Moran's I: Vary 'type' argument, controlling the connectivity network & weighting of points.
# 1. For each 'type', we do the following: Create a connectivity network, calculate morans I & accompanying metrics , and create a dataframe
# 2. Join all dataframes
# 3. Some exploratory work showing that different weighting of points doesn't influence the empirical Moran's I estimate
# 4. Create correlograms with both the raw size data & the residuals

## GLS models: 
# 1. Leave default arguments on the gls() function, letting gls() estimate the parameters (range & nugget) that define the correlation structure
# 2. Manually vary range & nugget values
# 3. Compare AIC values & parm estimates between models 

## Lessons learned: 
# The predominant confusion came as I was 1) examining spatial autocorrelation at too small of a spatial scale (tens or hundreds of km, instead of thousands), & then 2) using techniques meant for data with SAC when there was very little SAC in my data. These 2 issues resulted in highly unstable results (e.g., modeled parameter estimates & standard errors) that were confusing & challenging to interpret. The key resolution was 1) using default parameters for the variograms (which examine SAC at a large spatial scale), & 2) using the variograms to then set the initialization values for the gls() function (although the gls() function was probably not necessary at all). 

## Additional resources
# See extended script in the 'Extra_scripts' folder, or the Github commit on 12/4/24 
# Supporting Information of the Bergmann’s rule annual cycle manuscript 
# Diniz-Filho, J. A. F., L. M. Bini, and B. A. Hawkins (2003). Spatial autocorrelation and red herrings in geographical ecology. Global Ecology and Biogeography 12:53–64.

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(adespatial)
library(spdep)
library(geoR)
library(ggpubr)
library(sp)
library(cowplot)
library(AICcmodavg)
library(nlme)
library(broom.mixed)
ggplot2::theme_set(theme_cowplot())

load("Rdata/Capri_dfs_07.09.24.Rdata")
#load("Rdata/Spatial_autocorrelation_11.12.24.Rdata")
# vignette(package = "spdep", 'nb')


# Variograms -------------------------------------------------------------
# Remember, semivariance is a measure of how the difference in values of a variable changes as the distance between pairs of points increases. So if points close together have similar values, the semivariance will be small, indicating strong spatial correlation or similarity.

# NOTE:: Following along from the following website:
# https://www.paulamoraga.com/book-spatial/geostatistical-data-1.html

# >nlme::Variogram.gls ----------------------------------------------------
## Understand the different spatial correlation structures & how they are parameterized. This is also nice given the variogram is from nlme, so we can be sure that it is 'speaking the same language' as nlme::gls()
?Variogram.gls

# Example code from Zuur pg 167, can play around with parameters here
cor_str_long <- c("spherical", "exponential", "gaussian", "linear", "rational")
my.df <- data.frame(D = seq(0, 5, .1))
map(cor_str_long, \(cor_str){
  cor_fun <- corSpatial(value = c(1), type = cor_str, form = ~D, nugget = FALSE)
  init <- Initialize(cor_fun, my.df)
  varg <- Variogram(init)
  rbind(cor_fun, init) # Notice that the range & nuggets are exponentiated & then logged again
  plot(varg, smooth = FALSE, type = "l", col = 1) # The c.gaus parameters match the variogram plot better
})

# >Prep data ------------------------------------------------------------
# The overall top models, where Whip mass is with temperature
top.mods.chr <- map(.x = aictab_list3$Breeding, \(aic.tab){
  top.mods.chr <- aic.tab %>%
    slice_head() %>%
    pull(Modnames)
})

# Create lm objects
top.mods <- map2(top.mods.chr, njdf.l.br.age, \(mods, df){
  lm(as.formula(mods), data = df)
})

# Create an Sf object with jittered coordinates in meters (projected)
# Projections suggested from ChatGPT for distance calculation based on the extent of each data set
projections <- c(rep(5070, 2), rep(3035, 2), rep(5070, 2))

# NOTE:: The coordinates in the dataframe are still in decimal degrees, but the sf geometries have been transformed (into meters)
njdf_sf <- map2(njdf.l.br.age, projections, \(df, proj){
  # Prep spatial data
  sf_df <- st_as_sf(df,
    coords = c("B.Long", "B.Lat"),
    crs = 4326, remove = FALSE
  )
  # Jitter data as some birds were inevitably captured at the same location
  sf_df <- st_jitter(sf_df, amount = 0.001)
  sf_df <- sf_df %>% st_transform(proj) # Project so the units of the variogram are in meters
})

## Generate values for create_variog function
## Distances -- What is the appropriate spatial scale?  Ultimately, I think we are only worried about correlation in residuals at large spatial scales.. For example, food abundance -> body mass… If there is variation in food within a site, we would expect some birds to be larger than predicted by the model (residuals above the fitted line) & some birds to be smaller than predicted by the model (residuals below the fitted line). On the other hand, if there is variation in food between sites, then all birds at a given site would be above or below the expected value (the fitted model line), resulting in spatial autocorrelation. The Diniz-Filho article helped that thousands of kilometers is the appropriate scale for examining autocorrelation for these macrogeographical questions.

## DVs
DVs <- rep(c("Mass.comb", "Wing.comb"), 3)

## Define create_variog() function to cycle through distances
# NOTE:: DVs argument allows for testing the raw size data if desired
create_variog <- function(mods, dfs, DVs = NULL, method = "bin", band = NULL) {
  DVs <- if (is.null(DVs)) rep(list(NULL), length(mods)) else DVs
  pmap(.l = list(mods, dfs, DVs), \(mod, df, DV){
    tm.res <- rstandard(mod) # Extract standardized residuals from top model
    if (!is.null(DV)) {
      # For coding convenience, just overwrite tm.res even though name is no longer applicable
      tm.res <- df[[DV]]
    }

    # Prepare band for smooth option (ultiply by 1000 given units are in meters)
    band_value <- if (method == "smooth") {
      band * 1000
    } else {
      NULL
    }

    # Run variogram.
    geoR::variog(
      coords = st_coordinates(df), # Switch between variog & variog4
      data = tm.res,
      option = method,
      # pairs.min = 10,
      # breaks = seq(0, dist * 1000, l = 30),
      band = band_value, # ,
      max.dist = dist * 1000
    ) # km * 1000 m / km
    # tol=pi/8)
  })
}
# Modified function that runs with max.dist argument set to the default
create_variog_def_dist <- function(mods, dfs, DVs = NULL, dist = NULL, method = "bin", band = NULL) {
  DVs <- if (is.null(DVs)) rep(list(NULL), length(mods)) else DVs
  pmap(.l = list(mods, dfs, DVs), \(mod, df, DV){
    tm.res <- rstandard(mod) # Extract standardized residuals from top model
    if (!is.null(DV)) {
      # For coding convenience, just overwrite tm.res even though name is no longer applicable
      tm.res <- df[[DV]]
    }

    # Prepare band for smooth option (ultiply by 1000 given units are in meters)
    band_value <- if (method == "smooth") {
      band * 1000
    } else {
      NULL
    }

    # Run variogram.
    geoR::variog(
      coords = st_coordinates(df), # Switch between variog & variog4
      data = tm.res,
      option = method,
      # pairs.min = 10,
      # breaks = seq(0, dist * 1000, l = 30),
      band = band_value
    )
    # tol=pi/8)
  })
}

# >Calculate variograms ---------------------------------------------------
# Run with method = "bin & default parameters for range (need to remove max.dist argument in create_variog function)
variog_def_dist <- create_variog_def_dist(top.mods, njdf_sf, method = "bin") 
variog_def_dist_raw <- create_variog_def_dist(top.mods, njdf_sf, method = "bin", DVs = DVs) 

# >Identify starting values for gls() models -------------------------------------------------------------
# Use the variogram to identify starting values for gls() models to initialize at
variog_parms <- map(variog_def_dist, \(varg){
  data.frame(
    sv_y = varg$v[1], # sv_y = semi-variance y-intercept (first value, not techinically y-intercept)
    sill = max(varg$v)
  ) %>%
    mutate(nugget_start = sv_y / sill)
})

# >Plot variograms --------------------------------------------------------
# Plot raw data & with default distances
pdf(file = "Plots/Spatial_autocorrelation/Default_parms/gg_Variograms_bin_default_raw2.pdf", 
    width = 8.5, height = 11, bg = "white")
plot_list_raw <- imap(variog_def_dist_raw, \(varg, names) {
  varg_df <- data.frame(u = varg$u * 0.001, v = varg$v)
  
  # Generate the ggplot for each variogram
  ggplot(varg_df, aes(x = u, y = v)) +
    geom_point() + # Scatter plot of the variogram
    geom_smooth(se = FALSE, color = "red") + # Add loess smoother
    labs(title = names, x = "Distance (km)", y = "Semivariance") +
    ylim(0, max(varg$v, na.rm = TRUE))
})
ggarrange(plotlist = plot_list_raw, ncol = 2, nrow = 3)
dev.off()

# Plot residuals & with the default distances
pdf(file = "Plots/Spatial_autocorrelation/Default_parms/gg_Variograms_bin_default_resid2.pdf", 
    width = 8.5, height = 11, bg = "white")
plot_list_resid <- imap(variog_def_dist, \(varg, names) {
  varg_df <- data.frame(u = varg$u * 0.001, v = varg$v)

  # Generate the ggplot for each variogram
  ggplot(varg_df, aes(x = u, y = v)) +
    geom_point() + # Scatter plot of the variogram
    geom_smooth(se = FALSE, color = "red") + # Add loess smoother
    labs(title = names, x = "Distance (km)", y = "Semivariance") +
    ylim(0, max(varg$v, na.rm = TRUE))
})
ggarrange(plotlist = plot_list_resid, ncol = 2, nrow = 3)
dev.off()

# NOTE:: When using the raw size data (DVs = DVs in create_variog function), we might expect the variogram have a shape indicating autocorrelation.. However, the entire model only explains <40% of variation (see R^2 values), so B.Lat & B.Long explain a relatively small amount of variation in body size. This means that even raw size values are not THAT autocorrelated. 

# sp::bubble -------------------------------------------------------------
par(mfrow = c(2, 1)) # Not working, not sure why
pdf(
  file = "Plots/Spatial_autocorrelation/Bubble_plots.pdf",
  width = 8.5, height = 11, bg = "white"
)
pmap(list(njdf.l.br.age, top.mods, names(njdf.l.br.age)), \(df, mods, names){
  coordinates(df) <- c("B.Long", "B.Lat")
  if (names %in% c("Whip-poor-will_Mass.combBT", "Whip-poor-will_Wing.comb")) {
    df@coords <- jitter(df@coords, amount = 2)
  } else {
    df@coords <- jitter(df@coords, amount = 5)
  }
  df$resid <- rstandard(mods)
  bubble(df, zcol = "resid", fill = FALSE, main = names)
})
dev.off()

# Moran's I----------------------------------------------------------------
## Moran's I looks at spatial autocorrelation of model residuals -- when the morans I > 0.1 & the p-value is < 0.05 we may have issues of spatial autocorrelation that violate lm() model assumptions of independence
# There are two primary things we must define when calculating Moran's I: 1) who are 'neighbors' (the points that we look for autocorrelation in), & 2) how do we weight the points that we define as neighbors.
# NOTE:: In function chooseCN() the 'Type' argument defines both how neighbors are selected & their weights. If we want more direct control of how neighbors are weighted, we can set the chooseCN(result.type = "nb"), & then control weights using nb2listw(style = ). By default the weights style = "W" , where W is row standardised (sums over all links to n)

# First step for all types is to extract the coordinates
coords <- map(njdf_sf, \(df){
  df %>% st_coordinates()
})


# >Type = 5: Distance based nearest neighbor ------------------------------
# threshold = smallest distance (in m) to keep everything connected
threshold <- map_dbl(coords, \(xygd){
  give.thresh(dist(xygd)) # Huge thresholds, hundreds or thousands of km
})

# Calculate connectivity network
dnn_l <- map2(coords, threshold, \(xygd, thresh){
  chooseCN(xygd, type = 5, result.type = "listw", d1 = 0.01, d2 = thresh, edit.nb = F)
})

# Calculate morans I & relevant metrics
moran5 <- pmap(
  list(top.mods, dnn_l, threshold),
  \(tm, dnn, thresh) {
    mt <- lm.morantest(tm, dnn)
    data.frame(
      type = 5, def.parm = "d2", parm.val = thresh,
      avg.nn = round(mean(card(dnn$neighbours)), 1),
      obs.mor = mt$estimate[1], p = mt$p.value
    ) %>%
      mutate(across(where(is.numeric), ~ round(., 2))) # Correct rounding
  }
)
# Bind rows & remove rownames
moran_df_dnn <- bind_rows(moran5, .id = "SppDV") %>%
  rownames_to_column() %>%
  select(-rowname)
moran_df_dnn %>% filter(obs.mor > .1)

## See extended script 
# For code examining type 6 (k nearest neighbours), & type = 7 (inverse spatial)
# When type = 6, we conduct a K nearest neighbours analysis, where spatial autocorrelation is examined amongst the k closest points.
# When type = 7, the spatial weights are directly proportional to the inverse spatial distances. the argument "a" is the the exponent of the inverse distance matrix, which controls how quickly the influence of neighboring points diminishes with increasing distance.

# >Correlogram -------------------------------------------------------------
# What are lag distances (order argument)? Lag 1 refers to the direct neighbors, lag 2 includes neighbors of neighbors
# Ultimately, for the supporting information in the Berg's rule manuscript I ended up just showing that there was some spatial autocorrelation at low lags 

# >>Raw size ---------------------------------------------------------------
# Calculate nearest neighbor network
threshold / 1000 # Number of kilometers
dnn_nb <- map2(coords, threshold, \(xygd, thresh){
  chooseCN(xygd, type = 5, result.type = "nb", d1 = 10, d2 = thresh, edit.nb = F)
})

# Calculate correlograms, w/ default style argument (compared all styles & "W" is good)
correlograms_raw <- pmap(list(dnn_nb[c(1,2,5,6)], njdf_sf[c(1,2,5,6)], DVs[c(1,2,5,6)]), \(nb, df, dv){
  sp.correlogram(nb, var = df[[dv]], order = 5, method = "I", zero.policy = TRUE)
})

# Create dataframes of estimates & variances
dfs_raw <- map_depth(correlograms_raw, .depth = 1, \(cor){
  cor$res %>%
    data.frame() %>%
    rename_all(~ c("estimate", "expectation", "variance")) %>%
    rownames_to_column(var = "Lag")
})

# Examine max Moran's I estimate
max_m <- map_dbl(dfs_raw, \(df){
  max(df$estimate)
})
max_m
range(max_m)

# Plot correlogramms in base R
Spp_dv <- expand.grid(Long = c("Nighthawk", "European nightjar", "Whip-poor-will"), 
                      DV = c("Mass", "Wing")) %>% 
  mutate(Full = paste0(Long, DV)) %>% 
  arrange(Long)

pdf("Plots/Spatial_autocorrelation/Correlograms_raw.pdf", width = 6.14, height = 8, bg = "white")
par(mfrow = c(3, 2))
map2(correlograms_raw, Spp_dv$Full[c(1,2,5,6)], \(corr, names){
  plot(corr, main = names, ylim = c(-.5, .5))
})
dev.off()

# >>Residuals --------------------------------------------------------------
# lags <- c(5, 5, 9, 9) # Note sure it makes sense to change lag values for different species..
correlograms_resid <- pmap(
  list(dnn_nb[c(1,2,5,6)], top.mods[c(1,2,5,6)]), # [c(1,2,5,6)]
  \(nb, tm, lag){
    sp.correlogram(nb,
      var = resid(tm),
      order = 5, method = "I", zero.policy = TRUE
    )
  }
)

# Create dataframes of estimates & variances
dfs_l <- map_depth(correlograms_resid, .depth = 1, \(cor){
  cor$res %>%
    data.frame() %>%
    rename_all(~ c("estimate", "expectation", "variance")) %>%
    rownames_to_column(var = "Lag")
})

# Examine max Moran's I estimate
map_dbl(dfs_l, \(df){
  max(df$estimate)
})

# Plot correlograms
pdf("Plots/Spatial_autocorrelation/Correlograms_resid.pdf", width = 6.14, height = 8, bg = "white")
par(mfrow = c(3, 2)) # c(2, 2)
map2(correlograms_resid, Spp_dv$Full[c(1,2,5,6)], \(corr, names){
  plot(corr, main = names, ylim = c(-.5, .5))
})
dev.off()

# Examine cardinalities: the numbers correspond to the frequencies of observations with specific cardinalities for each lag order
map(correlograms_resid, \(corr){
  corr$cardnos
})

# For example, for Nighthawk mass & at the first lag order, 7 units have 7 neighbors, 1 unit has 8 neighbors, 7 units have 23 neighbors, & so on
correlograms_resid$Nighthawk_Mass.combBT$cardnos[[1]]

# Note at larger lag orders, not every spatial unit will have neighbors (especially if the spatial data is sparse or irregular). As a result, the cardinality distribution changes, and some cardinalities may no longer be observed, leading to shorter lists.
correlograms_resid$Nighthawk_Mass.combBT$cardnos

# Can examine this for each spp & lag order if desired
barplot(correlograms_resid$Nighthawk_Mass.combBT$cardnos[[1]],
  main = "Lag Order 1",
  xlab = "Cardinality (Number of Neighbors)", ylab = "Frequency"
)

# gls() mods -----------------------------------------------------------------
# If spatial autocorrelation is observed, GLS models would be an appealing option given that our sites are not readily grouped into discrete sites, which we would want if we were fitting a random intercept.
# NOTE:: But remember, the GLS models still fit the correlation structure at the same time as the fixed effects (not after the fixed effects), and so this can still soak up the variation of the spatial fixed effects & result in altered parameter estimates.

# Create projected coordinates columns
njdf.proj <- map2(njdf_sf, coords, \(df, proj_coords){
  df %>% mutate(
    B.Long.proj = proj_coords[, 1],
    B.Lat.proj = proj_coords[, 2]
  )
})

# >Explore parms ---------------------------------------------------------
## The default values for the gls() function don't result in meaningful ranges, likely because there is little SAC in the data (see extended script). However, we already estimated where the nugget value should initialize at from the empirical variogram, so let's use these initial values (variog_parms$nugget_start). There are two other parameters that we have control over in the gls function -- the correlation structure & the range - this is an attempt to understand how these values influence our results (parameter estimates & standard errors)

cor_str <- setNames(c("corGaus", "corRatio", "corSpher"), c("corGaus", "corRatio", "corSpher")) # "corLin", "corExp"

# Calculate the maximum distance between observations for each species
max.dists <- map_dbl(njdf_sf, \(df){
  max.dist <- max(st_distance(df))
  units(max.dist) <- "km"
  max.dist
})

# Set range in kilometers
range_val_proj <- map(max.dists, \(md){
  distances <- round(seq(from = 1000, to = md, length.out = 5), 0)
  distances <- setNames(distances * 1000, distances)
  distances
})

# Function to fit GLS models with various correlation structures and ranges
imp_gls_spp <- function(df, cor_str, range_vals, nugg, mod) {
  # Loop through each correlation structure type, fitting models with specified ranges
  gls_models <- map(cor_str, function(cor) {
    # Fit models with specified correlation structure and ranges
    models_with_correlation <- map(range_vals, function(range) {
      # Define the correlation structure for the current range
      cor_structure <- get(cor, envir = asNamespace("nlme"))(value = c(range, nugg), nugget = TRUE,
        form = ~ B.Long.proj + B.Lat.proj)

      # Fit the GLS model with given cor_structure
      tryCatch(
        {
          nlme::gls(
            model = as.formula(mod),
            method = "REML",
            correlation = cor_structure,
            data = df
          )
        },
        error = function(e) {
          message(
            "Error fitting model for correlation: ", cor,
            " with range: ", range,
            ". Error: ", e$message
          )
          NULL # Return NULL if an error occurs
        }
      )
    })

    # Add the no-correlation model for each correlation structure type
    no_corr_model <- nlme::gls(
      model = as.formula(mod),
      method = "REML",
      correlation = NULL,
      data = df
    )

    # Combine models with correlation and no-correlation model for the current structure
    c(models_with_correlation, list(No_corr = no_corr_model))
  })
}

# Map over datasets and models in parallel
# NOTE:: corExp is supposed to be the fastest
gls_proj_sp <- pmap(list(njdf.proj[c(1,2, 5,6)], top.mods.chr[c(1,2, 5,6)], 
                         variog_parms[c(1,2, 5,6)], range_val_proj[c(1,2, 5,6)]),
  \(df, mod, nugg, range_proj){
    imp_gls_spp(
      df = df, cor_str = cor_str, range_vals = range_proj,
      nugg = nugg$nugget_start, mod = mod
    ) # cor_str = cor_str
  }
)

# >AIC table -------------------------------------------------------------
## Create AIC table comparing correlation structures (including no correlation structure) & different ranges
# Flatten so all correlation structures are sublists under SppDV
gls_sp_flat <- map(gls_proj_sp, function(gls_mods) {
  list_flatten(gls_mods) # name_spec = "{outer}_{inner}",
})
# Create AIC table for each SppDV
# NOTE:: Warnings are OK given that we are maintaining fixed effects the same & only varying correlation structures
aic_cor_spp <- map(gls_sp_flat, ~ .x %>%
  aictab() %>%
  as.data.frame() %>%
  mutate(
    Corr_str = str_split_i(Modnames, "_", 1),
    Range = str_split_i(Modnames, "_", 2)
  ))

# Of supported models the spatial structure is generally >100 km
ranges_2 <- map(aic_cor_spp, ~ .x %>%
  filter(Delta_AICc < 2) %>%
  pull(Modnames))
table(unlist(ranges_2))

# The No correlation model performs much worse than the top model w/ a correlation structure
map(aic_cor_spp, \(aictab) {
  aictab %>%
    filter(str_detect(Modnames, "No_corr")) %>%
    pull(Delta_AICc) %>%
    unique()
})

# NOTE:: Despite the no correlation model performing significantly worse than the models with correlation structures, the Diniz-Filho et al (2003) article argues convincingly that that doesn't necessarily mean it is the correct model to make inference from

# >df of parm estimates ---------------------------------------------------
# Process the nested list, flattening to a dataframe of parameter estimates
parms.df <- map_df(gls_proj_sp, \(Spp_list) {
  map_df(Spp_list, \(model_list) {
    if (inherits(model_list, "gls")) {
      tidy(model_list)
    } else {
      map_df(model_list, tidy, .id = "range")
    }
  }, .id = "Corr_str")
}, .id = "SppDV")

## Work with parameter estimates to summarize meaningfully
# Create a filtered parm df with relevant parameters & ranges
parms.filt <- parms.df %>%
  filter(term %in% c("B.Elev", "B.Lat", "B.Long", "B.Tavg")) %>%
  # select(-Corr_str) %>% #Corr_str barely impacts parm estimates
  filter(range >= 100 | range == "No_corr") %>% # ranges
  mutate(
    Corr_str = ifelse(range == "No_corr", NA, Corr_str),
    corr.mod = ifelse(range == "No_corr", "No", "Yes")
  ) %>%
  distinct()

# Plot:: Visually inspect the differences between models with & without correlation structure 
# NOTE:: Regardless of model type or range, the estimates don't change very much!
parms.filt %>%
  filter(Corr_str %in% c("corGaus") & range == 1000 | is.na(Corr_str)) %>% 
  ggplot(aes(
    x = term, y = estimate,
    group = interaction(corr.mod, SppDV), color = SppDV
  )) + # Corr_str
  geom_point(aes(shape = corr.mod), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
    position = position_dodge(width = 0.5)
  ) +
  theme_minimal()

# Inspect results & notice that the correlation type doesn't influence results much
parms.filt %>%
  summarize(across(c(estimate, std.error, statistic, p.value), \(x) mean(x, na.rm = TRUE)),
    .by = c(SppDV, term, corr.mod, Corr_str)
  ) %>%
  arrange(SppDV, term, Corr_str)

# Summarize parameter estimates, averaging across the different correlation structures & examining the difference in magnitude in the p-value
parms.red <- parms.filt %>%
  summarize(across(c(estimate, std.error, statistic, p.value), mean, na.rm = TRUE),
    .by = c(SppDV, term, corr.mod)
  ) %>%
  mutate(
    sig = p.value < .05,
    p.value = formatC(p.value, format = "e", digits = 2), # Convert to scientific notation
    magnitude = as.numeric(str_split_i(p.value, "-", 2)),
    across(where(is.numeric), ~ round(., 3))
  ) %>%
  arrange(SppDV, term, desc(corr.mod))
parms.red %>% summarize(Diff.mag = diff(magnitude), .by = c(SppDV, term))

# Export Rdata ------------------------------------------------------------
#rm(list = ls()[!(ls() %in% c("aictab_list3", "njdf.l.br.age", "top.geo.mods.chr", "variog_tm_bin", "variog_tm_cloud", "variog_tm_smooth", "dnn_l", "knn_l", "inv_dist_l", "moran5", "moran6", "moran7", "moran_df", "gls_proj_sp", "range.def", "parms.def", "parms.df", "correlograms_raw", "correlograms_resid", "Spp_dv"))])
#save.image(paste0("Rdata/Spatial_autocorrelation_", format(Sys.Date(), "%m.%d.%y"), ".Rdata"))