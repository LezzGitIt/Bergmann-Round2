## Bergmann's rule across full-annual-cycle in Caprimulgids ##
## Data wrangling 01 -- Summarize environmental covariates 
# This script loads environmental datasets downloaded from Google Earth Engine, summarizes over the relevant months, and condenses so there is a single row per uniqID

#Contents
#1) Format elevation
#2) Average the environmental variables over the relevant months 
#3) Merge different data sets & seasons into a single data frame of EnviCovs

# Libraries & load data ---------------------------------------------------
library(tidyverse)
library(readxl)
library(chron)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(purrr::map)

capriBA <- read.csv(paste0("Intermediates/Capri_BA_", 
                            format(Sys.Date(), "%m.%d.%y"), ".csv")) %>%
  mutate(Banding.Time = chron(times = Banding.Time))

# NOTE:: THESE DATA SETS ARE NOT PROVIDED DUE TO GOOGLE'S STRICT DATA POLICIES. THE OUTPUT OF THIS SCRIPT IS PROVIDED ('DATA/Envi_Covs_07.29.25.xlsx') SO YOU CAN CONTINUE ON TO THE NEXT SCRIPT (02). 

if(file.exists("Data/Envi_Covs_07.29.25.xlsx")){
  print("The output file created from this script already exists. Please proceed to the next script.")
  stop()
}

## Bring in data sets
wc <- read.csv("Data/Envi_Covs_Raw/Wordclim-Buffer.csv") 
evi <- read.csv("Data/Envi_Covs_Raw/EVI-Buffer.csv")
tc <- read.csv("Data/Envi_Covs_Raw/Terraclim-buffer.csv")
elev <- read.csv("Data/Envi_Covs_Raw/DEM-Buffer.csv")

# Elevation format  -------------------------------------------------------
#NOTE: Elevation df is structured differently
elev2 <- elev %>% mutate(Season_abb = ifelse(season == "Breed", "B.", "W.")) %>% 
  dplyr::select(uniqID, Season_abb, mean) %>% 
  group_by(uniqID) %>% 
  pivot_wider(names_from = Season_abb, names_glue = "{Season_abb}Elev", values_from = mean) #Notice names glue so season starts the name

# Average environmental covs ----------------------------------------------
## Using for loop for each Spp, average the environmental variables over the relevant months (see Table 1 in manuscript) for each species
Spp <- c("CONI", "EWPW", "EUNI")
months <- list(months_br_end = c(8,9,8), months_wi_end = c(3,3,2))
br.wc <- wi.wc <- br.evi <- wi.evi <-  br.tc <- wi.tc <- list()
for(i in 1:length(Spp)){
  print(i)
  wc_df <- wc %>% filter(Species == Spp[i])
  evi_df <- evi %>% filter(Species == Spp[i])
  tc_df <- tc %>% filter(Species == Spp[i])
  
  #World Clim 
  br.wc[[i]] <- wc_df %>% group_by(uniqID) %>% 
    filter(season == "Breed", covmonth >= 5 & covmonth <= months[[1]][i]) %>% 
    summarize(B.Prec = mean(prec), 
              B.CVprec = raster::cv(prec), 
              B.Tavg = mean(tavg), 
              B.Tcv = ((sd(tavg/10))/(mean(tavg/10) + 273.15))*100)
  wi.wc[[i]] <- wc_df %>% group_by(uniqID) %>% 
    filter(season == "Winter", covmonth >= 11 | covmonth <= months[[2]][i]) %>% 
    summarize(W.Prec = mean(prec, na.rm = T), 
              W.CVprec = raster::cv(prec, na.rm = T), 
              W.Tavg = mean(tavg, na.rm = T), 
              W.Tcv = ((sd(tavg/10, na.rm = T))/(mean(tavg/10, na.rm = T) + 273.15))*100)
  
  #EVI
  br.evi[[i]] <- evi_df %>% group_by(uniqID) %>% 
    filter(season == "Breed", covmonth >= 5 & covmonth <= months[[1]][i]) %>% 
    summarize(B.EVI = mean(EVI, na.rm = T), B.EviCV = raster::cv(EVI, na.rm = T))
  wi.evi[[i]] <- evi_df %>% group_by(uniqID) %>% 
    filter(season == "Winter", covmonth >= 11 | covmonth <= months[[2]][i]) %>% 
    summarize(W.EVI = mean(EVI, na.rm = T), W.EviCV = raster::cv(EVI, na.rm = T))
  
  #Terraclim 
  br.tc[[i]] <- tc_df %>% group_by(uniqID) %>% 
    filter(season == "Breed", covmonth >= 5 & covmonth <= months[[1]][i]) %>% 
    summarize(B.Srad = mean(srad, na.rm = T))
  wi.tc[[i]] <- tc_df %>% group_by(uniqID) %>% 
    filter(season == "Winter", covmonth >= 11 | covmonth <= months[[2]][i]) %>% 
    summarize(W.Srad = mean(srad, na.rm = T)) 
}

# Create single data frame ------------------------------------------------
# Merge different data sets & seasons into a single data frame of EnviCovs
dfs <- list(br.wc, wi.wc, br.evi, wi.evi, br.tc, wi.tc, elev2)
names(dfs) <- c("br.wc", "wi.wc", "br.evi", "wi.evi", "br.tc", "wi.tc", "elev")
dfs <- lapply(dfs, bind_rows)
EnviCovs <- dfs %>% purrr::reduce(full_join, by = "uniqID") %>% 
  mutate(across(where(is.numeric), round, 3),
         uniqID = str_replace(uniqID, " ", ""))

#Join with Envi Covs df w/ a few key identifier variables from njdf
EnviCovs2 <- EnviCovs %>% 
  inner_join(capriBA[,c("uniqID", "Band.Number", "Banding.Date", "Species")], 
             by = "uniqID") 

# Export ------------------------------------------------------------------
nrow(EnviCovs2)
EnviCovs2 %>% as.data.frame() %>%
  write.csv(paste0("Intermediates/Envi_Covs_", format(Sys.Date(), "%m.%d.%y"), ".csv"), row.names = FALSE)
