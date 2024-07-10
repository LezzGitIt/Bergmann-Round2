##Bergmann's rule across full-annual-cycle in Caprimulgids##
##Data wrangling 01 -- Summarize environmental covariates 
#This script loads environmental datasets that Elly Knight downloaded from GEE, summarizes over the relevant months, and condenses so there is a single row per uniqID

#Contents
#1) Format elevation
#2) Average the environmental variables over the relevant months 
#3) Merge different data sets & seasons into a single data frame of EnviCovs

# Libraries & load data ---------------------------------------------------
library(tidyverse)
library(readxl)
library(xlsx)
library(chron)

#load("Data/EnviCovs_script.Rdata")

capriBA <- read_xlsx("Intermediate_products/Capri_BA_07.03.24.xlsx", trim_ws = TRUE) %>%
  mutate(Banding.Time = chron(times = Banding.Time))

##Bring in data sets (or load .Rdata file)
wc <- read.csv("Data/EK_Envi_Vars/Covs_12.24.23/Wordclim-Buffer.csv")
evi <- read.csv("Data/EK_Envi_Vars/Covs_12.24.23/EVI-Buffer.csv")
tc <- read.csv("Data/EK_Envi_Vars/Covs_12.24.23/Terraclim-buffer.csv")
elev <- read.csv("Data/EK_Envi_Vars/Covs_12.24.23/DEM-Buffer.csv")

# Elevation format  -------------------------------------------------------
#NOTE: Elevation df is structured differently
elev2 <- elev %>% mutate(Season_abb = ifelse(season == "Breed", "B.", "W.")) %>% 
  dplyr::select(uniqID, Season_abb, mean) %>% 
  group_by(uniqID) %>% 
  pivot_wider(names_from = Season_abb, names_glue = "{Season_abb}Elev", values_from = mean) #Notice names glue so season starts the name


# Average environmental covs ----------------------------------------------
##Using for loop for each Spp, average the environmental variables over the relevant months (see Table 1 in manuscript) for each species
Spp <- c("CONI", "EWPW", "EUNI")
months <- list(months_br_end = c(8,9,8), months_wi_end = c(3,3,2))
br.wc <- wi.wc <- br.evi <- wi.evi <-  br.tc <- wi.tc <- list()
for(i in 1:length(Spp)){
  print(i)
  wc_df <- wc %>% filter(Species == Spp[i])
  evi_df <- evi %>% filter(Species == Spp[i])
  tc_df <- tc %>% filter(Species == Spp[i])
  
  #World Clim - Previously had WiTmax = mean(tmax), WiTMAX = max(tmax), WiTmin = mean(tmin), WiTMIN = min(tmin), WiMDR = sum(tmax - tmin) / 5
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
  
  #Terraclim - Previously had Aet & Vap as well
  br.tc[[i]] <- tc_df %>% group_by(uniqID) %>% 
    filter(season == "Breed", covmonth >= 5 & covmonth <= months[[1]][i]) %>% 
    summarize(B.Srad = mean(srad, na.rm = T))
  wi.tc[[i]] <- tc_df %>% group_by(uniqID) %>% 
    filter(season == "Winter", covmonth >= 11 | covmonth <= months[[2]][i]) %>% 
    summarize(W.Srad = mean(srad, na.rm = T)) 
}


# Create single data frame ------------------------------------------------
#Merge different data sets & seasons into a single data frame of EnviCovs
dfs <- list(br.wc, wi.wc, br.evi, wi.evi, br.tc, wi.tc, elev2)
names(dfs) <- c("br.wc", "wi.wc", "br.evi", "wi.evi", "br.tc", "wi.tc", "elev")
#dfs <- lapply(dfs, function(x){setNames(x, Spp)}) Not necessary
dfs <- lapply(dfs, bind_rows)
EnviCovs <- dfs %>% purrr::reduce(full_join, by = "uniqID") %>% 
  mutate(across(where(is.numeric), round, 3))

#Join with Envi Covs df w/ a few key identifier variables from njdf & export
EnviCovs2 <- EnviCovs %>% inner_join(capriBA[,c("uniqID", "Band.Number", "Banding.Date", "Species")], by = "uniqID") 
nrow(EnviCovs2) 

# Export ------------------------------------------------------------------
#EnviCovs2 %>% as.data.frame() %>%
 # write.xlsx(paste0("Intermediate_products/Envi_Covs_", format(Sys.Date(), "%m.%d.%y"), ".xlsx"), row.names = F)


#Notes on MDR
#MDR = mean diurnal range, BIO2. Many of these are similar to the 19 worldclim vars (good comparison to make in manuscript, as people are familiar with these vars). Notice particularly that BrTcv is not equivalent to using the cv function (i.e compCV = cv(tavg, na.rm = T) produces slightly different results). 
#Should MDR go into seasonality or temp regulation... Elly & I think TR, b/c seasonality is ACROSS the entirety of the season, MDR is within each month.
