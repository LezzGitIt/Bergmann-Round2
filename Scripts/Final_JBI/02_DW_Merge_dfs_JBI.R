##Bergmann's rule across full-annual-cycle in Caprimulgids##
## Data wrangling 02 -- Merge njdf & environmental covs 
## Data wrangling script to merge nightjar df with environmental covariates df.

#Contents
#1. Remove additional rows in capriBA to create capri.fin (final)
#2. Link capri.fin with EnviCovs
#3. Remove superflous columns
#4. Export capri_analysis and capri_fac if necessary

library(tidyverse)
library(readxl)
library(chron)
library(gridExtra)
library(mvnormtest)
library(ggpubr)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(purrr::map)

# Read in CapriBA --------------------------------------------------------
capriBA <- read.csv(paste0("Intermediates/Capri_BA_", 
                             format(Sys.Date(), "%m.%d.%y"), ".csv")) %>%
  mutate(Banding.Time = chron(times = Banding.Time))

#1. Create capri.fin --------------------------------------------------------

# Remove repeat individuals, selecting adults when possible 
capriBAnr <- capriBA %>% group_by(Band.Number) %>% #nr = no repeats 
  arrange(is.na(W.Lat), Age) %>% 
  slice_head() %>% 
  ungroup()

### Data exploration indicated several issues that could bias size ~ covariate relationships. Filter the dataset to produce the final data set of caprimulgids (capri.fin)

## Important effects of date in European nightjars (e.g. migratory fattening, presence of migrants; see supporting information) 
# In creation of Band.md it assumes it is the present year, so this ensures our filter matches the Band.md column (i.e., Year = present year)
Year <- format(Sys.Date(), "%Y")
Year

# Filter based on date
euniRes <- capriBAnr %>% 
  filter(
    Species == "Nightjar" & Band.md > as.POSIXct(paste0(.env$Year,"-05-16")) & Band.md < as.POSIXct(paste0(.env$Year, "-08-01")) & is.na(W.Lat)
  )
# Do not remove individuals if they have winter data 
euniFAC <- capriBAnr %>% filter(Species == "Nightjar" & !is.na(W.Lat))
coni.ewpw <- capriBAnr %>% 
  filter(Species != "Nightjar" & Band.md > as.POSIXct(paste0(Year, "-04-30"))) 
# 'capriRes' object = Caprimulgid residents
capriRes <- rbind(euniRes, euniFAC, coni.ewpw)

# Sex is an important covariate in all models, and latitudes are sampled about evenly after 2010
capri.fin <- capriRes %>% filter(Sex != "U" & Year >= 2010)

#2. Merge -------------------------------------------------------------------
## Link capri.fin with EnviCovs

# Bring in EnviCovs2
EnviCovs2 <- read_excel("Data/Envi_Covs_07.29.25.xlsx")

EnviCovs3 <- EnviCovs2 %>% select(-c(Species, Banding.Date, Band.Number))
capriA <- capri.fin %>% left_join(EnviCovs3, by = "uniqID") #capri Analysis

#3. Remove superflous columns -----------------------------------------------
HypVars <- vector("list", length = 5)
HypVars[[1]] <- data.frame(Vars = c("Lat", "Long", "Elev")) #Geography
HypVars[[2]] <- data.frame(Vars = c("Srad", "Tavg")) #Temp Regulation (TR)
HypVars[[3]] <- data.frame(Vars = c("EVI", "Prec")) #Productivity
HypVars[[4]] <- data.frame(Vars = c("EviCV", "CVprec", "Tcv")) #Seasonality
HypVars[[5]] <- data.frame(Vars = c("Mig.dist"))
names(HypVars) <- c("Geo", "TR", "Prod", "Seas", "Mig.Dist") #"Post.Hoc"
HypVarsDf <- bind_rows(HypVars, .id = "Hypothesis")
HypVarsDf <- HypVarsDf %>% mutate(Full = c("Latitude", "Longitude", "Elevation", "Solar radiation", "Temperature", "EVI", "Precipitation", "EVI CV", "Precipitation CV", "Temperature CV", "Migratory distance")) #"Temperature", "Precipitation"

# Analysis file reduced to include only important variables
capriA.red <- capriA %>% select(c("uniqID","Band.Number", "Site.name", "Species","Age","Sex", "tsss.comb", "Wing.comb", "Mass.combBT", "Mass.comb", "Mig.dist", paste0("B.", HypVarsDf$Vars[1:10]), paste0("W.", HypVarsDf$Vars[1:10])))

## Final data frames 
# Remove a few coastal individuals (from breeding grounds) that have no environmental data 
capri_analysis <- capriA.red %>% drop_na(starts_with("B")) 

# Create full annual cycle dataset by filtering on wintering latitude
capri_fac <- capri_analysis %>% filter(!is.na(W.Lat))

# Export ------------------------------------------------------------------
# If desired can export these final 2 data frames, but this is not necessary as they are already provided in the 'Data/Analysis' folder 

if(!file.exists("Data/Analysis/capri_analysis07.28.25.xlsx") || !file.exists("Data/Analysis/capri_fac07.28.25.xlsx")){
  write.csv(capri_analysis, file = paste0("Data/Analysis/capri_analysis", format(Sys.Date(), "%m.%d.%y"), ".csv"), row.names = FALSE)
  write.csv(capri_fac, file = paste0("Data/Analysis/capri_fac", format(Sys.Date(), "%m.%d.%y"), ".csv"), row.names = FALSE)
} else{
  print("The output files created from this script already exist, thus proceed to the next script.")
}
