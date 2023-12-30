##Data wrangling 02 -- Merge njdf & environmental covs 
##Data wrangling script to merge nightjar df with environmental covariates df and ultimately create the environment needed to begin the analysis. This will happen in at least 6 steps
#1. Remove additional rows in capriBA to create capri.fin (final)
#2. Link capri.fin with EnviCovs
#3. Remove superflous columns
#4. Generates correlation plots from predictor variables
#5. Scale numeric variables 
#6&7. Generated subsetted dfs based on species & DV, and breeding / winter 
#8. Export files & save the .Rdata file

library(tidyverse)
library(readxl)
library(xlsx)
library(chron)

#1. Filter capriBA ---------------------------------------------------------------
#Remove additional rows in capriBA to create capri.fin (final)
#Remove repeat individuals, selecting adults when possible 
levels(capriBA$Age) #Should get NULL, this should NOT be a factor in order for arrange to work correctly
capriBAnr <- capriBA %>% group_by(Band.Number) %>% #nr = no repeats 
  arrange(is.na(W.Lat), Age) %>% 
  slice_head() %>% 
  ungroup()
nrow(capriBAnr)

#Filter based on date
euniRes <- capriBAnr %>% filter(Species == "EUNI" & Band.md > as.POSIXct("2023-05-16") & Band.md < as.POSIXct("2023-08-01") & is.na(W.Lat))
euniFAC <- capriBAnr %>% filter(Species == "EUNI" & !is.na(W.Lat))
coni.ewpw <- capriBAnr %>% filter(Species != "EUNI" & Band.md > as.POSIXct("2023-04-30")) 
capriRes <- rbind(euniRes, euniFAC, coni.ewpw)
nrow(capriRes) #capri residents

#Remove last few things
#Remove Greg's EDB project due to lack of precision in decimals, Unknown sex, and year < 2010. 
#Sex is an important covariate in all models, and latitudes are sampled about evenly after 2010
capri.fin <- capriRes %>% filter(Project != "Greg_EDB" & Sex != "U" & Year >= 2010)

#2. Merge -------------------------------------------------------------------
#Link capri.fin with EnviCovs
#Bring in EnviCovs2
EnviCovs2 <- read_xlsx()

EnviCovs3 <- EnviCovs2 %>% select(-c(Species, Banding.Date, Band.Number))
capriA <- capri.fin %>% left_join(EnviCovs3, by = "uniqID") #"uniqID" #capri Analysis
   
names(EnviCovs2)
#3. Remove superflous columns -----------------------------------------------
HypVars <- vector("list", length = 5)
HypVars[[1]] <- data.frame(Vars = c("Lat", "Long", "Elev")) #Geography
HypVars[[2]] <- data.frame(Vars = c("Srad", "Tavg")) #Temp Regulation (TR)
HypVars[[3]] <- data.frame(Vars = c("EVI", "Prec")) #Productivity
HypVars[[4]] <- data.frame(Vars = c("EviCV", "CVprec", "Tcv")) #Seasonality
HypVars[[5]] <- data.frame(Vars = c("Mig.dist"))
names(HypVars) <- c("Geo", "TR", "Prod", "Seas", "Mig.Dist")
HypVarsDf <- bind_rows(HypVars, .id = "Hypothesis")
HypVarsDf <- HypVarsDf %>% mutate(Full = c("Latitude", "Longitude", "Elevation", "Solar radiation", "Temperature", "EVI", "Precipitation", "EVI CV", "Precipitation CV", "Temperature CV", "Migratory distance"))

capriA.red <- capriA %>% 
  select(contains(c("uniqID","Band.Number","Project","Species","Age","Sex","comb", HypVarsDf$Vars))) %>%
  select(-c(WingFlat, Band.Age, Mass.comb))

#Remove a few coastal individuals (from breeding grounds) that have no environmental data 
capriA.red2 <- capriA.red %>% drop_na(starts_with("B")) 
nrow(capriA.red2) 

#4. Correlation plots ---------------------------------------------------
namesEC <- names(capriA.red2)
namesECb <- c(namesEC[str_detect(namesEC, "^B")], "Mig.dist")
namesECw <- c(namesEC[str_detect(namesEC, "^W")], "Mig.dist")
#Create df to loop through function
loop_spp_sea <- expand.grid(Species = c("Nighthawk", "Whip-poor-will", "Nightjar"), 
                            Season = c("Breeding", "Winter"),
                            stringsAsFactors = FALSE) %>% 
  arrange(Species, Season)

#Function selecting the species and season
select_spp_sea <- function(df, Spp, Season){
  df %>% filter(Species == Spp) %>%
    select(if (Season == "Breeding") {{ namesECb }} else {{ namesECw }})
}

panel <- paste0(LETTERS[1:6], ")")
Corr.plots <- list()
for(i in 1:nrow(loop_spp_sea)){
  df_spp_sea <- select_spp_sea(capriA.red2, loop_spp_sea[i,1], loop_spp_sea[i,2])
  Corr.plots[[i]] <- GGally::ggcorr(df_spp_sea, label = T, label_size = 3, hjust = 0.75, size = 3) + 
    ggtitle(paste(loop_spp_sea[i,1], loop_spp_sea[i,2], "Predictors"), subtitle = panel[i]) 
}

#Print plots to PDF
pdf("Plots/Corr_IVs2.pdf")
print(marrangeGrob(Corr.plots, ncol = 1, nrow = 1))
dev.off()

# 5.  Scale numeric variables  --------------------------------------------
capriA.s <- capriA.red2 %>% mutate(across(where(is.numeric), ~ scale(.)[,1])) 
#Now that all data wrangling is finished..
capri.fac <- capriA.s %>% filter(!is.na(W.Lat))
nrow(capri.fac)

#Of FAC data: 9 wing NAs, 2 massBT NAs, 3 Mig dist NAs
capri.fac %>% filter_all(any_vars(is.na(.)))

# 6.  Subset breeding dfs list --------------------------------------------------
#Subsetted dfs based on species & DV, these dfs will be used in future analysis script

subset.df <- function(df, spp, var){
  df %>% filter(Species == spp & !is.na(.data[[var]]))
}

loopSppDV <- expand.grid(Species = c("Nighthawk", "Whip-poor-will", "Nightjar"), 
                         DV = c("Wing.comb", "Mass.combBT"),
                         stringsAsFactors = FALSE) %>% 
  arrange(Species, DV)

njdf.list.br <- list()
for(i in 1:nrow(loopSppDV)){
  njdf.list.br[[i]] <- subset.df(capriA.s, spp = loopSppDV[i, "Species"], var = loopSppDV[i, "DV"])
}
names(njdf.list.br) <- paste0(loopSppDV[,"Species"], ".", loopSppDV[,"DV"])
lapply(njdf.list.br, nrow)

# 7. Subset winter dfs list --------------------------------------------------

njdf.list.fac <- list()
for(i in 1:nrow(loopSppDV)){
  njdf.list.fac[[i]] <- subset.df(capri.fac, spp = loopSppDV[i, "Species"], var = loopSppDV[i, "DV"])
}
names(njdf.list.fac) <- paste0(loopSppDV[,"Species"], ".", loopSppDV[,"DV"])
lapply(njdf.list.fac, nrow)

#DELETE once Elly gets back 
njdf.list.fac <- lapply(njdf.list.fac, function(x) {x %>% filter(!is.na(Mig.dist))})

# #8.  Export -------------------------------------------------------------

#save.image("Rdata/Capri_dfs.Rdata")
capriA.red %>% as.data.frame() %>%
  write.xlsx("Intermediate_products/capriA.red_12.21.23.xlsx", row.names = F)
