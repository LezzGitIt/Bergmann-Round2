##Data wrangling 01 -- Summarize environmental covariates 
#This script loads environmental datasets that Elly Knight downloaded from GEE, summarizes over the relevant months, and condenses so there is a single row per uniqID
library(tidyverse)
library(readxl)
library(xlsx)
library(chron)

capriBA <- read_xlsx("Intermediate_products/Capri_BA_12.19.23.xlsx", trim_ws = TRUE) %>%
  mutate(Banding.Time = chron(times = Banding.Time))

##Bring in data sets
#rm(list= ls()[!(ls() %in% c('wc', "evi", "tc", "elev"))])
wc <- read.csv("Data/EK_Envi_Vars/Covs_12.24.23/Wordclim-Buffer.csv")
evi <- read.csv("Data/EK_Envi_Vars/Covs_12.24.23/EVI-Buffer.csv")
tc <- read.csv("Data/EK_Envi_Vars/Covs_12.24.23/Terraclim-buffer.csv")
elev <- read.csv("Data/EK_Envi_Vars/Covs_12.24.23/DEM-Buffer.csv")

#Note elevation is structured differently
elev2 <- elev %>% mutate(Season_abb = ifelse(season == "Breed", "B.", "W.")) %>% 
  dplyr::select(uniqID, Season_abb, mean) %>% 
  group_by(uniqID) %>% 
  pivot_wider(names_from = Season_abb, names_glue = "{Season_abb}Elev", values_from = mean) #Notice names glue so season starts the name

##Average the environmental variables over the relevant months (Table 1) for each species
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

#Create a single data frame of EnviCovs
dfs <- list(br.wc, wi.wc, br.evi, wi.evi, br.tc, wi.tc, elev2)
names(dfs) <- c("br.wc", "wi.wc", "br.evi", "wi.evi", "br.tc", "wi.tc", "elev")
#dfs <- lapply(dfs, function(x){setNames(x, Spp)}) Not necessary
dfs <- lapply(dfs, bind_rows)
EnviCovs <- dfs %>% purrr::reduce(full_join, by = "uniqID") %>% 
  mutate(across(where(is.numeric), round, 3))


###B/c rowID got messed up (likely due to changing Marja's banding times), we have to overwrite some rowIDs so we get a good match 
CheckDFelly <- read.csv("/Users/aaronskinner/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Grad_School/MS/EWPW/Writing_Exit_Seminar/Bergs_Rule/CapDfElly.csv")
#capriBA <- capriBA %>% mutate(rowID = row_number())
#We can see that things got messed up
identical(capriBA[,c("rowID", "Band.Number")] , CheckDFelly[,c("rowID", "Band.Number")])

capriBA$Banding.Date <- as.character(capriBA$Banding.Date)
#Let's check for repeats in these key columns that we'll use to match. 
nrow(capriBA) #Should be 5927 (run through removal of age)
TF <- capriBA[,c("Band.Number", "Banding.Date", "W.Lat")] == CheckDFelly[,c( "Band.Number", "Banding.Date", "W.Lat")]
#NOT RELEVANT AT PRESENT::If you were to delete rows where rowID is same but band numbers differ, then things could be problematic when Band.Numbers are the same and something else differs. Biggest key is that this is not issue w/ W.Lat, but the Banding.Date is also not ideal (but it's OK b/c just 5 individuals). Probably best to just delete these 189 individuals
data.frame(TF) %>% filter(Band.Number != Banding.Date | Band.Number != W.Lat | W.Lat != Banding.Date)
TF <- capriBA[,c("Band.Number")] == CheckDFelly[,c( "Band.Number")]
table(TF)

#Merge on this and you'll see this won't be super straightforward b/c there are 199 instead of 189
test <- merge(capriBA[!TF,], CheckDFelly[,c("Band.Number",  "Banding.Date", "W.Lat", "rowID")], by = c("Band.Number",  "Banding.Date", "W.Lat"), all.x = T) %>% arrange(Band.Number)
test <- test %>% filter_all(any_vars(!is.na(.)))
nrow(test) #199 individuals, meaning there are problems with duplicates
test %>% filter(!duplicated(Band.Number) & !is.na(W.Lat)) 

##Loop to fix the problem
rowIDs <- vector()
for(i in 1:nrow(capriBA)){
  if(!is.na(CheckDFelly$Band.Number[i]) & !is.na(capriBA$Band.Number[i]) & capriBA$Band.Number[i] != CheckDFelly$Band.Number[i]){
    print(i)
    capriBA$rowID[i] <- CheckDFelly %>% 
      filter(Band.Number == capriBA$Band.Number[i] & Banding.Date == capriBA$Banding.Date[i] & is.na(W.Lat) == is.na(capriBA$W.Lat[i])) %>% 
      pull(rowID)
  }
}

CheckDFelly %>% 
  filter(Band.Number == capriBA$Band.Number[i] & is.na(W.Lat) == is.na(capriBA$W.Lat[i]))

##CHECK:: Can confirm loop worked by checking that these are the same individuals
i <- 5927
CheckDFelly[1154,]
capriBA[i,]

##FINAL CHECK, IMPORTANT:: Band numbers should be identical at this point upon merging
test <- merge(capriBA[,c("Band.Number", "rowID")], CheckDFelly[, c("Band.Number", "rowID")], by = "rowID") 
table(test$Band.Number.x == test$Band.Number.y)

#REMEMBER::Finally, b/c there were some individuals caught on same date and w/out wintering information, they were assigned the same rowID.
length(unique(capriBA$rowID)) #There are some individuals that have repeated rowID, but this does not matter for assignment of environmental covariates
capriBA %>% filter(duplicated(rowID))
##End of rowID fix 


#Link Environmental covs w/ capriBA
EnviCovs2 <- EnviCovs %>% inner_join(capriBA[,c("uniqID", "Band.Number", "Banding.Date", "Species")], by = "uniqID") #"B.Lat", "B.Long",  "W.Lat", "W.Long", "Mig.dist"
nrow(EnviCovs2) 
ggplot(data = EnviCovs2, aes(x = B.Lat, y = BrEVI, color = Species)) + 
  geom_jitter(alpha = .2) + 
  geom_smooth(method = "lm")
EnviCovs2 %>% group_by(Species) %>% 
  summarize(cor = cor(BrPrec, B.Long, use = "complete.obs"))
#ggsave("Plots/EVI_Blat_Rd2.png", bg = "white")

#Notes on MDR
#MDR = mean diurnal range, BIO2. Many of these are similar to the 19 worldclim vars (good comparison to make in manuscript, as people are familiar with these vars). Notice particularly that BrTcv is not equivalent to using the cv function (i.e compCV = cv(tavg, na.rm = T) produces slightly different results). 
#Should MDR go into seasonality or temp regulation... Elly & I think TR, b/c seasonality is ACROSS the entirety of the season, MDR is within each month.

#save.image("/Users/aaronskinner/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Grad_School/MS/EWPW/Writing:Exit_Seminar/Bergs_Rule/Data_share/Data/EK_Envi_Vars/Correlations/EC_workspace.Rdata")
