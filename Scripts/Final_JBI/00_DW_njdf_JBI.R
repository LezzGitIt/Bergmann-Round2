##Bergmann's rule across full-annual-cycle in Caprimulgids##

# General  -----------------------------------------------------------------
#SessionInfo:
#R version 4.3.2 (2023-10-31)
#Platform: x86_64-apple-darwin20 (64-bit)
#Running under: macOS Sonoma 14.4.1

# Instructions / suggestions for this repository
#All necessary data and code is provided to reproduce all analyses, figures, and tables in the main text. The file structure is provided, and if you open the .Rproj file everything will run as is. Otherwise, all that is necessary is that you set your working directory to the repository folder, like so: 

directory <- getwd()
setwd(directory) #Personalize, depending on where you've stored the repository folder

#Data wrangling 00 -- Create base night jar data frame (njdf)
##This script takes the combined data (unshown script joins all the raw data from providers into a single data frame) and formats, generates the time since sunset (tsss) variable, averages masses & wings from multiple captures by band age, & tidies and filters data to produce the capriBA data frame that will be used in downstream analyses

#Contents
# 1) Load & format data 
# 2) Condense age classes to Young, adult, or unknown
# 3) Calculate the difference between time since sunset & when an individual was captured
# 4) Calculate migratory distance 
# 5) Create data frames for export 

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(stringi)
library(naniar)
library(readxl)
library(chron)
library(lutz)
library(suncalc)
library(zoo) 
library(geosphere)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(purrr::map)

# Load and format data ----------------------------------------------------
capri.df <- read_excel("Data/Capri_df_combined_07.25.25.xlsx", sheet = "Data")

#Adjust time & date
capri.df <- capri.df %>% 
  replace_with_na_all(condition = ~.x %in% c(-99,-990, 9999, "<NA>", "-", ".", "na", 'NONABAND')) %>% 
  as.data.frame() %>% 
  filter(!is.na(B.Lat) & !is.na(Band.Number))
capri.df <- capri.df %>% 
  mutate(Species = ifelse(capri.df$Species == "Ceur" | capri.df$Species == "European Nightjar" | capri.df$Species == "European Nigthtjar", "EUNI", capri.df$Species), 
         Band.Number = stri_replace_all_regex(capri.df$Band.Number,
                                              pattern=c('-', ' '),
                                              replacement=c(''),
                                              vectorize=FALSE)) %>%
  mutate(uniqID = paste0(Band.Number, "_", Banding.Date, "_", !is.na(W.Lat))) %>% 
  mutate(Species = case_when(
    Species == "EUNI" ~ "Nightjar",
    Species == "EWPW" ~ "Whip-poor-will",
    Species == "CONI" ~ "Nighthawk",
  )) %>%
  arrange(is.na(W.Lat), Species) 

capri.df$Banding.Time <- str_pad(capri.df$Banding.Time, 4, pad = "0")
capri.df$Banding.Time <- sapply(str_split(parse_date_time(capri.df[,c("Banding.Time")], c("HMS"), truncated = 3), " "), function(x){x[2]})
capri.df$Banding.Time <- chron(times = capri.df$Banding.Time)
capri.df$Year <- str_pad(capri.df$Year, 3, pad = "0")
capri.df$Year <- str_pad(capri.df$Year, 4, pad = "2")

#Format times & dates, data types#
capri.df[,c("Wing.Chord","Mass", "B.Lat", "B.Long", "W.Lat", "W.Long", "Mig.dist", "Year")] <- lapply(capri.df[,c("Wing.Chord","Mass","B.Lat", "B.Long", "W.Lat", "W.Long", "Mig.dist", "Year")], as.numeric)
capri.df[c("Banding.Date", "B.dep", "W.arr")] <- lapply(capri.df[c("Banding.Date", "B.dep", "W.arr")], parse_date_time, c("mdy", "ymd"))

## Create month day for understanding timing
capri.df[,c("Band.md","Bdep.md", "Warr.md")] <- lapply(capri.df[,c("Banding.Date", "B.dep", "W.arr")], format, "%m/%d")
capri.df[,c("Band.md","Bdep.md", "Warr.md")] <- lapply(capri.df[,c("Band.md","Bdep.md", "Warr.md")], as.Date, "%m/%d")

# All years = present year
Year <- format(Sys.Date(), "%Y")
Year
# Add a year to maintain correct temporal distance (present year + 1)
capri.df$Warr.md <- zoo::as.Date(ifelse(capri.df$Warr.md < as.POSIXct(paste0(Year,"-04-01")), capri.df$Warr.md + lubridate::years(1), capri.df$Warr.md))

#Remove 13 individuals banded during the day (improbable capture time). 
capri.df2 <- capri.df %>% 
  filter(hms(Banding.Time) < hms("08:05:00") | 
           hms(Banding.Time) > hms("18:00:00") | 
           is.na(Banding.Time))

# Age classes condense ----------------------------------------------------
#Caprimulgids reach adult size (at least for feather length) after their second calendar year of life, so we need to condense age classes into Young, Adult, or Unknown

#Don't run this more than once
table(capri.df2$Age) #, capri.df$Species)
# #L & 1 are nestling / pullus, 3 = HY, 4 = AHY, 5 = SY, 6 = ASY, 8 = ASY
capri.df3 <- capri.df2 %>% filter(Age != "1" & Age != "L" & Age != "3") %>% 
  mutate(Age = case_when(Age == "4" ~ "Unk", #This should be adult? 
                         Age == "AHY" ~ "Unk",
                         Age == "ASY?" ~ "Unk",
                         Age == "5" ~ "Young",
                         Age == "SY" ~ "Young",
                         Age == "6" ~ "Adult",
                         Age == "4Y" ~ "Adult",
                         Age == "8" ~ "Adult",
                         Age == "A4Y" ~ "Adult",
                         Age == "A5Y" ~ "Adult",
                         Age == "ASY" ~ "Adult",
                         Age == "ATY" ~ "Adult",
                         Age == "TY" ~ "Adult"),
         Sex = trimws(Sex))

# Mass Correction ----------------------------------------------------------
#Nightjars will continue feeding throughout the night (dependent upon moonlight) & their mass will continue to increase. We capture individuals throughout the night which could introduce bias in body mass that is unrelated to the latitude of capture. Thus, we want to include the time since sunset as a covariate in body mass models to account for this bias.

#Add in row_number here to serve as primary key for linking back up with tsss
capri.df3 <- capri.df3 %>% mutate(rowID = row_number())

## Format for suncalc functions
df <- capri.df3 %>% summarize(date = as.Date(Banding.Date), 
                              lat = B.Lat, 
                              lon = B.Long,
                              DateTime = ymd_hms(paste(
                                as.character(Banding.Date), Banding.Time)
                                ),
                              Banding.Time = Banding.Time, 
                              rowID = rowID) 

#Extract time zone for each individual 
df <- df %>% mutate(tz=tz_lookup_coords(lat, lon, method="accurate"))

#Extract the time of sunset for each individual based on geographic coordinates & date
tzs <- unique(df$tz)
dftz <- list()
for(i in 1:length(tzs)){
  #Separate by timezone
  dftz[[i]] <- df %>% filter(tz == tzs[i])
  dftz[[i]]$sunset <- getSunlightTimes(data=dftz[[i]], keep="sunset", tz=tzs[i])$sunset
}

#Force the sunset column to be in UTC timezone. This allows for correct calculation of tsss b/c the banding times are assumed to be in UTC (we don't specify their actual tz)
dftz <- lapply(dftz, function(x){ force_tz(x, "UTC")})
sun_df <- bind_rows(dftz) 

#Subtract a day from sunset time if bird was caught early in the AM 
sun_df2 <- sun_df %>% mutate(sunset = as_datetime(ifelse(hms(Banding.Time) < hms("09:00:00"), 
                                                         ymd_hms(sunset) - lubridate::days(1), 
                                                         ymd_hms(sunset))))

#Calculate the tsss (time since sunset)
sun_df2$tsss <- as.numeric(difftime(sun_df2$DateTime, sun_df2$sunset), units="hours")
sun_df2 <- sun_df2 %>% mutate(across(where(is.numeric), round, 3))

capri.df4 <- merge(capri.df3, sun_df2[,c("rowID", "tsss", "sunset")], by = "rowID")

# Migratory dist -----------------------------------------------------------
#Calculate the great-circle distance from individuals breeding locations to winter location. NOTE:: These are not the actual migratory paths taken
capri.df4 <- transform(capri.df4, 
                       Str8line = distHaversine(
                         cbind(capri.df4$B.Long, capri.df4$B.Lat), 
                         cbind(capri.df4$W.Long, capri.df4$W.Lat)) / 1000)

# Band.Age -------------------------------------------
#Given that some individuals were captured multiple times we want to incorporate that information to have a better estimate of body size. Sizes differ by age so we want to combine morphological variables of individuals that are the same age.

#Create Band.Age variable, where individuals of the same age will have the same 'Band.Age', and use mutate to average morphologies (and time since sunset cov) by Band.Age
capri.df4 <- capri.df4 %>% mutate(Band.Age = paste0(Age, "_", Band.Number)) %>%
  group_by(Band.Age) %>% 
  mutate(Wing.comb = mean(Wing.Chord, na.rm = TRUE), 
         Mass.comb = mean(Mass, na.rm = TRUE)) 
capri.df5 <- capri.df4 %>% filter(!is.na(Banding.Time)) %>%
  mutate(Mass.combBT = mean(Mass, na.rm = T)) %>% 
  filter(!is.na(Mass)) %>% 
  mutate(tsss.comb = mean(tsss, na.rm = T)) %>% 
  rbind(capri.df4 %>% filter(is.na(Mass) | is.na(Banding.Time))) %>% 
  tidyr::fill(tsss.comb, Mass.combBT, .direction = "downup") 

# CapriBA --------------------------------------------------------
##CapriBA has only a single row for each Band & Age combo
capriBA <- capri.df5 %>% group_by(Band.Age) %>% #BA = band age
  arrange(is.na(W.Lat), Year, .by_group = TRUE) %>% 
  slice_head() %>% 
  relocate(uniqID, .after = rowID) %>%
  data.frame() 

# Export capriBA ---------------------------------------------------------
#Write capriBA.. Note uniqID is truly unique, but this df will still be further filtered
write.csv(capriBA, file = paste0("Intermediates/Capri_BA_", format(Sys.Date(), "%m.%d.%y"), ".csv"), row.names = FALSE)
