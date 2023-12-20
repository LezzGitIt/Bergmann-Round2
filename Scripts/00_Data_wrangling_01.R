##Bergmann's Rule in Caprimulgids##
##This script loads raw data from providers, combines data into single df, & tidies and filters data to produce various dataframes that will be used in downstream analyses


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(naniar)
library(readxl)
library(xlsx)
library(chron)
library(lutz)
library(suncalc)
library(stringi)

library(lme4)
library(lmerTest)
library(ggpubr)
library(viridis)
library(sp)
library(sf)
library(geosphere)
library(gtools)
library(cowplot)
library(AICcmodavg)
library(MuMIn)
library(gridExtra)
library(lmtest)
library(nlme)
library(ggrepel)
#library(conflicted)
select <- dplyr::select
ggplot2::theme_set(theme_cowplot())

#load(paste0(bs, "Desktop/Grad_School/R_Files/MS/BergAnalysis7.27.22.Rdata"))

# Load and format data ----------------------------------------------------
#load(paste0(bs,"Desktop/Grad_School/MS/EWPW/Data_Whips/Analysis/R_Files/spdf.all.tr.Rdata"))
EKconiFAC <- read.csv("Data/Raw_Data/EK_DataSheet_CONI_Feb16.csv") #Notice there are 7 or 8 rows of repeat individuals 2nd wintering locations. This is removed further down w/ band_age given that it's same band # and same deployment year
EKconiBr <- read.csv("Data/Raw_Data/EK_CONI_Breeding.csv")
ASewpw <- read.csv("Data/Raw_Data/Tonra_Ward_EWPW_Berg_Combined.csv")
AKewpw <-read.csv("Data/Raw_Data/AK_Bergs_19Sep2023.csv")
MBewpw <-read.csv( "Data/Raw_Data/MB_VItz_EWPW_Massachusetts.csv")
JHeunj <-read.csv( "Data/Raw_Data/JH_Finland_v6.csv")
LJeunj <- read.csv("Data/Raw_Data/Nightjar_data_Denmark_Thorup_Jacobsen/LJ_eunj_2023_DK.csv")##SEE NEW sheet w/ updated age code
GCeunj <- read.csv("Data/Raw_Data/GJC_Bergmann_UK_updated.csv")
GNeunj <- read.csv("Data/Raw_Data/GN_Sweden_2023-04-20.csv")
REeunj <- read.csv("Data/Raw_Data/Evens_EUNI/EvensLathouwers_data_sheet_corrected_stomach.csv")
ITeunj <- read.csv("Data/Raw_Data/Boano_Italy_2023-03-23.csv")
VariousGC <- read.csv("Data/Raw_Data/Various_EUNI_Breeding_only.csv")

#Make necessary changes for data formatting
cols <- c("Project", "Year", "Banding.Date", "Recap.", "Banding.Time", "Species", "Site.name", "Country", "Band.Number","TagID", "Age", "Sex", "CP","BP","Fat", "Wing.Chord", "WingFlat", "Tail.Length", "Tarsus","Mass","B.Lat","B.Long","W.Lat", "W.Long", "B.dep","W.arr", "Mig.dist", "Mig.n", "Temp.res.mig")

##Handle the pecularities of each data set to allow for combination and more efficient cleaning##
#Korpach's database
#Confirmed that all mig.dist.filt values are < mig.dist values (abline slope = 1), which makes sense
migs <- c("Mig.dist", "Mig.n", "Temp.res.mig")
AKewpw <- AKewpw %>% select(-migs) %>%
  rename_with(.cols = contains("filt"), .fn = ~str_remove(., pattern = ".filt")) %>% 
  mutate(Project = "Korpach")
#Ensure columns are in the correct order (names don't have to be identical)
AKnames <- names(AKewpw)[1:29]
data.frame(AKnames, cols, AKnames == cols)

#Greg Various data sets
table(VariousGC$Project)
Various <- VariousGC %>% mutate(Year = sapply(str_split(Banding.Date, "/"), function(x){x[3]}), 
                                Project = paste0("Greg_", Project),
                                Banding.Date = as.character(dmy(Banding.Date))) %>% 
  select(1,30, 2:29)

#Merge Elly's data to create a single data frame 
#EKconi2 <- smartbind(EKconiFAC, EKconiBr)
EKconiFAC <- EKconiFAC %>% filter(Wint.Loc == 1)
EKconiFAC$uniq.ID <- with(EKconiFAC, paste0(BandNumber, "_", TagID))
EKconiBr$uniq.ID <- with(EKconiBr, paste0(BandNumber, "_", TagID))
EKconi <- merge(EKconiBr, EKconiFAC[c("uniq.ID", "W.Lat", "W.Long", "Wint.Loc", "Dist.Fall", "N.Fall", "Dist.Spring", "N.Spring", "B.dep" , "W.arr", "W.dep", "B.arr")], by = "uniq.ID", all.x = T)
EKconi <- EKconi %>% 
  mutate(Temp.res.mig = ifelse(Year == 2017 | Year == 2018, 7, 10),
         WingChord = WingChord * 10,
         WingChord = ifelse(BandNumber == "1212-58492", 198, WingChord), #Fix typo
         TailLength = TailLength * 10, 
         Mass = ifelse(Mass == 0, NA, Mass),
         B.dep = as.character(as.Date(B.dep,    
                                      origin = as.Date(paste0(Year, "-01-01")))),
         W.arr = as.character(as.Date(W.arr,   
                                      origin = as.Date(paste0(Year, "-01-01")))),
         BandTime = ifelse(BandNumber =="1352-97758" | BandNumber == "1352-97768", "21:15:00", BandTime)) %>%
  #Remove 12 individuals that Elly says may not be unique captures
  filter(!(SiteName == "Alberta" & B.Lat > 57 & BandNumber == "Unbanded" & is.na(W.Lat))) %>%
  select(c(2:23,25,26,32,33,28,29,36,24)) %>%
  mutate(BandNumber = ifelse(BandNumber == "Unbanded" & !is.na(TagID), TagID, BandNumber)) %>%
  mutate(BandNumber = ifelse(BandNumber == "Unbanded" | BandNumber == "-99", paste0("Unbanded", c(1:nrow(EKconi))), BandNumber))
#CHECK::Ensure that the TagID replaced 'Unbanded' and that Unbandeds are now unique
EKconi %>% arrange(BandNumber) %>% select(BandNumber, TagID)


#Greg 9 FAC birds
GCeunj <- GCeunj %>% 
  mutate(Project = "UK_GJC", 
         Wing.Chord = WingFlatStraight,
         WingFlat = "Y", 
         W.Lat = ifelse(Band.Number == "RR05040", NA, W.Lat),
         W.Long = ifelse(Band.Number == "RR05040", NA, W.Long),
         W.arr = ifelse(Band.Number == "LA42163" | Band.Number == "RR05040", NA, W.arr), 
         W.arr = ifelse(Band.Number == "LE31898", "3/11/2021", W.arr)) %>% 
  select(-WingFlatStraight)
#Greg's is only dates in the wrong format
GCeunj[c("Banding.Date", "B.dep", "W.arr")] <- lapply(GCeunj[c("Banding.Date", "B.dep", "W.arr")], function(x){as.character(dmy(x))}) #1 failing to parse is OK 

#Lars 
LJeunj <- LJeunj %>% mutate(Species = "EUNI",
                            Project = "LJdenmark", 
                            B.dep = ifelse(B.dep == "1-sep", "9/1/10", B.dep))

# Create single df --------------------------------------------------------
#Create list, adjust columns and colnames, create single dataframe
njdfs_all <- list(EKconi, ASewpw, AKewpw, MBewpw, JHeunj, LJeunj, GCeunj, GNeunj, REeunj, ITeunj, Various)
names(njdfs_all) <- c("EKconi", "ASewpw", "AKewpw", "MBewpw", "JHeunj", "LJeunj", "GCeunj", "GNeunj", "REeunj", "ITeunj", "Various")
njdfs_all <- lapply(njdfs_all, function(x){x[,c(1:29)]})
njdfs_all <- lapply(njdfs_all, setNames, cols)
capri.df <- rbind(njdfs_all[[1]], njdfs_all[[2]], njdfs_all[[3]],njdfs_all[[4]], njdfs_all[[5]], njdfs_all[[6]], njdfs_all[[7]], njdfs_all[[8]], njdfs_all[[9]], njdfs_all[[10]], njdfs_all[[11]])
#lapply(njdfs_all, function(x){head(x$Banding.Date)})
#lapply(njdfs_all, function(x){parse_date_time(x[,c("Banding.Date")], c("mdy", "ymd"))})

#Adjust time & date
capri.df <- capri.df %>% #Few individuals removed w/ no B.Lat
  replace_with_na_all(condition = ~.x %in% c(-99,-990, 9999, "<NA>", "-", ".", "na", 'NONABAND')) %>% 
  as.data.frame() %>% 
  filter(!is.na(B.Lat) & !is.na(Band.Number))
nrow(capri.df) 
capri.df <- capri.df %>% 
  mutate(Species = ifelse(capri.df$Species == "Ceur" | capri.df$Species == "European Nightjar" | capri.df$Species == "European Nigthtjar", "EUNI", capri.df$Species), 
         Band.Number = stri_replace_all_regex(capri.df$Band.Number,
                                              pattern=c('-', ' '),
                                              replacement=c(''),
                                              vectorize=FALSE)) %>%
  mutate(uniqID = paste0(Band.Number, "_", Banding.Date, "_", !is.na(W.Lat))) %>% 
  arrange(is.na(W.Lat), Species, Project) #%>% 
#mutate(rowID = row_number()) #Remove nestlings first then create rowID

capri.df$Banding.Time <- str_pad(capri.df$Banding.Time, 4, pad = "0")
capri.df$Banding.Time <- sapply(str_split(parse_date_time(capri.df[,c("Banding.Time")], c("HMS"), truncated = 3), " "), function(x){x[2]})
capri.df$Banding.Time <- chron(times = capri.df$Banding.Time)
capri.df$Year <- str_pad(capri.df$Year, 3, pad = "0")
capri.df$Year <- str_pad(capri.df$Year, 4, pad = "2")
#CHECK::Ensure that this worked
capri.df %>% select(Project, Year) %>% 
  mutate(ncharYr = nchar(as.character(Year))) %>% 
  arrange(ncharYr, desc(Year))

#Format times & dates, data types#
capri.df[,c("Wing.Chord","Mass","W.Lat", "W.Long", "Mig.dist", "Year")] <- lapply(capri.df[,c("Wing.Chord","Mass","W.Lat", "W.Long", "Mig.dist", "Year")], as.numeric)
capri.df[c("Banding.Date", "B.dep", "W.arr")] <- lapply(capri.df[c("Banding.Date", "B.dep", "W.arr")], parse_date_time, c("mdy", "ymd")) 
#Create month day (all years = 2023) for understanding timing
capri.df[,c("Band.md","Bdep.md", "Warr.md")] <- lapply(capri.df[,c("Banding.Date", "B.dep", "W.arr")], format, "%m/%d")
capri.df[,c("Band.md","Bdep.md", "Warr.md")] <- lapply(capri.df[,c("Band.md","Bdep.md", "Warr.md")], as.Date, "%m/%d")
capri.df$Warr.md <- as.Date(ifelse(capri.df$Warr.md < as.POSIXct("2023-04-01"), capri.df$Warr.md + lubridate::years(1), capri.df$Warr.md)) 

#CHECK::Are there individuals banded during the day? Do they have wintering information? 
day.caps <- capri.df %>% 
  filter(hms(Banding.Time) > hms("06:00:00") & hms(Banding.Time) < hms("18:00:00")) %>%
  select(Project, Band.Number, Banding.Time, W.Lat, Age) %>% #tsss,
  arrange(Band.Number, Project, Banding.Time) 
day.caps %>% filter(!is.na(W.Lat))

#CHECK::Are there individuals banded during the day? Do they have wintering information? 
day.caps <- capri.df %>% 
  filter(hms(Banding.Time) > hms("06:00:00") & hms(Banding.Time) < hms("18:00:00")) %>%
  select(Project, Band.Number, Banding.Time, W.Lat, Age) %>% #tsss,
  arrange(Band.Number, Project, Banding.Time) 
day.caps %>% filter(!is.na(W.Lat))

#Remove 13 individuals banded during the day. DO NOT REMOVE FOR ROWID CHECKELLY DF 
capri.df2 <- capri.df %>% 
filter(hms(Banding.Time) < hms("08:05:00") | 
         hms(Banding.Time) > hms("18:00:00") | 
         is.na(Banding.Time))
nrow(capri.df2)

# Age classes condense ----------------------------------------------------
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
table(capri.df3$Age) 
nrow(capri.df3)

# MassCorrection ----------------------------------------------------------
##See other script 
capri.df3 <- capri.df3 %>% mutate(rowID = row_number()) #8.14.23 Adding row_number() here, but notice row_number was previously added above. NEED TO CONFIRM BEFORE MERGING ENVIRONMENTAL DATA W/ ELLY
#nrow(capri.df2) #5914 rows

#1. Read in data---
df <- capri.df3 %>% summarize(date = as.Date(Banding.Date), 
                              lat = B.Lat, 
                              lon = B.Long,
                              DateTime = ymd_hms(paste(as.character(Banding.Date), Banding.Time)),
                              Banding.Time = Banding.Time, 
                              Project = Project,
                              rowID = rowID) %>% 
  mutate(tz=tz_lookup_coords(lat, lon, method="accurate"))

tzs <- unique(df$tz)
dftz <- list()
for(i in 1:length(tzs)){
  dftz[[i]] <- df %>% filter(tz == tzs[i])
  #sunsets[[i]] <- getSunlightTimes(data = dftz)$sunset
  dftz[[i]]$sunset <- getSunlightTimes(data=dftz[[i]], keep="sunset", tz=tzs[i])$sunset
}

dftz <- lapply(dftz, function(x){ force_tz(x, "UTC")})
sun_df <- bind_rows(dftz) #Confirmed this is good 4.26.23
#CHECK:: Ensure sunset times are reasonable
sun_df %>% arrange(format(sun_df$sunset, format = "%H:%M:%S")) %>% pull(sunset)

#Subtract a day from sunset time if bird was caught early in the AM 
sun_df2 <- sun_df %>% mutate(sunset = as_datetime(ifelse(hms(Banding.Time) < hms("09:00:00"), 
                                     ymd_hms(sunset) - lubridate::days(1), 
                                     ymd_hms(sunset))))

#CHECK: Number of NAs in sunset should be equal to NAs in Banding time
sun_df2 %>% filter(is.na(sunset)) %>% nrow()
capri.df3 %>% filter(is.na(Banding.Date) | is.na(Banding.Time)) %>% nrow()

#Calculate the tsss (time since sunset)
sun_df2$tsss <- as.numeric(difftime(sun_df2$DateTime, sun_df2$sunset), units="hours")
sun_df2 <- sun_df2 %>% mutate(across(where(is.numeric), round, 3))

capri.df4 <- merge(capri.df3, sun_df2[,c("rowID", "tsss", "sunset")], by = "rowID")
#Visualize, note hour 0 is sunset 
hist(capri.df4$tsss)

# Str8line dist -----------------------------------------------------------
capri.df4 <- transform(capri.df4, Str8line = 
                         distHaversine(cbind(capri.df4$B.Long, capri.df4$B.Lat), 
                                       cbind(capri.df4$W.Long, capri.df4$W.Lat)) / 1000)

# Band.Age -------------------------------------------
#Create Band.Age variable, and use mutate to combine morphologies by Band.Age
capri.df4 <- capri.df4 %>% mutate(Band.Age = paste0(Age, "_", Band.Number)) %>%
  group_by(Band.Age) %>% 
  mutate(Wing.comb = mean(Wing.Chord, na.rm = TRUE), 
         Mass.comb = mean(Mass, na.rm = TRUE)) 
capri.df5 <- capri.df4 %>% filter(!is.na(Banding.Time)) %>%
  mutate(Mass.combBT = mean(Mass, na.rm = T)) %>% 
  filter(!is.na(Mass)) %>% 
  mutate(tsss.comb = mean(tsss, na.rm = T)) %>% 
  rbind(capri.df4 %>% filter(is.na(Mass) | is.na(Banding.Time))) %>% 
  tidyr::fill(tsss.comb, Mass.combBT, .direction = "downup") #Capri Mutated. Tail.comb = mean(Tail.Length, na.rm = TRUE)

#How many times was each individual caught on average
capri.df5 %>% group_by(Species) %>% 
  count(Band.Age) %>%
  summarize(mean = mean(n), sd = sd(n))
capri.df5 %>% ungroup() %>% 
  count(Band.Age) %>%
  summarize(mean = mean(n), sd = sd(n), min = min(n), max = max(n))

#CHECK::This individual of Alicia's shows it worked for mass 
capri.df5 %>% 
  filter(Band.Number == "135268737") %>% 
  select(Band.Age, Banding.Time, Wing.Chord, Mass, Mass.comb, tsss.comb, Mass.combBT, W.Lat) %>% 
  arrange(Mass.combBT)
#This individual shows it worked for tsss.comb 
capri.df5 %>% filter(Band.Number == "114281962") %>% 
  select(Band.Age, Banding.Time, Mass, Mass.comb, tsss, tsss.comb, Mass.combBT)
mean(c(.109,.698,.461,.33,1.56)) #tsss.comb equal to this
#Another example
capri.df5 %>% filter(Band.Number == "135268713") %>% 
  select(Band.Age, Banding.Time, Mass, Mass.comb, tsss, tsss.comb, Mass.combBT)
mean(c(57, 68, 61)) #Mass.combBT equal to this
mean(c(1.25, 2.68, .655)) #tsss.comb equal to this

# CapriBA --------------------------------------------------------
##CapriBA has only a single row for each Band & Age combo
capriBA <- capri.df5 %>% group_by(Band.Age) %>% #BA = band age
  arrange(is.na(W.Lat), Year, .by_group = TRUE) %>% 
  slice_head() %>% 
  relocate(uniqID, .after = rowID) %>%
  select(-rowID) %>%
  data.frame() 
#Ensure that each uniqID is truly unique
nrow(capriBA) 
length(unique(capriBA$uniqID))

#Write capriBA.. Note uniqID is truly unique, but this df will still be further filtered
capriBA %>% as.data.frame() %>%
  write.xlsx("Intermediate_products/Capri_BA_12.19.23.xlsx", row.names = F)


# capriBAnr ---------------------------------------------------------------
#Remove repeat individuals, selecting adults when possible 
levels(capriBA$Age) #Should get NULL, this should NOT be a factor in order for arrange to work correctly
capriBAnr <- capriBA %>% group_by(Band.Number) %>% 
  arrange(is.na(W.Lat), Age) %>% 
  slice_head() #nr = no repeats 
nrow(capriBAnr)

#Filter year
capriBAnr2 <- capriBAnr %>% filter(Year >= 2010)

#Filter based on date
nrow(capriBAnr2)
euniRes <- capriBAnr2 %>% filter(Species == "EUNI" & Band.md > as.POSIXct("2023-05-16") & Band.md < as.POSIXct("2023-08-01") & is.na(W.Lat))
euniFAC <- capriBAnr2 %>% filter(Species == "EUNI" & !is.na(W.Lat))
coni.ewpw <- capriBAnr2 %>% filter(Species != "EUNI" & Band.md > as.POSIXct("2023-04-30")) 
capriA <- rbind(euniRes, euniFAC, coni.ewpw)
nrow(capriA) #capri Analysis



