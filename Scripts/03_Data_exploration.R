###Data exploration###

library(ggrepel)

capri.df.int <- read_xlsx("Intermediate_products/capri.df.int_12.18.23.xlsx", trim_ws = TRUE) %>%
  mutate(Banding.Time = chron(times = Banding.Time))

# Random checks -----------------------------------------------------------
#Greg's precision w/ decimals, remove EDB in the end as it will add noise
Greg_dec <- str_sub(sapply(str_split(Various$B.Lat, "[.]"), function(x){x[2]}), start= -4)
data.frame(Various[c("Project")], Greg_dec) %>% filter(Project == "Greg_EDB")

##A number of random checks to better understand data and ensure it's good to go
#See if there are repeat individuals included across multiple data sets
multProj <- capri.df.int %>% count(Band.Number, Project) %>% 
  count(Band.Number, sort = T) %>% 
  filter(n>1)
#Identify where there could be additional overlap 
capri.df.int %>% count(Project, Country, Species) %>% 
  filter(Species == "EUNI") %>% 
  count(Country, sort = T)
#filter(n>1)

#This is showing no overlap in Band numbers outside of Stewart and IAN, and Greg. Both are OK 
capri.df.int %>% filter(Band.Number %in% multProj$Band.Number) %>% select(Project)
#Delete all of this 
JHeunj <- JHeunj %>% replace_with_na_all(condition = ~.x %in% c("na"))
jhband <- paste0("A", JHeunj$Band.Number)
table(jhband %in% Various[Various$Country == "Finland",]$Band.Number)
TF <- jhband %in% Various[Various$Country == "Finland",]$Band.Number
jhband[TF]
capri.df.int %>% mutate(Band.Number = ifelse(Project == "FMNH/Ceur", paste0("A", Band.Number), Band.Number)) %>% filter(Band.Number == "A770945" & Year == 2018) %>% select (Project, Year, Band.Number, B.Lat)

#Identify and remove(?) individuals w/ multiple sexes recorded
multSex <- capri.df.int %>% filter (Sex != "U") %>% 
  count(Band.Number, Sex) %>%
  count(Band.Number) %>% 
  filter(n>1) %>% pull(Band.Number)
capri.df.int %>% filter(Band.Number %in% multSex) %>% select(Band.Number, Sex, Project, W.Lat) %>% filter(Sex != "U")  %>% arrange(Band.Number) %>% count(Band.Number, Sex) %>% count(Band.Number) %>% filter(n == 2) %>% pull(Band.Number)
capri.df.int %>% filter(Band.Number == "LA42163")

#All coords should be 3+ digits, and breeding was not an issue but winter offenders wre Naskaswe & our own data! You can see now that there are no issues
decimals <- sapply(str_split(capri.df.int$W.Lat, "[.]"), function(x){x[2]})
#A few have 2 decimals due to rounding but this is OK 
data.frame(capri.df.int$Project, capri.df.int$Band.Number, decimals) %>% 
  filter(!is.na(decimals)) %>% 
  mutate(ncharW = nchar(decimals)) %>% 
  filter(ncharW < 4)

##Take a closer look at the birds with wintering data
facRed <- capri.df.int %>% filter(!is.na(W.Lat)) %>% 
  select(Project, Banding.Date, Banding.Time, Year, Band.Number, Age, 
         Sex, Mass, Wing.Chord, B.dep, W.arr, W.Lat, W.Long, Mig.dist) #fac reduced
lapply(facRed, function(x) {table(is.na(x))})
table(facRed$Age)
nrow(facRed)
facRed %>% filter(!duplicated(Band.Number)) %>% nrow() #Should end with 189 unique individuals!

#3 banding time NAs
facRed %>% filter(is.na(Banding.Time))
#Note 31 individuals w/out wing, mass or mig distance. None of these 
missDV <- facRed %>% filter(is.na(Wing.Chord) | is.na(Mass) | is.na(Mig.dist)) %>% pull(Band.Number)

facRed %>% filter(Band.Number %in% missDV) %>% View()
  count(Band.Number, sort = T) %>% 
  filter(n > 1)

capri.df.int %>% filter(Band.Number == 137257703)
capri.df.int %>% filter(Band.Number == 137228576)

# Data entry errors / outliers for Wing & mass ----------------------------
#Note there were some typos but these have been corrected now
capri.df.int %>% group_by(Species) %>% 
  summarize(mnWing = mean(Wing.Chord, na.rm = T), 
            mn.Mass = mean(Mass, na.rm = T))

capri.df.int %>% arrange(Wing.Chord) %>% 
  select(Project, Species, Age, Wing.Chord, Mass, Year, Band.Number) %>% 
  slice_head(n = 10)
capri.df.int %>% arrange(desc(Wing.Chord)) %>% 
  select(Project, Species, Age, Wing.Chord, Mass, Year, Band.Number) %>% 
  slice_head(n = 10)
capri.df.int %>% arrange(desc(Mass)) %>% 
  select(Project, Species, Age, Wing.Chord, Mass, Year, Band.Number) %>% 
  slice_head(n = 10)
capri.df.int %>% arrange(Mass) %>% 
  select(Project, Species, Age, Wing.Chord, Mass, Year, Band.Number) %>% 
  slice_head(n = 10)
capri.df.int %>% ggplot(aes(Mass, Species)) + geom_boxplot()
capri.df.int %>% ggplot(aes(Wing.Chord, Species)) + geom_boxplot()


# Correct mass ------------------------------------------------------------
hist(capri.df.int$tsss)

#1. Look at recaptures---- ###NEED to reorder this 
recap <- capri.df.int %>% 
  group_by(Band.Number, Year) %>% 
  summarize(captures = n()) %>% 
  # ungroup() %>% 
  dplyr::filter(captures > 1) %>% 
  left_join(capri.df.int)

#Visualize---
Sys.setenv(TZ='GMT')
capri.df.int %>% ggplot(aes(x = Banding.Time)) + 
  geom_histogram(color = "black") + 
  scale_x_chron(format="%H:%M")

#Can't figure out how to rename the facets of this plot, might be important
#Species <- c("Nighthawk", "Nightjar", "Whip-poor-will")
#names(Species) <-  c("Nighthawk", "Nightjar", "Whip-poor-will")
#, labeller = labeller(Species = Species) in facet_wrap call

capri.df.int %>% ggplot(aes(x=tsss, y = Mass)) +
  geom_jitter(alpha = .3) +
  geom_smooth() + 
  facet_wrap(~Species) + 
  labs(x = "Time since sunset (hrs)", y = "Mass")
ggsave("/Users/aaronskinner/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Grad_School/MS/EWPW/Writing_Exit_Seminar/Bergs_Rule/Results/Figures/Tsss_Mass_facet.png", bg = "white")

#Does latitude influence relationship between tsss & mass?
ggplot(capri.df.int) +
  geom_jitter(alpha = .3, aes(x=tsss, y=Mass, colour=factor(round(B.Lat, -1)))) +
  geom_smooth(aes(x=tsss, y=Mass, colour=factor(round(B.Lat, -1)))) +
  facet_wrap(~Species)

#Would be cool to show the average slope of all lines, but this is just the average slope of all points (can tell by commenting out geom_point & geom_line), so the overall line isn't really informative 
recap %>% filter(tsss > -1) %>% #Remove one outlier
ggplot(aes(x=tsss, y=Mass)) +
  #geom_point(aes(colour=Band.Number)) +
  #geom_line(aes(colour=Band.Number)) +
  geom_smooth() + 
  theme(legend.position = "none")


# >>Stomach -----------------------------------------------------------------
stomach <- read.csv("Data/Stomach/EvensLathouwers_data_sheet_corrected_stomach.csv") %>% 
  filter(!is.na(Stomach))
nrow(stomach)
stomach$Stomach

stomach$Banding.Time <- sapply(str_split(parse_date_time(stomach[,c("Banding.Time")], c("HMS"), truncated = 3), " "), function(x){x[2]})
stomach$Banding.Time <- chron(times = stomach$Banding.Time)
stomach <- stomach %>% mutate(Banding.Time = ifelse(Banding.Time < .5, Banding.Time + 1, Banding.Time))

#This is a pretty cool plot to include in manuscript as a supplementary figure, although would want to remember what the stomach score represents, and would want to convert to tsss 
Sys.setenv(TZ='GMT')
#tsss is better than Banding time as sunset is key for when foraging starts
stomach.tsss <- stomach %>% mutate(Band.Number = as.character(Band.Number)) %>% 
  left_join(capriBAnr[,c("Band.Number", "tsss")], by = "Band.Number")
stomach.tsss %>% ggplot(aes(x = tsss, y = Stomach)) + 
  geom_jitter(width = .05, alpha = .3) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se=T) +
  labs(x = "Time since sunset (hrs)", y = "Stomach fullness")
#scale_x_chron(format="%H:%M")

stomach.tsss2 <- stomach.tsss %>% mutate(Stomach = factor(Stomach, ordered = TRUE, levels = c(0:5))) %>% 
  filter(!is.na(tsss))
mod_stomach <- glm(Stomach ~ poly(tsss,2), data = stomach.tsss2, family = binomial)
summary(mod_stomach)

##Examine EUNI distribution of latitudes and how tsss impacts this 
capriBAnr %>% filter(!is.na(tsss) & Species == "EUNI") %>%
  #count(round(B.Lat, -1))
  ggplot(aes(x = B.Lat, fill = Country)) +
  geom_histogram()
names(capriBAnr)

capriBAnr %>% filter(Project == "EvensLathouwers" & is.na(tsss)) %>% 
  count(round(B.Lat, 0))
  

# Year effect -------------------------------------------------------------
#Visualize temporal spread of sampling by species
ggplot(data = capri.df.int, aes(x = Year, y = B.Lat)) + 
  geom_jitter(data = filter(capri.df.int, Species == "EUNI"), alpha = .2, 
              width = .4, height = .8, shape = 1, aes(color = Species)) +
  geom_jitter(data = filter(capri.df.int, Species %in% c("CONI","EWPW")), alpha = .4, 
              width = .4, height = .8, shape = 1, aes(color = Species)) +
  scale_color_hue(labels = c("Nighthawk", "Nightjar", "Whip-poor-will")) +
  labs(y = "Breeding latitude") + 
  geom_vline(aes(xintercept= 2010), color="black", linetype="dashed") 
#+ ggtitle("Temporal sampling resolution \n by breeding latitude")
ggsave("/Users/aaronskinner/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Grad_School/MS/EWPW/Writing_Exit_Seminar/Bergs_Rule/Results/Figures/Temporal_sampling_Blat.png", bg = "white")

euni.df <- capri.df.int %>% filter(Species == "EUNI" & Age != "Unk" & Year > 1989) 
nrow(euni.df)

#Visualize EUNI
ggplot(data = euni.df, aes(Year, B.Lat)) + 
  geom_hex()
ggplot(data = euni.df, aes(x = Year, y = B.Lat)) + 
  geom_smooth() + 
  geom_point(alpha = .3)
ggplot(data = euni.df, aes(x = Year, y = Mass)) + 
  geom_smooth(method = "lm") + 
  geom_point(alpha = .3)
ggplot(data = euni.df, aes(x = Year, y = Wing.Chord)) +
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm")

#Confirm effects are important even controlling for lat & age
summary(lm(Mass ~ scale(Year) + B.Lat + Age, data = euni.df)) 
summary(lm(Wing.Chord ~ scale(Year) + B.Lat + Age, data = euni.df)) 
summary(lm(B.Lat ~ scale(Year), data = euni.df)) 

#For now, I think >2010 is best option
capri.df.int <- capri.df.int %>% filter(Year >= 2010 )
euni.df <- euni.df %>% filter(Year >= 2010)

# Banding date effect -----------------------------------------------------
#Ensure banding dates are reasonable, notice all most extreme banding dates are EUNI
capri.df.int %>% filter(Species == "EUNI" ) %>% 
  dplyr::select(Species, Project, Country, Band.md, W.Lat, Band.Number) %>% 
  arrange(Band.md) %>% 
  slice_head(n = 10)
capri.df.int %>% filter(Species == "EUNI") %>% 
  dplyr::select(Species, Project, Country, Band.md, W.Lat, Band.Number) %>% 
  arrange(desc(Band.md)) %>% 
  slice_head(n = 10)

#Effect of band month statistically
summary(lm(Mass ~ Band.md + B.Lat + Age, data = euni.df)) #poly() gives error for band.md
#Suggests we're getting some individuals after they're molting and their wings are a bit longer?
summary(lm(Wing.Chord ~ Band.md + B.Lat + Age, data = euni.df)) 

#Plotting, Fat ~ band date
#Gabriel's EUNI data 
euni.df %>% filter(Project == "NASKASWE") %>%
  mutate(Fat = as.numeric(Fat)) %>% 
  ggplot(aes(x = Band.md, y = Fat)) + 
  geom_smooth() + 
  geom_point() + 
  ggtitle("Norevik EUNI data")
#ggsave('Fat_date_Norevik.png')

#Res = residents, 07-30 is final answer
euni.dfRes <- euni.df %>% 
  filter(Band.md > as.POSIXct("2023-05-16") & Band.md < as.POSIXct("2023-07-30"))

#Visualize, little to no relationship with band date at this point
ggplot(data = euni.dfRes, aes(x = Band.md, y = Wing.Chord)) + 
  geom_point(alpha = .3) + 
  geom_smooth()
ggplot(data = euni.dfRes, aes(x = Band.md, y = Mass)) + 
  geom_point(alpha = .3) + 
  geom_smooth() 

#Visualize EWPW & CONI species data
capri.df.int %>% filter(Species == "EWPW") %>% 
  ggplot(aes(x = Band.md, y = as.numeric(Fat))) + 
  geom_smooth() + 
  geom_point() + 
  ggtitle("EWPW fattening")
capri.df.int %>% filter(Species == "CONI") %>% 
  ggplot(aes(x = Band.md, y = as.numeric(Fat))) + 
  geom_smooth() + 
  geom_point() + 
  ggtitle("CONI fattening")

#Visualize EUNI
euni.df %>% filter(!is.na(Fat)) %>%
  count(Lat = round(B.Lat,0))
fat_p1 <- euni.df %>% ggplot(aes(x = Band.md, y = as.numeric(Fat))) + 
  geom_smooth() + 
  geom_point() + 
  labs(x = "Banding date", y = "Fat")
fat_p2 <- euni.dfRes %>% ggplot(aes(x = Band.md, y = as.numeric(Fat))) + 
  geom_smooth() + 
  geom_point() + 
  labs(x = "Banding date", y = "Fat")
ggpubr::ggarrange(fat_p1, fat_p2,
                  ncol = 2, nrow = 1,
                  align='hv', labels = c("A","B"),
                  legend = "none")
#ggsave("/Users/aaronskinner/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Grad_School/MS/EWPW/Writing_Exit_Seminar/Bergs_Rule/Results/Figures/EUNI fattening.png", bg = "white")

#Some way of seeing sequential variation explained? 
summary(lm(Mass ~ tsss + B.Lat + Band.md, data = euni.dfRes))
summary(lm(Mass ~ tsss + B.Lat + Band.md, data = ewpw))
summary(lm(Mass ~ tsss + B.Lat + Band.md, data = coni.m))
af <- anova(lm(Mass ~ tsss + B.Lat + Band.md, data = euni.dfRes))
af %>% mutate(pctExp = `Sum Sq` / sum(`Sum Sq`) * 100)

#Remove potential migrants, leaving just residents (res). 5974 individuals to 5593. Note there are some migrants that would be excluded based on these dates so have to do some filter & join gymnastics

#There are EUNI with wintering data that will be removed by the following dates so need to be added back into the data 
euni.df %>% filter(!is.na(W.Lat) & (Band.md < as.POSIXct("2023-05-16") | 
                                      Band.md > as.POSIXct("2023-07-30"))) %>% 
  count(Project)

# MOVE: Months Envi covs --------------------------------------------------
##Determine relevant months for environmental covariates. This ultimately used to create Table 1##
capri.df.int %>% filter(Warr.md < as.POSIXct("2024-03-01")) %>% dplyr::select(Species, Bdep.md, Warr.md) %>% group_by(Species) %>% summarize(N = n(), MeanDep = mean(Bdep.md, na.rm = T), sdDep = sd(Bdep.md, na.rm = T), MeanArr = mean(Warr.md, na.rm = T), sdArr = sd(Warr.md, na.rm = T))
#Determine dates when migrants may be present
capri.df.int %>% dplyr::select(Species, Project, B.Lat, Band.md, Bdep.md) %>% mutate(BlatR = round(B.Lat, -1)) %>% group_by(Species, BlatR) %>% summarize(n = n(), MinBand = min(Band.md, na.rm = T), MeanBand = mean(Band.md, na.rm = T), MaxBand = max(Band.md, na.rm = T), minDep = min(Bdep.md, na.rm = T), MeanDep = mean(Bdep.md, na.rm = T), sdDep = sd(Bdep.md, na.rm = T)) %>% mutate(Cutoff = MeanDep - sdDep, TF = Cutoff > MaxBand)
#Not an issue for CONI or EWPW. For EUNI though.. Birds leave as early as August 1 at 60+ degrees North. Average date of departure at 60N is 8-17, so let's subtract 1 sd and get cutoff date of 08-06. Alternatively, just remove everything banded after 08-01
##For spring migration for EUNI, both Evens (2017) and Norevik (2017) agree that birds depart wintering grounds in late February, but arrival date is "early May" for Evens (earliest April 28); and average is May 16 (earliest is May 5) for Norevik (2017), 12 birds in each paper but only 9 made it for spring migration in Evens. Given May 16 is so close to cut-off anyways, probably makes sense to include month of May. 
#For EWPW, cite my paper for spring departure as March 20th and arrival as April 17th, English (2017) departure = March 21, arrival = May 1. These 2 refs are in agreement. Korpach (2022) does not have spring migration dates.
#For CONI, Elly says to use Nov - March for winter and May - August for breeding.

# Multiple Wint locs  -------------------------------------------
##CapriBA (Band age) has only a single row for each individual & Age combo
capriBA <- capri.df.int %>% group_by(Band.Age) %>% 
  arrange(is.na(W.Lat), Year, .by_group = TRUE) %>% 
  slice_head() %>% 
  data.frame()
nrow(capriBA) 

#CHECK::Ensure that function is picking the correct row, if the Adult row has W.Lat then we're good
capriBA %>% filter(Band.Number == "135268737") %>% select(Band.Age, Banding.Time, Wing.Chord, Mass, Mass.comb, Wing.comb, tsss.comb, Mass.combBT, W.Lat) 

#Visualize morphological difference between captures
capriBA %>% select(Band.Number, Wing.Chord, Wing.comb) %>% 
  group_by(Band.Number) %>% 
  summarize(dif = Wing.Chord - Wing.comb) %>% 
  arrange(dif)
capri.df.int %>% select(Band.Number, Band.Age, Mass, Mass.comb, Species) %>% 
  group_by(Band.Age) %>% 
  summarize(difAvg = Mass - Mass.comb, difMax = max(Mass) - min(Mass), across()) %>% 
  arrange(desc(difMax)) #Difference of up to 28.1g in capture weights

#Visualize duplicate birds that have multiple years of winter data
DupBirds <- capri.df.int %>% group_by(Band.Age, Year) %>% 
  arrange(W.Lat, .by_group = TRUE) %>% 
  slice_head() %>% 
  data.frame() 
DupBirdsFac <- filter(DupBirds, !is.na(W.Lat))
dups <- data.frame(DupBirdsFac  %>% filter(duplicated(Band.Number) | duplicated(Band.Number, fromLast = TRUE)) %>% dplyr::select(Species, Project, Year, Band.Number, TagID, Age, W.Lat, W.Long, Wing.Chord, Mass) %>% arrange(Project, Band.Number))
#The W.Lat difference between years is never very large (see column "Wlat.diff")
df.Wlat.diff <- dups %>% group_by(Band.Number) %>% mutate(Wlat.diff = max(W.Lat) - min(W.Lat)) %>% arrange(Wlat.diff) %>% slice_head()

#Calculate distances between wintering locations
dups.sf <- st_as_sf(dups,
                    coords = c("W.Long", "W.Lat"),
                    crs = 4326)
uniqID <- unique(dups.sf$Band.Number)
IDs <- dist.all <- dist.all2 <- dist <- vector("list", length = length(uniqID))
for(i in 1:length(uniqID)){
  IDs[[i]] <- subset(dups.sf, Band.Number == uniqID[i])
  dist.all[[i]] <- max(st_distance(IDs[[i]]))
  #diag(dist.all[[i]]) <- NA
  #dist[i] <- max(round(colMeans(dist.all[[i]], na.rm = TRUE),2)) #In meters
}

dists <- data.frame(Band.Number = uniqID, Dists = round(unlist(dist.all) / 1000, 2)) 
df.dists <- merge(dists, df.Wlat.diff, by = "Band.Number") %>% arrange(Species, Dists)
df.dists %>% group_by(Species) %>% 
  summarize(mdn = median(Dists), mn = mean(Dists), sd = sd(Dists), N = n())

# Str8line vs mig distance ------------------------------------------------------
capriBAnr <- capriBA %>% group_by(Band.Number) %>% 
  arrange(is.na(W.Lat), Age) %>% 
  slice_head() #nr = no repeats 

#Create the FAC data frame with no repeat individuals 
capri.fac <- filter(capriBAnr, !is.na(W.Lat)) #FAC = full annual cycle
nrow(capri.fac) #189 unique birds..

TF <- capri.fac$Str8line < capri.fac$Mig.dist #13,14,19, 26
capri.fac[ !TF ,c("Str8line", "Mig.dist","Project", "Band.Number")] #Alicia used a different system to calculate migratory distance, so the numbers are slightly different. When she recalculated differences were <20km, so this is not very important on the scale of 4-5k km. 

#Overall, trying to determine whether we can use the actual migration distance or need to use the straight-line distance to account for differences in sampling resolution 
#Migration distance correlations by species
capri.fac %>% group_by(Species) %>% 
  summarize(cor = cor(Str8line, Mig.dist, use = "complete.obs"))

##Overall plot. The low correlation in EUNI means that str8line distance is not a good summary in EUNI, and that there may be mistakes in the actual migration distances. 
ggplot(data= capri.fac, aes(x = Str8line, y = Mig.dist)) + 
  geom_smooth(aes(color = Species), method = "lm", se = F, fullrange = T) + 
  geom_point(aes(color = Species), size = 3, position = "jitter", alpha=.3) + 
  geom_abline(slope= 1, linetype = "dashed", color="Black") + #+ xlab("Departure Date") + ylab("Migration\nRate (km / day)") 
  ggtitle("Mig dist calculations all species")

#Inspect EUNI more closely
capri.fac %>% filter(Species == "EUNI") %>% 
  ggplot(aes(x = Str8line, y = Mig.dist)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = Project)) +
  geom_abline(slope= 1, linetype = "dashed", color="Red") + 
  geom_text_repel(aes(label = TagID)) +
  ggtitle("EUNI Mig dist calculations")
#facet_wrap(~Species)

#Elly suggested plotting the geographic locations compared to migratory distances (size of circles). Idea was that this might help us locate any errors, but nothing really stands out, meaning it's likely OK to press forward assuming these migratory distance values are right. 
data.frame(capri.fac[,c("Species","Band.Number", "Mig.dist", "Project")], stack(capri.fac[,c("B.Lat", "W.Lat")]), stack(capri.fac[,c("B.Long", "W.Long")])) %>% 
  rename(Lat = values, Long = values.1) %>% 
  filter(Species == "EUNI") %>% 
  ggplot(aes(x = Long, y = Lat)) + geom_point(alpha  = .5,  aes(size = Mig.dist, color = Project)) + geom_line(alpha = .3, aes(group = Band.Number)) + scale_size(range = c(0, 8))
#ggsave("/Users/aaronskinner/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Grad_School/MS/EWPW/Writing_Exit_Seminar/Bergs_Rule/Results/Figures/EUNI_mig_dist_Lat_Long.png", bg = "white")

#Migration distance models. Does the variation in temporal schedule influence migration distance? Depends how you look at it.. Overall (model w/ RE), temp.res.mig is not significant once you include breeding latitude (which obviously influences migration distance).
#Model all species combined
summary(lmer(Mig.dist ~ Temp.res.mig + B.Lat + (1| Species), data = capri.fac))

#However, looking at each species separately you can see there is high sd in tag fix rate in EWPW, so let's look at EWPW more closely
capri.fac %>% group_by(Species) %>% 
  summarize(mean = mean(Temp.res.mig,na.rm = T),
            sd = sd(Temp.res.mig,na.rm = T),
            max = max(Temp.res.mig,na.rm = T),
            min = min(Temp.res.mig,na.rm = T))
summary(capri.fac$Temp.res.mig) #Ranges from 7 points a day to 1 point every 10 days

##IDEAL WOULD BE TO BRING IN ORIGINAL KORPACH DATA HERE AND SOMEHOW REJOIN WITH REST OF THE DATA TO BE ABLE TO CALCULATE BOTH SAMPLING RESOLUTIONS SIMULTANEOUSLY, BUT I SPENT GOOD AMOUNT OF TIME AND THIS ISN'T EASY (TRIED SMARTBIND, JOINS COULD WORK IF YOU FORMAT THE BAND NUMBERS THE SAME, ETC.). EASIEST IS JUST RECREATE CAPRI.DF FILE IN DATA_WRANGLING_01 SCRIPT AND NOT REPLACE THE ORIGINAL SAMPLING RESOLUTION DATA.

#By species, only EWPW is significantly impacted by temp sampling resolution.
Spp <- capri.fac %>% filter(Species == "EWPW") 
#This is complicated by the fact that breeding latitude is highly correlated with sampling resolution
cor(Spp$B.Lat, Spp$Temp.res.mig, use = "complete.obs") 
summary(lm(Mig.dist ~ B.Long + B.Lat + Temp.res.mig, data = Spp))

#Control for lat & long
af <- anova(lm(Mig.dist ~ B.Lat + B.Long + Temp.res.mig, data = Spp)) #Don't understand why PctExpplained decreases for B.Long when it becomes the first variable? 

#Overall, this suggests that temporal sampling resolution explains a small but significant amount of the variation in Mig.dist, even after controlling for lat & long. 
af %>% mutate(pctExp = `Sum Sq` / sum(`Sum Sq`) * 100)

#Could consider resampling to the lowest temporal resolution (5 days) just for EWPW?
Spp %>% group_by(Project) %>% 
  summarize(mean = mean(Temp.res.mig,na.rm = T),
            sd = sd(Temp.res.mig,na.rm = T),
            max = max(Temp.res.mig,na.rm = T),
            min = min(Temp.res.mig,na.rm = T))

# Leapfrog migration ------------------------------------------------------
#Evidence for leapfrog migration in EWPW and EUNI, but not CONI. Has this been reported previously in EUNI? 
ggplot(data= capri.fac, aes(x = B.Lat, y = W.Lat)) + 
  geom_smooth(method = "lm", se = TRUE, fullrange = F, aes(color = Species)) + 
  geom_point(aes(color = Species), size =3, position = "jitter", alpha=.65) 

#How would you test this statistically? 
summary(lmer(W.Lat ~ B.Lat + (1 | Species), data = capri.fac))
summary(lmer(W.Long ~ B.Long + (1 | Species), data = capri.fac))

#save.image(paste0(bs, "Desktop/Grad_School/R_Files/MS/BergAnalysis9.21.22.Rdata"))

