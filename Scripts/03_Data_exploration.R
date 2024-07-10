##Bergmann's rule across full-annual-cycle in Caprimulgids##
##Data exploration 03 -- 
##Code to confirm data accuracy, explore data, look for biases that could influence analysis, & produce figures in the supporting information 

#Contents: Miscellaneous, see table of contents

# Load libraries & dfs ----------------------------------------------------
library(tidyverse)
library(readxl)
library(xlsx)
library(chron)
library(ggrepel)
library(viridis)
library(lme4)
map <- purrr::map

#load("Rdata/Capri_dfs_3.14.24.Rdata")

capri.df.int <- read_xlsx("Intermediate_products/capri.df.int_07.03.24.xlsx", trim_ws = TRUE) %>%
  mutate(Banding.Time = chron(times = Banding.Time))
capriBA <- read_xlsx("Intermediate_products/Capri_BA_07.03.24.xlsx", trim_ws = TRUE) %>%
  mutate(Banding.Time = chron(times = Banding.Time))

#No repeats in 'nr' file
capriBAnr <- capriBA %>% group_by(Band.Number) %>% #nr = no repeats 
  arrange(is.na(W.Lat), Age) %>% 
  slice_head() %>% 
  ungroup()

#Create the FAC data frame with no repeat individuals 
capri.fac <- filter(capriBAnr, !is.na(W.Lat)) #FAC = full annual cycle
nrow(capri.fac) #189 unique birds..


# Explore SA / V & tsss^2 -------------------------------------------------
#Note that although SA / V is not significant in any case that could imply that is interesting in of itself as it says something about the rate that mass & wing are increasing together
map(njdf.list.br, function(df) {
  df$SA <- df$Wing.comb^2 / df$Mass.combBT
  summary(lm(SA ~ B.Lat, data = df))
})

#Could scaling data be problematic?
njdf.br.am %>% group_by(Species) %>%
  summarize(mean_Wing = mean(Wing.comb, na.rm = T), 
            mean_Mass = mean(Mass.combBT, na.rm = T),
            sd_Wing = sd(Wing.comb, na.rm = T), 
            sd_Mass = sd(Mass.combBT, na.rm = T), 
            ratio_wing = mean_Wing / sd_Wing,
            ratio_mass = mean_Mass / sd_Mass)

##What are the effects of including tsss^2?##
#To start, examine the correlation between Latitude, Mass & tsss & tsss^2 (1 & 2 in output)
map(njdf.list.br[c(1,3,5)], function(df) {
  df <- cbind(df, poly(df$tsss.comb, 2, raw = FALSE)) %>%  # poly(scale(df$tsss.comb, center = TRUE, scale = FALSE), 2, raw = TRUE))
    rename_with(.cols = c("Mass.combBT", "1","2"), ~c("Mass", "tsss", "tsss^2"))
  round(cor(df[,c("B.Lat", "Mass", "tsss", "tsss^2")]), 2)
})

map(njdf.list.br[c(1,3,5)], function(df) {
  round(cor(df[,c("B.Lat", "Mass.combBT", "tsss.comb")]), 2)
})

map(njdf.l.br.age[c(1,3,5)], function(df) {
  df <- cbind(df, poly(df$tsss.comb, 2, raw = FALSE)) %>% 
    rename_with(.cols = c("1","2"), ~c("tsss", "tsss2"))
  summary(lm(Mass.combBT ~ tsss, data = df))
  #ggplot(data = df, aes(x = tsss.comb, y = B.Lat)) + 
    #geom_point() +
    #stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
    #xlim(c(0, 8.5))
})

map(njdf.l.br.age[c(1,3,5)], function(df) {
  boxplot(log(df$tsss.comb))
})


#Examine the amount of variance explained by breeding latitude in models with lat + tsss^2 vs models with just latitude
#Models with tsss
map(njdf.list.br[c(1,3,5)], function(df) {
  mod <- lm(Mass.combBT ~ poly(tsss.comb,2) + B.Lat , data = df)
  af <- anova(mod)
  cbind(af, PctExp = af$"Sum Sq" / sum(af$"Sum Sq") * 100)
})
#Variance explained: 9, 5, 16

#Models with just latitude
map(njdf.list.br[c(1,3,5)], function(df) {
  mod <- lm(Mass.combBT ~ B.Lat, data = df)
  af <- anova(mod)
  cbind(af, PctExp = af$"Sum Sq" / sum(af$"Sum Sq") * 100)
})
#Variance explained: 23, 13, 19

#Compare the estimates of the slope for latitude and associated standard errors
#NOTE:: Final models have Sex (& Age in most cases) so coefficient estimates are different. But at the moment I'm unable to reproduce the same coefficient estimate for nighthawks for Breeding latitude.. Going to ignore for now but ultimately should be able to reproduce the parameter coefficients when you add in the appropriate nuisance variables
map(njdf.list.br[c(1,3,5)], function(df) {
  mod1 <- lm(Mass.combBT ~ poly(tsss.comb,2) + B.Lat, data = df)
  mod1_df <- data.frame(Model = "Mass ~ tsss^2 + Lat", tidy(mod1))
  mod2 <- lm(Mass.combBT ~ B.Lat , data = df)
  mod2_df <- data.frame(Model = "Mass ~ Lat", tidy(mod2))
  bind_rows(mod1_df, mod2_df) %>% mutate_if(is.numeric, round, 2) %>% 
    filter(term != "(Intercept)")
})

# Random checks -----------------------------------------------------------
#Greg's precision w/ decimals, remove EDB in the end as it will add noise
VariousGC <- capri.df.int %>% filter(str_detect(Project, "^Greg"))
Greg_dec <- str_sub(sapply(str_split(VariousGC$B.Lat, "[.]"), function(x){x[2]}), start= -4)
data.frame(VariousGC[c("Project")], Greg_dec) %>% filter(Project == "Greg_EDB")

##A number of random checks to better understand data and ensure it's good to go
#See if there are repeat individuals included across multiple data sets
multProj <- capri.df.int %>% count(Band.Number, Project) %>% 
  count(Band.Number, sort = T) %>% 
  filter(n>1)
#Identify where there could be additional overlap 
capri.df.int %>% count(Project, Country, Species) %>% 
  filter(Species == "Nightjar") %>% 
  count(Country, sort = T)
#filter(n>1)

#This is showing no overlap in Band numbers outside of Stewart and IAN, and Greg. Both are OK 
capri.df.int %>% filter(Band.Number %in% multProj$Band.Number) %>% select(Project)
#Delete all of this 
JHeunj <- capri.df.int %>% filter(Project == "FMNH/Ceur")
JHeunj <- JHeunj %>% replace_with_na_all(condition = ~.x %in% c("na"))
jhband <- paste0("A", JHeunj$Band.Number)
table(jhband %in% VariousGC[VariousGC$Country == "Finland",]$Band.Number)
TF <- jhband %in% VariousGC[VariousGC$Country == "Finland",]$Band.Number
jhband[TF]
capri.df.int %>% 
  mutate(Band.Number = ifelse(Project == "FMNH/Ceur", paste0("A", Band.Number), Band.Number)) %>%
  filter(Band.Number == "A770945" & Year == 2018) %>% select (Project, Year, Band.Number, B.Lat)

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
#Note 31 individuals w/out wing, mass or mig distance
missDV <- facRed %>% filter(is.na(Wing.Chord) | is.na(Mass) | is.na(Mig.dist)) %>% pull(Band.Number)

facRed %>% filter(Band.Number %in% missDV) %>%
  count(Band.Number, sort = T) %>% 
  filter(n > 1)

#Note some of these have information from other years so it's not a big deal in those cases
capri.df.int %>% filter(Band.Number == 137257703) %>% 
  select(Project, Age, Wing.Chord, Mass, Mig.dist)
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

#Generate plots showing quadratic relationship between tsss & Mass 
capri.df.int %>% ggplot(aes(x=tsss, y = Mass)) +
  geom_jitter(alpha = .3) +
  geom_smooth() + 
  facet_wrap(~Species) + 
  labs(x = "Time since sunset (hrs)", y = "Mass")
ggsave("Plots/Tsss_Mass_facet.png", bg = "white")

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

#Correct mass just in breeding data, or also in FAC data? 
capri.fac %>% ggplot(aes(x=tsss.comb, y = Mass.combBT)) +
  geom_jitter(alpha = .3) +
  #geom_smooth(se=FALSE, method = "glm", formula= y ~ poly(x,2)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Species) + 
  labs(x = "Time since sunset (hrs)", y = "Mass")
ggsave("Plots/FAC_tsss_mass_quadratic.png", bg = "white")

#If there is no relationship between B.Lat & tsss this would imply that the variation in tsss is fairly random and should not systematically affect size at any given latitude. This would suggest it is probably not super important to control for.. However, remember it is not just latitude that is relevant, but the linear combination of several variables.. Including tsss is also more generative I think, we're including a variable that helps explain some of the noise 
gg_plot(df = capriA.red2, x = B.Lat, y = tsss.comb) + #method = "lm"
  facet_wrap(facet = ~Species) +
  labs(x = "Breeding latitude")

capriA.red2 %>% group_by(Species) %>% 
  summarize(cor.lat.tsss = cor(B.Lat, tsss.comb, use = "complete.obs"), n = n())

# Stomach -----------------------------------------------------------------
stomach <- read.csv("Data/Stomach/EvensLathouwers_data_sheet_corrected_stomach.csv") %>% 
  filter(!is.na(Stomach))
nrow(stomach)
stomach$Stomach

stomach$Banding.Time <- sapply(str_split(parse_date_time(stomach[,c("Banding.Time")], c("HMS"), truncated = 3), " "), function(x){x[2]})
stomach$Banding.Time <- chron(times = stomach$Banding.Time)
stomach <- stomach %>% mutate(Banding.Time = ifelse(Banding.Time < .5, Banding.Time + 1, Banding.Time))

#This is a pretty cool plot to include in manuscript as a supplementary figure, although would want to remember what the stomach score represents
Sys.setenv(TZ='GMT')

stomach.tsss <- stomach %>% mutate(Band.Number = as.character(Band.Number)) %>% 
  left_join(capriBA[,c("Band.Number", "tsss")], by = "Band.Number")
stomach.tsss %>% filter(!is.na(tsss) & !is.na(Stomach)) %>% nrow()
#499 European nightjars
stomach.tsss %>% ggplot(aes(x = tsss, y = Stomach)) + 
  geom_jitter(width = .05, alpha = .3) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se=T) +
  labs(x = "Time since sunset (hrs)", y = "Stomach fullness")
ggsave("Plots/Supplementary/stomach_tsss.png", bg = "white")

#There is a highly significant relationship
stomach.tsss2 <- stomach.tsss %>% mutate(Stomach = factor(Stomach, ordered = TRUE, levels = c(0:5))) %>% 
  filter(!is.na(tsss))
mod_stomach <- glm(Stomach ~ poly(tsss,2), data = stomach.tsss2, family = binomial)
summary(mod_stomach)

##Examine EUNI distribution of latitudes and how tsss impacts this 
capriBAnr %>% filter(!is.na(tsss.comb) & Species == "Nightjar") %>%
  #count(round(B.Lat, -1))
  ggplot(aes(x = B.Lat, fill = Country)) +
  geom_histogram()

capriBAnr %>% filter(Project == "EvensLathouwers" & is.na(tsss)) %>% 
  count(round(B.Lat, 0)) 


# Age / Sex effect --------------------------------------------------------
#Examine how ages & sexes are distributed across latitudes, as well as the impact of age & sex on body size for each species

#First, Examine the raw counts of each species 
count.NV <- function(NV){
  capriA.red2 %>% group_by(Species) %>% 
    count(.data[[NV]])
}

NV <- c("Age", "Sex")
counts_NV <- lapply(NV, count.NV)

#Function to plot the distribution of nuisance variables by latitude
plot.NV.dist <- function(Spp, NV){  #NV.dis = Distribution of nuisance variables
  df <- capriA.red2 %>% filter(Species == Spp)
  #if(Spp == "Nighthawk"){df <- njdf.list.br[[paste0(Spp, ".", DV)]]}
  ggplot(data = df, aes(x = B.Lat, color = .data[[NV]])) +
    geom_density(alpha=0.5) + #, aes(y = after_stat(scaled))
    labs(x = "Breeding latitude", y = "Density", title = Spp) 
}
#Generate loop to cycle through this & next function 
loopNV <- data.frame(rbind(loopSppDV, loopSppDV), NV = c(rep("Age", 6), rep("Sex", 6))) %>% 
  arrange(Species)
#Plot distribution of nuisance variables by latitude
NV_list_dist <- apply(distinct(loopNV, Species, NV), 1, function(row) {
  plot.NV.dist(row["Species"], row["NV"])
})

#Plot the effect of the nuisance variables on body size, looking for additive effects vs interactions
plot.NV <- function(Spp, DV, NV){  #NV = nuisance variable
  df <- njdf.list.age[[paste0(Spp, ".", DV)]]
  if(Spp == "Nighthawk"){df <- njdf.list.br[[paste0(Spp, ".", DV)]]}
    ggplot(data = df, aes(x = B.Lat, y = .data[[DV]], color = .data[[NV]])) + 
    geom_point() +
    geom_smooth(method = "lm") + 
    labs(x = "Breeding latitude", y = DV, title = paste(Spp, DV))
}

NV_list <- apply(loopNV, 1, function(row) {
  plot.NV(row["Species"], row["DV"], row["NV"])
})
names(NV_list) <- with(loopNV, paste0(Species, DV, NV), sep = "_")

#Generate combined pdf
pdf(file = "Plots/Nuisance_variables2.pdf", width = 8.5, height = 11, bg = "white")
#print(marrangeGrob(grobs = counts_NV, ncol = 1, nrow = 2))
print(marrangeGrob(grobs = NV_list_dist, ncol = 2, nrow = 3, layout_matrix = matrix(1:6, 3, 2, TRUE)))
print(marrangeGrob(grobs = NV_list, ncol = 2, nrow = 2, layout_matrix = matrix(1:4, 2, 2, TRUE)))
dev.off()

# Year effect -------------------------------------------------------------
#Visualize temporal spread of sampling by species
ggplot(data = capri.df.int, aes(x = Year, y = B.Lat)) + 
  geom_jitter(data = filter(capri.df.int, Species == "Nightjar"), alpha = .2,
              width = .4, height = .8, shape = 1, aes(color = Species)) +
  geom_jitter(data = filter(capri.df.int, Species %in% c("Nighthawk","Whip-poor-will")), alpha = .4, 
              width = .4, height = .8, shape = 1, aes(color = Species)) +
  scale_color_manual(values= colorblind_pal()(8)[c(4,6,8)], 
                     guide = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(y = "Breeding latitude") + 
  geom_vline(aes(xintercept= 2010), color="black", linetype="dashed") 
#+ ggtitle("Temporal sampling resolution \n by breeding latitude")
ggsave("Plots/Supplementary/Temporal_sampling_Blat.png", bg = "white")

euni.df <- capri.df.int %>% filter(Species == "Nightjar" & Age != "Unk" & Sex != "U" & Year > 1989) 
euni.df.mass <- euni.df %>% filter(!is.na(tsss.comb))
nrow(euni.df)

#Visualize EUNI
ggplot(data = euni.df, aes(Year, B.Lat)) + 
  geom_hex()
ggplot(data = euni.df, aes(x = Year, y = B.Lat)) + 
  geom_smooth() + 
  geom_point(alpha = .3)
ggplot(data = euni.df, aes(x = Year, y = Mass)) + 
  geom_smooth() + #method = "lm"
  geom_point(alpha = .3)
ggplot(data = euni.df, aes(x = Year, y = Wing.Chord)) +
  geom_point(alpha = .3) + 
  geom_smooth() #method = "lm"



#Confirm effects are important even controlling for lat, sex, & age. Note increase in wing, but not in mass, suggesting shape-shifting through time
euni.yr.m <- lm(Mass ~ scale(Year) + B.Lat + Age + Sex + poly(tsss.comb, 2), data = euni.df.mass)
euni.yr.w <- lm(Wing.Chord ~ scale(Year) + B.Lat + Age + Sex, data = euni.df)
summary(lm(B.Lat ~ scale(Year), data = euni.df)) 

tidy(euni.yr.m, conf.int = TRUE, conf.level = .95) #Create tidy table of parameter estimates

#For now, I think >2010 is best option
capri.df.int <- capri.df.int %>% filter(Year >= 2010 )
euni.df <- euni.df %>% filter(Year >= 2010)

# Banding date effect -----------------------------------------------------
#Ensure banding dates are reasonable, notice all most extreme banding dates are EUNI
capri.df.int %>% filter(Species == "Nightjar" ) %>% 
  dplyr::select(Species, Project, Country, Band.md, W.Lat, Band.Number) %>% 
  arrange(Band.md) %>% 
  slice_head(n = 10)
capri.df.int %>% filter(Species == "Nightjar") %>% 
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
  filter(Band.md > as.POSIXct("2023-05-16") & Band.md < as.POSIXct("2023-08-01"))

#Visualize, little to no relationship with band date at this point
ggplot(data = euni.dfRes, aes(x = Band.md, y = Wing.Chord)) + 
  geom_point(alpha = .3) + 
  geom_smooth()
ggplot(data = euni.dfRes, aes(x = Band.md, y = Mass.comb)) + 
  geom_point(alpha = .3) + 
  geom_smooth() + 
  labs(x = "Band date", y = "Mass")
ggsave("Plots/Supplementary/Mass_band.date.png", bg = "white")

#Visualize EWPW & CONI species data
capri.df.int %>% filter(Species == "Whip-poor-will") %>% 
  ggplot(aes(x = Band.md, y = as.numeric(Fat))) + 
  geom_smooth() + 
  geom_point() + 
  ggtitle("EWPW fattening")
capri.df.int %>% filter(Species == "Nighthawk") %>% 
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
  labs(x = "Banding date", y = "Fat") +
  geom_vline(aes(xintercept = 2023-05-19), color="black", linetype="dashed") 
fat_p2 <- euni.dfRes %>% ggplot(aes(x = Band.md, y = as.numeric(Fat))) + 
  geom_smooth() + 
  geom_point() + 
  labs(x = "Banding date", y = "Fat") + 
  scale_x_datetime(limits = as.POSIXct(c('2023-05-05','2023-08-01'))) #+
  #geom_hline(aes(yintercept = 1.3))

ggpubr::ggarrange(fat_p1, fat_p2,
                  ncol = 2, nrow = 1,
                  align='hv', labels = c("A","B"),
                  legend = "none")
ggsave("Plots/EUNI_fattening.png", bg = "white")

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
capri.df.int %>% filter(Warr.md < as.POSIXct("2024-03-01")) %>% 
  dplyr::select(Species, Bdep.md, Warr.md) %>% group_by(Species) %>% 
  summarize(N = n(), MeanDep = mean(Bdep.md, na.rm = T), sdDep = sd(Bdep.md, na.rm = T), MeanArr = mean(Warr.md, na.rm = T), sdArr = sd(Warr.md, na.rm = T))
#Determine dates when migrants may be present
capri.df.int %>% dplyr::select(Species, Project, B.Lat, Band.md, Bdep.md) %>% 
  mutate(BlatR = round(B.Lat, -1)) %>% 
  group_by(Species, BlatR) %>% 
  summarize(n = n(), MinBand = min(Band.md, na.rm = T), MeanBand = mean(Band.md, na.rm = T), MaxBand = max(Band.md, na.rm = T), minDep = min(Bdep.md, na.rm = T), MeanDep = mean(Bdep.md, na.rm = T), sdDep = sd(Bdep.md, na.rm = T)) %>% mutate(Cutoff = MeanDep - sdDep, TF = Cutoff > MaxBand)
#Not an issue for CONI or EWPW. For EUNI though.. Birds leave as early as August 1 at 60+ degrees North. Average date of departure at 60N is 8-17, so let's subtract 1 sd and get cutoff date of 08-06. Alternatively, just remove everything banded after 08-01
##For spring migration for EUNI, both Evens (2017) and Norevik (2017) agree that birds depart wintering grounds in late February, but arrival date is "early May" for Evens (earliest April 28); and average is May 16 (earliest is May 5) for Norevik (2017), 12 birds in each paper but only 9 made it for spring migration in Evens. Given May 16 is so close to cut-off anyways, probably makes sense to include month of May. 
#For EWPW, cite my paper for spring departure as March 20th and arrival as April 17th, English (2017) departure = March 21, arrival = May 1. These 2 refs are in agreement. Korpach (2022) does not have spring migration dates.
#For CONI, Elly says to use Nov - March for winter and May - August for breeding.

# Multiple Wint locs  -------------------------------------------
#CHECK::Ensure that function is picking the correct row, if the Adult row has W.Lat then we're good
capriBA %>% filter(Band.Number == "135268737") %>% select(Band.Age, Banding.Time, Wing.Chord, Mass, Mass.comb, Wing.comb, tsss.comb, Mass.combBT, W.Lat) 

#Visualize morphological difference between captures
capriBA %>% select(Band.Number, Wing.Chord, Wing.comb) %>% 
  group_by(Band.Number) %>% 
  summarize(dif = Wing.Chord - Wing.comb) %>% 
  arrange(dif)
capriBA %>% select(Band.Number, Band.Age, Mass, Mass.comb, Species) %>% 
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
TF <- capri.fac$Str8line < capri.fac$Mig.dist #13,14,19, 26
capri.fac[ !TF ,c("Str8line", "Mig.dist","Project", "Band.Number")] #Alicia used a different system to calculate migratory distance, so the numbers are slightly different. When she recalculated differences were <20km, so this is not very important on the scale of 4-5k km. 

#Overall, trying to determine whether we can use the actual migration distance or need to use the straight-line distance to account for differences in sampling resolution 
#Migration distance correlations by species
capri.fac %>% group_by(Species) %>% 
  summarize(cor.lat = cor(B.Lat, Mig.dist, use = "complete.obs", method = "pearson"), n = n())
capriA.red2 %>% group_by(Species) %>% 
  summarize(cor = cor(B.Lat, B.Tcv, use = "complete.obs", method = "spearman"), n = n())

capriA.red2 %>% filter(Species == "Whip-poor-will") %>% 
  ggplot(aes(x = B.Lat, y = Mig.dist)) +
  geom_point() + 
  geom_smooth(method = "lm")

capriA.red2 %>% group_by(Species) %>% 
  summarize(cor = cor(B.Lat, B.Tcv, use = "complete.obs"))

##Overall plot. The low correlation in EUNI means that str8line distance is not a good summary in EUNI, and that there may be mistakes in the actual migration distances. 
ggplot(data= capri.fac, aes(x = Str8line, y = Mig.dist)) + 
  geom_smooth(aes(color = Species), method = "lm", se = F, fullrange = T) + 
  geom_point(aes(color = Species), size = 3, position = "jitter", alpha=.3) + 
  geom_abline(slope= 1, linetype = "dashed", color="Black") + #+ xlab("Departure Date") + ylab("Migration\nRate (km / day)") 
  ggtitle("Mig dist calculations all species")

#Inspect EUNI more closely
capri.fac %>% filter(Species == "Nightjar") %>% 
  ggplot(aes(x = Str8line, y = Mig.dist)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = Project)) +
  geom_abline(slope= 1, linetype = "dashed", color="Red") + 
  geom_text_repel(aes(label = TagID)) +
  ggtitle("EUNI Mig dist calculations")
#facet_wrap(~Species)


# >Temporal sampling resolution --------------------------------------------
##Migration distance models. Does the variation in temporal schedule influence migration## distance? Depends how you look at it.. Overall (model w/ RE), temp.res.mig is not significant once you include breeding latitude (which obviously influences migration distance).
#Model all species combined
summary(lmer(Mig.dist ~ Temp.res.mig + B.Lat + (1| Species), data = capri.fac))

#However, looking at each species separately you can see there is high sd in tag fix rate in EWPW, so let's look at EWPW more closely
capri.fac %>% group_by(Species) %>% 
  summarize(mean = mean(Temp.res.mig,na.rm = T),
            sd = sd(Temp.res.mig,na.rm = T),
            max = max(Temp.res.mig,na.rm = T),
            min = min(Temp.res.mig,na.rm = T))
summary(capri.fac$Temp.res.mig) #Ranges from 7 points a day to 1 point every 10 days
capri.fac %>% filter(Project == "Korpach") %>% 
  summarize(mean = mean(Temp.res.mig,na.rm = T),
            sd = sd(Temp.res.mig,na.rm = T))

##IDEAL WOULD BE TO BRING IN ORIGINAL KORPACH DATA HERE AND SOMEHOW REJOIN WITH REST OF THE DATA TO BE ABLE TO CALCULATE BOTH SAMPLING RESOLUTIONS SIMULTANEOUSLY, BUT I SPENT GOOD AMOUNT OF TIME AND THIS ISN'T EASY (TRIED SMARTBIND, JOINS COULD WORK IF YOU FORMAT THE BAND NUMBERS THE SAME, ETC.). EASIEST IS JUST RECREATE CAPRI.DF FILE IN DW_njdf SCRIPT AND NOT REPLACE THE ORIGINAL SAMPLING RESOLUTION DATA.

#By species, only EWPW is significantly impacted by temp sampling resolution.
Spp <- capri.fac %>% filter(Species == "Nightjar") #Nighthawk
#This is complicated by the fact that breeding latitude is highly correlated with sampling resolution
cor(Spp$B.Lat, Spp$Temp.res.mig, use = "complete.obs") 
summary(lm(Mig.dist ~ B.Lat + B.Long + Temp.res.mig, data = Spp))

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

# Allometric scaling -----------------------------------------------------
###SMI -- scaled mass index 
#SMA is standard in studies of allometry. This youtube video explains what's going on very well, and the difference between OLS regression. Go to SMA vs OLS section of video, particularly around 9:05 and see subsequent plots. In SMA trying to minimize the distance to the line in both the X and Y axes (instead of just the Y axis). This different method to estimate residuals results in different lines of best fit. https://www.youtube.com/watch?v=dvXEcYYnask&ab_channel=MethodsinExperimentalEcologyI
library(smatr)
#vignette(package = "smatr") 

Spp <- c("Nighthawk", "Nightjar", "Whip-poor-will")
#Generate Wing mass models
gen.modWM <- function(Spp){
  df <- njdf.br.am %>% filter(Species == Spp)
  #Remove log = "xy" to extract the beta estimate on the unlogged scale (see supp materials). NOTE:: The d
  sma(Wing.comb ~ Mass.comb, data = df, slope.test = 1)  #log = "xy"
}
modWM <- lapply(Spp, gen.modWM)
names(modWM) <- Spp
#Remember implications of sample sizes when interpreting these results.. CONI has relatively small sample size, so it has less power to detect a difference in slope than (especially) EUNI & EWPW. This means the results are more likely to say CONI exhibits isometric scaling.
modWM

#NOTE that ggplot has no default method for plotting SMA (e.g., method = "sma"). Could try the function lmodel2() which can do SMA, but things get confusing when you log() the morphological variables. It is easiest to just adapt the base R plotting
library(scales)
plot_sma <- function(Spp){
  labs <- c(xlab = "Mass [log scale]", ylab = "Wing [log scale]")
  if(Spp == "Nightjar"){
    plot(modWM[[Spp]], main = Spp, xlab = labs[1], ylab = labs[2], col = alpha("blue", 0.1), type = "p")
  } else{
    plot(modWM[[Spp]],  main = Spp, xlab = labs[1], ylab = labs[2], col = alpha("blue", 0.3), type = "p")
  }
  plot(modWM[[Spp]], main = Spp, xlab = labs[1], ylab = labs[2], type = "l", lwd = 2, add = T)
  mtext(paste("b =", round(coef(modWM[[Spp]])[2],2)), side = 3, adj = 1, cex = .7)
}
#Check assumptions 
par(mfrow=c(1,3))
lapply(Spp, function(Spp){
  plot(modWM[[Spp]], which = "res", main = Spp, col = alpha("blue", 0.5))
})


#Print final plots
#pdf(file = "Plots/Allometry_untransformed.pdf", width = 8.5, height = 5, bg = "white")
par(mfrow=c(1,3))
lapply(Spp, plot_sma)
#do.call("grid.arrange", c(scaling.plots[c(1,3,5)], ncol = 2))
dev.off()


##Additional functionality
#Generate logged & non-logged models w/ all species for plotting lines & extracting slope coefficients on the non-logged scale
mod.comb.log <- sma(Wing.comb ~ Mass.comb * Species, data = njdf.br.am, log = "xy", robust = TRUE)
mod.comb <- sma(Wing.comb ~ Mass.comb * Species, data = njdf.br.am, robust = TRUE)

#Generate distinct models for CONI vs EUNI & EWPW for additional plotting flexibility
ewpw.coni.allometry <- njdf.br.am %>% filter(Species != "Nightjar")
euni.allometry <- njdf.br.am %>% filter(Species == "Nightjar")
mod.ewpw.coni <- sma(Wing.comb ~ Mass.comb * Species, data = ewpw.coni.allometry, log = "xy", robust = TRUE)
mod.euni <- sma(Wing.comb ~ Mass.comb , data = euni.allometry, log = "xy", robust = TRUE)
#Would be cool to interpret the intercepts and slopes.. The intercepts don't really make sense (a bird of weight 0 has wing of e.g., 94, 140, 126?), & intercept seems like species that are bigger have larger intercepts. On the other hand, slopes are helpful & suggest that nightjars increase in wing length more slowly than nighthawks or whip-poor-wills 
summary(mod.comb)

#Plot#
#Green = nighthawk, blue = nightjar, pink = whip-poor-will
#Start EWPW & Nighthawk points darker 
plot(mod.ewpw.coni, col = alpha(colorblind_pal()(8)[c(4,8)], .3), 
     type = "p", xlab = "Mass [log scale]", ylab = "Wing [log scale]") 
#Add in EUNI points lighter
plot(mod.euni, col = alpha(colorblind_pal()(8)[c(6)], .1), 
     type = "p", xlab = "Mass [log scale]", ylab = "Wing [log scale]", add = T) #[log scale]
#Plot lines on top
plot(mod.comb.log, type = "l", lwd = 2, col = colorblind_pal()(8)[c(4,6,8)], add = T)

#Create labels for top of plot
top.labs <- paste(paste(row.names(coef(mod.comb.log)), "b =", round(coef(mod.comb.log)$slope, 2)), collapse = "; ")
mtext(top.labs, side = 3, adj = 0, cex = 1.3)

#Tried several other options that I couldn't get to work, see Warton et al. 2012 
#shift = TRUE, multcomp = TRUE, multcompmethod = "adjust"

#Can use to determine if mass and wing scale isometrically. If log = "xy" then I think slope.test should be equal to 0.33. See Warton 2012 paper for more information
sma(Wing.comb ~ Mass.comb, data = df, slope.test = .33) 

#Additional things to play with if you want.. method = "OLS", '* Age' or '+ Age'
###Example from smatr package documentation 
data(leaflife)

# Leapfrog migration ------------------------------------------------------
#Evidence for leapfrog migration in EWPW and EUNI, but not CONI. Has this been reported previously in EUNI? 
ggplot(data= capri.fac, aes(x = B.Lat, y = W.Lat)) + 
  geom_smooth(method = "lm", se = TRUE, fullrange = F, aes(color = Species)) + 
  geom_point(aes(color = Species), size =3, position = "jitter", alpha=.65) 

#How would you test this statistically? 
summary(lmer(W.Lat ~ B.Lat + (1 | Species), data = capri.fac))
summary(lmer(W.Long ~ B.Long + (1 | Species), data = capri.fac))

#save.image(paste0(bs, "Desktop/Grad_School/R_Files/MS/BergAnalysis9.21.22.Rdata"))

