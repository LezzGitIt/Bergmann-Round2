##Bergmann's rule across full-annual-cycle in Caprimulgids##
##Data wrangling 02 -- Merge njdf & environmental covs 
##Data wrangling script to merge nightjar df with environmental covariates df and ultimately create the environment needed to begin the analysis. This will happen in at least 6 steps

#Contents
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
library(gridExtra)
library(mvnormtest)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(purrr::map)

#1. Filter capriBA --------------------------------------------------------

capriBA <- read_xlsx("Intermediate_products/Capri_BA_07.25.25.xlsx", 
                     trim_ws = TRUE) %>%
  mutate(Banding.Time = chron(times = Banding.Time))
nrow(capriBA)

#Remove additional rows in capriBA to create capri.fin (final)
#Remove repeat individuals, selecting adults when possible 
levels(capriBA$Age) #Should get NULL, this should NOT be a factor in order for arrange to work correctly
capriBAnr <- capriBA %>% group_by(Band.Number) %>% #nr = no repeats 
  arrange(is.na(W.Lat), Age) %>% 
  slice_head() %>% 
  ungroup()

## Filter European nightjar based on date
# All years = present year
Year <- format(Sys.Date(), "%Y")

# Res = residents
euniRes <- capriBAnr %>% 
  filter(Species == "Nightjar" & Band.md > as.POSIXct(paste0(.env$Year,"-05-16")) & Band.md < as.POSIXct(paste0(.env$Year, "-08-01")) & is.na(W.Lat))
# Maintain all birds with winter latitude
euniFAC <- capriBAnr %>% filter(Species == "Nightjar" & !is.na(W.Lat))
coni.ewpw <- capriBAnr %>% filter(Species != "Nightjar" & Band.md > as.POSIXct(paste0(Year, "-04-30"))) 

# Bind together data frames
capriRes <- rbind(euniRes, euniFAC, coni.ewpw)
nrow(capriRes)

#Remove last few things#
#Remove Greg's EDB project due to lack of precision in decimals, Unknown sex, and year < 2010.
#Sex is an important covariate in all models, and latitudes are sampled about evenly after 2010
capri.fin <- capriRes %>% filter(Project != "Greg_EDB" & Sex != "U" & Year >= 2010)
nrow(capri.fin)
euni.fin <- capri.fin %>% filter(Species == "Nightjar")

#2. Merge -------------------------------------------------------------------
#Link capri.fin with EnviCovs
#Bring in EnviCovs2
EnviCovs2 <- read_xlsx("Intermediate_products/Envi_Covs_07.25.25.xlsx")

EnviCovs3 <- EnviCovs2 %>% select(-c(Species, Banding.Date, Band.Number))
capriA <- capri.fin %>% left_join(EnviCovs3, by = "uniqID") #capri Analysis

#3. Remove superflous columns -----------------------------------------------
#Can add in 'Post.Hoc' model (Temp + Prec) here 
HypVars <- vector("list", length = 5)
HypVars[[1]] <- data.frame(Vars = c("Lat", "Long", "Elev")) #Geography
HypVars[[2]] <- data.frame(Vars = c("Srad", "Tavg")) #Temp Regulation (TR)
HypVars[[3]] <- data.frame(Vars = c("EVI", "Prec")) #Productivity
HypVars[[4]] <- data.frame(Vars = c("EviCV", "CVprec", "Tcv")) #Seasonality
HypVars[[5]] <- data.frame(Vars = c("Mig.dist"))
#HypVars[[6]] <- data.frame(Vars = c("Prec", "Tavg"))
names(HypVars) <- c("Geo", "TR", "Prod", "Seas", "Mig.Dist") #"Post.Hoc"
HypVarsDf <- bind_rows(HypVars, .id = "Hypothesis")
HypVarsDf <- HypVarsDf %>% mutate(Full = c("Latitude", "Longitude", "Elevation", "Solar radiation", "Temperature", "EVI", "Precipitation", "EVI CV", "Precipitation CV", "Temperature CV", "Migratory distance")) #"Temperature", "Precipitation"

# Analysis file reduced to include only important variables
capriA.red <- capriA %>% select(c("uniqID","Band.Number","Project", "Site.name", "Species","Age","Sex", "tsss.comb", "Wing.comb", "Mass.combBT", "Mass.comb", "Mig.dist", paste0("B.", HypVarsDf$Vars[1:10]), paste0("W.", HypVarsDf$Vars[1:10])))

# Remove a few coastal individuals (from breeding grounds) that have no environmental data 
capri_analysis <- capriA.red %>% drop_na(starts_with("B")) 

# Create full annual cycle dataset by filtering on wintering latitude
capri.fac <- capri_analysis %>% filter(!is.na(W.Lat))
nrow(capri.fac)
# NOTE: Of FAC data there are 9 wing NAs, 2 massBT NAs, 3 Mig dist NAs
capri.fac %>% filter_all(any_vars(is.na(.)))

# 4. Subset breeding dfs list -----------------------------------------------
#Subsetted dfs based on species & DV, these dfs will be used in future analysis script
subset.df <- function(df, spp, var){
  df %>% filter(Species == spp & !is.na(.data[[var]]))
}

loopSppDV <- expand.grid(Species = c("Nighthawk", "Whip-poor-will", "Nightjar"), 
                         DV = c("Wing.comb", "Mass.combBT"),
                         stringsAsFactors = FALSE) %>% 
  arrange(Species, DV)

njdf.list.br.ns <- list()
for(i in 1:nrow(loopSppDV)){
  njdf.list.br.ns[[i]] <- subset.df(capri_analysis, spp = loopSppDV[i, "Species"], var = loopSppDV[i, "DV"])
}
names(njdf.list.br.ns) <- paste0(loopSppDV[,"Species"], "_", loopSppDV[,"DV"])

# 5. Handle age in breeding data ---------------------------------------------
#Remove 'Unk' aged birds from those where Species != Nighthawk, only in breeding b/c Age is not included as a covariate in FAC data 
njdf.l.br.age <- lapply(njdf.list.br.ns, function(x) {
  if (unique(x$Species) != "Nighthawk") {
    x %>% filter(Age != "Unk")
  } else {
    x  # If Species is "CONI", return the original data frame
  }
})
lapply(njdf.l.br.age, nrow)
lapply(njdf.list.br.ns, nrow)

#Create a single data frame with all adult male nightjars (and unknown for CONI)
njdf.br.am <- njdf.l.br.age[str_detect(names(njdf.l.br.age), "Mass.combBT")] %>% 
  bind_rows() %>% 
  filter(Sex == "M")

# 6.  Scale numeric variables  --------------------------------------------
scale.num <- function(df){
  df %>% group_by(Species) %>%
    mutate(across(where(is.numeric), ~ scale(.)[,1])) %>% 
    ungroup()
}

njdf.list.br <- lapply(njdf.l.br.age, scale.num)
lapply(njdf.list.br, nrow)

# 7. Subset winter dfs list --------------------------------------------------
njdf.list.fac <- lapply(njdf.list.br, function(df){
  df %>% filter(!is.na(W.Lat))
})
lapply(njdf.list.fac, nrow)

#Remove one individual CONI that has no fall or spring migration distances
njdf.list.fac <- lapply(njdf.list.fac, function(x) {x %>% filter(!is.na(Mig.dist))})

# 8. Bivariate normality of pred vars ----------------------------------------
#Important to assess bivariate normality b/c we want to determine if it makes sense to use spearman or pearson correlation when examining collinearity of predictor vars. This is also important in model construction (analysis script). 
#If 2 predictor variables are NOT normally distributed it is better to use the spearman correlation, which is a rank order correlation, instead of Pearson which assumes bivariate normality
#H0 is that data is normally distributed, so a p < 0.05 means we reject the null & can say these two vars are not bivariate normal. 
#Just examined Prod & TR hypotheses (2 vars), otherwise need to alter code to look at all pairwise subsets. However it seems clear that these vars generally do not exhibit bivariate normality
map(.x = njdf.l.br.age, .f = \(x) 
    x %>% dplyr::select(matches(paste0("B.", HypVars$Prod[[1]]), ignore.case = F)) %>% #TR
      drop_na() %>%
      t() %>% 
      mshapiro.test())

#Just proof that it is possible to observe bivariate normality (ie, obtain p > 0.05)
test <- data.frame(var1 = rnorm(100, mean = 2, sd = 2), var2 = rnorm(100, mean = -2, sd = .2))
mshapiro.test(t(test))

#Conclusion is that we should use the spearman correlation

# 9. Correlation plots ---------------------------------------------------
namesEC <- names(capri_analysis)
namesECb <- c(namesEC[str_detect(namesEC, "^B\\.")], "Mig.dist")
namesECw <- c(namesEC[str_detect(namesEC, "^W\\.")], "Mig.dist")
#Create df to loop through function
loop_spp_sea <- expand.grid(Species = c("Nighthawk", "Whip-poor-will", "Nightjar"),
                            Season = c("Breeding", "Winter"),
                            stringsAsFactors = FALSE) %>% 
  mutate(Full.name = rep(c("Nighthawk", "Whip-poor-will", "European nightjar"), 2)) %>%
  arrange(Species, Season)

#Function selecting the species and season
select_spp_sea <- function(df, Spp, Season){
  df <- df %>% filter(Species == Spp)
  if(Season == "Winter"){
    df <- df %>% filter(!is.na(W.Lat))
  }
  df %>% select(if (Season == "Breeding") {{ namesECb }} else {{ namesECw }})
}

panel <- paste0(LETTERS[1:6], ")")
Corr.plots <- list()
for(i in 1:nrow(loop_spp_sea)){
  df_spp_sea <- select_spp_sea(capri_analysis, loop_spp_sea[i,1], loop_spp_sea[i,2])
  cor.mat <- round(cor(df_spp_sea, use = "complete.obs", method = "spearman"), 2)
  Corr.plots[[i]] <- GGally::ggcorr(data = NULL, cor_matrix = cor.mat, label = T, 
                                    label_size = 2, label_round = 2, hjust = 0.75, size = 3,
                                    layout.exp = 1.01) +
    labs(title = paste0(loop_spp_sea[i, 3], "\n", loop_spp_sea[i, 2]), #, "Predictors"
         subtitle = panel[i], caption = if (loop_spp_sea[i, 2] == "Winter") {
           paste("N =", nrow(df_spp_sea))
         } else {
           paste0("N = ", nrow(df_spp_sea), "; Migration distance N = ", 
                  nrow(filter(df_spp_sea, !is.na(Mig.dist))))
         }
    )
}

#Print plots to PDF
#pdf(paste0("Plots/Corr_IVs_spearman", format(Sys.Date(), "%m.%d.%y"), ".pdf"))
#print(marrangeGrob(Corr.plots, ncol = 1, nrow = 1))
#dev.off()

# Create a function to generate the plot and save as a .png for easy insertion into Word doc
plot_and_save <- function(index) {
  plot_corr <- ggarrange(Corr.plots[[index]], Corr.plots[[index + 1]], common.legend = TRUE, legend = "right")
  ggsave(paste0("Plots/Corr_IVs/Corr_IVs_", loop_spp_sea[index, 1], ".png"), plot_corr, bg = "white", height = 4)
}

# Use map to apply the function to each value of i
i <- c(1,3,5)
#map(i, plot_and_save)

# 10. Export -------------------------------------------------------------
if(FALSE){
  # Breeding dataset
  capri_analysis %>% select(-c(uniqID, Mass.comb)) %>%
    as.data.frame() %>%
    write.xlsx(paste0("Intermediate_products/capri_analysis", format(Sys.Date(), "%m.%d.%y"), ".xlsx"), row.names = F, showNA = FALSE)
  
  # Full annual cycle dataset
  capri.fac %>% select(-c(uniqID, Mass.comb)) %>%
    as.data.frame() %>%
    write.xlsx(paste0("Intermediate_products/capri_fac", format(Sys.Date(), "%m.%d.%y"), ".xlsx"), row.names = F, showNA = FALSE) 
}

# Export Rdata file
#rm(list= ls()[!(ls() %in% c("njdf.list.br", "njdf.list.br.ns", "njdf.list.fac", "njdf.l.br.age", "njdf.br.am", "capri_analysis", 'capri.fac', "HypVars", "HypVarsDf", "loopSppDV"))])
#save.image(paste0("Rdata/Capri_dfs_", format(Sys.Date(), "%m.%d.%y"), ".Rdata"))
