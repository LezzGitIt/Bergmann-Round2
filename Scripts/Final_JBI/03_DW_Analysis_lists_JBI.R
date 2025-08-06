## Bergmann's rule across full-annual-cycle in Caprimulgids ##
## Data wrangling 03 -- Create the nightjar dataframe ('njdf') lists that will be used in analysis 
## Data wrangling script that generates lists of data frames based on Species x size (mass or wing) combinations, and breeding vs annual cycle. Ultimately, create the Global Evironment needed to begin the analysis.

## Contents: 
# Create breeding dfs list
# Generates correlation plots from predictor variables
# Scale numeric variables 
# Generate subsetted dfs based on species & DV, and breeding / winter 
# Create correlation plots 
# Save the .Rdata file

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(readxl)
library(ggpubr)

# Load in merged data frames -------------------------------------------------

capri_analysis <- read_excel("Data/Analysis/Capri_analysis07.28.25.xlsx", 
                      trim_ws = TRUE, 
                      guess_max = 3000)

# Subset breeding dfs list ----------------------------------------------
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

# Handle age in breeding data ---------------------------------------------
#Remove 'Unk' aged birds from those where Species != Nighthawk, only in breeding b/c Age is not included as a covariate in FAC data 
njdf.l.br.age <- lapply(njdf.list.br.ns, function(x) {
  if (unique(x$Species) != "Nighthawk") {
    x %>% filter(Age != "Unk")
  } else {
    x  # If Species is "CONI", return the original data frame
  }
})

#Create a single data frame with all adult male nightjars (and unknown for CONI)
njdf.br.am <- njdf.l.br.age[str_detect(names(njdf.l.br.age), "Mass.combBT")] %>% 
  bind_rows() %>% 
  filter(Sex == "M")

# Scale numeric variables  --------------------------------------------
scale.num <- function(df){
  df %>% group_by(Species) %>%
    mutate(across(where(is.numeric), ~ scale(.)[,1])) %>% 
    ungroup()
}

njdf.list.br <- lapply(njdf.l.br.age, scale.num)

# Subset winter dfs list --------------------------------------------------
njdf.list.fac <- lapply(njdf.list.br, function(df){
  df %>% filter(!is.na(W.Lat))
})

#Remove one individual CONI that has no fall or spring migration distances
njdf.list.fac <- lapply(njdf.list.fac, function(x) {x %>% filter(!is.na(Mig.dist))})

# Correlation plots ---------------------------------------------------
namesEC <- names(capri_analysis)
namesECb <- c(namesEC[str_detect(namesEC, "^B\\.")], "Mig.dist")
namesECw <- c(namesEC[str_detect(namesEC, "^W\\.")], "Mig.dist")
# Create df to loop through function
loop_spp_sea <- expand.grid(
  Species = c("Nighthawk", "Whip-poor-will", "Nightjar"),
  Season = c("Breeding", "Winter"), stringsAsFactors = FALSE) %>% 
  mutate(
    Full.name = rep(c("Nighthawk", "Whip-poor-will", "European nightjar"), 2)
  ) %>% arrange(Species, Season)

# Function selecting the species and season
select_spp_sea <- function(df, Spp, Season){
  df <- df %>% filter(Species == Spp)
  if(Season == "Winter"){
    df <- df %>% filter(!is.na(W.Lat))
  }
  df %>% select(if (Season == "Breeding") {{ namesECb }} else {{ namesECw }})
}

#Save plots to Corr.plots object
panel <- paste0(LETTERS[1:6], ")")
Corr.plots <- list()
for(i in 1:nrow(loop_spp_sea)){
  df_spp_sea <- select_spp_sea(capri_analysis, loop_spp_sea[i,1], loop_spp_sea[i,2])
  cor.mat <- round(cor(df_spp_sea, use = "complete.obs", method = "spearman"), 2)
  Corr.plots[[i]] <- GGally::ggcorr(
    data = NULL, cor_matrix = cor.mat, label = T, label_size = 2, label_round = 2, hjust = 0.75, size = 3, layout.exp = 1.01
  ) + labs(title = paste0(loop_spp_sea[i, 3], "\n", loop_spp_sea[i, 2]), #, "Predictors"
           subtitle = panel[i], caption = if (loop_spp_sea[i, 2] == "Winter") {
             paste("N =", nrow(df_spp_sea))
           } else {
             paste0("N = ", nrow(df_spp_sea), "; Migration distance N = ", 
                    nrow(filter(df_spp_sea, !is.na(Mig.dist))))
           }
  )
}

# Custom function to generate the plot and save as a .png for easy insertion into Word doc
plot_and_save <- function(index) {
  plot_corr <- ggarrange(Corr.plots[[index]], Corr.plots[[index + 1]], common.legend = TRUE, legend = "right")
  ggsave(paste0("Plots/Corr_IVs/Corr_IVs_", loop_spp_sea[index, 1], ".png"), plot_corr, bg = "white", height = 4)
}

# Use map to apply the function to each value of i
i <- c(1,3,5)
# Select 1 (Yes) to create new directory 
map(i, plot_and_save)

# Export -------------------------------------------------------------
rm(list= ls()[!(ls() %in% c("njdf.list.br", "njdf.list.br.ns", "njdf.list.fac", "njdf.l.br.age", "njdf.br.am", "capri_analysis", "HypVars", "HypVarsDf", "loopSppDV"))])
save.image(paste0("Rdata/Capri_lists_", format(Sys.Date(), "%m.%d.%y"), ".Rdata"))