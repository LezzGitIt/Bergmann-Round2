##Bergmann's rule across full-annual-cycle in Caprimulgids##
##Analysis 04 -- 
#This script conducts all analyses for main text of manuscript 

#Contents: 
# 1) Examine nuisance variables
# 2) Conduct model selection, generate & format AIC tables
# 3) Examine residuals of top models to ensure model assumptions are met
# 4) Conduct model averaging and extract parameter coefficients from important models

#NOTE:: There are two places where one must manually change the input depending on if using just adult males or the full annual cycle data set (see manuscript for more details). These are indicated with the following text: 'CHANGE FAC DATA SET HERE IF DOING ADULT MALES'

# Libraries & load key dfs ------------------------------------------------
library(AICcmodavg)
library(lme4)
library(MuMIn)
library(tidyverse)
library(naniar)
library(readxl)
library(stringi)
library(sjPlot)
library(broom)
library(mapview)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(purrr::map)

# Setup env ---------------------------------------------------------------
# Load in .Rdata created in previous script
load(paste0("Rdata/Capri_lists_", format(Sys.Date(), "%m.%d.%y"), ".Rdata"))

# Recreate the Hypothesis variables (HypVars) list & dataframe
HypVars <- vector("list", length = 5)
HypVars[[1]] <- data.frame(Vars = c("Lat", "Long", "Elev")) #Geography
HypVars[[2]] <- data.frame(Vars = c("Srad", "Tavg")) #Temp Regulation (TR)
HypVars[[3]] <- data.frame(Vars = c("EVI", "Prec")) #Productivity
HypVars[[4]] <- data.frame(Vars = c("EviCV", "CVprec", "Tcv")) #Seasonality
HypVars[[5]] <- data.frame(Vars = c("Mig.dist"))
names(HypVars) <- c("Geo", "TR", "Prod", "Seas", "Mig.Dist") #"Post.Hoc"
HypVarsDf <- bind_rows(HypVars, .id = "Hypothesis")
HypVarsDf <- HypVarsDf %>% mutate(Full = c("Latitude", "Longitude", "Elevation", "Solar radiation", "Temperature", "EVI", "Precipitation", "EVI CV", "Precipitation CV", "Temperature CV", "Migratory distance")) 

# Nuis vars mod selection -----------------------------------------------
#Goal is to examine the impact of nuisance variables Age, sex, and time since sunset (tsss; only on Mass) for each Spp * DV combination. We include tsss as a quadratic variable based on visualization & rationale laid out in paper

#Custom function to build global model & dredge all submodels
get.nuis.mods <- function(df, IVs, DV){
  globNV <- lm(as.formula(paste0({{ DV }}, "~" , {{ IVs }})), 
               data = df, na.action = "na.fail")
  drgNuis <- dredge(globNV, m.lim = c(0,6)) #, fixed = "B.Lat"
  candNuis <- get.models(object = drgNuis, subset = T)
  return(candNuis)
}

#Define DVs & IVs, where Age is not included for CONI & and tsss is only included for mass
IVs <- c("B.Lat + Sex + tsss.comb + poly(tsss.comb, 2)", "B.Lat + Sex", #CONI
         rep(c("B.Lat + Age + Sex + tsss.comb + poly(tsss.comb, 2)", "B.Lat + Age + Sex"), 2))
DV <- c(rep(c("Mass.combBT", "Wing.comb"), 3))

#Map custom function over list to generate candidate models
candNuis <- pmap(.l = list(njdf.list.br, IVs, DV), .f = get.nuis.mods)
NamesNuis <- map(candNuis, ~sapply(.x, function(x) paste(x$call)[2]))

#Generate AIC table & remove duplicates (NOTE: AICwt column messed up due to duplicates)
aictabNuis <- map2(.x = candNuis, .y = NamesNuis, 
                   .f = \(x , y) aictab(cand.set = x, modnames = y, sort = TRUE))
aictabNuis2 <- map(aictabNuis, ~distinct(.x, K, AICc, AICcWt, LL, .keep_all = TRUE))

#Examine impact of age in the ~50 nighthawks that were aged
coni.age <- capri_analysis %>% filter(Species == "Nighthawk" & Age != "Unk")
summary(lm(Mass.combBT ~ B.Lat + Age, coni.age))
summary(lm(Wing.comb ~ B.Lat + Age, coni.age))

# NuisVars FAC databases --------------------------------------------------
##Full FAC database (all ages & sexes)
njdf.fac.mass <- njdf.list.fac[str_detect(names(njdf.list.fac), "Mass.combBT")]

#Map custom function over list to generate candidate models
candNuisFAC <- pmap(.l = list(njdf.fac.mass, IVs = "B.Lat + tsss.comb + poly(tsss.comb, 2)", 
                           DV = "Mass.combBT"), .f = get.nuis.mods)
NamesNuisFAC <- map(candNuisFAC, ~sapply(.x, function(x) paste(x$call)[2]))

#Generate AIC table & remove duplicates
aictabNuisFAC <- map2(.x = candNuisFAC, .y = NamesNuisFAC, 
                   .f = \(x , y) aictab(cand.set = x, modnames = y, sort = TRUE))
aictabNuis2FAC <- map(aictabNuisFAC, ~distinct(.x, K, AICc, AICcWt, LL, .keep_all = TRUE))

#Manually record top models for each Spp * DV combo
fac.NV.all <- c("", "", rep(c("+ tsss.comb", ""),2))
names(fac.NV.all) <- names(njdf.list.fac)

##Adult male FAC database##
#Create list with just FAC adult males
njdf.fac.am <- lapply(njdf.list.fac, function(df) {
  if (unique(df$Species) != "Nighthawk") {
    df %>% filter(Age != "Unk" & Age != "Young" & Sex != "F")
  } else {
    df %>% filter(Age != "Young" & Sex != "F") #If Species is "CONI", just remove Young & F
  }
})
lapply(njdf.fac.am, nrow)
njdf.fac.mass.am <- njdf.fac.am[str_detect(names(njdf.fac.am), "Mass.combBT")]

#Map custom function over list to generate candidate models
candNuisFACam <- pmap(.l = list(njdf.fac.mass.am, 
                                IVs = "B.Lat + tsss.comb + poly(tsss.comb, 2)", 
                                DV = "Mass.combBT"), .f = get.nuis.mods)
NamesNuisFACam <- map(candNuisFACam, ~sapply(.x, function(x) paste(x$call)[2]))

#Generate AIC table & remove duplicates
aictabNuisFACam <- map2(.x = candNuisFACam, .y = NamesNuisFACam, 
                        .f = \(x , y) aictab(cand.set = x, modnames = y, sort = TRUE))
aictabNuis2FACam <- map(aictabNuisFACam, ~distinct(.x, K, AICc, AICcWt, LL, .keep_all = TRUE))

fac.NV.am <- c("", "", rep(c("+ tsss.comb", ""), 2))
names(fac.NV.am) <- names(njdf.list.fac)

# Format & export Nuis vars --------------------------------------------------------
#Rename Modnames so they only display the IVs & are not Season specific
only.IVs <- function(list) {
  lapply(list, function(df) {
    df <- df %>% mutate(Modnames = sapply(strsplit(Modnames, split = "~ "), function(x) {x[2]})) %>% 
      mutate(Modnames = str_remove_all(Modnames, "B.|W.")) %>%
      return(df)
  })
}
Nuis_dfs <- only.IVs(aictabNuis2)
Nuis_dfs_fac <- only.IVs(aictabNuis2FAC)

#Extract the minimum AIC to append as footnote to AIC table
extr.minAIC <- function(aictab_list){
  lapply(aictab_list, function(df){
    df %>% slice_min(AICc) %>% 
      pull(AICc) %>% round(2)
  })
}
minAIC_Nuis <- extr.minAIC(Nuis_dfs)
minAIC_Nuis_fac <- extr.minAIC(Nuis_dfs_fac)

#Remove columns and rename the remaining columns 
cols.format <- function(list, cols.rm, reloc = FALSE, vars.reloc, var.after) {
  lapply(list, function(df) {
    df <- df %>% as.data.frame() %>% 
      select(-{{ cols.rm }}) %>% #-Cum.Wt
      rename_at(c(1:5), ~ c("Model", "K", "ΔAICc", "wt", "-2 Log-likelihood"))
    if (reloc) {
      df <- df %>% relocate({{ vars.reloc }}, .after = {{ var.after }})
    }
    return(df)
  })
}
Nuis_dfs2 <- cols.format(Nuis_dfs, cols.rm = c(AICc, ModelLik, Cum.Wt))
Nuis_dfs_fac2 <- cols.format(Nuis_dfs_fac, cols.rm = c(AICc, ModelLik, Cum.Wt))

#Remove Nuisance variables from lists for printing 
Rm.NV <- function(aic_list, string){
  lapply(aic_list, function(df){
    df %>% mutate(Model = str_remove_all(Model, string))
  })
}

Nuis_dfs3 <- Rm.NV(Nuis_dfs2, string = " \\+ 1") # \\+ Lat
Nuis_dfs_fac3 <- Rm.NV(Nuis_dfs_fac2, string = " \\+ 1")

#Replace variables in the Null hypothesis with "1", and then replace "1" with "Null"
replace.null <- function(df, Null_vars) {
  df %>% mutate(Model = ifelse(Model %in% Null_vars, "1", Model)) %>% 
    mutate(Model = ifelse(Model == "1", "Null", Model))
}
Nuis_dfs4 <- lapply(Nuis_dfs3, replace.null, Null_vars = NULL)
Nuis_dfs_fac4 <- lapply(Nuis_dfs_fac3, replace.null, Null_vars = NULL)

#Final formatting
Nuis_dfs5 <- lapply(Nuis_dfs4, function(df){
  df %>% mutate(Model = str_replace_all(
    Model, pattern = c("poly\\(tsss\\.comb, 2\\)" = "Time since sunset²",
                       "tsss\\.comb" = "Time since sunset")
    ))
})
Nuis_dfs_fac5 <- lapply(Nuis_dfs_fac4, function(df){
  df %>% mutate(Model = str_replace_all(
    Model, pattern = c("poly\\(tsss\\.comb, 2\\)" = "Time since sunset²",
                       "tsss\\.comb" = "Time since sunset")
  ))
})

Nuis_dfs_print <- lapply(Nuis_dfs5, function(df){
  df %>% mutate(across(where(is.numeric), \(x) round(x, 3)))
})
Nuis_dfs_fac_print <- lapply(Nuis_dfs_fac5, function(df){
  df %>% mutate(across(where(is.numeric), \(x) round(x, 3)))
})

#Extract the top model of nuisance variables for each Spp*DV combination. Notice all species are the same, Age (not included in Nighthawk models) + Sex for wing chord, and Age + Sex + tsss for mass
Br.NV <- sapply(Nuis_dfs4, function(df){df$Model[1]})
Br.NV <- sapply(Br.NV, str_remove_all, " \\+ Lat|Lat \\+ ")

#Export Nuisance variable table
#Names to print
names.p <- loopSppDV %>% 
  mutate(DV = str_replace_all(DV, 
                              pattern = c("Mass.combBT" = "Mass", "Wing.comb" = "Wing"))) %>% 
  mutate(name = paste(Species, DV)) %>% 
  pull(name)

names.fac.p <- c("Nighthawk Mass", "European nightjar Mass", "Whip-poor-will Mass")

#Breeding dataset
tab_dfs(x = Nuis_dfs_print, 
        titles = names.p,
        footnotes = paste(rep("AICc value of best model =", 6), unlist(minAIC_Nuis)),
        show.footnote = TRUE,
        alternate.rows = TRUE,
        file = paste0("Tables/Nuisance_Vars_AIC", 
                      format(Sys.Date(), "%m.%d.%y"), ".doc"))

#Full annual cycle dataset
tab_dfs(x = Nuis_dfs_fac_print, 
        titles = names.fac.p,
        footnotes = paste(rep("AICc value of best model =", 3), unlist(minAIC_Nuis_fac)),
        show.footnote = TRUE,
        alternate.rows = TRUE,
        file = paste0("Tables/Nuisance_Vars_fac_AIC", 
                      format(Sys.Date(), "%m.%d.%y"), ".doc"))

# Loop dfs ----------------------------------------------------------------
#Dfs to cycle through loop
loop.fac <- expand.grid(Species = c("Nighthawk", "Whip-poor-will", "Nightjar"),
                        DV = c("Wing.comb", "Mass.combBT"),
                        Season = c("Breed", "Winter"),
                        Hypothesis = c("Geo", "TR", "Prod", "Seas"),
                        stringsAsFactors = FALSE)
loop.fac2 <- expand.grid(Species = c("Nighthawk", "Whip-poor-will", "Nightjar"),
                         DV = c("Wing.comb", "Mass.combBT"),
                         Season = "NA",
                         Hypothesis = "Mig.Dist",
                         stringsAsFactors = FALSE)
loop.fac <- rbind(loop.fac, loop.fac2) %>%
  arrange(Species, DV, Hypothesis)

loop.br <- expand.grid(Species = c("Nighthawk", "Whip-poor-will", "Nightjar"),
                       DV = c("Wing.comb", "Mass.combBT"),
                       Hypothesis = c("Geo", "TR", "Prod", "Seas"), #"Post.Hoc"
                       stringsAsFactors = FALSE) %>%
  arrange(Species, DV, Hypothesis)

# for loop wrapped in function --------------------------------------------
#This is definitely not the elegant way to do this! But working with the for loop that I already had..
cormat <- TF.cormat <- globHyp <- drgHyp <- cand.mods <- list()
#Create function to extract the model names 
get.modNames <- function(njdf.list, loop.df, NuisVars, Season = FALSE, FAC = FALSE){
  modNames <- list()
  
  #df <- njdf.list.fac[paste0(loop.fac[i,1], "_", loop.fac[i,2])][[1]]

  #Run for loop inside of function
  for(i in 1:nrow(loop.df)){ 
    print(paste("i =", i))
    df <- njdf.list[paste0(loop.df[i,1], "_", loop.df[i,2])][[1]]
    xcols <- df %>% select(contains("comb"), "Age", "Sex", "Mig.dist") #X for extra
    #CHANGE
    if(Season == TRUE){
      df <- df[, substr(names(df), 1, 1) == substr(loop.df[i,"Season"], 1, 1)] #Select columns for the season in question (winter or breeding)
    } else{
      df <- df %>% select(starts_with("B"))
    }
    df <- cbind(df, xcols)
    Vars <- dplyr::select(df, matches(HypVars[[loop.df[i,"Hypothesis"]]][[1]], ignore.case = F)) 
    VarsVect <- names(Vars)
    
    #Data are not multivariate normal, so better to use spearman correlations
    cormat[[i]] <- round(cor(Vars, use = "complete.obs", method = "spearman"), 2)
    TF.cormat[[i]] <- apply(cormat[[i]], 2, function(x){ifelse(x < .7 & x > -.7, T, F)})
    TF.cormat[[i]][upper.tri(TF.cormat[[i]], diag = T)] <- NA
    if(loop.df[i, "Hypothesis"] == "Mig.Dist"){
      TF.cormat[[i]] <- matrix(data = TRUE)
      rownames(TF.cormat[[i]]) <- "Mig.dist"
      colnames(TF.cormat[[i]]) <- "Mig.dist"
    }
    
    predictors <- paste(VarsVect, collapse = " + ")
    if(FAC == TRUE){ #Add limit of 3 vars per model due to small sample size
      globHyp[[i]] <- lm(as.formula(paste(loop.df[i,2], "~", predictors, NuisVars[paste0(loop.df[i,1], "_", loop.df[i,2])])), na.action = "na.fail", data = df)
      drgHyp[[i]] <- dredge(globHyp[[i]], subset = TF.cormat[[i]], evaluate = FALSE, m.lim = c(0,3))
      #Sloppy code, but idea is only include fixed argument if nchar() in NuisVars > 0
      if(nchar(NuisVars[paste0(loop.df[i,1], "_", loop.df[i,2])]) > 0){
        drgHyp[[i]] <- dredge(globHyp[[i]], subset = TF.cormat[[i]], evaluate = FALSE, m.lim = c(0,3), fixed = as.formula(paste("~", NuisVars[paste0(loop.df[i,1], "_", loop.df[i,2])])))
      }
    } else{ #Add fixed variable structure to dredge
      globHyp[[i]] <- lm(as.formula(paste(loop.br[i,2], "~", predictors, "+", NuisVars[paste0(loop.df[i,1], "_", loop.df[i,2])])), na.action = "na.fail", data = df)
      drgHyp[[i]] <- dredge(globHyp[[i]], subset = TF.cormat[[i]], fixed = as.formula(paste("~", NuisVars[paste0(loop.df[i,1], "_", loop.df[i,2])])), evaluate = FALSE)
    }
    modNames[[i]] <- lapply(drgHyp[[i]], formula)
  }
  return(modNames)
}

#br = 3 spp * 2 DVs * 4 hypotheses = 24
#fac = 3 spp * 2 DVs * 4 hypotheses * 2 seasons + 6 Mig.dist = 54
modNames.br <- get.modNames(njdf.list = njdf.list.br, loop.df = loop.br, NuisVars = Br.NV, Season = FALSE, FAC = FALSE)
modNames.fac <- get.modNames(njdf.list = njdf.list.fac, loop.df = loop.fac, NuisVars = fac.NV.all, Season = TRUE, FAC = TRUE)
modNames.fac.am <- get.modNames(njdf.list = njdf.fac.am, loop.df = loop.fac, NuisVars = fac.NV.am, Season = TRUE, FAC = TRUE)

names(modNames.br) <- apply(loop.br[,c("Species", "DV")], 1 , paste , collapse = "_" )

#Function to extract model names and create a data frame with every model for each Spp*DV combination
extract.ModN <- function(modNames, loop.df){ #Extract modNames
  modNames2 <- lapply(modNames, sapply, deparse, width.cutoff = 100L)
  modNames2 <- lapply(modNames2, str_remove, pattern = '\\+ 1 | \\+ 1|1 \\+ ')
  modNames2 <- lapply(modNames2, function(x){data.frame(modName = x)})
  names(modNames2) <- apply(loop.df[,c("Species", "DV")], 1 , paste , collapse = "_" )
  modNamesDf <- bind_rows(modNames2, .id = "SpeciesDV") 
  modNamesDf2 <- modNamesDf %>% 
    mutate(Species = sapply(str_split(SpeciesDV, "_"), function(x){x[1]}),
           DV = sapply(str_split(SpeciesDV, "_"), function(x){x[2]})) %>% 
    select(-SpeciesDV) %>% 
    relocate(modName, .after = "DV") %>% 
    distinct()
  return(modNamesDf2)
}

#Create lists for the extract modNames function# 
#CHANGE FAC DATA SET HERE IF DOING ADULT MALES
modNamesL <- list(modNames.br = modNames.br, modNames.fac = modNames.fac) #modNames.fac.am

loop.dfL <- list(loop.br = loop.br, loop.fac = loop.fac)
#Extract model names using custom function
modNamesDfs <- Map(extract.ModN, modNamesL, loop.df = loop.dfL)

# Strings -> formulas -> lms -> aictab -----------------------------------------------------
##Mods list with formulas as strings, split by SPP*DV
mods.Lstr_list <- lapply(modNamesDfs, function(df) {
  df_l <- df %>% group_split(Species, DV)
  names(df_l) <- sapply(df_l, function(df){
    unique(paste0(df$Species, "_", df$DV))
  })
  return(df_l)
  #names(df_l) <- paste0(df$Species, "_", df$DV) #Previously had this & it stopped working 
})

##Function to convert character strings to formulas
convert_to_formula <- function(df, col_name) {
  lapply(df[[col_name]], formula)
}

#Apply convert_to_formula function to each data frame in the list to generate formulas
mods.Lfrm_list <- lapply(mods.Lstr_list, function(mods.Lstr) {
  lapply(mods.Lstr, convert_to_formula, "modName")
})

##Take list of formulas and turn them into linear models
frm_to_lm <- function(list, df){
  lapply(list, lm, na.action = "na.fail", data = df)
}

##Create a list of lms by mapping custom frm_to_lm function to the list of formulas
#CHANGE FAC DATA SET HERE IF DOING ADULT MALES
njdf.list.all <- list(njdf.list.br = njdf.list.br, njdf.list.fac = njdf.list.fac) #njdf.fac.am
mods.Llm_list <- lapply(seq_along(mods.Lfrm_list), function(i) {
  Map(frm_to_lm, mods.Lfrm_list[[i]], df = njdf.list.all[[i]])
})
names(mods.Llm_list) <- c("Breeding", "FAC")

##Ensure Variance Inflation Factor reasonable in models##
#Extract max VIF function
ext.max.vif <- function(mod){
  if(length(mod$coefficients) > 2) {
    #vif() gives different output structure depending on whether polynomial of tsss is present
    #if(any(str_detect(names(mod$coefficients), "tsss"))){
     # max(car::vif(mod)[,"GVIF"])
    #} else {
      max(car::vif(mod))
    #}
  }
}


#Apply ext.max.vif function to generate list of lists
mods.vif <- lapply(mods.Llm_list, function(mods.list) {
 lapply(mods.list, lapply, ext.max.vif)
})
vif.unlist <- lapply(mods.vif, lapply, unlist)
#Determine which is model with VIF > 3
lapply(vif.unlist, lapply, max)
lapply(mods.Llm_list$Breeding[[1]], ext.max.vif)
mods.Llm_list$Breeding[[1]][[8]] #Nighthawk Mass


# AIC tables ---------------------------------------------------------------
#Create function to generate and name AIC tables
gen_aictab <- function(mods.list, names.list){
  aictab_list <- lapply(seq_along(mods.list), function(i) {
    aictab(cand.set = mods.list[[i]], modnames = names.list[[i]]$modName, sort = FALSE)
  })
  names(aictab_list) <- apply(loopSppDV[,c("Species", "DV")], 1 , paste , collapse = "-" )
  return(aictab_list)
}
#Use 
aictab_list <- Map(gen_aictab, mods.list = mods.Llm_list, names.list = mods.Lstr_list)
names(aictab_list) <- c("Breeding", "FAC")

gen.R2 <- function(mods.list){
  lapply(seq_along(mods.list), function(i){
    sapply(mods.list[[i]], function(mod){
      round(broom::glance(mod)$r.squared, 2)
      })
  })
}

R2_list <- lapply(mods.Llm_list, gen.R2)

#Add the R2 value to each df & then sort by DeltaAIC. Can also do same thing using Map()
add.R2 <- function(aictab_l, R2_l){
  lapply(seq_along(aictab_l), function(i){
    data.frame(aictab_l[[i]]) %>% mutate(R2 = R2_l[[i]]) %>% 
      arrange(Delta_AICc)
  })
}

aictab_list2 <- Map(add.R2, aictab_l = aictab_list, R2_l = R2_list)
aictab_list2 <- lapply(aictab_list2, setNames, paste0(loopSppDV$Species, "_", loopSppDV$DV))

#Create 'generate columns' function, in this case season, hypothesis & Null columns 
gen.cols <- function(list){
  lapply(list, function(df){
    df <- data.frame(df) %>% 
      mutate(Term = sapply(strsplit(Modnames, split = "+ "), function(x){x[3]})) %>%
      mutate(Season = case_when(
        substr(Term, 1, 1) == "B" ~ "Breed", 
        substr(Term, 1, 1) == "W" ~ "Winter",
        TRUE ~ NA
      ),
      Null = ifelse(K == min(K), "Yes", "No")) %>%
      mutate(Term = str_remove(Term, "B.|W.")) %>% 
      left_join(HypVarsDf[,c("Hypothesis", "Vars")], by = join_by("Term" == "Vars")) %>% 
      select(-Term) %>% 
      mutate(across(where(is.numeric), round, 2))
    return(df)
  })
}

#Apply custom function to aictab_list and remove Season from the Breeding season list
aictab_list3 <- lapply(aictab_list2, gen.cols)
aictab_list3$Breeding <- lapply(aictab_list3$Breeding, function(df) {df %>% select(-Season)})

#Export as it's needed to generate the top geographic model for figure 3
save(aictab_list3, file = "Rdata/aictab_list3.Rdata")

# Format tables -----------------------------------------------------------
#Apply custom only.IVs function to display only the IVs & make them not Season specific
aictab_list4 <- lapply(aictab_list3, only.IVs)

#Rename abbreviated variables to their full names using HypVarsDf
full.names <- function(list) {
  lapply(list, function(df) {
    df %>% mutate(
      Modnames = str_replace_all(Modnames, 
                                 pattern = setNames(HypVarsDf$Full, HypVarsDf$Vars)),
                  Hypothesis = str_replace_all(Hypothesis, pattern = c("Geo" = "Geographic pattern", "Mig.Dist" = "Migratory Distance", "TR" = "Temperature Regulation", "Seas" = "Seasonality", "Prod" = "Productivity")))
  })
}
aictab_list5 <- lapply(aictab_list4, full.names)

#Extract the minimum AIC to include as footnote in formatted table
minAIC_l <- lapply(aictab_list5, extr.minAIC)

#Format columns with custom cols.format function
aictab_list6 <- lapply(aictab_list5, cols.format, cols.rm = c(AICc, ModelLik), 
                       reloc = TRUE, vars.reloc = R2, var.after = Hypothesis)

#Apply custom made replace.null function
NV.l <- list(Breeding = Br.NV, FAC = str_replace_all(fac.NV.all,"\\+ ", ""))
aictab_list7 <- lapply(names(aictab_list6), function(nm){
  map2(.x = aictab_list6[[nm]], .y = NV.l[nm], .f = replace.null)
})
names(aictab_list7) <- names(aictab_list6)

#Filter the AIC tables based on the Null model & ΔAICc values
filt.tab <- function(list, AIC) {
  lapply(list, function(df) {
   df <- df %>% filter(row_number() <= which(Null == "Yes")[1]) %>%
      filter({{ AIC }} < 4) %>% 
      select(-Null)
  })
}

#Create reduced list of AIC tables
aictab_list_red <- lapply(aictab_list7, filt.tab, AIC = ΔAICc)
aictab_list_red_print <- lapply(aictab_list_red, Rm.NV, " \\+ Age| \\+ poly\\(tsss\\.comb, 2\\)| \\+ tsss\\.comb| \\+ Sex")

#Use custom Rm.NV function to remove Nuisance variables from lists for printing 
aictab_list_print <- lapply(aictab_list7, Rm.NV, " \\+ Age| \\+ poly\\(tsss\\.comb, 2\\)| \\+ tsss\\.comb| \\+ Sex")
aictab_list_print <- lapply(aictab_list_print, lapply, select, -Null)

#Select a single model for each hypothesis
aictab_hyp <- lapply(aictab_list_print, map, .f = \(x) {
  x %>% group_by(Hypothesis) %>% 
    slice_min(ΔAICc) %>% 
    arrange(ΔAICc)
})

#Create a combined table from the breeding grounds showing a model for each hypothesis & plus all models within 4 AIC points of the top (unless Null comes first)
aictab_hyp_u4 <- list()
aictab_hyp_u4[[1]] <- map2(.x = aictab_hyp[[1]], .y = aictab_list_red_print[[1]], .f = \(x, y)
                           bind_rows(x, y) %>% 
                             distinct() %>% 
                             arrange(ΔAICc))
aictab_hyp_u4[[2]] <- data.frame(Var1 = "This allows the export.aictab function to work")

# Export tables -----------------------------------------------------------
format.NV <- function(str){
  NV.form <- sapply(str, str_replace_all,
                    pattern = c("poly\\(tsss\\.comb, 2\\)" = "Time since sunset²",
                                "poly\\(tsss\\.comb, 2\\)" = "Time since sunset²",
                                "tsss\\.comb" = "Time since sunset"))
  ifelse(NV.form == "", "None", NV.form)
}

NV.l.unform <- list(Br.NV, fac.NV.all, fac.NV.am)
NV.l <- map(.x = NV.l.unform, .f = format.NV)
NV.l[c(2,3)] <- map(.x = NV.l[c(2,3)], .f = \(x) str_replace(x, "\\+ ", ""))


#Create titles for the tables
titles.fac.all <- data.frame(Breeding = paste(names.p, "–", "Fixed variables:", NV.l[[1]]),
                             FAC = paste(names.p, "–", "Fixed variables:", NV.l[[2]]))
titles.fac.am <- data.frame(Breeding = paste(names.p, "–", "Fixed variables:", NV.l[[1]]),
                             FAC = paste(names.p, "–", "Fixed variables:", NV.l[[3]]))

##Function to export tables 
export.aictab <- function(aictab_list, dir, fac.am = FALSE, i = NULL) {
#Replaces indices with the i argument in the function (if provided), to avoid exporting extra tables which could be confusing
  if(is.null(i)) {
    indices <- seq_along(aictab_list)
  } else {
    indices <- i
  }
  
  lapply(indices, function(idx) {
    col.header <- names(aictab_list[[idx]][[1]])
    col.header[length(col.header)] <- "R²" 
    
    tab_dfs(
      x = aictab_list[[idx]],
      titles = if (fac.am) titles.fac.am[[idx]] else titles.fac.all[[idx]],
      col.header = col.header,
      footnotes = paste(rep("AICc value of best model =", 6), minAIC_l[[idx]]),
      show.footnote = TRUE,
      alternate.rows = TRUE, 
      file = paste0(dir, names(aictab_list)[idx], format(Sys.Date(), "%m.%d.%y"), ".doc")
    )
  })
}

#Call custom function to export the AIC tables
#All models for supplementary materials
export.aictab(aictab_list_print, dir = paste0("Tables/AIC_tab_all_models_"))
#To print FAC models from just adult males 
export.aictab(aictab_list_red_print, dir = "Tables/AIC_tab_am_", fac.am = TRUE, i = 2)

#Just models within 4 AIC points of the top model
export.aictab(aictab_list_red_print, dir = "Tables/AIC_tab_red_", i = 2)
export.aictab(aictab_hyp_u4, dir = "Tables/AIC_tab_hyp_u4_", i = 1)

# Residual plots ----------------------------------------------------------
#Generate strings of important model names
aictab3_red <- lapply(aictab_list3, filt.tab, AIC = Delta_AICc)
imp.strs <- lapply(aictab3_red, lapply, function(df){
  df["Modnames"]
})
#Generate formulas of important model names
imp.frm <- lapply(imp.strs, function(mods.Lstr) {
  lapply(mods.Lstr, convert_to_formula, "Modnames")
})
#Generate lm calls of important models
imp.mods <- lapply(seq_along(imp.frm), function(i) {
  Map(frm_to_lm, imp.frm[[i]], df = njdf.list.all[[i]])
})
names(imp.mods) <- c("Breeding", "FAC")


#Create function to examine residuals 
plot.resid <- function(mods.list, Season = c("Breeding", "FAC")){
  lapply(seq_along(mods.list), function(i) {
    Spp <- names(mods.list)[i]
    formula.L <- lapply(mods.list[[i]], formula)
    name.L <- lapply(formula.L, as.character)
    Map(f = plot, mods.list[[i]], which = 1, 
        main = lapply(name.L, function(x) paste(Season, "\n", Spp, x[1], x[3])))
  })
}

#Print residual plots
plot.resid(imp.mods$Breeding, Season = "Breeding")
plot.resid(imp.mods$FAC, Season = "FAC")



# ModAvg ------------------------------------------------------------------
#Create function to generate and name your AIC tables
ext.modavgs <- function(aictab_red, mods.list, var, names.list, ...){
  #Extract the important models from reduced AIC Table
  imp.vars <- lapply(aictab_red, function(df){
    df[,"Modnames"]
  })
  #Conduct model averaging just using important variables from reduced AIC table
modavg_list <- lapply(seq_along(mods.list), function(i) {
  if(any(str_detect(imp.vars[[i]], var))){
    modavg(cand.set = mods.list[[i]], modnames = names.list[[i]]$modName, 
           parm = var, ...) 
    }
  })
  names(modavg_list) <- apply(loopSppDV[,c("Species", "DV")], 1 , paste , collapse = "_" )
  return(modavg_list)
}

parm.estB <- parm.estFAC <- list()
HypVarsDf2 <- rbind(HypVarsDf, HypVarsDf) %>% 
  mutate(Season = c(rep("B.", 10), "", rep("W.", 10), ""))

#Each slot in the list is a parameter, and each slot has 6 lists (one for each Spp * DV combo)
for(i in 1:nrow(HypVarsDf2)){ 
  print(i)
  if(i %in% 1:11 & HypVarsDf[i,]$Vars != "Mig.dist"){
    parm.estB[[i]] <- ext.modavgs(aictab_red = aictab3_red$Breeding, mods.list = mods.Llm_list$Breeding, var = paste0("B.", HypVarsDf2[i,]$Vars), names.list = mods.Lstr_list$modNames.br, conf.level = .95)
  } 
  parm.estFAC[[i]] <- ext.modavgs(aictab_red = aictab3_red$FAC, mods.list = mods.Llm_list$FAC, var = paste0(HypVarsDf2[i,]$Season, HypVarsDf2[i,]$Vars), names.list = mods.Lstr_list$modNames.fac, conf.level = .95)
}
names(parm.estB) <- HypVarsDf2$Full[1:10]

#Create function to generate data frame with parameter estimates
create.parmdf <- function(parm.estL, fac = FALSE){
  parm.dfs <- lapply(parm.estL, bind_rows, .id = "SppDV")
  parm.df <- bind_rows(parm.dfs)
  parm.df2 <- parm.df %>%
    mutate(Species = sapply(strsplit(SppDV, "_"), function(x){x[1]}),
           DV = sapply(strsplit(SppDV, "_"), function(x){x[2]}),
           Season = sapply(strsplit(Parameter, "\\."), function(x){x[1]}),
           Parameter = sapply(strsplit(Parameter, "\\."), function(x){x[2]})) %>% 
    select(-c(SppDV, Mod.avg.table, Conf.level)) %>% 
    distinct() 
  return(parm.df2)
}

#Generate 85% confidence intervals.. The formula is mean +- z score (1.44 for 85%) * the standard deviation or standard error. I confirmed that this is equivalent to calculating from the modavg() function
gen.CI85 <- function(df, ds){
  df %>% create.parmdf() %>% 
    mutate(Data.set = ds, 
           LCI85 = Mod.avg.beta - 1.44*Uncond.SE,
           UCI85 = Mod.avg.beta + 1.44*Uncond.SE) %>% 
    mutate(Important85 = ifelse(LCI85 < 0 & UCI85 > 0, "N", "Y"))
}

#Create and combine data frames
parm.df.br <- gen.CI85(parm.estB, ds = "Br") 
parm.df.fac <- gen.CI85(parm.estFAC, ds = "Fac")
parm.df <- rbind(parm.df.br, parm.df.fac)

#Format parm.df
parm.df.form <- parm.df %>% 
  rename(Beta = Mod.avg.beta, LCI95 = Lower.CL, UCI95 = Upper.CL) %>% 
  select(Species, DV, Data.set, Season, Parameter, Beta:Important85) %>%
  mutate(across(where(is.numeric), round, 3),
         Season = case_when(
           Season == "B" ~ "Breeding", 
           Season == "W" ~ "Winter",
           Season == "Mig" ~ "Migratory"
         ), 
         DV = case_when(
           DV == "Wing.comb" ~ "Wing", 
           DV == "Mass.combBT" ~ "Mass"
         ),
         Parameter = str_replace_all(Parameter, pattern = setNames(HypVarsDf$Full, HypVarsDf$Vars)),
         Parameter = ifelse(Parameter == "dist", "Migratory distance", Parameter)) %>% 
  left_join(HypVarsDf[,c("Hypothesis", "Full")], by = join_by("Parameter" == "Full")) %>% 
  mutate(Parameter = factor(Parameter, levels = unique(Parameter[order(Hypothesis)])))

#Export parm.df.form dataframe
write.csv(parm.df.form, file = paste0("Intermediates/parm.df.form_", format(Sys.Date(), "%m.%d.%y"), ".csv"), row.names = FALSE)

# Export parm estimates table ---------------------------------------------
parm.tbl.print <- parm.df.form %>% 
  select(-c(Season, Important85, LCI95, UCI95)) %>% 
  arrange(Species, DV, desc(Beta)) %>%
  mutate(Hypothesis = str_replace_all(Hypothesis, pattern = c("Geo" = "Geographic pattern", "Mig.Dist" = "Migratory Distance", "TR" = "Temperature Regulation", "Seas" = "Seasonality", "Prod" = "Productivity")),
         Data.set = case_when(
           Data.set == "Br" ~ "Breeding", 
           Data.set == "Fac" ~ "Full-annual-cycle"
         )) %>%
  rename_at(vars(DV, Beta:UCI85), ~ c("Size measure","Beta estimate", "Standard error", "Lower 85% CI", "Upper 85% CI")) %>% 
  group_by(Species, Data.set)
parm.tbl.print.l <- parm.tbl.print %>% group_split(.keep = FALSE) 
names(parm.tbl.print.l) <- parm.tbl.print %>% group_keys() %>% 
  mutate(name = paste0(Species, " ", Data.set)) %>% 
  pull(name)

#Export table
parm.tbl.print.l %>%
    tab_dfs(titles = names(parm.tbl.print.l),
        alternate.rows = TRUE,
        file = paste0("Tables/Parm_est_", format(Sys.Date(), "%m.%d.%y"), ".doc"))
