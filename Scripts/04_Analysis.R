##Analysis -- 
#This script conducts all primary analyses, including: 1) examining nuisance variables, 2) conducting model selection & generating AIC tables, and 3) extracting parameter coefficients from important models

# Libraries & load key dfs ------------------------------------------------
library(AICcmodavg)
library(MuMIn)
library(tidyverse)
library(naniar)
library(readxl)
library(xlsx)
library(stringi)
library(sjPlot)

load("Data/Capri_dfs_12.30.23.Rdata")

# Nuisance variables by spp -----------------------------------------------
#Goal is to examine the impact of nuisance variables Age, sex, and time since sunset (tsss; only on Mass) for each Spp * DV combination. We include tsss as a quadratic variable based on visualization (see Data exploration.R script)
njdf.list.age <- lapply(njdf.list.br, function(x){x[x$Age != "Unk",]})
lapply(njdf.list.age, nrow) #Nighthawk only has 50 individuals that are aged. Age is not in top model for Wing or Mass (w/ the 50 bird df), so let's leave in all individuals and remove Age from model

globNuis <- drgNuis <- candNuis <- aictabNuis <- sumTM <- TM <- resid.plots <- list()
for(i in 1:nrow(loopSppDV)){
  print(paste("i =", i))
  if(loopSppDV[i,1] == "Nightjar" | loopSppDV[i,1] == "Whip-poor-will"){
    df <- njdf.list.age[paste0(loopSppDV[i,1], ".", loopSppDV[i,2])][[1]]
    globNuis[[i]] <- lm(as.formula(paste(loopSppDV[i,2],"~", c("B.Lat + Age + Sex"))), na.action =
                          "na.fail", data = df)
    }
if(loopSppDV[i,1] == "Nighthawk"){#Overwrite Nuis global; include Unk age birds for Nighthawk (via njdf.list)
  df <- njdf.list.br[paste0(loopSppDV[i,1], ".", loopSppDV[i,2])][[1]]
#Remove age from model
  globNuis[[i]] <- lm(as.formula(paste(loopSppDV[i,2],"~", c("B.Lat + Sex"))), na.action = "na.fail",
                      data = df)
  }
if(loopSppDV[i,2] == "Mass.combBT"){
  df <- df %>% filter(!is.na(tsss.comb))
  globNuis[[i]] <- update(globNuis[[i]], ~. + poly(tsss.comb,2)) #Add tsss.comb to model
  }
  drgNuis[[i]] <- dredge(globNuis[[i]], m.lim = c(0,6))
  candNuis[[i]] <- get.models(object = drgNuis[[i]], subset = T)
  NamesNuis <- sapply(candNuis[[i]], function(x){paste(x$call)}[2]) #Why +1?
  aictabNuis[[i]] <- aictab(cand.set = candNuis[[i]], modnames = NamesNuis, sort = TRUE)
  TM[[i]] <- lm(as.formula(aictabNuis[[i]]$Modnames[1]), na.action = "na.fail", data = df) #Top model
  sumTM[[i]] <- summary(TM[[i]])  #Summary of the top model
#resid.plots[[i]] <- plot(TM[[i]], which = 1, main = paste(loopSppDV[i,1], loopSppDV[i,2]))
##Residual plots by Age & Sex, they all look reasonable
  boxplot(resid(TM[[i]]) ~ Sex, data = df, main = paste(loopSppDV[i,1], loopSppDV[i,2],
                                                        "Heterogeneity Sex"), ylab = "Residuals")
  if(loopSppDV[i,1] == "Whip-poor-will" | loopSppDV[i,1] == "Nightjar"){
    boxplot(resid(TM[[i]]) ~ Age, data = df, main = paste(loopSppDV[i,1], loopSppDV[i,2],
                                                          "Heterogeneity Age"), ylab = "Residuals")
  }
  }
names(aictabNuis) <- paste0(loopSppDV[,1], "_", loopSppDV[,2])

#Notice all species are the same, Age (not included in Nighthawk models) + Sex for wing chord, and Age + Sex + tsss for mass
lapply(aictabNuis, slice_head, n = 5)
#Examine impact of age in nighthawks
summary(lm(Mass.combBT ~ B.Lat + Age, njdf.list.age$Nighthawk.Mass.combBT))
names(sumTM) <- paste0(loopSppDV[,1], "_", loopSppDV[,2])

#Create combined data frame of nuisance variable model selection results for sharing
NuisVarsModSelect <- bind_rows(lapply(aictabNuis, slice_head, n = 5))

#Extract the top model of nuisance variables for each Spp*DV combination
Nuis_mods <- lapply(aictabNuis, function(x){x$Modnames[1]})
#Strip the top model down to just the predictor vars without B.Lat or +1
Nuis_mods2 <- lapply(Nuis_mods, str_remove, pattern = '\\+ 1 | \\+ 1|1 \\+ ')
Nuis_mods3 <- lapply(Nuis_mods2, str_remove, pattern = "B.Lat \\+ |B.Lat")
NuisVars <- lapply(Nuis_mods3, function(x){str_split(x, "~ ")[[1]][2]})
NuisVars

#Remove 'Unk' aged birds from those where Species != Nighthawk, only in breeding b/c Age is not included as a covariate in FAC data 
njdf.list.br <- lapply(njdf.list.br, function(x) {
  if (unique(x$Species) != "Nighthawk") {
    x %>% filter(Age != "Unk")
    } else {
      x  # If Species is "CONI", return the original data frame
      }
  })
lapply(njdf.list.br, nrow)

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
                       Hypothesis = c("Geo", "TR", "Prod", "Seas"),
                       stringsAsFactors = FALSE) %>%
  arrange(Species, DV, Hypothesis)

# for loop wrapped in function --------------------------------------------
cormat <- TF.cormat <- globHyp <- drgHyp <- cand.mods <- list()
#Create function to extract the model names 
get.modNames <- function(njdf.list, loop.df, Season = FALSE, FAC = FALSE){
  modNames <- list()
  
  #Run for loop inside of function
  for(i in 1:nrow(loop.df)){ 
    print(paste("i =", i))
    df <- njdf.list[paste0(loop.df[i,1], ".", loop.df[i,2])][[1]]
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
    
    cormat[[i]] <- round(cor(Vars, use = "complete.obs", method = "spearman"), 2)
    TF.cormat[[i]] <- apply(cormat[[i]], 2, function(x){ifelse(x < .7 & x > -.7, T, F)})
    TF.cormat[[i]][upper.tri(TF.cormat[[i]], diag = T)] <- NA
    if(loop.df[i, "Hypothesis"] == "Mig.Dist"){
      TF.cormat[[i]] <- matrix(data = TRUE)
      rownames(TF.cormat[[i]]) <- "Mig.dist"
      colnames(TF.cormat[[i]]) <- "Mig.dist"
    }
    
    predictors <- paste(VarsVect, collapse = "+")
    if(FAC == TRUE){ #Add limit of 3 vars per model due to small sample size
      globHyp[[i]] <- lm(as.formula(paste(loop.df[i,2], "~", predictors)), na.action = "na.fail", data = df)
      drgHyp[[i]] <- dredge(globHyp[[i]], subset = TF.cormat[[i]], evaluate = FALSE, m.lim = c(0,3))
    } else{ #Add fixed variable structure to dredge
      NV <- NuisVars[paste0(loop.df[i,1], "_", loop.df[i,2])][[1]]
      globHyp[[i]] <- lm(as.formula(paste(loop.df[i,2], "~", predictors, "+", NV)), na.action = "na.fail", data = df)
      drgHyp[[i]] <- dredge(globHyp[[i]], subset = TF.cormat[[i]], fixed = as.formula(paste("~", NV)), evaluate = FALSE)
    }
    modNames[[i]] <- lapply(drgHyp[[i]], formula)
  }
  return(modNames)
}

#br = 3 spp * 2 DVs * 4 hypotheses , fac = 3 spp * 2 DVs * 4 hypotheses * 2 seasons + 6 Mig.dist
modNames.br <- get.modNames(njdf.list = njdf.list.br, loop.df = loop.br, Season = FALSE, FAC = FALSE)
modNames.fac <- get.modNames(njdf.list = njdf.list.fac, loop.df = loop.fac, Season = TRUE, FAC = TRUE)
warnings() #NEED TO INVESTIGATE!!
names(modNames.br) <- apply(loop.br[,c("Species", "DV")], 1 , paste , collapse = "_" )


#Function to extract model names and create a dataframe with every model for each Spp*DV combination
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

#Create lists for the extract modNames function
modNamesL <- list(modNames.br = modNames.br, modNames.fac = modNames.fac)
loop.dfL <- list(loop.br = loop.br, loop.fac = loop.fac)
#Extract model names using custom function
modNamesDfs <- Map(extract.ModN, modNamesL, loop.df = loop.dfL)

# Strings -> formulas -> lms -> aictab -----------------------------------------------------
##Mods list with formulas as strings, split by SPP*DV
mods.Lstr_list <- lapply(modNamesDfs, function(df) {
  df %>% group_split(Species, DV)
})

##Function to convert character strings to formulas
convert_to_formula <- function(df, col_name) {
  lapply(df[[col_name]], formula)
}
#Apply the function to each data frame in the list
mods.Lfrm_list <- lapply(mods.Lstr_list, function(mods.Lstr) {
  lapply(mods.Lstr, convert_to_formula, "modName")
})

##Create a list of lms by mapping custom function to the list of formulas
frm_to_lm <- function(list, df){
  lapply(list, lm, na.action = "na.fail", data = df)
}

njdf.list.all <- list(njdf.list.br = njdf.list.br, njdf.list.fac = njdf.list.fac)
mods.Llm_list <- lapply(seq_along(mods.Lfrm_list), function(i) {
  Map(frm_to_lm, mods.Lfrm_list[[i]], df = njdf.list.all[[i]])
})
names(mods.Llm_list) <- c("Breeding", "FAC")

#Create function to generate and name AIC tables
gen_aictab <- function(mods.list, names.list){
  aictab_list <- lapply(seq_along(mods.list), function(i) {
    aictab(cand.set = mods.list[[i]], modnames = names.list[[i]]$modName, sort = TRUE)
  })
  names(aictab_list) <- apply(loopSppDV[,c("Species", "DV")], 1 , paste , collapse = "-" )
  return(aictab_list)
}
aictab_list <- Map(gen_aictab, mods.Llm_list, mods.Lstr_list)
names(aictab_list) <- c("Breeding", "FAC")

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
aictab_list2 <- lapply(aictab_list, gen.cols)
aictab_list2$Breeding <- lapply(aictab_list2$Breeding, function(df) {df %>% select(-Season)})

# Format tables -----------------------------------------------------------
#Rename Modnames so they only display the IVs & are not Season specific
only.IVs <- function(list) {
  lapply(list, function(df) {
    df<- df %>% mutate(Modnames = sapply(strsplit(Modnames, split = "~ "), function(x) {x[2]})) %>% 
      mutate(Modnames = str_remove_all(Modnames, "B.|W.")) %>%
      return(df)
  })
}
aictab_list3 <- lapply(aictab_list2, only.IVs)

#Rename abbreviated variables to their full names using HypVarsDf
full.names <- function(list) {
  lapply(list, function(df) {
    df %>% mutate(Modnames = str_replace_all(Modnames, pattern = setNames(HypVarsDf$Full, HypVarsDf$Vars)),
                  Hypothesis = str_replace_all(Hypothesis, pattern = c("Geo" = "Geography", "Mig.Dist" = "Migratory Distance", "TR" = "Temperature Regulation", "Seas" = "Seasonality", "Prod" = "Productivity")))
  })
}
aictab_list4 <- lapply(aictab_list3, full.names)

#Remove columns and rename the remaining columns 
cols.format <- function(list) {
  lapply(list, function(df) {
    df %>% select(-c(AICc, ModelLik, Cum.Wt)) %>% 
      rename_at(vars(1:5), ~ c("Model", "K", "ΔAICc", "wt", "-2 Log-likelihood"))
  })
}
aictab_list5 <- lapply(aictab_list4, cols.format)

#Filter the AIC tables based on the Null model & ΔAICc values
filt.tab <- function(list, AIC) {
  lapply(list, function(df) {
    df %>% filter(row_number() <= which(Null == "Yes")[1]) %>%
      filter(.data[[AIC]] < 4)
  })
}

#Create reduced list of AIC tables
aictab_list_red <- lapply(aictab_list5, filt.tab, AIC = "ΔAICc")
aictab_list_red <- lapply(aictab_list_red, lapply, select, -Null)

# Export tables -----------------------------------------------------------
lapply(seq_along(aictab_list_red), function(i) {
  tab_dfs(x = aictab_list_red[[i]],  # Use [[i]] to subset the list correctly
          titles = names(aictab_list_red[[i]]),
          footnotes = rep("AICc value of best model =", 6),
          show.footnote = TRUE,
          alternate.rows = TRUE, 
          file = paste0("Tables/AIC_tab_", names(aictab_list_red)[i], 
                        format(Sys.Date(), "%m.%d.%y"), ".doc")
  )
})

# ModAvg ------------------------------------------------------------------
#Create function to generate and name your AIC tables
ext.modavgs <- function(aictab_red, mods.list, var, names.list){
  imp.vars <- lapply(aictab_red, function(df){
    df[,"Modnames"]
  })
  modavg_list <- lapply(seq_along(mods.list), function(i) {
    if(any(str_detect(imp.vars[[i]], var))){
      modavg(cand.set = mods.list[[i]], modnames = names.list[[i]]$modName, parm = var, conf.level = .95)
    }
  })
  names(modavg_list) <- apply(loopSppDV[,c("Species", "DV")], 1 , paste , collapse = "_" )
  return(modavg_list)
}

parm.estB <- parm.estFAC <- list()
HypVarsDf2 <- rbind(HypVarsDf, HypVarsDf) %>% 
  mutate(Season = c(rep("B.", 10), "", rep("W.", 10), ""))
aictab2_red <- lapply(aictab_list2, filt.tab, AIC = "Delta_AICc")

#Each slot in the list is a parameter, and each slot has 6 lists (one for each Spp * DV combo)
for(i in 1:nrow(HypVarsDf2)){ # c(1:6,8:17,19:22)
  print(i)
  if(i %in% 1:11 & HypVarsDf[i,]$Vars != "Mig.dist"){
    parm.estB[[i]] <- ext.modavgs(aictab_red = aictab2_red$Breeding, mods.list = mods.Llm_list$Breeding, var = paste0("B.", HypVarsDf2[i,]$Vars), names.list = mods.Lstr_list$modNames.br)
  } 
  parm.estFAC[[i]] <- ext.modavgs(aictab_red = aictab2_red$FAC, mods.list = mods.Llm_list$FAC, var = paste0(HypVarsDf2[i,]$Season, HypVarsDf2[i,]$Vars), names.list = mods.Lstr_list$modNames.fac)
}
names(parm.estB) <- HypVarsDf2$Full[1:10]

#Create function to generate data frame with parameter estimates
create.parmdf <- function(parm.estL, fac = FALSE){
  parm.dfs <- lapply(parm.estL, bind_rows, .id = "SppDV")
  parm.df <- bind_rows(parm.dfs)
  parm.df2 <- parm.df %>%
    mutate(Species = sapply(strsplit(SppDV, "_"), function(x){x[1]}),
           DV = sapply(strsplit(SppDV, "_"), function(x){x[2]}),
           Important = ifelse(Lower.CL < 0 & Upper.CL > 0, "N", "Y"),
           Season = sapply(strsplit(Parameter, "\\."), function(x){x[1]}),
           Parameter = sapply(strsplit(Parameter, "\\."), function(x){x[2]})) %>% 
    select(-c(SppDV, Mod.avg.table, Conf.level)) %>% 
    distinct() 
  return(parm.df2)
}

parm.df.br <- create.parmdf(parm.estB) %>% mutate(Data.set = "Br")
parm.df.fac <- create.parmdf(parm.estFAC) %>% mutate(Data.set = "Fac")
parm.df <- rbind(parm.df.br, parm.df.fac)

#Format parm.df
parm.df.form <- parm.df %>% rename(Beta = Mod.avg.beta, LCI = Lower.CL, UCI = Upper.CL) %>% 
  select(Species, DV, Data.set, Season, Parameter, Beta:Important) %>%
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

parm.df.form %>% filter(Data.set == "Fac") %>%
  count(Season)
parm.df.form %>% filter(Data.set == "Br") %>%
  count(Hypothesis)

#Export parm.df.form dataframe
parm.df.form %>% as.data.frame() %>%
  write.xlsx("Intermediate_products/parm.df.form.xlsx", row.names = FALSE)
