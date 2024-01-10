##Analysis -- 
#This script conducts all primary analyses, including: 1) examining nuisance variables, 2) conducting model selection & generating AIC tables, and 3) extracting parameter coefficients from important models

#TO DO: 
#Look more closely at residual plots and determine if you should transform any of the IVs? 
#Anything counting the number of times a hypothesis appears is probably not really valid 
#Want to do fewest parms within delta AIC < 2? Don't think it is necessary
#Likelihood ratio test to examine goodness of fit?
#Add age in for EWPW FAC models? Or just run adult males? 
#Format parm estimates plot , consider FAC to supp material
#Plot IVs & DVs of important models (Or at least latitude!) using predict? 
#Change from mass.comb to mass for FAC models? Or include tsss as predictor in FAC models as well?
#Look back at draft w/ comments & paste over a few relevant comments
#Look at feedback from the UBC Biol 347 class
#See 3 specific predictions at end of the intro to see if these can be easily tested 

#USE THE DRY ERASE BOARD, HAVE YOUR LITTLE NOTEBOOK WITH YOU, DRAW IT OUT, TRY TO GET OUT OF JUST BANGING YOUR HEAD AGAINST THE WALL IN R. IF YOU FEEL YOURSELF SAYING WOW IM DOING A LOT OF WORK TO DO SOMETHING THAT SEEMS PRETTY EASY.. PROBABLY TIME TO TAKE A STEP BACK. 

gg_plot(df = njdf.list.fac$Nightjar.Wing.comb, x = "B.Prec", y = "Wing.comb", method = "lm") 
hist(njdf.list.fac$Nightjar.Wing.comb$B.Prec)

#Examining IVs
lapply(njdf.list.br, function(df){
  df <- df %>% select(where(is.numeric))
  lapply(df, hist)
})

# Libraries & load key dfs ------------------------------------------------
library(AICcmodavg)
library(MuMIn)
library(tidyverse)
library(naniar)
library(readxl)
library(xlsx)
library(stringi)
library(sjPlot)

load("Rdata/Capri_dfs_12.30.23.Rdata")

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
    
    #PEARSON CORR COEFFICIENT
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
#Apply the function to each data frame in the list to generate formulas
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

#Ensure Variance Inflation Factor reasonable in models
#Extract max VIF function
ext.max.vif <- function(mod){
  if(length(mod$coefficients) > 2) {
    #vif() gives different output structure depending on whether polynomial of tsss is present
    if(any(str_detect(names(mod$coefficients), "tsss"))){
      max(car::vif(mod)[,"GVIF"])
    } else {
      max(car::vif(mod))
    }
  }
}

#Apply function to generate list of lists
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

#Add the R2 value to each df & then sort by DeltaAIC
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

# Format tables -----------------------------------------------------------
#Rename Modnames so they only display the IVs & are not Season specific
only.IVs <- function(list) {
  lapply(list, function(df) {
    df<- df %>% mutate(Modnames = sapply(strsplit(Modnames, split = "~ "), function(x) {x[2]})) %>% 
      mutate(Modnames = str_remove_all(Modnames, "B.|W.")) %>%
      return(df)
  })
}
aictab_list4 <- lapply(aictab_list3, only.IVs)

#Rename abbreviated variables to their full names using HypVarsDf
full.names <- function(list) {
  lapply(list, function(df) {
    df %>% mutate(Modnames = str_replace_all(Modnames, 
                                             pattern = setNames(HypVarsDf$Full, HypVarsDf$Vars)),
                  Hypothesis = str_replace_all(Hypothesis, pattern = c("Geo" = "Geography", "Mig.Dist" = "Migratory Distance", "TR" = "Temperature Regulation", "Seas" = "Seasonality", "Prod" = "Productivity")))
  })
}
aictab_list5 <- lapply(aictab_list4, full.names)

#Remove columns and rename the remaining columns 
cols.format <- function(list) {
  lapply(list, function(df) {
    df %>% select(-c(AICc, ModelLik)) %>% #-Cum.Wt
      rename_at(vars(1:5), ~ c("Model", "K", "ΔAICc", "wt", "-2 Log-likelihood")) %>% 
      relocate(R2, .after = Hypothesis)
  })
}
aictab_list6 <- lapply(aictab_list5, cols.format)

#DELETE, Change to list6 if need to continue 
aictab_list5$Breeding$`Nightjar-Mass.combBT`
aictab_list_red

summary(lm(Mass.combBT ~ B.Lat + B.Elev + Age + poly(tsss.comb, 2) + Sex, data = njdf.list.br$Nightjar.Mass.combBT))
summary(lm(Mass.combBT ~ B.EviCV + B.Tcv + Age + poly(tsss.comb, 2) + Sex, data = njdf.list.br$Nightjar.Mass.combBT))
summary(lm(Mass.combBT ~ B.Srad + B.Tavg + Age + poly(tsss.comb, 2) + Sex, data = njdf.list.br$Nightjar.Mass.combBT))

num.hyp <- lapply(aictab_list2, lapply, function(x){
  x %>% count(Hypothesis)
})
bind_rows(num.hyp, .id = "Spp_DV")
#THROUGH HERE

#Filter the AIC tables based on the Null model & ΔAICc values
filt.tab <- function(list, AIC) {
  lapply(list, function(df) {
    df %>% filter(row_number() <= which(Null == "Yes")[1]) %>%
      filter(.data[[AIC]] < 4)
  })
}

#Create reduced list of AIC tables
aictab_list_red <- lapply(aictab_list6, filt.tab, AIC = "ΔAICc")
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

# Residual plots ----------------------------------------------------------
#Generate strings of important model names
aictab3_red <- lapply(aictab_list3, filt.tab, AIC = "Delta_AICc")
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

#Examine summary() of linear models, notice Latitude usally has highest t values
lapply(imp.mods, lapply, lapply, summary)$Breeding

#This is basically a for loop as far as I can tell
plot.resid <- function(mods.list, Season = c("Breeding", "FAC")){
  lapply(seq_along(mods.list), function(i) {
    Spp <- names(mods.list)[i]
    formula.L <- lapply(mods.list[[i]], formula)
    name.L <- lapply(formula.L, as.character)
    Map(f = plot, mods.list[[i]], which = 1, 
        main = lapply(name.L, function(x) paste(Season, "\n", Spp, x[1], x[3])))
  })
}

#Figure out how to add names to this later, may be problem with CONI Wing models
pdf(file = "Plots/Model_assumption_plots/plot1_fitted_resid2.pdf", 
    width = 8.5, height = 5, bg = "white")
plot.resid(imp.mods$Breeding, Season = "Breeding")
plot.resid(imp.mods$FAC, Season = "FAC")
dev.off()


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
           Important95 = ifelse(Lower.CL < 0 & Upper.CL > 0, "N", "Y"),
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
           UCI85 = Mod.avg.beta + 1.44*Uncond.SE)
}

#Create and combine data frames
parm.df.br <- gen.CI85(parm.estB, ds = "Br") 
parm.df.fac <- gen.CI85(parm.estFAC, ds = "Fac")
parm.df <- rbind(parm.df.br, parm.df.fac)

#Format parm.df
parm.df.form <- parm.df %>% rename(Beta = Mod.avg.beta, LCI95 = Lower.CL, UCI95 = Upper.CL) %>% 
  select(Species, DV, Data.set, Season, Parameter, Beta:UCI85) %>%
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


#Examine the number of models compared in each dataset, and the number of breeding vs. winter models in the FAC dataset
lapply(aictab_list5, lapply, nrow)
40 + 42 + 38 + 38 + 36 + 36
21 + 21 + 17 + 17 + 20 + 20
#Wintering with 10 more models overall
lapply(aictab_list5$FAC,function(x){table(x$Season)})
2 + 2 + 4 + 4 + 0 - 2 

#Examine support for Breeding vs winter hypotheses
aictab_list_red$FAC
parm.df.form %>% filter(Data.set == "Fac") %>%
  count(Season)
#Examine support for different hypotheses on the breeding grounds
parm.df.form %>% filter(Data.set == "Br") %>%
  count(Hypothesis)
parm.df.form %>% filter(Data.set == "Br" & Important == "Y") %>% 
  arrange(desc(abs(Beta))) #%>% filter(Parameter == "Elevation")

#Export parm.df.form dataframe
parm.df.form %>% as.data.frame() %>%
  write.xlsx("Intermediate_products/parm.df.form.xlsx", row.names = FALSE)
