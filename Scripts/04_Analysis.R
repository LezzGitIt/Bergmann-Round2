##Analysis -- 
#This script conducts all primary analyses, including: 1) examining nuisance variables, 2) conducting model selection & generating AIC tables, and 3) extracting parameter coefficients from important models

#DIFFERENCES W/ FAC PREVIOUS APPROACH:
#NOTICE SAMPLE SIZES DIFFERENT (ESPECIALLY EUNI)! NO AGE IN EWPW, MODEL SELECTION PROCEDURE (NOT JUST MOST PARSIMONIOUS MODEL FOR EACH SPP*DV COMBO), 3 CONI FEWER DUE TO MIG.DIST, MIG.DIST IS ACTUAL DISTANCE, MASS.COMBBT VS MASS ALSO GOING TO CHANGE DATA SET SLIGHTLY (NOT ONLY WHO'S PRESENT, BUT THE MASS.COMB VAR IS AVERAGED FROM SLIGHTLY DIFFERENT NUMBERS), MAY BE A FEW MINOR ERRORS CORRECTED THROUGH THAT TIME
#EWPW IS THE SAME, CONI IS THE SAME, EUNI 

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
      globHyp[[i]] <- lm(as.formula(paste(loop.df[i,2], "~", predictors)), na.action = "na.fail", data = df, m.lim = c(0,3))
      drgHyp[[i]] <- dredge(globHyp[[i]], subset = TF.cormat[[i]], evaluate = FALSE)
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

mods.Lstr_list$modNames.br[[3]] %>% View()

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

#CHECK::Confirmed that modNames linked up correctly. DELETE
aictab_list$FAC[[6]]
aictab_list2$FAC[[6]]
mods.Llm_list$FAC[[6]][[4]]

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

# Export tables -----------------------------------------------------------
library(sjPlot)
lapply(seq_along(aictab_list_red), function(i) {
  tab_dfs(x = aictab_list_red[[i]],  # Use [[i]] to subset the list correctly
          titles = names(aictab_list_red[[i]]),
          footnotes = rep("AICc value of best model =", 6),
          show.footnote = TRUE,
          alternate.rows = TRUE, 
          file = paste0("Tables/AIC_tab_", names(aictab_list_red)[i], ".doc")
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


imp.vars <- lapply(aictab2_red$Breeding, function(df){
  df[,"Modnames"]
})

i <- 2
var <- "Long"

ext.modavgs(aictab_red = aictab2_red$Breeding, mods.list = mods.Llm_list$Breeding, var = paste0("B.", HypVarsDf2[i,]$Vars), names.list = mods.Lstr_list$modNames.br)
any(str_detect(imp.br.vars$`Whip-poor-will-Wing.comb`, HypVarsDf2[1,"Vars"]))

#Each slot in the list is a parameter, and each slot has 6 lists (one for each Spp * DV combo)
for(i in 1:nrow(HypVarsDf2)){ # c(1:6,8:17,19:22)
  print(i)
  if(i %in% 1:11 & HypVarsDf[i,]$Vars != "Mig.dist"){
    parm.estB[[i]] <- ext.modavgs(aictab_red = aictab2_red$Breeding, mods.list = mods.Llm_list$Breeding, var = paste0("B.", HypVarsDf2[i,]$Vars), names.list = mods.Lstr_list$modNames.br)
  } 
  parm.estFAC[[i]] <- ext.modavgs(aictab_red = aictab2_red$FAC, mods.list = mods.Llm_list$FAC, var = paste0(HypVarsDf2[i,]$Season, HypVarsDf2[i,]$Vars), names.list = mods.Lstr_list$modNames.fac)
}
names(parm.estB) <- HypVarsDf2$Full[1:10]

#DELETE
parm.estFAC$Longitude$Nighthawk_Mass.combBT
aictab_list$Breeding$`Nightjar-Mass.combBT`
aictab_list6$Breeding$`Nighthawk-Mass.combBT`

aictab_list6$FAC
aictab2_red$FAC

summary(lm(Mass.combBT ~ B.Elev + B.Lat + B.Long + poly(tsss.comb, 2) + Sex, data = njdf.list.br$Nighthawk.Mass.combBT))


#Create the data frame with parameter estimates.
create.parmdf <- function(parm.estL, fac = FALSE){
  parm.dfs <- lapply(parm.estL, bind_rows, .id = "SppDV")
  parm.df <- bind_rows(parm.dfs)
  parm.df2 <- parm.df %>% mutate(Species = sapply(strsplit(SppDV, "_"), function(x){x[1]}),
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

# Find error --------------------------------------------------------------

##CHECK THESE PLOTS AGAINST 1) SOME LMS FROM THE RAW DATA, 2) THE MODEL SELECTION TABLES. ULTIMATELY CAN GO BACK & CHECK THINGS BIT BY BIT UNTIL FIND THE ERROR, AND MAY HAVE TO ASK ELLY & ALICIA TO HELP SCOUR THE CODE TO FIND THE ERROR 
mods.Llm_list 
NuisVars

#This is model that it's basing Long parameter estimate off of (estimate = +1.01)
summary(lm(Mass.combBT ~ B.Long + B.Elev + Age + poly(tsss.comb, 2) + Sex, data = njdf.list.br$Nightjar.Mass.combBT))
#Note reversal of direction when using B.Lat instead of B.Elev
summary(lm(Mass.combBT ~ B.Long + B.Lat + Age + poly(tsss.comb, 2) + Sex, data = njdf.list.br$Nightjar.Mass.combBT))


# DELETE  -----------------------------------------------------------------
for(i in 1:nrow(loop)){ #
  print(paste("i =", i))
  df <- njdf.list.br[paste0(loop[i,1], ".", loop[i,2])][[1]]
  xcols <- df %>% select(contains("comb"), "Age", "Sex", "Mig.dist") #X for extra
  #CHANGE
  df <- df[, substr(names(df), 1, 1) == substr(loop[i,3], 1, 1)] #Select columns for the season in question (winter or breeding)
  df <- cbind(df, xcols)
  Vars <- dplyr::select(df, matches(HypVars[[loop[i,4]]][[1]], ignore.case = F)) 
  VarsVect <- names(Vars)
  
  cormat[[i]] <- round(cor(Vars, use = "complete.obs", method = "spearman"), 2)
  TF.cormat[[i]] <- apply(cormat[[i]], 2, function(x){ifelse(x < .7 & x > -.7, T, F)})
  TF.cormat[[i]][upper.tri(TF.cormat[[i]], diag = T)] <- NA
  if(loop[i,4] == "Mig.Dist"){
    TF.cormat[[i]] <- matrix(data = TRUE)
    rownames(TF.cormat[[i]]) <- "Mig.dist"
    colnames(TF.cormat[[i]]) <- "Mig.dist"
  }
  
  predictors <- paste(VarsVect, collapse = "+")
  NV <- NuisVars[paste0(loop[i,1], "_", loop[i,2])][[1]]
  globHyp[[i]] <- lm(as.formula(paste(loop[i,2], "~", predictors, "+", NV)), na.action = "na.fail", data = df)
  #Removed m.lim = c(0,3) as with breeding data this won't be an issue
  drgHyp[[i]] <- dredge(globHyp[[i]], subset = TF.cormat[[i]], fixed = as.formula(paste("~", NV)), evaluate = FALSE)
  modNames[[i]] <- lapply(cand.mods[[i]], formula)
}
