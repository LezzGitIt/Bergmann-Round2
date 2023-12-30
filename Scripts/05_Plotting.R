##Plotting -- 
#This script conducts all plotting, including the capture locs & migratory paths figure, and the parameter estimates figure across DV & data set (Breeding, FAC) combinations
library(tidyverse)
library(rnaturalearthdata)
library(smoothr)
library(rnaturalearth)
library(sf)
library(ggpubr)

# Mig paths ---------------------------------------------------------------
sf_use_s2(FALSE)
#Download continents where bird data was collected
study.sites <- ne_countries(scale = 50, returnclass = "sf") %>%
  filter(continent %in% c("North America", "South America", "Europe", "Africa") & name_sort != "Russian Federation" & name_sort != "Equatorial Guinea") %>% 
  st_set_precision(1e6) %>%
  st_union() %>%
  st_geometry() 

#Generate sf objects from njdf
capriA.sf <- capriA.red2 %>% st_as_sf(coords = c("B.Long", "B.Lat"), crs = 4326, remove = FALSE)
capri.fac.sf <- capriA.sf %>% filter(!is.na(W.Lat))

#Create linestrings from breeding & winter locations
BrWiConn <- vector("list", length(capri.fac.sf$uniqID))
for(i in 1:length(capri.fac.sf$uniqID)){
  BrWiConn[[i]] <- st_linestring(matrix(c(cbind(capri.fac.sf$B.Long[i], capri.fac.sf$W.Long[i]), cbind(capri.fac.sf$B.Lat[i], capri.fac.sf$W.Lat[i])), nrow = 2, ncol = 2))
}

#Segmentize lines 
BrWiConn <- st_sfc(BrWiConn, crs = 4326) %>% 
  st_segmentize(units::set_units(100, km)) 
BrWiConnDf <- st_sf(uniqID = capri.fac.sf$uniqID, Band.Number = capri.fac.sf$Band.Number, Species = capri.fac.sf$Species, geometry = BrWiConn)

#How does mig distance differ by species?
capri.fac.sf %>% st_drop_geometry() %>% 
  group_by(Species) %>% 
  summarize(Mn.dist = mean(Mig.dist, na.rm = T), N = n())

#Create background plot p, bounding box, & legend
p <- ggplot(data = study.sites) + geom_sf() + 
  theme(legend.position = "none") + 
  theme(axis.title = element_blank()) 
bbox <- st_bbox(BrWiConnDf)
legend <- get_legend(p + geom_sf(data = BrWiConnDf, aes(color = Species)) +
                       scale_color_manual(values= colorblind_pal()(8)[c(4,6,8)]) +
                       theme(legend.position = "top"))

#Create plot with capture locations and migratory paths
BrWiMap <- p + geom_sf(data = BrWiConnDf, alpha = 0.2, aes(color = Species)) + 
  geom_sf(data = capriA.sf, aes(color = Species), alpha = 0.2, shape = 1) + 
  coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) + 
  scale_color_manual(values= colorblind_pal()(8)[c(4,6,8)]) 

ggarrange(BrWiMap, legend.grob = legend)
ggsave("Plots/BrWiBergMap.png", bg = "white")

# Parm estimates plot -----------------------------------------------------

#parm.df.imp <- parm.df.form %>% filter(Important == "Y")
parms.impL <- parm.df.form %>% filter(Important == "Y") %>%
  group_split(DV, Data.set)

#Create function to plot parameters
plot.parms <- function(parms.df){
  ggplot(data = parms.df, aes(x = Parameter, y = Beta, 
                              group = interaction(Species, Season), color = Species, shape = Season)) + 
    geom_point(size = 4, position = position_dodge(width = 0.75)) + 
    geom_errorbar(aes(ymin=LCI, ymax=UCI), 
                  alpha = .6, width=0, size = 1.5, position = position_dodge(width = 0.75)) +
    scale_y_continuous(name = expression(paste("Parameter estimate (", beta, ")"))) + 
    scale_color_manual(values= colorblind_pal()(8)[c(4,6,8)], 
                       breaks = c("Nighthawk", "Nightjar", "Whip-poor-will")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12), 
          axis.text.x = element_text(size = 12, vjust = .58, angle = 60), 
          legend.title = element_text(size=14), legend.text = element_text(size=12), 
          plot.margin = margin(10,10,8,25)) + 
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) + 
    ggtitle(paste(parms.df$DV, parms.df$Data.set))
}

#Apply custom function across list
plots.list <- lapply(parms.impL, plot.parms)
ggarrange(plots.list[[1]], plots.list[[2]], plots.list[[3]], plots.list[[4]], common.legend = TRUE)

ggsave("Plots/Parm_ests_imp_aic_less4_12.29.23.png", bg = "white", width = 9)
