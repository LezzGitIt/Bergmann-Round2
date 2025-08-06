## Bergmann's rule across full-annual-cycle in Caprimulgids##
## Plotting 05  
# This script conducts all plotting for the main text of the manuscript, outside of Fig 1 (produced in Adobe Illustrator) 

#Contents
#Fig2: Joint data figure including the capture locs & migratory paths figure
#Fig3: Isoclines (map) predicting size based on top geography models
#Fig4: Parameter estimates from breeding & FAC data

# Load libraries & data ---------------------------------------------------
# Load libraries
library(tidyverse)
library(readxl)
library(rnaturalearthdata)
library(rnaturalearth)
library(rnaturalearthhires)
library(smoothr)
library(sf)
library(ggpubr)
library(ggrepel)
library(ggthemes)
library(viridis)
library(terra)
library(metR)
library(raster)
library(cowplot)
library(conflicted)
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(purrr::map)

# Load data
load(paste0("Rdata/Capri_lists_", format(Sys.Date(), "%m.%d.%y"), ".Rdata"))
capri_fac <- read_excel(paste0("Data/Analysis/Capri_fac07.28.25", ".xlsx"), trim_ws = TRUE)

#Fig2: Breeding map ----------------------------------------------------------
#Function to round a number to the nearest number that the user specifies (with base)
mround <- function(x,base){
  base*round(x/base)
}

##Summarize breeding locations spatially##
capriA.sf <- capri_analysis %>% 
  mutate(B.Lat = mround(B.Lat, .25), B.Long = mround(B.Long, .5)) %>%
  count(Species, B.Lat, B.Long, sort = T) %>% 
  st_as_sf(coords = c("B.Long", "B.Lat"), crs = 4326, remove = FALSE)

#Create base map to plot points on top of#
sf_use_s2(FALSE)
world <- map_data("world")
lake <- map_data("lakes")
states <- map_data("state")
provinces <- ne_states(returnclass = "sf", country = "canada")

#Set theme
map.theme <- theme(text=element_text(size=12, family="Arial"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  theme(plot.margin = margin(t = 0, r = .1, b = 0, l = .1))

p <- ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group), 
                             fill = "gray80", colour = "black", size = 0.3) +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), 
               fill = "gray80", colour = "black", size = 0.2) + 
  geom_sf(data = provinces, fill = "gray80", colour = "black", size = 0.3) +
  geom_polygon(data = lake, aes(x = long, y = lat, group = group), 
               fill = "white", colour = "black", size = 0.3) +
  map.theme

#Create df with species short & species long names 
Spp <- data.frame(Short = c("Nighthawk", "Nightjar", "Whip-poor-will"), 
                  Long = c("Nighthawk", "European \nnightjar", "Whip-poor-will"))

#Create plot with capture locations 
plot.map <- function(Species1, Species2 = NULL, Species3 = NULL, leg.pos = "none", ...){
  capriA.sf.red <- capriA.sf %>% filter(Species %in% c(Species1, Species2, Species3))
  bbox <- st_bbox(capriA.sf.red)
  #colors.i <- c(4,6,8)[seq_along(c(Species1, Species2, Species3))]
  p + geom_sf(data =  capriA.sf.red, aes(color = Species, size = n), alpha = 1, shape = 1) + 
    coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax), expand = TRUE) + 
    scale_color_manual(values= colorblind_pal()(8)[c(4,6,8)], 
                       breaks = Spp$Short) +
    scale_size_continuous(limits = c(1, 500), range = c(2,7), 
                          breaks = c(20, 50, 100, 200, 350, 500), 
                          name = "Sample size") +
    scale_x_continuous(breaks = seq(mround(max(capriA.sf$B.Long), 5), 
                                    mround(min(capriA.sf$B.Long), 5), by = -15)) + 
    theme(legend.position = leg.pos, 
          plot.margin = margin(0,0,0,0)) + 
    guides(...)
}

#Alternatively, if doing each species individually can use lapply
maps <- lapply(Spp$Short, function(species) {
  plot.map(species, legend.position = "none")
})

#Create function to plot boxplots
plot.boxplot <- function(df, Species1, var, flip = FALSE){
  df %>% filter(Species %in% c(Species1)) %>% 
    ggplot(aes(x = Species, y = .data[[var]], color = Species)) + 
    geom_boxplot(outlier.alpha = .6) +
    #geom_point(position=position_jitterdodge(jitter.width = .4), alpha = 0.2) +
    scale_color_manual(values= colorblind_pal()(8)[c(4,6,8)], breaks = Spp$Short) +
    theme_void() +
    theme(legend.position = "none") +
    labs(x = NULL, y = NULL) +
    if(flip == TRUE) coord_flip()
}

#Latitude boxplot
boxp.Lat <- lapply(Spp$Short, function(species) {
  plot.boxplot(df = capri_analysis, Species1 = species, var = "B.Lat") +
    theme(plot.margin = margin(0,0,0,0))
})

#Longitude boxplot
boxp.Long <- lapply(Spp$Short, function(species) {
  plot.boxplot(df = capri_analysis, Species1 = species, var = "B.Long", flip = TRUE) +
    theme(plot.margin = margin(0,0,0,0))
})

#Additional NULLs allow for control of spacing between boxplots & map
maps.boxp <- lapply(seq_along(Spp$Short), function(i){
  #Notice entire row & column of NULLs to control spacing 
  ggarrange(boxp.Long[[i]], NULL, NULL, 
            NULL, NULL, NULL,
            maps[[i]], NULL, boxp.Lat[[i]], 
            align = "hv", widths = c(4, -.24, 1), heights = c(1, -.24, 4),
            ncol = 3, nrow = 3)
})

#Print to final plot
ggarrange(maps.boxp[[1]], maps.boxp[[2]], maps.boxp[[3]], nrow = 1) 
ggsave(paste0("Plots/Data_fig2/Br_Berg_Boxplots_Map_", format(Sys.Date(), "%m.%d.%y"), ".png"), bg = "white", height = 2)

# Fig2: Mig paths ---------------------------------------------------------------
sf_use_s2(FALSE)
Continents <- ne_countries(scale = 50, returnclass = "sf") %>%
  filter(continent %in% c("North America", "South America","Europe", "Africa") & name_sort != "Russian Federation" & name_sort != "Equatorial Guinea") %>% #"North America", "South America"
  st_set_precision(1e6) %>%
  st_union() %>%
  st_geometry()

Lakes <- ne_download(scale = 50, type = 'lakes', 
                     category = 'physical', returnclass = "sf") %>%
  st_set_precision(1e6) %>%
  st_union() %>%
  st_geometry()

migDplot <- list()
for(i in 1:nrow(Spp)){
  #Use stack to create data frame 2x as long where breeding & winter Lats and Longs are split into two rows
  stack.df.spp <- data.frame(capri_fac[,c("Species","Band.Number", "Mig.dist")], stack(capri_fac[,c("B.Lat", "W.Lat")]), stack(capri_fac[,c("B.Long", "W.Long")])) %>%
    rename(Lat = values, Long = values.1) %>% 
    st_as_sf(coords = c("Long", "Lat"), crs = 4326, remove = FALSE) %>%
    filter(Species %in% c(Spp$Short[i])) 
  bbox <- st_bbox(stack.df.spp)
  Lines.sf <- stack.df.spp %>% group_by(Band.Number) %>% 
    summarize(do_union = F) %>% 
    st_cast("LINESTRING")
  
  migDplot[[i]] <- ggplot(data = Continents) + geom_sf() + 
    geom_sf(data = Lakes, fill = "white") +
    geom_sf(data = stack.df.spp, alpha = .5, aes(color = Mig.dist), shape = 1, size = 3) +
    geom_sf(data = Lines.sf, alpha = .03) + 
    coord_sf(xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
    geom_hline(yintercept = c(23.5, -23.5), linetype = "dashed", alpha = .3) + 
    geom_hline(yintercept = 0, alpha = .5) +
    scale_color_viridis_c(option = "magma", 
                        breaks=c(min(stack.df.spp$Mig.dist) + 500, 
                                 max(stack.df.spp$Mig.dist) - 500), 
                        labels=c("Short", "Long")) +
    theme_void() #+ 
    #labs(title = Spp$Long[i])
}
fac.maps <- ggarrange(migDplot[[1]],migDplot[[2]],migDplot[[3]], nrow = 1, ncol = 3, 
                      legend = "none", widths = c(1, 1, 1)) #common.legend = TRUE

#Create density plot 
dens.migD <- capri_fac %>% ggplot(aes(x = Mig.dist, fill = Species)) + 
  geom_density(alpha=0.2, aes(y = after_stat(scaled))) + 
  scale_fill_manual(values= colorblind_pal()(8)[c(4,6,8)], labels = Spp$Long) +
  labs(x = "Migratory distance", y = "Density")
dens_leg <- ggpubr::get_legend(dens.migD)
dens.migD <- dens.migD + theme(legend.position = "none")
#Control width of legend
dens.migD.f <- ggarrange(NULL, dens.migD, dens_leg, widths = c(.1, 3, 1), nrow = 1) #f for formatted

#Produce & save plot
ggarrange(fac.maps, dens.migD.f, nrow = 2, heights = c(2, 1))
ggsave("Plots/Data_fig2/Mig.dist.png", bg = "white")
#Manually arrange in Ppt for final figure

#Fig3: Map - Isoclines predicting size -----------------------------------------------------
#Extract the top 'Geographic pattern' model for each hypothesis & run a linear model that will be used for predictions#
load( "Rdata/aictab_list3.Rdata")
top.lms <- map2(.x = aictab_list3$Breeding[c(1,3,5)], .y = njdf.l.br.age[c(1,3,5)], 
                .f = \(aic.tab, df){
                  top.geo <- aic.tab %>% filter(Hypothesis == "Geo") %>% 
                    slice_head() %>% 
                    pull(Modnames) 
                  lm(as.formula(top.geo), data = df)
                })

#Rasterize ranges#
ranges <- read_sf("Spatial_files/Ranges/BOTW_9-3_Caprimulgids.shp")
sci.names <- c("Chordeiles minor", "Caprimulgus europaeus", "Caprimulgus vociferus")

#Create SVs
SVs <- map(.x = sci.names,.f = \(spp){
  #Create SpatVectors
  SV <- ranges %>% #Create SpatVectors
    dplyr::filter(SCINAME== spp,
                  SEASONA==2) %>% 
    vect()
})

#Create rasters from SVs
Rs <- map(.x = SVs, \(SV)
          rast(ext(SV), resolution=0.1, crs=crs(SV))
)

#Rasterize & turn into dataframe
R.dfs <- map2(.x = SVs, .y = Rs, \(SV, r)
              rasterize(x=SV, y=r, field = "SCINAME") %>% 
                as.data.frame(xy=TRUE)
)

#Set Ages, Sex to M, & the tsss to the mean for eacah species
#Obtain mean tsss for each species
mean.tsss <- map_dbl(njdf.l.br.age[c(1,3,5)], \(df)
                     mean(df$tsss.comb)
)

Ages <- c("Unk", "Adult", "Adult")
spp.dfs <- pmap(.l = list(R.dfs, Ages, mean.tsss), .f = \(Rdf, Age, tsss)
                Rdf %>% mutate(Age = Age, 
                               Sex = "M",
                               tsss.comb = tsss) %>% 
                  rename(B.Lat = y, B.Long = x)
)

#Extract elevations from CONI raster
elev_coni <- rast('Spatial_files/ETOPO1_coni.tif')
elevs <- terra::extract(x = elev_coni, y = spp.dfs[[1]][,c("B.Long", "B.Lat")], 
                          ID = FALSE) %>% 
  pull()
spp.dfs[[1]] <- spp.dfs[[1]] %>% mutate(B.Elev = as.numeric(elevs))


#Make predictions#
pred.dfs <- map2(.x = top.lms, .y = spp.dfs, .f = \(lm, df)
                 df %>% mutate(pred.mass = predict(lm, newdata = df))
)

# >Make map ---------------------------------------------------------------
pred.map <- list()
pred.map <- map(seq_along(Spp$Long), \(i){
  bbox <- st_bbox(c(xmin = min(pred.dfs[[i]]$B.Long), xmax = max(pred.dfs[[i]]$B.Long), 
                    ymin = min(pred.dfs[[i]]$B.Lat), ymax = max(pred.dfs[[i]]$B.Lat))) 
  ggplot(data = pred.dfs[[i]]) +
    geom_polygon(data=world, aes(x=long, y=lat, group=group), 
                 fill="gray80", colour = "gray80", size=0.3) +
    geom_polygon(data=lake, aes(x=long, y=lat, group=group), 
                 fill="white", colour = "white", size=0.3) +
    geom_tile(aes(x = B.Long, y = B.Lat, fill = pred.mass)) +
    geom_contour(aes(x = B.Long, y = B.Lat, z = pred.mass),
                 color = "white", size = .4, alpha = .7) + #, bins = 6
    geom_text_contour(aes(x = B.Long, y = B.Lat, z = round(pred.mass, 0)), 
                      stroke = 0.1, check_overlap = TRUE, 
                      min.size = 60, skip = 0, rotate = FALSE,
                      label.placer = label_placer_fraction(frac = 0.5)) +
    coord_sf(xlim = bbox[c(1,3)], ylim = bbox[c(2,4)], expand = FALSE, crs=4326) +
    scale_fill_viridis_c(option = "magma", 
                         name = "Mass (g)") + 
    map.theme +
    #labs(title = Spp$Long[i])
    theme(legend.position = "none")
})

#Note NULLs and widths of 0 don't alter this figure at all, but this is how you would adjust the spacing between plots if needed
isocline.fig <- ggarrange(pred.map[[1]], NULL, pred.map[[2]], NULL, pred.map[[3]], 
                          nrow = 1, widths = c(1, 0, 1, 0, 1))
ggsave(paste0("Plots/isoclines_fig", format(Sys.Date(), "%m.%d.%y"), ".png"), 
              bg = "white", height = 2.93, width = 8.5)

# Fig4: Parm estimates plot -----------------------------------------------------
#For whatever reason the order of explanatory vars is changed when excel is uploaded here
parm.df.form <- read.csv(paste0("Intermediates/parm.df.form_", format(Sys.Date(), "%m.%d.%y"), ".csv"))

parms.group <- parm.df.form %>% #filter(Important == "Y") %>% 
  group_by(DV, Data.set)
parms.impL <- parms.group %>% group_split()
names(parms.impL) <- parms.group %>% group_keys() %>%
  mutate(name = paste0(DV, "_", Data.set)) %>% 
  pull(name)

#Create function to plot parameters
plot.parms <- function(parms.df){
  ggplot(data = parms.df, aes(x = Parameter, y = Beta, 
                              group = interaction(Species, Season), color = Species)) +
    geom_errorbar(aes(ymin=LCI95, ymax=UCI95), 
                  alpha = .3, width=0, size = 1.5, position = position_dodge(width = 0.75)) +
    geom_errorbar(aes(ymin=LCI85, ymax=UCI85), 
                  alpha = .6, width=0, size = 3, position = position_dodge(width = 0.75)) +
    geom_point(size = 4, position = position_dodge(width = 0.75), aes(shape = Season)) + 
    scale_y_continuous(name = expression(paste("Parameter estimate (", beta, ")"))) + 
    scale_color_manual(values= colorblind_pal()(8)[c(4,6,8)], 
                       breaks = c("Nighthawk", "Nightjar", "Whip-poor-will")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12), 
          axis.text.x = element_text(size = 12, vjust = .58, angle = 60), 
          legend.title = element_text(size=14), legend.text = element_text(size=12),
          legend.position = "top",
          plot.margin = margin(10,10,8,25)) + 
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) + 
    ggtitle(paste(parms.df$DV)) #+ 
    #guides(color = "none", shape = "none")
}

#Apply custom function across list
#Breeding grounds data set plots.list[[1]] & [[3]]
plots.list <- lapply(parms.impL, plot.parms)
ggarrange(plots.list[[1]], plots.list[[3]], common.legend = TRUE, labels = "AUTO")
ggsave(paste0("Plots/Parm_ests_all_aic_less4_85CI_Br", format(Sys.Date(), "%m.%d.%y"), ".png"), bg = "white", width = 9)
