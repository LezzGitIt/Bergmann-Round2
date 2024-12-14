# Evaluation of Bergmann's Rule in Migratory Nightjars

This repository contains updated and additional code associated with the paper *"Environmental pressures on the breeding grounds drive Bergmannian clines in nightjars* (in review at the Journal of Biogeography). The study investigates whether three migratory nightjar species conform to Bergmann's rule, and explores the environmental and geographical factors influencing body size across their breeding and wintering grounds. The study uses GPS tracking data to assess competing hypotheses explaining variation in body size for each species, based on their breeding and wintering locations.

## Overview

This repository contains the full code for the project, which contains additional scripts that were used for data exploration and model validation. If interested in reproducing the analysis associated with *"Environmental pressures on the breeding grounds drive Bergmannian clines in nightjars"*, see this Dryad link [LINK AVAILABLE UPON PUBLICATION].

**Research Questions**:
- Do migratory nightjar species adhere to Bergmann's rule across their migratory ranges?
- How do environmental factors at breeding and wintering sites influence body size in these species?
- What are the key mechanistic hypotheses that explain Bergmannian patterns in body size?

## Methods

- **Data**: Morphometric data from >3,000 nightjars, and GPS tracking data from 189 individual birds
- **Analysis**: Linear models, AIC model selection, model averaging 

## Usage

### Requirements

- R (version 4.0 or higher)
- Required R packages: R packages are called at the top of every script 

### How to Run

To reproduce the analysis associated with *"Environmental pressures on the breeding grounds drive Bergmannian clines in nightjars"*

1. Download the data and scripts from Dryad [LINK AVAILABLE UPON PUBLICATION]

2. Install required R packages:
    ```R
  packages <- c(
  "tidyverse", "readxl", "xlsx", "chron", "gridExtra", "mvnormtest", "janitor", 
  "broom.mixed", "ggrepel", "ggthemes", "viridis", "AICcmodavg", "lme4", "nlme", 
  "sf", "smatr", "scales", "geoR", "cowplot", "adespatial", "spdep", "ggpubr", 
  "sp", "naniar", "stringi", "sjPlot", "broom", "mapview", "conflicted", 
  "rnaturalearthdata", "rnaturalearth", "rnaturalearthhires", "smoothr", "terra", 
  "metR", "raster", "MuMIn"
)
install.packages("requiRements")
requiRements::install(packages)
    ```
    
3. Load the scripts and run to reproduce the analysis and all figures in the main text:
    ```R
script_files <- list.files(path = "Scripts", pattern = "\\.R$", full.names = TRUE)
map(script_files, \(files){
 source(files)
 })
    ```

## Data

The data and code needed to reproduce all analyses, figures, and tables in the main text of the manuscript are available on Dryad (LINK AVAILABLE UPON PUBLICATION). The collation of this rich morphological dataset with information across the annual cycle is the most complete for Caprimulgids to date. We would be interested in future collaborations using these data. 

## License

Feel free to use, modify, and distribute this code however you wish. No restrictions apply.

## Acknowledgments

Special thanks to the many researchers who provided morphological and tracking data on nightjars, the many field assistants and their tremendous effort capturing nightjars, and to Elly Knight and Alicia Korpach for their integral part in these analyses and the resulting manuscript. 

## Contact

For more information or questions, please contact skinnerayayron93 [at] gmail [dot] com
