# clear things ####

rm(list = ls())

# packages ####

pacman::p_load(tidyverse,
               metafor,
               pander,
               stringr,
               ape,
               kableExtra,
               patchwork,
               here,
               lme4,
               readxl,
               metaAidR
)

# loading data ####

dat <- read_csv(here("data", "cleaned.csv"), na = c("", "NA")) 

names(dat) 

#  getting effect size #### 


