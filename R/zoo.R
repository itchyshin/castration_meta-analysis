# quick analysis of zoo data

# loading
pacman::p_load(tidyverse,
               metafor,
               pander,
               stringr,
               ape,
               here,
               kableExtra,
               patchwork,
               lme4,
               readxl,
               metaAidR,
               rotl,
               orchaRd,
               emmeans,
               clubSandwich,
               png,
               grid,
               here
)

# data
dat <- read_csv(here("data", "zoo.csv"), na = c("", "NA"))
glimpse(dat)

# filtering and getting SD
dat %>% mutate()
  
  
# getting effect size

dat <- escalc("ROM", 
              m1i = ,
              m2i = ,
              sd1i = ,
              sd2i = ,
              n1i = ,
              n2i = ,
              data = dat,
              )