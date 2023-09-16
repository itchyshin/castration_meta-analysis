# zoo 4 - creating sterilzed or not. 

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
               #metaAidR,
               rotl,
               orchaRd,
               emmeans,
               clubSandwich,
               png,
               grid,
               here
)

library(apextra)

########################
# data: life expectancy
########################

dat_full <- read_csv(here("data", "zoo.csv"), na = c("", "NA"))

# trimming data

dim(dat_full)



# extra
tax <- read.csv(here("data", "vertlife_taxonomy_translation_table.csv"))

# animal order and tree relationships 
#dford <- readRDS(here("Rdata", "dford.RDS"))

dat_full %>% left_join(tax, by = c("species" = "zims.species")) -> dat_full

dat_full %>% filter(species != "Chrysocyon brachyurus" &
                      species != "Crocuta crocuta" &
                      species != "Neofelis nebulosa" &
                      species != "Panthera uncia" &
                      species != "Pseudocheirus peregrinus") %>% 
  mutate(phylogeny = gsub(" ", "_", vertlife.species)) -> dat

# getting a tree

tree <- read.tree( here("data", "tree_zoo4.tre"))

# life span data 
to_drop <-
  tree$tip.label[which(!(tree$tip.label %in% unique(dat$phylogeny)))]

tree <- drop.tip(tree, to_drop)
length(tree$tip.label)

tree <- as.ultrametric(tree)

#tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)

# F normal vs M normal

