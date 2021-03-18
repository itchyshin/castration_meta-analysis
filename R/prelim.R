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
               metaAidR,
               rotl
)

# loading data ####

dat <- read_csv(here("data", "cleaned_ML.csv"), na = c("", "NA")) 

names(dat) 
str(dat) 

# phylo tree
myspecies <- as.character(unique(dat$Species_Latin)) #get list of species
length(myspecies) #14 species
str_sort(myspecies) #visual check
taxa <- tnrs_match_names(names = myspecies) # using *rotl* package to retrieve synthetic species tree from Open Tree of Life
dim(taxa) #14 species - all matched
table(taxa$approximate_match) #0 approximate matches
tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name")  #retrieve the tree
plot(tree, cex=.6, label.offset =.1, no.margin = TRUE) #good
is.binary(tree) #no polytomies
str(tree) 
tree$tip.label <- gsub("_"," ", tree$tip.label) #get rid of the underscores from tree tip labels
intersect(myspecies, tree$tip.label) #14 - perfect overlap and differences with myspecies list
write.tree(tree, file="./data/tree_rotl.tre") #save the tree
# tree <- read.tree(file="./data/tree_rotl.tre") #if you need to read in the tree

#  getting effect size #### 


