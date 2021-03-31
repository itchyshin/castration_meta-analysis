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

# function for getting d for propotional data

smdp <- function(m1, m2, n1, n2) {
  mean_diff <- (car::logit(m1) - car::logit(m2))
  j_cor <- 1 - (3 /(4*(n1 + n2) -9))
  smd <- mean_diff/(pi/sqrt(3))*j_cor
  var <- ((n1 + n2)/(n1*n2) + smd^2/(2*(n1 + n2)))*(j_cor^2)
  invisible(data.frame(yi = smd , vi = var))
}

d_to_lnor <- function(smd){
  lnor <- smd * (sqrt(3)/pi)
  return(lnor)
} 

# invisible

# loading data ####

dat <- read_csv(here("data", "cleaned_ML.csv"), na = c("", "NA")) 

names(dat) 
str(dat) 

# deleting unusable rows #####

dat %>% filter(is.na(Treatment_lifespan_variable) == FALSE) -> dat1

# separating two kinds

groups <- str_detect(dat1$Lifespan_parameter, "Me")

# longevity
# this needs more work here
# SEM, Probable error???, SD, CI, Interquartile range
# pobably easier to calcuate by hand
dat_long <- dat1[groups == TRUE, ]

# I think we assume - binomial error
# then this will be all fine...
dat_surv <- as.data.frame(dat1[groups == FALSE, ])

######################################################################################
######################################################################################
# phylo tree ####
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
######################################################################################
######################################################################################

#  getting effect size #### 

effect <- smdp(dat_surv$Treatment_lifespan_variable, dat_surv$Control_lifespan_variable, dat_surv$Sample_size_sterilization, dat_surv$Sample_size_control)

effect2 <- smdp(dat_surv$Treatment_lifespan_variable, dat_surv$Opposite_sex_lifespan_variable, dat_surv$Sample_size_sterilization, dat_surv$Sample_size_opposite_sex)

names(effect2) <- c("yi2", "vi2")

dat_surv <-cbind(dat_surv, effect, effect2)
dat_surv$Effect_ID <- factor(1:dim(dat_surv)[1])

mod_surv <- rma.mv(yi, V = vi, random = list(~1|Study, ~1|Effect_ID), data = dat_surv)
summary(mod_surv) 


# opposite sex - it is still reducing death rate
mod_opp <- rma.mv(yi2, V = vi2, random = list(~1|Study, ~1|Effect_ID), data = dat_surv)
summary(mod_opp) 
  
  
mod_surv1 <- rma.mv(yi, V = vi, mod = ~ Type_of_sterilization - 1, random = list(~1|Study, ~1|Effect_ID), data = dat_surv)
summary(mod_surv1) 

mod_surv2 <- rma.mv(yi, V = vi, mod = ~ Sex - 1, random = list(~1|Study, ~1|Effect_ID), data = dat_surv)
summary(mod_surv2) 



