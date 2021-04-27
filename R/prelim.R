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
               rotl,
               orchaRd,
               emmeans,
               clubSandwich
)

# function for getting d for proportional data

smdp <- function(m1, m2, n1, n2) {
  mean_diff <- (car::logit(m1) - car::logit(m2))
  j_cor <- 1 - (3 /(4*(n1 + n2) -9))
  smd <- (mean_diff/(pi/sqrt(3)))*j_cor
  var <- ((n1 + n2)/(n1*n2) + (smd^2)/(2*(n1 + n2)))*(j_cor^2)
  invisible(data.frame(yi = smd , vi = var))
}

d_to_lnor <- function(smd){
  lnor <- smd * (sqrt(3)/pi)
  return(lnor)
} 

# funciont for getting d for traditional data

smdm <- function(m1, m2, n1, n2, sd1, sd2) {
  mean_diff <- m1 - m2
  sd_pool <- sqrt( ((n1 - 1)*sd1^2 + (n2 - 1)*sd2) / (n1 + n2 - 2) )
  j_cor <- 1 - (3 /(4*(n1 + n2) -9))
  smd <- (mean_diff/sd_pool)*j_cor
  var <- ((n1 + n2)/(n1*n2) + (smd^2)/(2*(n1 + n2)))*(j_cor^2)
  invisible(data.frame(yi = smd , vi = var))
}

d_to_lnor <- function(smd){
  lnor <- smd * (sqrt(3)/pi)
  return(lnor)
} 

lnor_to_d <- function(lnor){
  smd <- lnor*(pi/sqrt(3))
}
# invisible

# loading data ####

dat_full <- read_csv(here("data", "dat_09042020.csv"), na = c("", "NA")) 

#names(dat_full) 
#str(dat_full) 

# deleting unusable rows #####

dat_full %>% filter(is.na(Treatment_lifespan_variable) == FALSE) -> dat

dim(dat)
dim(dat_full)
# separating two kinds

effect_type <- str_detect(dat$Lifespan_parameter, "Me")

# creating effect sizes
dat$yi <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, dat$Error_experimental_SD, dat$Error_control_SD)[[1]], smdp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[1]])

dat$vi <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, dat$Error_experimental_SD, dat$Error_control_SD)[[2]], smdp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[2]])

dat$yi2 <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[1]], smdp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[1]])

dat$vi2 <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[2]], smdp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[2]])

dat$yi3 <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[1]], smdp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[1]])

dat$vi3 <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[2]], smdp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[2]])

# effect-level ID

dat$Effect_ID <- 1:nrow(dat)
dat$Phylogeny <- sub(" ", "_",  dat$Species_Latin,)
dat$Effect_type <- effect_type

# shared control 
# this does not seem to work
#V_matrix <- make_VCV_matrix(dat, V= "vi", cluster = "Shared_control", obs = "Effect_ID")

V_matrix <- impute_covariance_matrix(vi = dat$vi, cluster = dat$Shared_control, r = 0.5)

# phylogeny

tree <- read.tree(here("data/tree_rotl.tre"))

tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)
# meta-analysis basics
# phylo model
mod <-  rma.mv(yi, V = V_matrix, mod = ~ 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod) 
robust(mod, cluster = dat$Species_Latin)
round(i2_ml(mod)*100, 2)

funnel(mod)
funnel (mod, yaxis="seinv")

# alternative
mod0 <-  rma.mv(yi, V = vi, mod = ~ 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod0)
robust(mod0, cluster = dat$Study)

mod01 <-  rma.mv(yi, V = vi, mod = ~ 1, random = list( ~1|Species_Latin, ~1|Study, ~1|Effect_ID), data = dat, test = "t")
summary(mod01)
coef_test(mod01, vcov = "CR2", cluster = dat$Species_Latin)

mod02 <-  rma.mv(yi, V = vi, mod = ~ 1, random = list(~1|Study, ~1|Effect_ID), data = dat, test = "t")
summary(mod02)
coef_test(mod02, vcov = "CR2", cluster = dat$Study)

# sex effect
mod1 <-  rma.mv(yi, V = V_matrix, mod = ~ Sex -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod1) 

r2_ml(mod1)

res <- mod_results(mod1, "Sex")

orchard_plot(mod1, mod = "Sex", xlab = "Standardised mean difference")

# just at means of all continiosu variables (if we do not set anything)
res <- qdrg(object = mod1, data = dat)

# marginal means ; all groups are proportionally weighted
emmeans(res, specs = ~1, df = mod1$dfs, weights = "prop") 
emmeans(res, specs = "Sex", df = mod1$dfs, weights = "prop") 

#
mod2 <-  rma.mv(yi, V = V_matrix, mod = ~ Sex + Gonads_removed + Effect_type, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod2) 

res2 <- qdrg(object = mod2, data = dat)

# marginal means ; all groups are proportionally weighted
emmeans(res2, specs = ~1, df = mod1$dfs, weights = "prop") 
emmeans(res2, specs = "Sex", df = mod1$dfs, weights = "prop") 

# Sex difference: new set of analyses #### 

## recalculating - effect size

# reudcting data

dat <- dat[(is.na(dat$vi2) == FALSE), ]


V_matrix <- impute_covariance_matrix(vi = dat$vi, cluster = dat$Shared_control, r = 0.5)

# phylogeny
tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% dat$Phylogeny)])

tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)

###

model <-  rma.mv(yi, V = V_matrix, mod = ~ 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(model) 
round(i2_ml(model)*100,2)

model2 <-  rma.mv(yi, V = V_matrix, mod = ~ Sex -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(model2) 

orchard_plot(model2, mod = "Sex", xlab = "Standardised mean difference")

######## DO NOT RUN
# we can run - some heteroscad models
# this does not improve model
mod_h <-  rma.mv(yi, V = vi, mod = ~ Sex-1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~Sex|Effect_ID), rho = 0, struct = "HCS", R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod_h) 

...

#########################################

# longevity
# this needs more work here
# SEM, Probable error???, SD, CI, Interquartile range
# pobably easier to calcuate by hand
dat_long <- dat[effect_type == TRUE, ]

# I think we assume - binomial error
# then this will be all fine...
dat_surv <- as.data.frame(dat[effect_type == FALSE, ])


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
