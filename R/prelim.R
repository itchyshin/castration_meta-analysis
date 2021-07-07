# clear things ####

rm(list = ls())

# packages ####

#ochaRd

install.packages("devtools")
install.packages("tidyverse")
#install.packages("metafor")
install.packages("patchwork")
install.packages("R.rsp")

devtools::install_github("daniel1noble/orchaRd", force = TRUE, build_vignettes = TRUE)
remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_opts = "") 

#emmeans
remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_opts = "")
# metafor
install.packages("remotes")
remotes::install_github("wviechtb/metafor")

# loading
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





# loading data ####

dat_full <- read_csv(here("data", "dat_07072021.csv"), na = c("", "NA")) 

source(here("R","function.R"), chdir = TRUE)

#names(dat_full) 
#str(dat_full) 

# deleting unusable rows #####
# excluding Vasectomy 
dat_full %>% filter(is.na(Treatment_lifespan_variable) == FALSE) %>% filter(Type_of_sterilization != "Vasectomy") -> dat


dim(dat)
dim(dat_full)
# separating two kinds

effect_type <- str_detect(dat$Lifespan_parameter, "Me")

# creating effect sizes

# yi = Treatment - Control

#lnRR
dat$yi <- ifelse(effect_type == TRUE, lnrrm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, dat$Error_experimental_SD, dat$Error_control_SD)[[1]], lnrrp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[1]])

dat$vi <- ifelse(effect_type == TRUE, lnrrm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, dat$Error_experimental_SD, dat$Error_control_SD)[[2]], lnrrp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[2]])


# SMD
# dat$yi <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, dat$Error_experimental_SD, dat$Error_control_SD)[[1]], smdp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[1]])
# 
# dat$vi <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, dat$Error_experimental_SD, dat$Error_control_SD)[[2]], smdp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[2]])

# lnRR

dat$yi2 <- ifelse(effect_type == TRUE, lnrrm(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control, dat$Error_opposite_sex_SD, dat$Error_control_SD)[[1]], lnrrp(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control)[[1]])

dat$vi2 <- ifelse(effect_type == TRUE, lnrrm(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control, dat$Error_opposite_sex_SD, dat$Error_control_SD)[[2]], lnrrp(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control)[[2]])
# SMD
# yi2 = Oppsite  - Control (male - female or male - female)
# dat$yi2 <- ifelse(effect_type == TRUE, smdm(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control, dat$Error_opposite_sex_SD, dat$Error_control_SD)[[1]], smdp(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control)[[1]])
# 
# dat$vi2 <- ifelse(effect_type == TRUE, smdm(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control, dat$Error_opposite_sex_SD, dat$Error_control_SD)[[2]], smdp(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control)[[2]])

# 
# dat$yi3 <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[1]], smdp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[1]])
# 
# dat$vi3 <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[2]], smdp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[2]])

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

orchard_plot(mod1, mod = "Sex", xlab = "log response ratio (lnRR)")

mod1b <-  rma.mv(yi, V = V_matrix, mod = ~ Sex, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod1b) 

# exuding inflectional study

dat2 <- dat[dat$Order_extracted != 87 & dat$Order_extracted != 88, ]
V_matrix2 <- impute_covariance_matrix(vi = dat2$vi, cluster = dat2$Shared_control, r = 0.5)


mod12 <-  rma.mv(yi, V = V_matrix2, mod = ~ Sex -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat2, test = "t")
summary(mod12) 

r2_ml(mod12)

res <- mod_results(mod12, "Sex")

orchard_plot(mod12, mod = "Sex", xlab = "log response ratio (lnRR)")

# mod = age
model4 <-  rma.mv(yi, V = V_matrix, mod = ~ 1 + Maturity_at_treatment_ordinal, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(model4) 
regplot(model4)

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
emmeans(res2, specs = ~1, df = mod1$ddf, weights = "prop") 
emmeans(res2, specs = "Sex", df = mod1$ddf, weights = "prop") 


#########################
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

# Mod = sex
model1 <-  rma.mv(yi, V = V_matrix, mod = ~ Sex, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(model1) 

model2 <-  rma.mv(yi, V = V_matrix, mod = ~ Sex -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(model2) 

orchard_plot(model2, mod = "Sex", xlab = "Standardised mean difference")

# type sterilzation

model3 <-  rma.mv(yi, V = V_matrix, mod = ~ Type_of_sterilization -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(model3) 


# mod = age
model4 <-  rma.mv(yi, V = V_matrix, mod = ~ 1 + Maturity_at_treatment_ordinal, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(model4) 
regplot(model4)

# sex d
# male
mdat <- subset(dat, Sex == "Male")

V_matrix2 <- impute_covariance_matrix(vi = mdat$vi, cluster = mdat$Shared_control, r = 0.5)

sub1 <-  rma.mv(yi, V = V_matrix2, mod = ~ 1 + yi2, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = mdat, test = "t")
summary(sub1) 
regplot(sub1)
robust(sub1, cluster = mdat$Species_Latin)

qplot(yi2, yi, data = mdat)

# female
fdat <- subset(dat, Sex == "Female")

V_matrix3 <- impute_covariance_matrix(vi = fdat$vi, cluster = fdat$Shared_control, r = 0.5)

sub2 <-  rma.mv(yi, V = V_matrix3, mod = ~ 1 + yi2, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = fdat, test = "t")
summary(sub2) 
regplot(sub2)
robust(sub2, cluster = fdat$Species_Latin)
qplot(yi2, yi, data = fdat)




######## DO NOT RUN
# we can run - some heteroscad models
# this does not improve model
mod_h <-  rma.mv(yi, V = vi, mod = ~ Sex-1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~Sex|Effect_ID), rho = 0, struct = "HCS", R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod_h) 

...


###########################################
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
