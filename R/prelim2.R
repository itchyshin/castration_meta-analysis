# M/F vs CM/CF
# clear things ####

rm(list = ls())

# packages ####

#ochaRd
# 
# install.packages("devtools")
# install.packages("tidyverse")
# #install.packages("metafor")
# install.packages("patchwork")
# install.packages("R.rsp")
# 
# devtools::install_github("daniel1noble/orchaRd", force = TRUE, build_vignettes = TRUE)
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_opts = "") 
# 
# #emmeans
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_opts = "")
# # metafor
# install.packages("remotes")
# remotes::install_github("wviechtb/metafor")

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


# effect-level ID

dat$Effect_ID <- 1:nrow(dat)
dat$Phylogeny <- sub(" ", "_",  dat$Species_Latin,)
dat$Effect_type <- effect_type


# we create a longer data format

dat1 <- dat
dat2 <- dat
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

# here we create the ratio of M/F or F/M 
dat1$yi <- ifelse(effect_type == TRUE, lnrrm(dat$Control_lifespan_variable,dat$Opposite_sex_lifespan_variable,  dat$Sample_size_control, dat$Sample_size_opposite_sex, datl$Error_control_SD, dat$Error_opposite_sex_SD)[[1]], lnrrp(dat$Control_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_control,  dat$Sample_size_opposite_sex)[[1]])

dat1$vi <-ifelse(effect_type == TRUE, lnrrm(dat$Control_lifespan_variable,dat$Opposite_sex_lifespan_variable,  dat$Sample_size_control, dat$Sample_size_opposite_sex, datl$Error_control_SD, dat$Error_opposite_sex_SD)[[2]], lnrrp(dat$Control_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_control,  dat$Sample_size_opposite_sex)[[2]])

# here we create CM/F or CF/M
dat2$yi <- ifelse(effect_type == TRUE, lnrrm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[1]], lnrrp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[1]])

dat2$vi <- ifelse(effect_type == TRUE, lnrrm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[2]], lnrrp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[2]])
# putting two data frames
dat_long <- rbind(dat1, dat2)

# putt 2 new column

dat_long$Obs <- factor(1:dim(dat_long)[[1]])
dat_long$Comp_type <- as.factor(rep(c("both_normal", "one_neutured"), each = dim(dat_long)[[1]]/2))




# SMD
# yi2 = Oppsite  - Control (male - female or male - female)
# dat$yi2 <- ifelse(effect_type == TRUE, smdm(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control, dat$Error_opposite_sex_SD, dat$Error_control_SD)[[1]], smdp(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control)[[1]])
# 
# dat$vi2 <- ifelse(effect_type == TRUE, smdm(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control, dat$Error_opposite_sex_SD, dat$Error_control_SD)[[2]], smdp(dat$Opposite_sex_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_opposite_sex, dat$Sample_size_control)[[2]])

# 
# dat$yi3 <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[1]], smdp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[1]])
# 
# dat$vi3 <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[2]], smdp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[2]])



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

orchard_plot(mod, mod = "Int", xlab = "log response ratio (lnRR)")