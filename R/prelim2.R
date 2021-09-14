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

dat <- read_csv(here("data", "dat2_16082021.csv"), na = c("", "NA")) 



source(here("R","function.R"), chdir = TRUE)

# looking at data

dat
dim(dat)
names(dat)

effect_type <- str_detect(dat$Lifespan_parameter, "Me")


# effect-level ID

dat$Effect_ID <- 1:nrow(dat)
dat$Phylogeny <- sub(" ", "_",  dat$Species_Latin,)
dat$Effect_type <- effect_type
# yi = Treatment - Control

# we create a longer data format

dat1 <- dat
dat2 <- dat
# lnRR

# here we create the ratio of F/M 
dat1$yi <- ifelse(effect_type == TRUE, 
                  lnrrm(dat$Female_control_lifespan_variable, dat$Male_control_lifespan_variable,
                        dat$Sample_size_female_control, dat$Sample_size_male_control, 
                        dat$Error_female_control_SD, dat$Error_male_control_SD)[[1]], 
                  lnrrp(dat$Female_control_lifespan_variable, dat$Male_control_lifespan_variable, 
                        dat$Sample_size_female_control,  dat$Sample_size_male_control)[[1]])

dat1$vi <- ifelse(effect_type == TRUE, 
                  lnrrm(dat$Female_control_lifespan_variable, dat$Male_control_lifespan_variable,
                        dat$Sample_size_female_control, dat$Sample_size_male_control, 
                        dat$Error_female_control_SD, dat$Error_male_control_SD)[[2]], 
                  lnrrp(dat$Female_control_lifespan_variable, dat$Male_control_lifespan_variable, 
                        dat$Sample_size_female_control,  dat$Sample_size_male_control)[[2]])

# here we create CF/CM
dat2$yi <- ifelse(effect_type == TRUE, 
                  lnrrm(dat$Female_sterilization_lifespan_variable, dat$Male_sterilization_lifespan_variable,
                        dat$Sample_size_female_sterilization, dat$Sample_size_male_sterilization, 
                        dat$Error_female_sterilization_SD, dat$Error_male_sterilization_SD)[[1]], 
                  lnrrp(dat$Female_sterilization_lifespan_variable, dat$Male_sterilization_lifespan_variable, 
                        dat$Sample_size_female_sterilization,  dat$Sample_size_male_sterilization)[[1]])

dat2$vi <-  ifelse(effect_type == TRUE, 
                   lnrrm(dat$Female_sterilization_lifespan_variable, dat$Male_sterilization_lifespan_variable,
                         dat$Sample_size_female_sterilization, dat$Sample_size_male_sterilization, 
                         dat$Error_female_sterilization_SD, dat$Error_male_sterilization_SD)[[2]], 
                   lnrrp(dat$Female_sterilization_lifespan_variable, dat$Male_sterilization_lifespan_variable, 
                         dat$Sample_size_female_sterilization,  dat$Sample_size_male_sterilization)[[2]])

# merging data frames
dat_long <- rbind(dat1, dat2)

# putt 2 new column

dat_long$Obs <- factor(1:dim(dat_long)[[1]])
dat_long$Comp_type <- as.factor(rep(c("both_normal", "both_neutured"), each = dim(dat_long)[[1]]/2))


# shared control 
# this does not seem to work
#V_matrix <- make_VCV_matrix(dat, V= "vi", cluster = "Shared_control", obs = "Effect_ID")
# TODO we need to do for the longer data
V_matrix_long <- impute_covariance_matrix(vi = dat_long$vi, cluster = dat_long$Shared_control, r = 0.5)

# phylogeny

tree <- read.tree(here("data/tree_rotl.tre"))

tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)
# meta-analysis basics
### longer data


# TODO we will do hetero model
# we can run - some heteroscad models
# this does not improve model
mod_n <-  rma.mv(yi, V = V_matrix_long, mod = ~ Comp_type - 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID, ~1|Obs ), R = list(Phylogeny = cor_tree), data = dat_long, test = "t")
summary(mod_n) 

orchard_plot(mod_n, mod = "Comp_type", xlab = "log response ratio (lnRR)")

mod_h <-  rma.mv(yi, V = vi, mod = ~ Comp_type - 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~Comp_type|Study, ~1|Effect_ID, ~Comp_type|Obs), rho = 0, struct = "HCS", R = list(Phylogeny = cor_tree), data = dat_long, test = "t")
summary(mod_h) 

mod_nb <-  rma.mv(yi, V = vi, mod = ~ Comp_type, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID, ~1|Obs ), R = list(Phylogeny = cor_tree), data = dat_long, test = "t")
summary(mod_nb) 

# doing both group sepratly 
