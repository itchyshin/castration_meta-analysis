
#TODO 

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

#dat_full <- read_csv(here("data", "dat_07072021.csv"), na = c("", "NA")) 
dat_full <- read_csv(here("data", "data_14092021.csv"), na = c("", "NA"))
glimpse(dat_full)

source(here("R","function.R"), chdir = TRUE)

#names(dat_full) 
#str(dat_full) 

# deleting unusable rows #####
# excluding Vasectomy 
# TODO this is important!!!!!
dat_full %>% filter(is.na(Treatment_lifespan_variable) == FALSE) %>% 
  filter(Type_of_sterilization != "Vasectomy") %>% 
  mutate_if(is.character, as.factor) -> dat


dim(dat)
dim(dat_full)
# separating two kinds

effect_type <- ifelse(str_detect(dat$Lifespan_parameter, "Me"), "longevity", "mortality")


# effect-level ID

dat$Effect_ID <- 1:nrow(dat)
dat$Phylogeny <- sub(" ", "_",  dat$Species_Latin)
dat$Effect_type <- effect_type

# key variables
names(dat)

# Sex
# Wild_or_semi_wild
# Maturity_at_treatment_ordinal (missing values + contionous variable)
# Gonads_removed (only applies to female)
# Controlled_treatments
# Effect_typle

# let's get CVs

dat %>% group_by(Study) %>% summarise(cv2_cont = mean((Error_control_SD/Control_lifespan_variable)^2, na.rm = T), cv2_trt = mean((Error_experimental_SD/Treatment_lifespan_variable)^2, na.rm = T), cv2_opst = mean((Error_opposite_sex_SD/Opposite_sex_lifespan_variable)^2, na.rm = T), n_cont = mean(Sample_size_control, na.rm = T), n_trt =  mean(Sample_size_sterilization, na.rm = T), n_opst =  mean(Sample_size_opposite_sex, na.rm = T)) %>% 
  ungroup() %>% 
  summarise(cv2_cont = weighted.mean(cv2_cont, n_cont, na.rm = T), cv2_trt = weighted.mean(cv2_trt, n_trt, na.rm = T), cv2_opst = weighted.mean(cv2_opst, n_opst, na.rm = T)) -> cvs

# we create a longer data format

dat1 <- dat
dat2 <- dat
# creating effect sizes

# yi = Treatment - Control

# #lnRR
# # older version
# dat$yi <- ifelse(effect_type == TRUE, lnrrm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, dat$Error_experimental_SD, dat$Error_control_SD)[[1]], lnrrp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[1]])
# 
# dat$vi <- ifelse(effect_type == TRUE, lnrrm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, dat$Error_experimental_SD, dat$Error_control_SD)[[2]], lnrrp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[2]])

# lnRR
# using CV
dat$yi <- ifelse(effect_type == TRUE, lnrrm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, cvs[["cv2_trt"]],cvs[["cv2_cont"]])[[1]], lnrrp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[1]])

dat$vi <- ifelse(effect_type == TRUE, lnrrm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, cvs[["cv2_trt"]],cvs[["cv2_cont"]])[[2]], lnrrp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[2]])

# SMD
# dat$yi <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, dat$Error_experimental_SD, dat$Error_control_SD)[[1]], smdp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[1]])
# 
# dat$vi <- ifelse(effect_type == TRUE, smdm(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control, dat$Error_experimental_SD, dat$Error_control_SD)[[2]], smdp(dat$Treatment_lifespan_variable, dat$Control_lifespan_variable, dat$Sample_size_sterilization, dat$Sample_size_control)[[2]])



#lnRR

# M/F or F/M - control/control

# here we create the ratio of M/F or F/M 
# old 
# dat1$yi <- ifelse(effect_type == TRUE, 
#                   lnrrm(dat$Control_lifespan_variable,dat$Opposite_sex_lifespan_variable,  
#                         dat$Sample_size_control, dat$Sample_size_opposite_sex, 
#                         dat$Error_control_SD, dat$Error_opposite_sex_SD)[[1]], 
#                   lnrrp(dat$Control_lifespan_variable, dat$Opposite_sex_lifespan_variable, 
#                         dat$Sample_size_control,  dat$Sample_size_opposite_sex)[[1]])
# 
# dat1$vi <-ifelse(effect_type == TRUE, 
#                  lnrrm(dat$Control_lifespan_variable,dat$Opposite_sex_lifespan_variable,  
#                        dat$Sample_size_control, dat$Sample_size_opposite_sex, 
#                        dat$Error_control_SD, dat$Error_opposite_sex_SD)[[2]], 
#                  lnrrp(dat$Control_lifespan_variable, dat$Opposite_sex_lifespan_variable, 
#                        dat$Sample_size_control,  dat$Sample_size_opposite_sex)[[2]])
# 
# # here we create CM/F or CF/M
# dat2$yi <- ifelse(effect_type == TRUE, 
#                   lnrrm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable,
#                         dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, 
#                         dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[1]], 
#                   lnrrp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, 
#                         dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[1]])
# 
# dat2$vi <- ifelse(effect_type == TRUE, 
#                   lnrrm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable,
#                         dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, 
#                         dat$Error_experimental_SD, dat$Error_opposite_sex_SD)[[2]], 
#                   lnrrp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, 
#                         dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[2]])


dat1$yi <- ifelse(effect_type == TRUE, 
                  lnrrm(dat$Control_lifespan_variable,dat$Opposite_sex_lifespan_variable,  
                        dat$Sample_size_control, dat$Sample_size_opposite_sex, 
                        cvs[["cv2_cont"]],cvs[["cv2_opst"]])[[1]], 
                  lnrrp(dat$Control_lifespan_variable, dat$Opposite_sex_lifespan_variable, 
                        dat$Sample_size_control,  dat$Sample_size_opposite_sex)[[1]])

dat1$vi <-ifelse(effect_type == TRUE, 
                 lnrrm(dat$Control_lifespan_variable,dat$Opposite_sex_lifespan_variable,  
                       dat$Sample_size_control, dat$Sample_size_opposite_sex, 
                       cvs[["cv2_cont"]],cvs[["cv2_opst"]])[[2]], 
                 lnrrp(dat$Control_lifespan_variable, dat$Opposite_sex_lifespan_variable, 
                       dat$Sample_size_control,  dat$Sample_size_opposite_sex)[[2]])

# here we create CM/F or CF/M
dat2$yi <- ifelse(effect_type == TRUE, 
                  lnrrm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable,
                        dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, 
                        cvs[["cv2_trt"]],cvs[["cv2_opst"]])[[1]], 
                  lnrrp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, 
                        dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[1]])

dat2$vi <- ifelse(effect_type == TRUE, 
                  lnrrm(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable,
                        dat$Sample_size_sterilization, dat$Sample_size_opposite_sex, 
                        cvs[["cv2_trt"]],cvs[["cv2_opst"]])[[2]], 
                  lnrrp(dat$Treatment_lifespan_variable, dat$Opposite_sex_lifespan_variable, 
                        dat$Sample_size_sterilization, dat$Sample_size_opposite_sex)[[2]])




# putting two data frames
dat_long <- rbind(dat1, dat2)

# putt 2 new column

dat_long$Obs <- factor(1:dim(dat_long)[[1]])
dat_long$Comp_type <- as.factor(rep(c("both_normal", "one_neutured"), each = dim(dat_long)[[1]]/2))


dim(dat_long)

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


## species-wise funnel plot 



# alternative
# mod0 <-  rma.mv(yi, V = vi, mod = ~ 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
# summary(mod0)
# robust(mod0, cluster = dat$Study)
# 
# mod01 <-  rma.mv(yi, V = vi, mod = ~ 1, random = list( ~1|Species_Latin, ~1|Study, ~1|Effect_ID), data = dat, test = "t")
# summary(mod01)
# coef_test(mod01, vcov = "CR2", cluster = dat$Species_Latin)
# 
# mod02 <-  rma.mv(yi, V = vi, mod = ~ 1, random = list(~1|Study, ~1|Effect_ID), data = dat, test = "t")
# summary(mod02)
# coef_test(mod02, vcov = "CR2", cluster = dat$Study)

# missing

summary(dat$Sex)

# sex effect
mod1 <-  rma.mv(yi, V = V_matrix, mod = ~ Sex -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod1) 

r2_ml(mod1)

res <- mod_results(mod1, "Sex")

orchard_plot(mod1, mod = "Sex", xlab = "log response ratio (lnRR)")

mod1b <-  rma.mv(yi, V = V_matrix, mod = ~ Sex, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod1b) 


# excluding influential study

# dat2 <- dat[dat$Order_extracted != 87 & dat$Order_extracted != 88, ]
# V_matrix2 <- impute_covariance_matrix(vi = dat2$vi, cluster = dat2$Shared_control, r = 0.5)
# 
# 
# mod12 <-  rma.mv(yi, V = V_matrix2, mod = ~ Sex -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat2, test = "t")
# summary(mod12) 
# 
# r2_ml(mod12)
# 
# res <- mod_results(mod12, "Sex")
# 
# orchard_plot(mod12, mod = "Sex", xlab = "log response ratio (lnRR)")

# missing 

summary(dat$Wild_or_semi_wild)

# environment - there seems to be too many categories 
# what to put together - can farm and outdoor encloure can be put togehter???

mod3 <-  rma.mv(yi, V = V_matrix, mod = ~ Wild_or_semi_wild -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod3) 

r2_ml(mod3)

#res <- mod_results(mod1, "Environment")

orchard_plot(mod3, mod = "Environment", xlab = "log response ratio (lnRR)")



# mod = age

# missing 

summary(dat$Maturity_at_treatment_ordinal)

model4 <-  rma.mv(yi, V = V_matrix, mod = ~ 1 + Maturity_at_treatment_ordinal, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(model4) 
regplot(model4,  col = dat$Sex)


# sex mat interaction

model4b <-  rma.mv(yi, V = V_matrix, mod = ~ Sex*Maturity_at_treatment_ordinal, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(model4b) 

model4c <-  rma.mv(yi, V = V_matrix, mod = ~ Sex*Maturity_at_treatment_ordinal - 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(model4c) 


# drawing ggplots

# drawing sex - matruality interaciton 

pred_sex_matuarity <-predict.rma(model4c) 

# plotting

plot <-  dat %>% 
  filter(!is.na(Maturity_at_treatment_ordinal))  %>% # getting ride of NA values
  mutate(ymin = pred_sex_matuarity$ci.lb, 
         ymax = pred_sex_matuarity$ci.ub, 
         ymin2 = pred_sex_matuarity$cr.lb,
         ymax2 = pred_sex_matuarity$cr.ub,
         pred = pred_sex_matuarity$pred) %>% 
  ggplot(aes(x = Maturity_at_treatment_ordinal, y = yi, size = sqrt(1/vi), group = Sex)) +
  geom_point(aes( col = Sex)) +
  geom_abline(slope = (0.0274 -0.1194),
              intercept = 0.4547, col = "#00BFC4") +
  geom_abline(slope = (0.0274),
              intercept = 0.1396, col = "#F8766D") +
  
  geom_smooth(aes(y = ymin, col = Sex), method = "loess", alpha = 0.2, lty = "dotted", lwd = 0.25,se = F) + 
  
  
  
  # geom_smooth(aes(y = ymin2), method =  "loess", se = FALSE, lty =  "dotted", lwd = 0.25, colour = "#0072B2") +
  # geom_smooth(aes(y = ymax2), method =  "loess", se = FALSE, lty = "dotted", lwd = 0.25, colour = "#0072B2") +
  # geom_smooth(aes(y = ymin), method =  "loess", se = FALSE,lty = "dotted", lwd = 0.25, colour ="#D55E00") +
  # geom_smooth(aes(y = ymax), method =  "loess", se = FALSE, lty ="dotted", lwd = 0.25, colour ="#D55E00") + 
  # geom_smooth(aes(y = pred), method =  "loess", se = FALSE, lty ="dashed", lwd = 0.5, colour ="black") +  
 # ylim(-1, 2) + xlim(0, 1.5) +
  #geom_abline(intercept = mr_host_range_link_ratio$beta[[1]], slope = mr_host_range_link_ratio$beta[[2]], alpha = 0.7, linetype = "dashed", size = 0.5) +
  labs(x = "Maturity", y = "lnRR (effect isze)", size = "Precision") +
  #guides(fill = "none", colour = "none") +
  # themses
  theme_bw() +
  theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
  theme(legend.direction="horizontal") +
  #theme(legend.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.background = element_blank()) +
  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, angle = 90)) 

plot


# gonad
# Gonads_removed
summary(dat$Gonads_removed)
# TODO sex specific analysis - Gondas_removed

mod5 <-  rma.mv(yi, V = V_matrix, mod = ~ Gonads_removed -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod5) 

r2_ml(mod5)

mod5b <-  rma.mv(yi, V = V_matrix, mod = ~ Gonads_removed, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod5b) 

orchard_plot(mod5, mod = "Gonads_removed", xlab = "log response ratio (lnRR)")


# creating Sex + Gonad

dat$Sex_Gonads <- paste(dat$Sex, dat$Gonads_removed, sep = "_")

mod6 <-  rma.mv(yi, V = V_matrix, mod = ~ Sex_Gonads -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod6) 
r2_ml(mod6)


orchard_plot(mod6, mod = "Sex_Gonads", xlab = "log response ratio (lnRR)")



# quality
summary(dat$Controlled_treatments)


mod7 <-  rma.mv(yi, V = V_matrix, mod = ~ Controlled_treatments -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod7) 

orchard_plot(mod7, mod = "Controlled_treatments", xlab = "log response ratio (lnRR)")

mod7b <-  rma.mv(yi, V = V_matrix, mod = ~ Controlled_treatments, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod7b) 

# effect type 
summary(dat$Effect_type)


mod8 <-  rma.mv(yi, V = V_matrix, mod = ~ Effect_type -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), 
                data = dat, 
                test = "t",
                control = list(optimizer="optim", optmethod="BFGS"))
                #control=list(optimizer="optim", optmethod="Nelder-Mead"))
summary(mod8) 

orchard_plot(mod8, mod = "Controlled_treatments", xlab = "log response ratio (lnRR)")

mod8b <-  rma.mv(yi, V = V_matrix, mod = ~ Effect_type, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod8b) 


# sex vs enviroment

dat$Sex_Env <- paste(dat$Sex, dat$Wild_or_semi_wild, sep = "_")

mod9 <-  rma.mv(yi, V = V_matrix, mod = ~ Sex_Env -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod9) 
r2_ml(mod9)

orchard_plot(mod9, mod = "Sex_Env", xlab = "log response ratio (lnRR)")

# Evn vs Controlled

dat$Env_Controlled <- paste(dat$Wild_or_semi_wild, dat$Controlled_treatments, sep = "_")

mod10 <-  rma.mv(yi, V = V_matrix, mod = ~ Env_Controlled -1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod10) 
r2_ml(mod10)

orchard_plot(mod10, mod = "Env_Controlled", xlab = "log response ratio (lnRR)")




#### mulitvaraite ####




##### Condtional analysis #######

# just at means of all continiosu variables (if we do not set anything)
res <- qdrg(object = mod1, data = dat)

# marginal means ; all groups are proportionally weighted
emmeans(res, specs = ~1, df = mod1$ddfs, weights = "prop") 
emmeans(res, specs = "Sex", df = mod1$ddfs, weights = "prop") 

#
dat$Maturity_scaled <- scale(dat$Maturity_at_treatment_ordinal)

mod2 <-  rma.mv(yi, V = V_matrix, mod = ~ Wild_or_semi_wild -1 + Sex + Gonads_removed + Effect_type + Maturity_scaled, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = dat, test = "t")
summary(mod2) 

res2 <- qdrg(object = mod2, data = dat)
# 
#marginal means ; all groups are proportionally weighted
# emmeans(res2, specs = ~1, df = mod1$ddf, weights = "prop") 
emmeans(res2, specs = "Gonads_removed", df = mod2$ddf[[1]], weights = "prop") 
emmeans(res2, specs = "Sex", df = mod2$ddf[[1]], weights = "prop") 


# Gonads_removed emmean     SE df lower.CL upper.CL
# No              0.175 0.0643 75   0.0470    0.303
# Yes             0.076 0.0488 75  -0.0213    0.173


# Sex    emmean     SE df  lower.CL upper.CL
# Female 0.0942 0.0487 75 -0.002708    0.191
# Male   0.0982 0.0496 75 -0.000586    0.197




### longer data####


# TODO we will do hetero model

V_matrix_long <- impute_covariance_matrix(vi = dat_long$vi, cluster = dat_long$Shared_control, r = 0.5)

# we can run - some heteroscad models
# this does not improve model
mod_n <-  rma.mv(yi, V = V_matrix_long, mod = ~ Comp_type - 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID, ~1|Obs ), R = list(Phylogeny = cor_tree), data = dat_long, test = "t")
summary(mod_n) 

orchard_plot(mod_n, mod = "Comp_type", xlab = "log response ratio (lnRR)")

# mod_h <-  rma.mv(yi, V = V_matrix_long, mod = ~ Comp_type - 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~Comp_type|Study, ~1|Effect_ID, ~Comp_type|Obs), rho = 0, struct = "HCS", R = list(Phylogeny = cor_tree), data = dat_long, test = "t")
# summary(mod_h) 

# mod_nb <-  rma.mv(yi, V = V_matrix_long, mod = ~ Comp_type, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID, ~1|Obs ), R = list(Phylogeny = cor_tree), data = dat_long, test = "t")
# summary(mod_nb) 

# TODO - using V_matrix_long

# modeling both type and sex 

mod_n2 <-  rma.mv(yi, V = vi, mod = ~ Comp_type*Sex , random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID, ~1|Obs ), R = list(Phylogeny = cor_tree), data = dat_long, test = "t")
summary(mod_n2) 


# mod_h2 <-  rma.mv(yi, V = vi, mod = ~Comp_type*Sex, random = list(~1|Phylogeny, ~1|Species_Latin, ~Comp_type_Sex|Study, ~1|Effect_ID, ~Comp_type_Sex|Obs), rho = 0, struct = "HCS", R = list(Phylogeny = cor_tree), data = dat_long, test = "t")
# summary(mod_h2) 


# as interaction

dat_long$Comp_type_Sex <- paste(dat_long$Comp_type, dat_long$Sex, sep = "_")

mod_n3 <-  rma.mv(yi, V = vi, mod = ~ Comp_type_Sex -1 , random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID, ~1|Obs ), R = list(Phylogeny = cor_tree), data = dat_long, test = "t")
summary(mod_n3) 

orchard_plot(mod_n3, mod = "Comp_type_Sex", xlab = "log response ratio (lnRR)")

mod_h3 <-  rma.mv(yi, V = vi, mod = ~Comp_type_Sex -1 , random = list(~1|Phylogeny, ~1|Species_Latin, ~Comp_type_Sex|Study, ~1|Effect_ID, ~Comp_type_Sex|Obs), rho = 0, struct = "HCS", R = list(Phylogeny = cor_tree), data = dat_long, test = "t")
summary(mod_h3) 

# both group separate 

mod_neut <-  rma.mv(yi, V = vi, 
                   mod = ~ 1, 
                   random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), 
                   R = list(Phylogeny = cor_tree), 
                   test = "t",
                   data = subset(dat_long, Comp_type=="one_neutured"))


mod_norm <-  rma.mv(yi, V = vi, 
                   mod = ~ 1, 
                   random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), 
                   R = list(Phylogeny = cor_tree), 
                   test = "t",
                   data = subset(dat_long, Comp_type== "both_normal"))


summary(mod_neut) 
summary(mod_norm) 

round(i2_ml(mod_neut)*100,2)
round(i2_ml(mod_norm)*100, 2)

# for both sexes running separate
# female nuked

mod_neut_f <-  rma.mv(yi, V = vi, 
                    mod = ~ 1, 
                    random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), 
                    R = list(Phylogeny = cor_tree), 
                    test = "t",
                    data = subset(dat_long, Comp_type=="one_neutured" & Sex== "Female"))


mod_norm_f <-  rma.mv(yi, V = vi, 
                    mod = ~ 1, 
                    random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), 
                    R = list(Phylogeny = cor_tree), 
                    test = "t",
                    data = subset(dat_long, Comp_type== "both_normal" & Sex== "Female"))


summary(mod_neut_f) 
summary(mod_norm_f) 

round(i2_ml(mod_neut_f)*100,2)
round(i2_ml(mod_norm_f)*100, 2)

# make nuked

mod_neut_m <-  rma.mv(yi, V = vi, 
                      mod = ~ 1, 
                      random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), 
                      R = list(Phylogeny = cor_tree), 
                      test = "t",
                      data = subset(dat_long, Comp_type=="one_neutured" & Sex== "Male"))


mod_norm_m <-  rma.mv(yi, V = vi, 
                      mod = ~ 1, 
                      random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), 
                      R = list(Phylogeny = cor_tree), 
                      test = "t",
                      data = subset(dat_long, Comp_type== "both_normal" & Sex== "Male"))


summary(mod_neut_m) 
summary(mod_norm_m) 

round(i2_ml(mod_neut_m)*100,2)
round(i2_ml(mod_norm_m)*100, 2)


# #female nuked
# mod_nf <-  rma.mv(yi, V = vi, mod = ~ Comp_type - 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID, ~1|Obs ), R = list(Phylogeny = cor_tree), data = subset(dat_long, Sex== "Female"), test = "t")
# summary(mod_nf) 
# 
# mod_hf <-  rma.mv(yi, V = vi, mod = ~ Comp_type - 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID, ~Comp_type|Obs), rho = 0, struct = "HCS", R = list(Phylogeny = cor_tree), data = subset(dat_long, Sex== "Female"), test = "t")
# summary(mod_hf) 
# 
# #male nuked
# mod_nm <-  rma.mv(yi, V = vi, mod = ~ Comp_type - 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID, ~1|Obs ), R = list(Phylogeny = cor_tree), data = subset(dat_long, Sex== "Male"), test = "t")
# summary(mod_nm) 
# 
# mod_hm <-  rma.mv(yi, V = vi, mod = ~ Comp_type - 1, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID, ~Comp_type|Obs), rho = 0, struct = "HCS", R = list(Phylogeny = cor_tree), data = subset(dat_long, Sex== "Male"), test = "t")
# summary(mod_hm) 
# 
# 
# orchard_plot(mod_n, mod = "Comp_type", xlab = "log response ratio (lnRR)")
# 
# orchard_plot(mod_nf, mod = "Comp_type", xlab = "log response ratio (lnRR)")
# 
# orchard_plot(mod_nm, mod = "Comp_type", xlab = "log response ratio (lnRR)")
# 
# # separately 
# mod_h1 <-  rma.mv(yi, V = vi,  random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = subset(dat_long, Comp_type== "one_neutured"),  test = "t")
# summary(mod_h1) 
# 
# round(i2_ml(mod_h1)*100,2)
# 
# mod_h2 <-  rma.mv(yi, V = vi, random = list(~1|Phylogeny, ~1|Species_Latin, ~1|Study, ~1|Effect_ID), R = list(Phylogeny = cor_tree), data = subset(dat_long, Comp_type== "both_normal"), test = "t")
# summary(mod_h2) 
# 
# round(i2_ml(mod_h2)*100,2)
##########################
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
