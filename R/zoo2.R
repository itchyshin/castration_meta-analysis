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
               here,
               cowplot
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



# female normal vs. male surgical

dat_fm_ns <- dat %>% filter(is.na(Female_None_Mean) == FALSE, 
                            is.na(Male_Surgical_Mean) == FALSE) %>% 
  transmute(F_normal_m = Male_None_Mean,
         F_normal_sd = sqrt(Female.None)*Female_None_SE,
         F_normal_n = Female.None,
         M_surgical_m = Male_Surgical_Mean,
         M_surgical_sd = sqrt(Male.Surgical)*Male_Surgical_SE,
         M_surgical_n = Male.Surgical,
         species = species,
         vertlife.species = vertlife.species,
         phylogeny = phylogeny,
         category = "F normal/M surgical"
         )

# getting effect size

dat_fm_ns <- escalc("ROM", 
              m1i = F_normal_m,
              m2i = M_surgical_m,
              sd1i = F_normal_sd,
              sd2i = M_surgical_sd,
              n1i = F_normal_n,
              n2i = M_surgical_n,
              data = dat_fm_ns,
              )

# modeling
# 
# mod_fm_ns <- rma.mv(yi, V = vi, 
#                      random = list(
#                        ~1|species,
#                        ~1|phylogeny), 
#                      R = list(phylogeny = cor_tree), 
#                      data = dat_fm_ns)
# summary(mod_fm_ns)
# i2_ml(mod_fm_ns)
# 
# orchard_plot(mod_fm_ns, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

# female surgical vs. male normal

dat_fm_sn <- dat %>% filter(is.na(Female_Surgical_Mean) == FALSE,
                            is.na(Male_None_Mean) == FALSE) %>%
  transmute(F_surgical_m = Female_Surgical_Mean,
         F_surgical_sd = sqrt(Female.Surgical*Female_Surgical_SE),
         F_surgical_n = Female.Surgical,
         M_normal_m = Male_None_Mean,
         M_normal_sd = sqrt(Male.None)*Male_None_SE,
         M_normal_n = Male.None,
         species = species,
         vertlife.species = vertlife.species,
         phylogeny = phylogeny,
         category = "F hormonal/M normal"
         )

# getting effect size         
dat_fm_sn <- escalc("ROM", 
              m1i = F_surgical_m,
              m2i = M_normal_m,
              sd1i = F_surgical_sd,
              sd2i = M_normal_sd,
              n1i = F_surgical_n,
              n2i = M_normal_n,
              data = dat_fm_sn,
              )

# modeling

# mod_fm_sn <- rma.mv(yi, V = vi, 
#                      random = list(
#                        ~1|species,
#                        ~1|phylogeny), 
#                      R = list(phylogeny = cor_tree), 
#                      data = dat_fm_sn)
# summary(mod_fm_sn)
# i2_ml(mod_fm_sn)
# 
# orchard_plot(mod_fm_sn, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

# female hormonal vs. male normal

dat_fm_hn <- dat %>% filter(is.na(Female_Hormonal_Mean) == FALSE, 
                            is.na(Male_None_Mean) == FALSE) %>% 
  transmute(F_hormonal_m = Female_Hormonal_Mean,
         F_hormonal_sd = sqrt(Female.Hormonal*Female_Hormonal_SE),
         F_hormonal_n = Female.Hormonal,
         M_normal_m = Male_None_Mean,
         M_normal_sd = sqrt(Male.None)*Male_None_SE,
         M_normal_n = Male.None,
         species = species,
         vertlife.species = vertlife.species,
         phylogeny = phylogeny,
         category = "F hormonal/M normal"
         )

# getting effect size         
dat_fm_hn <- escalc("ROM", 
              m1i = F_hormonal_m,
              m2i = M_normal_m,
              sd1i = F_hormonal_sd,
              sd2i = M_normal_sd,
              n1i = F_hormonal_n,
              n2i = M_normal_m,
              data = dat_fm_hn,
              )

# modeling
# 
# mod_fm_hn <- rma.mv(yi, V = vi, 
#                      random = list(
#                        ~1|species,
#                        ~1|phylogeny), 
#                      R = list(phylogeny = cor_tree), 
#                      data = dat_fm_hn)
# summary(mod_fm_hn)
# i2_ml(mod_fm_hn)
# 
# orchard_plot(mod_fm_hn, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

# female surgical vs. male normal

dat_fm_sn <- dat %>% filter(is.na(Female_Surgical_Mean) == FALSE,
                            is.na(Male_None_Mean) == FALSE) %>%
  transmute(F_surgical_m = Female_Surgical_Mean,
         F_surgical_sd = sqrt(Female.Surgical*Female_Surgical_SE),
         F_surgical_n = Female.Surgical,
         M_normal_m = Male_None_Mean,
         M_normal_sd = sqrt(Male.None)*Male_None_SE,
         M_normal_n = Male.None,
         species = species,
         vertlife.species = vertlife.species,
         phylogeny = phylogeny,
         category = "F surgical/M normal"
         )

# getting effect size         
dat_fm_sn <- escalc("ROM", 
              m1i = F_surgical_m,
              m2i = M_normal_m,
              sd1i = F_surgical_sd,
              sd2i = M_normal_sd,
              n1i = F_surgical_n,
              n2i = M_normal_m,
              data = dat_fm_sn,
              )

# modeling
# 
# mod_fm_sn <- rma.mv(yi, V = vi, 
#                      random = list(
#                        ~1|species,
#                        ~1|phylogeny), 
#                      R = list(phylogeny = cor_tree), 
#                      data = dat_fm_sn)
# summary(mod_fm_sn)
# i2_ml(mod_fm_sn)
# 
# orchard_plot(mod_fm_sn, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

# female hormonal vs. male surgical

dat_fm_hs <- dat %>% filter(is.na(Female_Hormonal_Mean) == FALSE, 
                            is.na(Male_Surgical_Mean) == FALSE) %>% 
  transmute(F_hormonal_m = Female_Hormonal_Mean,
         F_hormonal_sd = sqrt(Female.Hormonal*Female_Hormonal_SE),
         F_hormonal_n = Female.Hormonal,
         M_surgical_m = Male_Surgical_Mean,
         M_surgical_sd = sqrt(Male.Surgical)*Male_Surgical_SE,
         M_surgical_n = Male.Surgical,
         species = species,
         vertlife.species = vertlife.species,
         phylogeny = phylogeny,
         category = "F hormonal/M surgical"
         )

# getting effect size         
dat_fm_hs <- escalc("ROM", 
              m1i = F_hormonal_m,
              m2i = M_surgical_m,
              sd1i = F_hormonal_sd,
              sd2i = M_surgical_sd,
              n1i = F_hormonal_n,
              n2i = M_surgical_n,
              data = dat_fm_hs,
              )

# modeling

# mod_fm_hs <- rma.mv(yi, V = vi, 
#                      random = list(
#                        ~1|species,
#                        ~1|phylogeny), 
#                      R = list(phylogeny = cor_tree), 
#                      data = dat_fm_hs)
# summary(mod_fm_hs)
# i2_ml(mod_fm_hs)
# 
# orchard_plot(mod_fm_hs, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

# female surgical vs. male surgical

dat_fm_ss <- dat %>% filter(is.na(Female_Surgical_Mean) == FALSE, 
                                                        is.na(Male_Surgical_Mean) == FALSE) %>% 
    transmute(F_surgical_m = Female_Surgical_Mean,
                 F_surgical_sd = sqrt(Female.Surgical*Female_Surgical_SE),
                 F_surgical_n = Female.Surgical,
                 M_surgical_m = Male_Surgical_Mean,
                 M_surgical_sd = sqrt(Male.Surgical)*Male_Surgical_SE,
                 M_surgical_n = Male.Surgical,
                 species = species,
                 vertlife.species = vertlife.species,
                 phylogeny = phylogeny,
                 category = "F surgical/M surgical"
                 )

# getting effect size         
dat_fm_ss <- escalc("ROM", 
                            m1i = F_surgical_m,
                            m2i = M_surgical_m,
                            sd1i = F_surgical_sd,
                            sd2i = M_surgical_sd,
                            n1i = F_surgical_n,
                            n2i = M_surgical_n,
                            data = dat_fm_ss,
                            )

# # modeling
# 
# mod_fm_ss <- rma.mv(yi, V = vi, 
#                                          random = list(
#                                              ~1|species,
#                                              ~1|phylogeny), 
#                                          R = list(phylogeny = cor_tree), 
#                                          data = dat_fm_ss)
# summary(mod_fm_ss)
# i2_ml(mod_fm_ss)
# 
# orchard_plot(mod_fm_ss, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

# female normal vs. male normal

dat_fm_nn <- dat %>% filter(is.na(Female_None_Mean) == FALSE, 
                                                        is.na(Male_None_Mean) == FALSE) %>% 
    transmute(F_normal_m = Female_None_Mean,
                 F_normal_sd = sqrt(Female.None*Female_None_SE),
                 F_normal_n = Female.None,
                 M_normal_m = Male_None_Mean,
                 M_normal_sd = sqrt(Male.None*Male_None_SE),
                 M_normal_n = Male.None,
                 species = species,
                 vertlife.species = vertlife.species,
                 phylogeny = phylogeny,
                 category = "F normal/M normal"
                 )

# getting effect size         
dat_fm_nn <- escalc("ROM", 
                            m1i = F_normal_m,
                            m2i = M_normal_m,
                            sd1i = F_normal_sd,
                            sd2i = M_normal_sd,
                            n1i = F_normal_n,
                            n2i = M_normal_n,
                            data = dat_fm_nn,
                            )

# modeling

# mod_fm_nn <- rma.mv(yi, V = vi, 
#                                          random = list(
#                                              ~1|species,
#                                              ~1|phylogeny), 
#                                          R = list(phylogeny = cor_tree), 
#                                          data = dat_fm_nn)
# summary(mod_fm_nn)
# i2_ml(mod_fm_nn)
# 
# orchard_plot(mod_fm_nn, xlab = "lnRR (male surgical)", group = "species", g = FALSE)
# 
# 
# # merging all data
# 
# 
# 
# rbind(
#   dat_fm_ns[ , 7:12], # 1
#   dat_fm_hn[ , 7:12], # 2 
#   dat_fm_sn[ , 7:12], # 3
#   dat_fm_hs[ , 7:12], # 4
#   dat_fm_ss[ , 7:12], # 5
#   dat_fm_nn[ , 7:12]  # 6
# ) -> dat_comb
# 
# 
# dat_comb$obs_id <- factor(1:nrow(dat_comb))
# 
# # reordering cateogry
# 
# dat_comb$category <- factor(dat_comb$category, levels = c( "F surgical/M surgical",
#                                                            "F hormonal/M surgical",  
#                                                            "F surgical/M normal",
#                                                            "F hormonal/M normal",
#                                                            "F normal/M surgical",
#                                                            "F normal/M normal"),
#                             )
# 
# 
# 
# VCV <- vcalc(vi, vertlife.species, rho = 0.5, data = dat_comb)
# 
# 
# mod_comb <- rma.mv(yi, V = VCV,
#                    mods = ~category - 1,
#                      random = list(
#                        ~1|vertlife.species,
#                        ~1|phylogeny,
#                        ~1|obs_id),
#                      R = list(phylogeny = cor_tree),
#                      data = dat_comb,
#                      control = list(optimizer = "Nelder-Mead"))
# summary(mod_comb)
# r2_ml(mod_comb)
# 
# robust(mod_comb, cluster = vertlife.species)  
# 
# orchard_plot(mod_comb, mod = "category",
#             xlab = "lnRR (all)", group = "vertlife.species", 
#             g = FALSE, angle = 45)
# 
# res <- all_models(mod_comb, mod = "category")  

# creating a table to show contrasts between categories

############
# wriring a function to get the contrasts
############

# custom functions

#' Title: Contrast name generator
#'
#' @param name: a vector of character strings
cont_gen <- function(name) {
  combination <- combn(name, 2)
  name_dat <- t(combination)
  names <- paste(name_dat[, 1], name_dat[, 2], sep = "-")
  return(names)
}

#' @title get_pred1: intercept-less model
#' @description Function to get CIs (confidence intervals) and PIs (prediction intervals) from rma objects (metafor)
#' @param model: rma.mv object 
#' @param mod: the name of a moderator 
get_pred1 <- function (model, mod = " ") {
  name <- firstup(as.character(stringr::str_replace(row.names(model$beta), mod, "")))
  len <- length(name)
  
   if (len != 1) {
        newdata <- diag(len)
        pred <- metafor::predict.rma(model, 
                                     newmods = newdata,
                                     tau2.levels = 1:len)
    }
    else {
        pred <- metafor::predict.rma(model)
  }
  estimate <- pred$pred
  lowerCL <- pred$ci.lb
  upperCL <- pred$ci.ub 
  lowerPR <- pred$cr.lb
  upperPR <- pred$cr.ub 
  
  table <- tibble(name = factor(name, levels = name, labels = name), estimate = estimate,
                  lowerCL = lowerCL, upperCL = upperCL,
                  pval = model$pval,
                  lowerPR = lowerPR, upperPR = upperPR)
}

#' @title get_pred2: normal model
#' @description Function to get CIs (confidence intervals) and PIs (prediction intervals) from rma objects (metafor)
#' @param model: rma.mv object 
#' @param mod: the name of a moderator 
get_pred2 <- function (model, mod = " ") {
  name <- as.factor(str_replace(row.names(model$beta), 
                                paste0("relevel", "\\(", mod,", ref = name","\\)"),""))
  len <- length(name)
  
  if(len != 1){
  newdata <- diag(len)
  pred <- predict.rma(model, intercept = FALSE, newmods = newdata[ ,-1])
  }
  else {
    pred <- predict.rma(model)
  }
  estimate <- pred$pred
  lowerCL <- pred$ci.lb
  upperCL <- pred$ci.ub 
  lowerPR <- pred$cr.lb
  upperPR <- pred$cr.ub 
  
  table <- tibble(name = factor(name, levels = name, labels = name), estimate = estimate,
                  lowerCL = lowerCL, upperCL = upperCL,
                  pval = model$pval,
                  lowerPR = lowerPR, upperPR = upperPR)
}

#' @title mr_results
#' @description Function to put results of meta-regression and its contrasts
#' @param res1: data frame 1
#' @param res1: data frame 2
mr_results <- function(res1, res2) {
  restuls <-tibble(
    `Fixed effect` = c(as.character(res1$name), cont_gen(res1$name)),
    Estimate = c(res1$estimate, res2$estimate),
    `Lower CI [0.025]` = c(res1$lowerCL, res2$lowerCL),
    `Upper CI  [0.975]` = c(res1$upperCL, res2$upperCL),
    `P value` = c(res1$pval, res2$pval),
    `Lower PI [0.025]` = c(res1$lowerPR, res2$lowerPR),
    `Upper PI  [0.975]` = c(res1$upperPR, res2$upperPR),
  )
}


#' @title all_models
#' @description Function to take all possible models and get their results
#' @param model: intercept-less model
#' @param mod: the name of a moderator 

all_models <- function(model, mod = " ", type = "homo") {
  
  # getting the level names out
  level_names <- levels(factor(model$data[[mod]]))
  dat2 <- model$data
  mod <- mod


  #model$data[[mod]] <- factor(model$data[[mod]], ordered = FALSE)
  # meta-regression: contrasts 
  # helper function to run metafor meta-regression
  run_rma1 <- function(name) {
      VCV1 <- vcalc(vi = dat2$vi,
             cluster = dat2$vertlife.species,
             rho = 0.5)
    rma.mv(yi, V = VCV1,
                   mods = ~relevel(dat2[[mod]], ref = name),
                     random = list(
                       ~1|vertlife.species,
                       ~1|phylogeny,
                       ~1|obs_id),
                     R = list(phylogeny = cor_tree),
                     data = dat2,
                     control = list(optimizer = "Nelder-Mead"))
   }

    run_rma2 <- function(name) {
    
            VCVa <- vcalc(abs_vi, vertlife.species, 
                    rho = 0.5, data = dat2)
               
               rma.mv(abs_yi, V = VCVa,
               mods = ~relevel(dat2[[mod]], ref = name),
                     random = list(
                       ~1|vertlife.species,
                       ~1|phylogeny,
                       ~1|obs_id),
                     R = list(phylogeny = cor_tree),
                     data = dat2,
                     control = list(optimizer = "Nelder-Mead"))
   }

# results of meta-regression including all contrast results; taking the last level out ([-length(level_names)])
# this does not work for hetero model?
if (type == "homo"){

    model_all <- purrr::map(level_names[-length(level_names)], run_rma1)

  } else {
  model_all <- purrr::map(level_names[-length(level_names)], run_rma2)
  }
  
  # getting estimates from intercept-less models (means for all the groups)
  res1 <- get_pred1(model, mod = mod)
  
  # getting estiamtes from all contrast models
  res2_pre <- purrr::map(model_all, ~ get_pred2(.x, mod = mod))
  
  # a list of the numbers to take out unnecessary contrasts
  contra_list <- Map(seq, from=1, to=1:(length(level_names) - 1))
  res2 <- purrr::map2_dfr(res2_pre, contra_list, ~.x[-(.y), ]) 
  # creating a table
  res_tab <- mr_results(res1, res2) %>% 
  kable("html",  digits = 3) %>%
  kable_styling("striped", position = "left") %>%
  scroll_box(width = "100%")
  
  # results
  res_tab

}


##########
# functions for absolute values


# folded mean
folded_mu <-function(mean, variance){
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_mu
} 

# folded variance
folded_v <-function(mean, variance){
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_se <- sqrt(mu^2 + sigma^2 - fold_mu^2)
  # adding se to make bigger mean
  fold_v <-fold_se^2
  fold_v
} 

# Aboslute analyses
# 
# dat_comb <- dat_comb %>% mutate(
#                       abs_yi = abs(yi),
#                       abs_yi2 = folded_mu(yi, vi), 
#                       abs_vi = folded_v(yi, vi))
# 
# dat_comb[which(dat_comb$abs_yi == max(dat_comb$abs_yi)), ]

#hist(log(dat_comb$abs_vi))


# modeling

# VCVa <- vcalc(abs_vi, vertlife.species, rho = 0.5, data = dat_comb)
# 
# mod_comb_a <- rma.mv(abs_yi, V = VCVa,
#                    mods = ~category - 1,
#                      random = list(
#                        ~1|vertlife.species,
#                        ~1|phylogeny,
#                        ~1|obs_id),
#                      R = list(phylogeny = cor_tree),
#                      data = dat_comb,
#                      control = list(optimizer = "Nelder-Mead"))
# summary(mod_comb_a)
# r2_ml(mod_comb_a)
# res2 <- all_models(mod_comb_a, mod = "category", type = "abs")
# 
# robust(mod_comb_a, cluster = vertlife.species)  
# 
# orchard_plot(mod_comb_a, mod = "category",
#             xlab = "lnRR (all)", group = "vertlife.species", 
#             g = FALSE, angle = 45)

##########################
# adding more conditions 
###########################

# female immunological vs. male normal

dat_fm_in <- dat %>% filter(is.na(Female_Immunological_Mean) == FALSE,
                            is.na(Male_None_Mean) == FALSE) %>%
  transmute(F_immunological_m = Female_Immunological_Mean,
            F_immunological_sd = sqrt(Female.Immunological*Female_Immunological_SE),
            F_immunological_n = Female.Immunological,
            M_normal_m = Male_None_Mean,
            M_normal_sd = sqrt(Male.None)*Male_None_SE,
            M_normal_n = Male.None,
            species = species,
            vertlife.species = vertlife.species,
            phylogeny = phylogeny,
            category = "F immunological/M normal"
  )

# getting effect size         
dat_fm_in <- escalc("ROM", 
                    m1i = F_immunological_m,
                    m2i = M_normal_m,
                    sd1i = F_immunological_sd,
                    sd2i = M_normal_sd,
                    n1i = F_immunological_n,
                    n2i = M_normal_n,
                    data = dat_fm_in,
)


# female normal vs. male hormonal

dat_fm_nh <- dat %>% filter(is.na(Female_None_Mean) == FALSE,
                            is.na(Male_Hormonal_Mean) == FALSE) %>%
  transmute(F_normal_m = Female_None_Mean,
            F_normal_sd = sqrt(Female.None*Female_None_SE),
            F_normal_n = Female.None,
            M_hornomal_m = Male_Hormonal_Mean,
            M_hornomal_sd = sqrt(Male.Hormonal)*Male_Hormonal_SE,
            M_hornomal_n = Male.Hormonal,
            species = species,
            vertlife.species = vertlife.species,
            phylogeny = phylogeny,
            category = "F normal/M hornomal"
  )

# getting effect size         
dat_fm_nh <- escalc("ROM", 
                    m1i = F_normal_m,
                    m2i = M_hornomal_m,
                    sd1i = F_normal_sd,
                    sd2i = M_hornomal_sd,
                    n1i = F_normal_n,
                    n2i = M_hornomal_n,
                    data = dat_fm_nh,
)

# female normal vs. male immunological

dat_fm_ni <- dat %>% filter(is.na(Female_None_Mean) == FALSE,
                            is.na(Male_Immunological_Mean) == FALSE) %>%
  transmute(F_normal_m = Female_None_Mean,
            F_normal_sd = sqrt(Female.None*Female_None_SE),
            F_normal_n = Female.None,
            M_immunological_m = Male_Immunological_Mean,
            M_immunological_sd = sqrt(Male.Immunological)*Male_Immunological_SE,
            M_immunological_n = Male.Immunological,
            species = species,
            vertlife.species = vertlife.species,
            phylogeny = phylogeny,
            category = "F normal/M immunological"
  )

# getting effect size         
dat_fm_ni <- escalc("ROM", 
                    m1i = F_normal_m,
                    m2i = M_immunological_m,
                    sd1i = F_normal_sd,
                    sd2i = M_immunological_sd,
                    n1i = F_normal_n,
                    n2i = M_immunological_n,
                    data = dat_fm_ni,
)

# female hormonal vs. male hormonal

dat_fm_hh <- dat %>% filter(is.na(Female_Hormonal_Mean) == FALSE,
                            is.na(Male_Hormonal_Mean) == FALSE) %>%
  transmute(F_hornomal_m = Female_Hormonal_Mean,
            F_hornomal_sd = sqrt(Female.Hormonal*Female_Hormonal_SE),
            F_hornomal_n = Female.Hormonal,
            M_hornomal_m = Male_Hormonal_Mean,
            M_hornomal_sd = sqrt(Male.Hormonal)*Male_Hormonal_SE,
            M_hornomal_n = Male.Hormonal,
            species = species,
            vertlife.species = vertlife.species,
            phylogeny = phylogeny,
            category = "F hornomal/M hornomal"
  )

# getting effect size         
dat_fm_hh <- escalc("ROM", 
                    m1i = F_hornomal_m,
                    m2i = M_hornomal_m,
                    sd1i = F_hornomal_sd,
                    sd2i = M_hornomal_sd,
                    n1i = F_hornomal_n,
                    n2i = M_hornomal_n,
                    data = dat_fm_hh,
)


# female immunological vs. male immunological

dat_fm_ii <- dat %>% filter(is.na(Female_Immunological_Mean) == FALSE,
                    is.na(Male_Immunological_Mean) == FALSE) %>%
  transmute(F_immunological_m = Female_Immunological_Mean,
            F_immunological_sd = sqrt(Female.Immunological*Female_Immunological_SE),
            F_immunological_n = Female.None,
            M_immunological_m = Male_Immunological_Mean,
            M_immunological_sd = sqrt(Male.Immunological)*Male_Immunological_SE,
            M_immunological_n = Male.Immunological,
            species = species,
            vertlife.species = vertlife.species,
            phylogeny = phylogeny,
            category = "F immunological/M immunological"
  )

# dat_fm_ii effect size         
dat_fm_ii <- escalc("ROM", 
                    m1i = F_immunological_m,
                    m2i = M_immunological_m,
                    sd1i = F_immunological_sd,
                    sd2i = M_immunological_sd,
                    n1i = F_immunological_n,
                    n2i = M_immunological_n,
                    data = dat_fm_ii,
)


#

rbind(
  dat_fm_ns[ , 7:12], # 1
  dat_fm_hn[ , 7:12], # 2 
  dat_fm_sn[ , 7:12], # 3
#  dat_fm_hs[ , 7:12], # 4
  dat_fm_ss[ , 7:12], # 4
  dat_fm_nn[ , 7:12],  # 5
  dat_fm_nh[ , 7:12],  # 6
  dat_fm_ni[ , 7:12],  # 7
  dat_fm_hh[ , 7:12],  # 8
  dat_fm_ii[ , 7:12]  # 9
) -> dat_comb


dat_comb$obs_id <- factor(1:nrow(dat_comb))

# reordering cateogry

dat_comb$category <- factor(dat_comb$category, levels = c("F immunological/M immunological",
                                                           "F hornomal/M hornomal",
                                                           "F surgical/M surgical",
                                                           "F immunological/M normal",
                                                           "F hormonal/M normal",  
                                                           "F surgical/M normal",
                                                           "F normal/M immunological",
                                                           "F normal/M hormonal",
                                                           "F normal/M surgical",
                                                           "F normal/M normal"),
                                                labels = c("F sterlized/M sterlized",
                                                           "F sterlized/M sterlized",
                                                           "F sterlized/M sterlized",
                                                           "F sterlized/M normal",
                                                           "F sterlized/M normal",  
                                                           "F sterlized/M normal",
                                                           "F normal/M sterlized",
                                                           "F normal/M sterlized",
                                                           "F normal/M sterlized",
                                                           "F normal/M normal"),
                             
)


VCV <- vcalc(vi, vertlife.species, rho = 0.5, data = dat_comb)


mod_comb <- rma.mv(yi, V = VCV,
                   mods = ~category - 1,
                   random = list(
                     ~1|vertlife.species,
                     ~1|phylogeny,
                     ~1|obs_id),
                   R = list(phylogeny = cor_tree),
                   data = dat_comb,
                   control = list(optimizer = "Nelder-Mead"))
summary(mod_comb)
r2_ml(mod_comb)

robust(mod_comb, cluster = vertlife.species)  

p_sexdiff <- orchard_plot(mod_comb, mod = "category",
             xlab = "log response ratio (lnRR)", group = "vertlife.species", 
             g = FALSE, angle = 45)

res <- all_models(mod_comb, mod = "category")  


####### 
# absolute analyses

# Aboslute analyses

dat_comb <- dat_comb %>% mutate(
  abs_yi = abs(yi),
  abs_yi2 = folded_mu(yi, vi), 
  abs_vi = folded_v(yi, vi))

dat_comb[which(dat_comb$abs_yi == max(dat_comb$abs_yi)), ]

#hist(log(dat_comb$abs_vi))


# modeling

VCVa <- vcalc(abs_vi, vertlife.species, rho = 0.5, data = dat_comb)

mod_comb_a <- rma.mv(abs_yi, V = VCVa,
                     mods = ~category - 1,
                     random = list(
                       ~1|vertlife.species,
                       ~1|phylogeny,
                       ~1|obs_id),
                     R = list(phylogeny = cor_tree),
                     data = dat_comb,
                     control = list(optimizer = "Nelder-Mead"))
summary(mod_comb_a)
r2_ml(mod_comb_a)
#res2 <- all_models(mod_comb_a, mod = "category", type = "abs")

robust(mod_comb_a, cluster = vertlife.species)  

p_abs <- orchard_plot(mod_comb_a, mod = "category",
             xlab = "absolute log response ratio (lnRR)", group = "vertlife.species", 
             g = FALSE, angle = 45)

#res2

# putting figs together


###################
## Figure 4
###################

library(cowplot)
plot_grid(p_phylo2) / (p_sexdiff + p_abs) + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(heights = c(1.5,1.0))


