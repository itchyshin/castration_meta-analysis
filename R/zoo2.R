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

mod_fm_ns <- rma.mv(yi, V = vi, 
                     random = list(
                       ~1|species,
                       ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat_fm_ns)
summary(mod_fm_ns)
i2_ml(mod_fm_ns)

orchard_plot(mod_fm_ns, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

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

mod_fm_sn <- rma.mv(yi, V = vi, 
                     random = list(
                       ~1|species,
                       ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat_fm_sn)
summary(mod_fm_sn)
i2_ml(mod_fm_sn)

orchard_plot(mod_fm_sn, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

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

mod_fm_hn <- rma.mv(yi, V = vi, 
                     random = list(
                       ~1|species,
                       ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat_fm_hn)
summary(mod_fm_hn)
i2_ml(mod_fm_hn)

orchard_plot(mod_fm_hn, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

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
              n2i = M_normal_m,
              data = dat_fm_sn,
              )

# modeling

mod_fm_sn <- rma.mv(yi, V = vi, 
                     random = list(
                       ~1|species,
                       ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat_fm_sn)
summary(mod_fm_sn)
i2_ml(mod_fm_sn)

orchard_plot(mod_fm_sn, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

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

mod_fm_hs <- rma.mv(yi, V = vi, 
                     random = list(
                       ~1|species,
                       ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat_fm_hs)
summary(mod_fm_hs)
i2_ml(mod_fm_hs)

orchard_plot(mod_fm_hs, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

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

# modeling

mod_fm_ss <- rma.mv(yi, V = vi, 
                                         random = list(
                                             ~1|species,
                                             ~1|phylogeny), 
                                         R = list(phylogeny = cor_tree), 
                                         data = dat_fm_ss)
summary(mod_fm_ss)
i2_ml(mod_fm_ss)

orchard_plot(mod_fm_ss, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

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

mod_fm_nn <- rma.mv(yi, V = vi, 
                                         random = list(
                                             ~1|species,
                                             ~1|phylogeny), 
                                         R = list(phylogeny = cor_tree), 
                                         data = dat_fm_nn)
summary(mod_fm_nn)
i2_ml(mod_fm_nn)

orchard_plot(mod_fm_nn, xlab = "lnRR (male surgical)", group = "species", g = FALSE)


# merging all data



rbind(
  dat_fm_ns[ , 7:12], # 1
  dat_fm_sn[ , 7:12], # 2
  dat_fm_hn[ , 7:12], # 3 
  dat_fm_hs[ , 7:12], # 4
  dat_fm_ss[ , 7:12], # 5
  dat_fm_nn[ , 7:12]  # 6
) -> dat_comb


dat_comb$obs_id <- factor(1:nrow(dat_comb))

# reordering cateogry

dat_comb$category <- factor(dat_comb$category, levels = c( "F surgical/M surgical",
                                                           "F hormonal/M surgical",  
                                                           "F surgical/M normal",
                                                           "F hormonal/M normal",
                                                           "F normal/M surgical",
                                                           "F normal/M normal"))

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

orchard_plot(mod_comb, mod = "category",
            xlab = "lnRR (all)", group = "vertlife.species", 
            g = FALSE, angle = 45)

# creating a table to show contrasts between categories

############
# wriring a function to get the contrasts
############



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

dat_comb <- dat_comb %>% mutate(
                      abs_yi = abs(yi),
                      abs_yi2 = folded_mu(yi, vi), 
                      abs_vi = folded_v(yi, vi))

dat_comb[which(dat_comb$abs_yi == max(dat_comb$abs_yi)), ]

hist(log(dat_comb$abs_vi))

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

robust(mod_comb_a, cluster = vertlife.species)  

orchard_plot(mod_comb_a, mod = "category",
            xlab = "lnRR (all)", group = "vertlife.species", 
            g = FALSE, angle = 45)
