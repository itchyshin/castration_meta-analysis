#zoo 3

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
               apextra,
               cowplot
)

########################
# data: life expectancy
########################

dat_pre <-  read.csv(here("data", "hormTabPre05.csv"))
dat_post <-  read.csv(here("data", "hormTabPost05.csv"))

# 
dim(dat_pre)
dim(dat_post)

# extra
tax <- read.csv(here("data", "vertlife_taxonomy_translation_table.csv"))

# animal order and tree relationships 
#dford <- readRDS(here("Rdata", "dford.RDS"))

dat_pre %>% left_join(tax, by = c("species" = "zims.species")) -> dat_pre
dat_post %>% left_join(tax, by = c("species" = "zims.species")) -> dat_post

dat_pre  %>% 
  transmute(F_normal_m = LifeExpNocontMean,
            F_normal_sd = sqrt(Nnocon)*LifeExpNocontSE,
            F_normal_n = Nnocon,
            F_hormonal_m = LifeExpHormMean,
            F_hormonal_sd = sqrt(Nhorm)*LifeExpHormSE,
            F_hormonal_n = Nhorm,
            species = species,
            species_tree = vertlife.species,
            phylogeny = gsub(" ", "_", vertlife.species)
  ) -> dat_pre

dat_post %>% 
  transmute(F_normal_m = LifeExpNocontMean,
          F_normal_sd = sqrt(Nnocon)*LifeExpNocontSE,
          F_normal_n = Nnocon,
          F_hormonal_m = LifeExpHormMean,
          F_hormonal_sd = sqrt(Nhorm)*LifeExpHormSE,
          F_hormonal_n = Nhorm,
          species = species,
          species_tree = vertlife.species,
          phylogeny = gsub(" ", "_", vertlife.species)
) -> dat_post

# getting effect size

dat_pre <- escalc("ROM", 
                    m2i = F_normal_m,
                    m1i = F_hormonal_m,
                    sd2i = F_normal_sd,
                    sd1i = F_hormonal_sd,
                    n2i = F_normal_n,
                    n1i = F_hormonal_n,
                    data = dat_pre,
)

dat_post <- escalc("ROM", 
                  m2i = F_normal_m,
                  m1i = F_hormonal_m,
                  sd2i = F_normal_sd,
                  sd1i = F_hormonal_sd,
                  n2i = F_normal_n,
                  n1i = F_hormonal_n,
                  data = dat_post,
)

# getting a tree

tree <- read.tree( here("data", "tree_zoo4.tre"))

# life span data 
#to_drop <-
#  tree$tip.label[which(!(tree$tip.label %in% unique(dat$phylogeny)))]

#tree <- drop.tip(tree, to_drop)
length(tree$tip.label)

tree <- as.ultrametric(tree)

#tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)

# modeling21

mod_pre <- rma.mv(yi, V = vi, 
                    random = list(
                      ~1|species,
                      ~1|phylogeny), 
                    R = list(phylogeny = cor_tree), 
                    data = dat_pre)
summary(mod_pre)
i2_ml(mod_pre)

p1 <- orchard_plot(mod_pre, xlab = "lnRR (Pre-2005)", group = "species", g = FALSE) + ylim(-0.5, 1.4)


mod_post <- rma.mv(yi, V = vi, 
                  random = list(
                    ~1|species,
                    ~1|phylogeny), 
                  R = list(phylogeny = cor_tree), 
                  data = dat_post)
summary(mod_post)
i2_ml(mod_post)

p2 <- orchard_plot(mod_post, xlab = "lnRR (Post-2005)", group = "species", g = FALSE)  + ylim(-0.5, 1.4)

p1/p2 +  plot_annotation(tag_levels = 'A')


