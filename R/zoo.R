# quick analysis of zoo data

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
               here,
               cowplot
)

########################
# data: life expectancy
########################

dat <- read_csv(here("data", "zoo.csv"), na = c("", "NA"))

#glimpse(dat)
#names(dat)

dim(dat)

# getting a tree

#dat$Phylogeny <- gsub("Perca_fluviatilis", "Lamperta_fluviatilis", dat$Phylogeny) #replace with the original name
#tree <- read.tree(here("data/tree_zoo.tre"))
#tree <- read.tree( here("data", "tree_zoo_all.tre"))
tree <- read.tree( here("data", "tree_zoo4.tre"))

# life span data 
#to_drop <-
#  tree$tip.label[which(!(tree$tip.label %in% unique(lifespan$phylogeny)))]

#tree <- drop.tip(tree, to_drop)

#fortify(tree)

#tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)

# extra
tax <- read.csv(here("data", "vertlife_taxonomy_translation_table.csv"))

setdiff( dat$species, tax$vertlife.species)
setdiff(tree$tip.label, gsub(" ","_", tax$vertlife.species))
match(tree$tip.label, gsub(" ","_", tax$vertlife.species))

# checking naming consistency between 
length(unique(dat$species))
setdiff(gsub(" ","_", dat$species), tree$tip.label)
#  "Aonyx_cinereus"    "Bubalus_bubalis"   "Cervus_canadensis" "Equus_asinus" 

setdiff(tree$tip.label, gsub(" ","_", dat$species))
# "Bubalus_arnee"   "Equus_africanus" "Aonyx_cinerea"  

# what to change for spp

# Aonyx_cinereus -> Aonyx_cinerea
# Bubalus_bubalis -> Bubalus_arnee
# Equus_asinus -> Equus_africanus

#setdiff(tree$tip.label, unique(dat_all$phylogeny))
# [1] "Chrysocyon_brachyurus" - maned wolf
# "Crocuta_crocuta" - spotted hyena
# "Panthera_uncia" - snow lepard
# "Neofelis_nebulosa" - 

# creating spp

dat$spp <- dat$species

dat$spp[dat$spp == "Aonyx cinereus"] <- "Aonyx cinerea"
dat$spp[dat$spp == "Bubalus bubalis"] <- "Bubalus arnee"
dat$spp[dat$spp == "Equus asinus"] <- "Equus_africanus"

dat <- dat %>% filter(species != "Pseudocheirus peregrinus")

dim(dat)

# filtering and getting SD

################
# male surgical
################

dat_m_surg <- dat %>% filter(is.na(Male_Surgical_Mean) == FALSE, is.na(Male_None_Mean) == FALSE) %>% 
  mutate(M_control_m = Male_None_Mean,
         M_control_sd = sqrt(Male.None)*Male_None_SE,
         M_control_n = Male.None,
         M_surgical_m = Male_Surgical_Mean,
         M_surgical_sd = sqrt(Male.Surgical)*Male_Surgical_SE,
         M_surgical_n = Male.Surgical,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "male",
         type = "surgical")

  
# getting effect size

dat_m_surg <- escalc("ROM", 
              m1i = M_surgical_m,
              m2i = M_control_m,
              sd1i = M_surgical_sd,
              sd2i = M_control_sd,
              n1i = M_surgical_n,
              n2i = M_control_n,
              data = dat_m_surg,
              )

# checking the match
matched <- match((dat_m_surg$phylogeny), colnames(cor_tree))

which(is.na(matched))

#grep("Equus",colnames(cor_tree))

dat_m_surg$phylogeny[8] <- "Aonyx_cinerea"
dat_m_surg$phylogeny[18] <- "Bubalus_arnee"
#dat_m_surg$phylogeny[31] <- "Cervus_elaphus"
dat_m_surg$phylogeny[45] <- "Equus_africanus"
# meta-analysis
#dat_f_horm$phylogeny[3] <- "Aonyx_cinerea"
#dat_f_horm$phylogeny[14] <- "Cervus_elaphus"


#[1]  8 18 31 45

# adjusting phylogeny

mod_m_surg <- rma.mv(yi, V = vi, 
                     random = list(
                       ~1|species,
                       ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat_m_surg)
summary(mod_m_surg)
i2_ml(mod_m_surg)

orchard_plot(mod_m_surg, xlab = "lnRR (male surgical)", group = "species", g = FALSE)

# there is an outliner
dat_m_surg[which(dat_m_surg$yi == max(dat_m_surg$yi)), "species"]
#[1] "Pseudocheirus peregrinus"


# exclusing outlier

dat_m_surg2 <- dat_m_surg[-which(dat_m_surg$yi == max(dat_m_surg$yi)), ]

mod_m_surg2 <- rma.mv(yi, V = vi, 
                     random = list(~1|species),
                     data = dat_m_surg2)
summary(mod_m_surg2)
i2_ml(mod_m_surg2)

orchard_plot(mod_m_surg2,xlab = "lnRR (male surgical)", group = "species",  g = FALSE)


#######################
# female hormonal data 
########################
dat_f_horm <- dat %>% filter(is.na(Female_Hormonal_Mean) == FALSE, is.na(Female_None_Mean) == FALSE) %>% 
  mutate(F_control_m = Female_None_Mean,
         F_control_sd = sqrt(Female.None)*Female_None_SE,
         F_control_n = Female.None,
         F_hormonal_m = Female_Hormonal_Mean,
         F_hormonal_sd = sqrt(Female.Hormonal)*Female_Hormonal_SE,
         F_hormonal_n = Female.Hormonal,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "female",
         type = "hormonal")

# getting effect size

dat_f_horm <- escalc("ROM", 
                     m1i = F_hormonal_m,
                     m2i = F_control_m,
                     sd1i = F_hormonal_sd,
                     sd2i = F_control_sd,
                     n1i = F_hormonal_n,
                     n2i = F_control_n,
                     data = dat_f_horm,
)

# matching tree names
matched <- match((dat_f_horm$phylogeny), colnames(cor_tree))

which(is.na(matched))

# # meta-analysis
dat_f_horm$phylogeny[3] <- "Aonyx_cinerea"
#dat_f_horm$phylogeny[14] <- "Cervus_elaphus"

mod_f_horm <- rma.mv(yi, V = vi, 
                     random = list(
                       ~1|species,
                       ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat_f_horm)
summary(mod_f_horm)
i2_ml(mod_f_horm)

orchard_plot(mod_f_horm,xlab = "lnRR (female hormonal)", group = "species", g = FALSE)


######################
# female surgical data
######################
dat_f_surg<- dat %>% filter(is.na(Female_Surgical_Mean) == FALSE, is.na(Female_None_Mean) == FALSE) %>% 
  mutate(F_control_m = Female_None_Mean,
         F_control_sd = sqrt(Female.None)*Female_None_SE,
         F_control_n = Female.None,
         F_surgical_m = Female_Surgical_Mean,
         F_surgical_sd = sqrt(Female.Surgical)*Female_Surgical_SE,
         F_surgical_n = Female.Surgical,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "female",
         type = "surgical")


# getting effect size

dat_f_surg <- escalc("ROM", 
                     m1i = F_surgical_m,
                     m2i = F_control_m,
                     sd1i = F_surgical_sd,
                     sd2i = F_control_sd,
                     n1i = F_surgical_n,
                     n2i = F_control_n,
                     data = dat_f_surg,
)

# matching tree data_f_surg
matched <- match((dat_f_horm$phylogeny), colnames(cor_tree))

which(is.na(matched))


mod_f_surg <- rma.mv(yi, V = vi, 
                     random = list(
                       ~1|species,
                       ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat_f_surg)
summary(mod_f_surg)
i2_ml(mod_f_surg)

robust(mod_f_surg, cluster = species)

orchard_plot(mod_f_surg, xlab = "lnRR (female surgical)", group = "species", g = FALSE)

########################
# putting data together
########################
# getting other dat

#  male hormonal
dat_m_horm <- dat %>% filter(is.na(Male_Hormonal_Mean) == FALSE, is.na(Male_None_Mean) == FALSE) %>% 
  mutate(M_control_m = Male_None_Mean,
         M_control_sd = sqrt(Male.None)*Male_None_SE,
         M_control_n = Male.None,
         M_hormonal_m = Male_Hormonal_Mean,
         M_hormonal_sd = sqrt(Male.Hormonal)*Male_Hormonal_SE,
         M_hormonal_n = Male.Hormonal,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "male",
         type = "hormonal")

dat_m_horm <- escalc("ROM", 
                     m1i = M_hormonal_m,
                     m2i = M_control_m,
                     sd1i = M_hormonal_sd,
                     sd2i = M_control_sd,
                     n1i = M_hormonal_n,
                     n2i = M_control_n,
                     data = dat_m_horm,
)

dim(dat_m_horm)

# checking the match
matched <- match((dat_m_horm$phylogeny), colnames(cor_tree))

which(is.na(matched))

dat_m_horm$phylogeny[1] <- "Aonyx_cinerea"

#grep("Equus",colnames(cor_tree))

# male immunological
dat_m_immu <- dat %>% filter(is.na(Male_Immunological_Mean) == FALSE, 
    is.na(Male_None_Mean) == FALSE) %>%
  mutate(M_control_m = Male_None_Mean,
         M_control_sd = sqrt(Male.None)*Male_None_SE,
         M_control_n = Male.None,
         M_immunol_m = Male_Immunological_Mean,
         M_immunol_sd = sqrt(Male.Immunological)*Male_Immunological_SE,
         M_immunol_n = Male.Immunological,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "male",
         type = "immunological")

dat_m_immu <- escalc("ROM", 
                     m1i = M_immunol_m,
                     m2i = M_control_m,
                     sd1i = M_immunol_sd,
                     sd2i = M_control_sd,
                     n1i = M_immunol_n,
                     n2i = M_control_n,
                     data = dat_m_immu,
)

dim(dat_m_immu)

# checking the match
matched <- match((dat_m_immu$phylogeny), colnames(cor_tree))

which(is.na(matched))

# female immunological
dat_f_immu <- dat %>% filter(is.na(Female_Immunological_Mean) == FALSE, 
    is.na(Female_None_Mean) == FALSE) %>%
  mutate(F_control_m = Female_None_Mean,
         F_control_sd = sqrt(Female.None)*Female_None_SE,
         F_control_n = Female.None,
         F_immunol_m = Female_Immunological_Mean,
         F_immunol_sd = sqrt(Female.Immunological)*Female_Immunological_SE,
         F_immunol_n = Female.Immunological,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "female",
         type = "immunological")

dat_f_immu <- escalc("ROM", 
                     m1i = F_immunol_m,
                     m2i = F_control_m,
                     sd1i = F_immunol_sd,
                     sd2i = F_control_sd,
                     n1i = F_immunol_n,
                     n2i = F_control_n,
                     data = dat_f_immu,
)

dim(dat_f_immu)

# checking the match
matched <- match((dat_f_immu$phylogeny), colnames(cor_tree))

which(is.na(matched))

rbind(
dat_m_horm[ , c(1, 59:64)], # 1
dat_m_surg[ , c(1, 59:64)], # 2
dat_f_horm[ , c(1, 59:64)], # 3 
dat_f_surg[ , c(1, 59:64)], # 4
dat_m_immu[ , c(1, 59:64)], # 5
dat_f_immu[ , c(1, 59:64)] # 6
) -> dat_all

dim(dat_all)

# checking the match
matched <- match((dat_all$phylogeny), colnames(cor_tree))

which(is.na(matched))

#dat_all$phylogeny[1] <- "Aonyx_cinerea"

dat_all$obs_id <- factor(1:nrow(dat_all))

dat_all %>% mutate(sex_type = paste(sex, type, sep = "_")) -> dat_all

dat_all <- dat_all %>% filter(!species == "Pseudocheirus peregrinus")

setdiff(tree$tip.label, unique(dat_all$phylogeny))

#to_drop <-
#  tree$tip.label[which(!(tree$tip.label %in% unique(dat_all$phylogeny)))]

#tree3 <- drop.tip(tree, to_drop)

#write.tree(tree3, here("data", "tree_zoo3.tre"))
#setdiff(tree$tip.label, unique(dat_all$phylogeny))
# [1] "Chrysocyon_brachyurus" - maned wolf
# "Crocuta_crocuta" - spotted hyena
# "Panthera_uncia" - snow lepard
# "Neofelis_nebulosa" - clouded leopard
# "Pseudocheirus_peregrinus" - posumm

# saving data_all

#dat_all %>% mutate(effect = "lifespan") -> dat_all

#saveRDS(dat_all, here("Rdata", "lifespan_all.RDS"))

dat_all <- readRDS(here("Rdata", "lifespan_all.RDS"))

setdiff(unique(dat_all$phylogeny), tree$tip.label)

VCV <- vcalc(vi, species, rho = 0.5, data = dat_all)


mod_all <- rma.mv(yi, V = VCV,
                     random = list(
                       ~1|species,
                       ~1|phylogeny,
                       ~1|obs_id),
                     R = list(phylogeny = cor_tree),
                     data = dat_all)
summary(mod_all)
i2_ml(mod_all)

robust(mod_all, cluster = species)  

orchard_plot(mod_all,xlab = "lnRR (all)", group = "species", g = FALSE)

# sex

mod_all1 <- rma.mv(yi, V = VCV,
                   mod = ~ sex,
                  random = list(
                    ~1|species,
                    ~1|phylogeny,
                    ~1|obs_id),
                  R = list(phylogeny = cor_tree),
                  data = dat_all)
summary(mod_all1)


mod_all1b <- rma.mv(yi, V = VCV,
                   mod = ~ sex - 1,
                   random = list(
                     ~1|species,
                     ~1|phylogeny,
                     ~1|obs_id),
                   R = list(phylogeny = cor_tree),
                   data = dat_all)
summary(mod_all1b)

r2_ml(mod_all1)

orchard_plot(mod_all1, mod = "sex",
             xlab = "lnRR (all)", group = "species", g = FALSE)

# TODO - contrast...

# types

dat_all$type <- factor(dat_all$type, 
                   levels = rev(c("surgical", "hormonal", "immunological")))



mod_all2 <- rma.mv(yi, V = VCV,
                   mod = ~ type,
                   random = list(
                     ~1|species,
                     ~1|phylogeny,
                     ~1|obs_id),
                   R = list(phylogeny = cor_tree),
                   data = dat_all)
summary(mod_all2)


mod_all2b <- rma.mv(yi, V = VCV,
                 mod = ~ type-1,
                 random = list(
                   ~1|species,
                   ~1|phylogeny,
                   ~1|obs_id),
                 R = list(phylogeny = cor_tree),
                 data = dat_all)

summary(mod_all2b)

r2_ml(mod_all2)


orchard_plot(mod_all2, mod = "type",
             xlab = "lnRR (all)", group = "species", g = FALSE, angle = 90)

# interaction

mod_all3 <- rma.mv(yi, V = VCV,
                   mod = ~ sex_type-1,
                   random = list(
                     ~1|species,
                     ~1|phylogeny,
                     ~1|obs_id),
                   R = list(phylogeny = cor_tree),
                   data = dat_all)
summary(mod_all3)

r2_ml(mod_all3)

orchard_plot(mod_all3, mod = "sex_type",
             xlab = "lnRR (all)", group = "species", g = FALSE, angle =45)


# Creating Fig 1

main <- mod_results(mod_all, group = "species")
sex_diff <- mod_results(mod_all1, mod = "sex", group = "species")

combined <- submerge(sex_diff, main)

combo <- orchard_plot(combined,
                      group = "species", 
                      xlab = "log response ratio (lnRR)") +
  scale_colour_manual(values = rev(c("#999999", "#88CCEE", "#CC6677"))) +
  scale_fill_manual(values = rev(c("#999999", "#88CCEE", "#CC6677")))

#combo


type_diff <- orchard_plot(mod_all2, mod = "type",
             xlab = "log response ratio (lnRR)", group = "species", g = FALSE, angle = 90)+
  scale_colour_manual(values = rev(c("#117733",  "#332288", "#DDCC77"))) +
  scale_fill_manual(values = rev(c("#117733",  "#332288", "#DDCC77")))

#type_diff
# "#117733",  "#332288", "#DDCC77" "#AA4499"

###################
## Figure 1
###################

library(cowplot)
plot_grid(p_phylo) / (combo + type_diff) + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(heights = c(1.5,1.0))

#############
#############
# NOT RUN
#############
#############

###################################

##############################
# data: gamma - mortality risk
##############################

# TODO - adding phylogeny

dat2 <- read.csv(here("data", "gammaBaSTA.csv"))

#dat2 <- tibble(dat2)
glimpse(dat2)
names(dat2)

# creating spp

dat2$spp <- dat$species

dat2$spp[dat$spp == "Aonyx cinereus"] <- "Aonyx cinerea"
dat2$spp[dat$spp == "Bubalus bubalis"] <- "Bubalus arnee"
dat2$spp[dat$spp == "Equus asinus"] <- "Equus_africanus"

dim(dat2)



# filtering and getting SD

################
# male surgical
################

dat2_m_surg <- dat2 %>% filter(is.na(Male_Surgical_Mean) == FALSE) %>% 
  mutate(yi = Male_Surgical_Mean ,
         vi = Male_Surgical_SE^2,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "male",
         type = "surgical")


matched <- match((dat2_m_surg$phylogeny), colnames(cor_tree))

which(is.na(matched))

dat2_m_surg$phylogeny[8] <- "Aonyx_cinerea"
dat2_m_surg$phylogeny[18] <- "Bubalus_arnee"
#dat2_m_surg$phylogeny[31] <- "Cervus_elaphus"
dat2_m_surg$phylogeny[45] <- "Equus_africanus"
# meta-analysis
#dat2_f_horm$phylogeny[3] <- "Aonyx_cinerea"
#dat2_f_horm$phylogeny[14] <- "Cervus_elaphus"


hist(dat2_m_surg$yi)

mod2_m_surg <- rma.mv(yi, V = vi, 
                     random = list(~1|species,
                     ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat2_m_surg)
summary(mod2_m_surg)
i2_ml(mod2_m_surg)

orchard_plot(mod2_m_surg,xlab = "lnHR (male surgical)", group = "species",g = FALSE)

# there is an outliner
dat2_m_surg[which(dat2_m_surg$yi == min(dat2_m_surg$yi)), "species"]
#[1] "Pseudocheirus peregrinus"

#######################
# female hormonal data 
########################

dat2_f_horm <- dat2 %>% filter(is.na(Female_Hormonal_Mean) == FALSE) %>% 
  mutate(yi = Female_Hormonal_Mean ,
         vi = Female_Hormonal_SE^2,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "female",
         type = "hormonal")


matched <- match((dat2_f_horm$phylogeny), colnames(cor_tree))

which(is.na(matched))

# meta-analysis
dat2_f_horm$phylogeny[3] <- "Aonyx_cinerea"
#dat2_f_horm$phylogeny[14] <- "Cervus_elaphus"

hist(dat2_f_horm$yi)


mod2_f_horm <- rma.mv(yi, V = vi, 
                     random = list(~1|species,
                                   ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat2_f_horm)
summary(mod2_f_horm)
i2_ml(mod2_f_horm)

orchard_plot(mod2_f_horm,xlab = "lnHR (female hormonal)", group = "species",  g = FALSE)

######################
# female surgical data
######################
dat2_f_surg<- dat2 %>% filter(is.na(Female_Surgical_Mean) == FALSE) %>% 
  mutate(yi = Female_Surgical_Mean ,
         vi = Female_Surgical_SE^2,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "female",
         type = "surgical")

hist(dat2_f_surg$yi)


mod2_f_surg <- rma.mv(yi, V = vi, 
                     random = list(~1|species,
                                   ~1|phylogeny), 
                     R = list(phylogeny = cor_tree), 
                     data = dat2_f_surg)
summary(mod2_f_surg)
i2_ml(mod2_f_surg)

orchard_plot(mod2_f_surg, xlab = "lnHR (female surgical)", group = "species", g = FALSE)

##########
# all data 
###########

# M horm

dat2_m_horm <- dat2 %>% filter(is.na(Male_Hormonal_Mean) == FALSE) %>% 
  mutate(yi = Male_Hormonal_Mean ,
         vi = Male_Hormonal_SE^2,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "male",
         type = "surgical")


matched <- match((dat2_m_horm$phylogeny), colnames(cor_tree))

which(is.na(matched))

dim(dat2_m_horm)

dat2_m_horm$phylogeny[1] <- "Aonyx_cinerea"

# M immu
dat2_m_immu <- dat2 %>% filter(is.na(Male_Immunological_Mean) == FALSE) %>% 
  mutate(yi = Male_Immunological_Mean ,
         vi = Male_Immunological_SE^2,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "male",
         type = "immunological")

dim(dat2_m_immu)

# F immu
dat2_f_immu <- dat2 %>% filter(is.na(Female_Immunological_Mean) == FALSE) %>% 
  mutate(yi = Female_Immunological_Mean ,
         vi = Female_Immunological_SE^2,
         species = species,
         species_tree = spp,
         phylogeny = gsub(" ","_", species),
         sex = "female",
         type = "immunological")

dim(dat2_f_immu)


#####
######

rbind(
  dat2_m_horm[ , c(1, 35:40)], # 1
  dat2_m_surg[ , c(1, 35:40)], # 2
  dat2_f_horm[ , c(1, 35:40)], # 3 
  dat2_f_surg[ , c(1, 35:40)], # 4
  dat2_m_immu[ , c(1, 35:40)], # 5
  dat2_f_immu[ , c(1, 35:40)]  # 6
) -> dat2_all

dim(dat2_all)

# analysis

dat2_all$obs_id <- factor(1:nrow(dat2_all))

dat2_all %>% mutate(sex_type = paste(sex, type, sep = "_")) -> dat2_all


dat2_all <- dat2_all %>% filter(!species == "Pseudocheirus peregrinus")

dat2_all %>% select("species", "species_tree", "phylogeny",
                    "sex", "type", "yi", "vi", "obs_id", "sex_type") -> dat2_all

dat2_all %>% mutate(effect = "risk") -> dat2_all

# saving data_all

saveRDS(dat2_all, here("Rdata", "risk_all.RDS"))

mod2_all <- rma.mv(yi, V = vi,
                  random = list(
                    ~1|species,
                    ~1|phylogeny,
                    ~1|obs_id),
                  R = list(phylogeny = cor_tree),
                  data = dat2_all)
summary(mod2_all)
i2_ml(mod2_all)

orchard_plot(mod2_all,xlab = "lnHR (all)", group = "species", g = FALSE)

# sex

mod2_all1 <- rma.mv(yi, V = vi,
                   mod = ~ sex,
                   random = list(
                     ~1|species,
                     ~1|phylogeny,
                     ~1|obs_id),
                   R = list(phylogeny = cor_tree),
                   data = dat2_all)
summary(mod2_all1)

orchard_plot(mod2_all1, mod = "sex",
             xlab = "lnRR (all)", group = "species", g = FALSE)

# types

mod2_all2 <- rma.mv(yi, V = vi,
                   mod = ~ type-1,
                   random = list(
                     ~1|species,
                     ~1|phylogeny,
                     ~1|obs_id),
                   R = list(phylogeny = cor_tree),
                   data = dat2_all)
summary(mod2_all2)

orchard_plot(mod2_all2, mod = "type",
             xlab = "lnRR (all)", group = "species", g = FALSE, angle = 90)

# interaction

mod2_all3 <- rma.mv(yi, V = vi,
                   mod = ~ sex_type-1,
                   random = list(
                     ~1|species,
                     ~1|phylogeny,
                     ~1|obs_id),
                   R = list(phylogeny = cor_tree),
                   data = dat2_all)
summary(mod2_all3)

orchard_plot(mod2_all3, mod = "sex_type",
             xlab = "lnRR (all)", group = "species", g = FALSE, angle =45)


# looking correlation between lifespan and risk of death 

tib1 <- tibble(lifespan = dat_all$yi, risk = dat2_all$yi)

ggplot(data = tib1, aes(y = lifespan, x = risk)) + 
  geom_point() + 
  geom_smooth(se = FALSE, method = lm)
