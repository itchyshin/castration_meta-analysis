# figs for zoo data


# pacakges

library(tidyverse)

# Phylogenetic tree:
library(tidytree)
library(ape)
library(ggtree)
library(ggimage)
library(ggstance)
library(here)

# Plot layout:
library(patchwork)


dat <- read_csv(here("data", "zoo.csv"), na = c("", "NA"))

tax <- read.csv(here("data", "vertlife_taxonomy_translation_table.csv"))

glimpse(dat)
names(dat)

#dat$Phylogeny <- gsub("Perca_fluviatilis", "Lamperta_fluviatilis", dat$Phylogeny) #replace with the original name
tree <- read.tree(here("data/tree_zoo.tre"))

tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)

# getting relevant data

lifespan <- readRDS(here("Rdata", "lifespan_all.RDS"))


p <-
  ggtree(
    tree,
    layout = "rectangular",
    ladderize = TRUE,
    size = 0.2,
    color = "#454747"
  )
for (i in 1:nrow(dford)) {
  p <-
    p + geom_cladelabel(
      node = dford[i, "node"],
      label = NA, 
      #label = dford[i, "order"], # uncomment to plot order names
      hjust = 0.4,
      align = T,
      offset.text = 10,
      barsize = 1,
      fontsize = 2.5,
      fill = "grey20",
      color = "grey50"
    )
}
p
# p1 <- ggtree(tree, branch.length = "branch.length") 

# p1





# note needed 

load(
  "./Rdata/maxCredTreeMammals.RData"
)
tree2 <- maxCred


# example



#####

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Project: Life expectancy differences in contracepted vs non-contracepted
# mammals
# Version: 01
# Author: Johanna Staerk
# Description: Phylogenetic tree Figure
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

rm(list = ls())

#install.packages(c("tidyverse" , "ape", "tidytree", "ggtree", "ggimage", "patchwork))

library(tidyverse)

# Phylogenetic tree:
library(tidytree)
library(ape)
library(ggtree)
library(ggimage)

# Plot layout:
library(patchwork)





# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# READ DATA ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# life expectancy values with and without contraception for ZIMS species
dif <- read.csv("DATA/life_expectancies/lifeExpBaSTA.csv")

# read phylogenetic tree mammals (Vertlife.org)
load(
  "./Rdata/maxCredTreeMammals.RData"
)
tree2 <- maxCred

# Taxonomic names from vertlife.org
tax <-
  read.csv(
    "DATA/vertlife/vertlife_taxonomy_translation_table.csv"
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PREP DATA ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# subset relevant columns in ZIMS data
dif <- dif %>% select(
  species,
  Female_None_Mean,
  Female_Surgical_Mean,
  Male_None_Mean,
  Male_Surgical_Mean,
  Female_Hormonal_Mean,
  Male_Hormonal_Mean
)

# calculate life expectancy differences for each sex:

# surgical contraception
dif$dif_surg_female <-
  (dif$Female_None_Mean - dif$Female_Surgical_Mean) / dif$Female_None_Mean
dif$dif_surg_male <-
  (dif$Male_None_Mean - dif$Male_Surgical_Mean) / dif$Male_None_Mean

# hormonal contraception
dif$dif_horm_female <-
  (dif$Female_None_Mean - dif$Female_Hormonal_Mean) / dif$Female_None_Mean
dif$dif_horm_male <-
  (dif$Male_None_Mean - dif$Male_Hormonal_Mean) / dif$Male_None_Mean

# remove species that have no surgical or hormonal contraceptiond data
dif <-
  dif %>% filter(!(is.na(dif_surg_male)) |
                   !(is.na(dif_surg_female)) |
                   !(is.na(dif_horm_female)) |
                   !(is.na(dif_horm_male)))

# match species names to tree labels
tree.label <- str_replace(tree$tip.label, "_", " ")
dif <- dif %>% left_join(tax, by = c("species" = "zims.species"))

# remove Cervus canadiensis as it is synynym to Cervus elaphus in tree taxonomy
dif <- dif[which(!(dif$species == "Cervus canadensis")), ]

# Remove extreme outlier in life exp differences (check?)
dif <- dif %>% filter(!species == "Pseudocheirus peregrinus")

# add species names as row names so tree can match
rownames(dif) <- gsub(" ", "_", dif$vertlife.species)

#  drop any species from tree that are not in final ZIMS dataset
keep <- dif$vertlife.species
keep <-  gsub(" ", "_", keep)
to_drop <-
  tree$tip.label[which(!(tree$tip.label %in% keep))]
tree_subset <- drop.tip(tree, to_drop)

# use getMRCA() to obtain common ancestor nodes to position the order silhouttes
tree.tibble <- tidytree::as_tibble(tree_subset)
ord <- unique(tax$order)

dford <- data.frame(order = ord, node = NA)
for (i in ord) {
  tip <-  dif[which(dif$order == i), "vertlife.species"]
  tip <-  gsub(" ", "_", tip)
  if (length(tip) > 1) {
    dford[which(dford$order == i), "node"] <-
      getMRCA(tree_subset, tip = tip) 
  } else {
    dford[which(dford$order == i), "node"] <-
      tree.tibble[which(tree.tibble$label == tip), "node"]
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Basic rectangular tree
p <-
  ggtree(
    tree_subset,
    layout = "rectangular",
    ladderize = TRUE,
    size = 0.2,
    color = "#454747"
  )
for (i in 1:nrow(dford)) {
  p <-
    p + geom_cladelabel(
      node = dford[i, "node"],
      label = NA, 
      #label = dford[i, "order"], # uncomment to plot order names
      hjust = 0.4,
      align = T,
      offset.text = 10,
      barsize = 1,
      fontsize = 2.5,
      fill = "grey20",
      color = "grey50"
    )
}
p

# Add order silhouettes
# (takes a while to display)

# Image directory
imgdir <- "DATA/phylopics/"

p1 = p +
  geom_image(
    x = 170,
    y = 148,
    image = paste0(
      imgdir,
      "Artiodactyla_PhyloPic.8b567be8.An-Ignorant-Atheist.Antilopinae.png"
    ),
    size = 0.04
  ) +
  geom_image(
    x = 170,
    y = 116,
    image = paste0(
      imgdir,
      "Perissodactlyla_PhyloPic.071ee517.David-Orr.Equus-ferus-caballus.png"
    ),
    size = 0.06
  ) +
  geom_image(
    x = 170,
    y = 100,
    image = paste0(
      imgdir,
      "Carnivora_PhyloPic.34e482b4.An-Ignorant-Atheist.Panthera.png"
    ),
    size = 0.06
  ) +
  geom_image(
    x = 170,
    y = 80,
    image = paste0(
      imgdir,
      "Chrioptera_PhyloPic.e7da460a.Margot-Michaud.Chiroptera_Eptesicus_Eptesicus-fuscus_Vespertilio-Noctilio_Vespertilionidae_Vespertilioniformes_Vespertilioninae_Vespertilionoidea.png"
    ),
    size = 0.06
  ) +
  geom_image(
    x = 170,
    y = 48,
    image = paste0(
      imgdir,
      "Primates_PhyloPic.071db0d0.Margot-Michaud.Papio_Papio-anubis.png"
    ),
    size = 0.06
  ) +
  geom_image(
    x = 170,
    y = 19,
    image = paste0(
      imgdir,
      "Rodentia_PhyloPic.570c7d9e.Alexandra-van-der-Geer.Rattus_Rattus-exulans.png"
    ),
    size = 0.07
  ) +
  geom_image(
    x = 170,
    y = 7,
    image = paste0(
      imgdir,
      "PhyloPic.b62bab6e.An-Ignorant-Atheist.Macropus-Macropus.png"
    ),
    size = 0.06
  )

p1

# Add barplot with life expectancy differences

# reorder dataframe to match the tree
d = fortify(tree_subset)
dd = subset(d, isTip)
ordered <- dd$label[order(dd$y, decreasing=TRUE)]
dif <- dif[match(ordered, row.names(dif)),]

# reformat data to long format
d1 <- dif %>% select(vertlife.species, starts_with("dif_"))
d1$order <- seq(1:nrow(d1))
d1 <- d1 %>% pivot_longer(cols = c(dif_surg_female, dif_surg_male,
                                   dif_horm_female, dif_horm_male))
d1 <- d1[which(!(is.na(d1$value))), ]
d1$sex <- ifelse(d1$name == "dif_surg_female" | d1$name == "dif_horm_female", "Female", "Male")
d1$method <- ifelse(d1$name == "dif_surg_female" | d1$name == "dif_surg_male", "Surgical", "Hormonal")
d1$bias <- ifelse(d1$value <= 0, "Contraception", "No contraception")

p2 <- ggplot(d1, aes(x=0, xend=value, y=reorder(vertlife.species, -order), yend= reorder(vertlife.species, -order), color = bias,
                     text = paste("species:", reorder(vertlife.species, -order)))) +
  geom_segment(linewidth = 1) +
  facet_wrap(method~sex, ncol = 4) +
  scale_color_manual(values = c("#336699", "#FF9933", NA), name = "") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.ticks  = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        panel.grid.minor.x = element_blank()) +
  xlab(bquote(~italic((l[n] - l[c]) / l[n]))) +
  ylab("")

# Plot both trees (takes a while)
p1 + p2
#ggsave("RESULTS/PLOTS/PlotRect.pdf", width = 25, height = 18, units = "cm")
dev.off()


#########
##########

# get class data
spp.class <- dplyr::distinct(dat, Phy, .keep_all = TRUE) %>% dplyr::select(Class, Phy)
names(spp.class)[2] <- "Species"

spp.class2 <- left_join(tip.label,spp.class,by = "Species")

# Actinopterygii and Osteichthyes all belong to fish, although there are differnt. In our dataset, only Actinopterygii only has one estimate. To make the figure easy, we change Actinopterygii to Osteichthyes 
spp.class2$Class[spp.class2$Class == "Actinopterygii"] <- "Osteichthyes"

# version 2
p1 <- ggtree(my.tr, branch.length = "branch.length") 
# https://yulab-smu.top/treedata-book/chapter7.html  
p2 <- p1 %<+% spp.class2 + geom_tiplab(aes(color  = Class),size=3,fontface = "italic") + geom_tippoint(aes(color=Class)) + xlim_expand(c(0,3), panel = "Tree") + guides(color="none")  # xlim issue https://stackoverflow.com/questions/48928250/ggtreefacet-plot-second-panel-uses-xlim-parameter-from-first-one
p3 <- facet_plot(p2, panel = 'N observations', data = re.spp.lnRR2, 
                 geom = geom_barh, 
                 mapping = aes(x = N_obs,fill=Class,color=Class), alpha = 0.4, stat='identity') +  guides(fill="none",color="none") + # #0072B2 https://yulab-smu.top/treedata-book/chapter12.html,
  geom_facet(panel = "Estimates of mean ratio (lnRR)", data = re.spp.lnRR2,
             geom = ggstance::geom_pointrangeh, # geom_pointrange
             mapping = aes(x = Mean, xmin = Lower_bound, xmax = Upper_bound, color = Class)) + # #CC79A7
  # add overall effect MLMA_lnRR$b
  geom_facet(panel = "Estimates of mean ratio (lnRR)", data = re.spp.lnRR2,
             geom = geom_vline, 
             mapping = aes(xintercept = -0.45),linetype="dashed",color="grey") +  # add multiple panels or multiple layers on the same panels (Figure 13.1) # https://yulab-smu.top/treedata-book/chapter7.html - at Chapter 13
  geom_facet(panel = "Estimates of CV ratio (lnCVR)", data = re.spp.lnCVR2,
             geom = ggstance::geom_pointrangeh, 
             mapping = aes(x = Mean, xmin = Lower_bound, xmax = Upper_bound, color = Class)) + 
  # add overall effect MLMA_lnCVR$b
  geom_facet(panel = "Estimates of CV ratio (lnCVR)", data = re.spp.lnRR2,
             geom = geom_vline, 
             mapping = aes(xintercept = 0.21), linetype="dashed",color="grey") + 
  theme_tree2() + theme(strip.background = element_rect(fill = "white"))

# lbs <- c(Tree = "Phylogenetic tree") # change facet name # https://guangchuangyu.github.io/cn/2018/09/facet_labeller/
# p3 <- ggtree::facet_labeller(p3,lbs)
p4 <- facet_widths(p3, c(Tree=0.7,`N observations`=0.2,`Estimates of mean ratio (lnRR)`=0.4,`Estimates of CV ratio (lnCVR)`=0.4)) # set the widths of specific panels (it also supports using a name vector)

png(filename = "Figure S7 (Phy est).png", width = 8, height = 6, units = "in", type = "windows", res = 400)
p4
dev.off()


# Figure S7
# getting blups (best linear predictors from the model) for lnRR
blups <- ranef(MLMA_lnRR) 
t.spp <- blups$Spp # this needs to be added
t.phy <- blups$Phy
t.study <- blups$StudyID # study
t.spp <- rownames_to_column(t.spp, var = "Species")
t.phy <- rownames_to_column(t.phy, var = "Species")
t.spp$Species <- t.phy$Species # uppercase
t.study <- rownames_to_column(t.study, var = "Study") # study

# average species
colnames(t.spp) <- c("Species", "Deviation", "SE", "Lower_bound", "Upper_bound")
colnames(t.phy) <- c("Species", "Deviation", "SE", "Lower_bound", "Upper_bound")
colnames(t.study) <- c("Study", "Deviation", "SE", "Lower_bound", "Upper_bound")

# sorting match between study and species
dat %>% group_by(StudyID) %>% summarise(Species = unique(Phy)) -> dat_match # warning information: https://stackoverflow.com/questions/62140483/how-to-interpret-dplyr-message-summarise-regrouping-output-by-x-override

index <- match(dat_match$Species, t.spp$Species)

spp.mean <- rep(MLMA_lnRR$b, dim(t.spp)[1]) + t.spp$Deviation + t.phy$Deviation
spp.se <- sqrt(MLMA_lnRR$se^2 + t.spp$SE^2 +  t.phy$SE^2) 
spp.lb <- spp.mean - spp.se * qnorm(0.975) 
spp.ub <- spp.mean + spp.se * qnorm(0.975)
re.spp <- tibble(Species = t.spp$Species, Mean = spp.mean, SE = spp.se, Lower_bound = spp.lb, Upper_bound = spp.ub) %>% arrange(Species)
N_obs <- dat %>% group_by(Phy) %>% summarise(N_obs=n()) # calculate sample size
names(N_obs)[1] <- "Species"
re.spp.lnRR <- left_join(re.spp,N_obs,by = "Species")
tip.label <- data.frame(Species = my.tr$tip.label) # extract tip label
re.spp.lnRR2 <- left_join(tip.label,re.spp.lnRR,by = "Species") 


# getting blups (best linear predictors from the model) for lnCVR
blups <- ranef(MLMA_lnCVR) 
t.spp <- blups$Spp # this needs to be added
t.phy <- blups$Phy
t.study <- blups$StudyID # study
t.spp <- rownames_to_column(t.spp, var = "Species")
t.phy <- rownames_to_column(t.phy, var = "Species")
t.spp$Species <- t.phy$Species # uppercase
t.study <- rownames_to_column(t.study, var = "Study") # study

# average species
colnames(t.spp) <- c("Species", "Deviation", "SE", "Lower_bound", "Upper_bound")
colnames(t.phy) <- c("Species", "Deviation", "SE", "Lower_bound", "Upper_bound")
colnames(t.study) <- c("Study", "Deviation", "SE", "Lower_bound", "Upper_bound")

# sorting match between study and species
dat %>% group_by(StudyID) %>% summarise(Species = unique(Phy)) -> dat_match # warning information: https://stackoverflow.com/questions/62140483/how-to-interpret-dplyr-message-summarise-regrouping-output-by-x-override

index <- match(dat_match$Species, t.spp$Species)

spp.mean <- rep(MLMA_lnCVR$b, dim(t.spp)[1]) + t.spp$Deviation + t.phy$Deviation
spp.se <- sqrt(MLMA_lnCVR$se^2 + t.spp$SE^2 +  t.phy$SE^2) 
spp.lb <- spp.mean - spp.se * qnorm(0.975) 
spp.ub <- spp.mean + spp.se * qnorm(0.975)
re.spp <- tibble(Species = t.spp$Species, Mean = spp.mean, SE = spp.se, Lower_bound = spp.lb, Upper_bound = spp.ub) %>% arrange(Species)
N_obs <- dat %>% group_by(Phy) %>% summarise(N_obs=n()) # calculate sample size
names(N_obs)[1] <- "Species"
re.spp.lnCVR <- left_join(re.spp,N_obs,by = "Species")
tip.label <- data.frame(Species = my.tr$tip.label) # extract tip label
re.spp.lnCVR2 <- left_join(tip.label,re.spp.lnCVR,by = "Species") 



# get class data
spp.class <- dplyr::distinct(dat, Phy, .keep_all = TRUE) %>% dplyr::select(Class, Phy)
names(spp.class)[2] <- "Species"

spp.class2 <- left_join(tip.label,spp.class,by = "Species")

# Actinopterygii and Osteichthyes all belong to fish, although there are differnt. In our dataset, only Actinopterygii only has one estimate. To make the figure easy, we change Actinopterygii to Osteichthyes 
spp.class2$Class[spp.class2$Class == "Actinopterygii"] <- "Osteichthyes"

# version 2
p1 <- ggtree(my.tr, branch.length = "branch.length") 
# https://yulab-smu.top/treedata-book/chapter7.html  
p2 <- p1 %<+% spp.class2 + geom_tiplab(aes(color  = Class),size=3,fontface = "italic") + geom_tippoint(aes(color=Class)) + xlim_expand(c(0,3), panel = "Tree") + guides(color="none")  # xlim issue https://stackoverflow.com/questions/48928250/ggtreefacet-plot-second-panel-uses-xlim-parameter-from-first-one
p3 <- facet_plot(p2, panel = 'N observations', data = re.spp.lnRR2, 
                 geom = geom_barh, 
                 mapping = aes(x = N_obs,fill=Class,color=Class), alpha = 0.4, stat='identity') +  guides(fill="none",color="none") + # #0072B2 https://yulab-smu.top/treedata-book/chapter12.html,
  geom_facet(panel = "Estimates of mean ratio (lnRR)", data = re.spp.lnRR2,
             geom = ggstance::geom_pointrangeh, # geom_pointrange
             mapping = aes(x = Mean, xmin = Lower_bound, xmax = Upper_bound, color = Class)) + # #CC79A7
  # add overall effect MLMA_lnRR$b
  geom_facet(panel = "Estimates of mean ratio (lnRR)", data = re.spp.lnRR2,
             geom = geom_vline, 
             mapping = aes(xintercept = -0.45),linetype="dashed",color="grey") +  # add multiple panels or multiple layers on the same panels (Figure 13.1) # https://yulab-smu.top/treedata-book/chapter7.html - at Chapter 13
  geom_facet(panel = "Estimates of CV ratio (lnCVR)", data = re.spp.lnCVR2,
             geom = ggstance::geom_pointrangeh, 
             mapping = aes(x = Mean, xmin = Lower_bound, xmax = Upper_bound, color = Class)) + 
  # add overall effect MLMA_lnCVR$b
  geom_facet(panel = "Estimates of CV ratio (lnCVR)", data = re.spp.lnRR2,
             geom = geom_vline, 
             mapping = aes(xintercept = 0.21), linetype="dashed",color="grey") + 
  theme_tree2() + theme(strip.background = element_rect(fill = "white"))

# lbs <- c(Tree = "Phylogenetic tree") # change facet name # https://guangchuangyu.github.io/cn/2018/09/facet_labeller/
# p3 <- ggtree::facet_labeller(p3,lbs)
p4 <- facet_widths(p3, c(Tree=0.7,`N observations`=0.2,`Estimates of mean ratio (lnRR)`=0.4,`Estimates of CV ratio (lnCVR)`=0.4)) # set the widths of specific panels (it also supports using a name vector)

png(filename = "Figure S7 (Phy est).png", width = 8, height = 6, units = "in", type = "windows", res = 400)
p4
dev.off()







# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Project: Life expectancy differences in contracepted vs non-contracepted
# mammals
# Version: 01
# Author: Johanna Staerk
# Description: Phylogenetic tree Figure
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

rm(list = ls())

#install.packages(c("tidyverse" , "ape", "tidytree", "ggtree", "ggimage", "patchwork))

library(tidyverse)

# Phylogenetic tree:
library(tidytree)
library(ape)
library(ggtree)
library(ggimage)

# Plot layout:
library(patchwork)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# READ DATA ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# life expectancy values with and without contraception for ZIMS species
dif <- read.csv("DATA/life_expectancies/lifeExpBaSTA.csv")

# read phylogenetic tree mammals (Vertlife.org)
load(
  "DATA/vertlife/maxCredTreeMammals.RData"
)
tree <- maxCred

# Taxonomic names from vertlife.org
tax <-
  read.csv(
    "DATA/vertlife/vertlife_taxonomy_translation_table.csv"
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PREP DATA ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# subset relevant columns in ZIMS data
dif <- dif %>% select(
  species,
  Female_None_Mean,
  Female_Surgical_Mean,
  Male_None_Mean,
  Male_Surgical_Mean,
  Female_Hormonal_Mean,
  Male_Hormonal_Mean
)

# calculate life expectancy differences for each sex:

# surgical contraception
dif$dif_surg_female <-
  (dif$Female_None_Mean - dif$Female_Surgical_Mean) / dif$Female_None_Mean
dif$dif_surg_male <-
  (dif$Male_None_Mean - dif$Male_Surgical_Mean) / dif$Male_None_Mean

# hormonal contraception
dif$dif_horm_female <-
  (dif$Female_None_Mean - dif$Female_Hormonal_Mean) / dif$Female_None_Mean
dif$dif_horm_male <-
  (dif$Male_None_Mean - dif$Male_Hormonal_Mean) / dif$Male_None_Mean

# remove species that have no surgical or hormonal contraceptiond data
dif <-
  dif %>% filter(!(is.na(dif_surg_male)) |
                   !(is.na(dif_surg_female)) |
                   !(is.na(dif_horm_female)) |
                   !(is.na(dif_horm_male)))

# match species names to tree labels
tree.label <- str_replace(tree$tip.label, "_", " ")
dif <- dif %>% left_join(tax, by = c("species" = "zims.species"))

# remove Cervus canadiensis as it is synynym to Cervus elaphus in tree taxonomy
dif <- dif[which(!(dif$species == "Cervus canadensis")), ]

# Remove extreme outlier in life exp differences (check?)
dif <- dif %>% filter(!species == "Pseudocheirus peregrinus")

# add species names as row names so tree can match
rownames(dif) <- gsub(" ", "_", dif$vertlife.species)

#  drop any species from tree that are not in final ZIMS dataset
keep <- dif$vertlife.species
keep <-  gsub(" ", "_", keep)
to_drop <-
  tree$tip.label[which(!(tree$tip.label %in% keep))]
tree_subset <- drop.tip(tree, to_drop)

# use getMRCA() to obtain common ancestor nodes to position the order silhouttes
tree.tibble <- tidytree::as_tibble(tree_subset)
ord <- unique(tax$order)

dford <- data.frame(order = ord, node = NA)
for (i in ord) {
  tip <-  dif[which(dif$order == i), "vertlife.species"]
  tip <-  gsub(" ", "_", tip)
  if (length(tip) > 1) {
    dford[which(dford$order == i), "node"] <-
      getMRCA(tree_subset, tip = tip) 
  } else {
    dford[which(dford$order == i), "node"] <-
      tree.tibble[which(tree.tibble$label == tip), "node"]
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOT ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Basic rectangular tree
p <-
  ggtree(
    tree_subset,
    layout = "rectangular",
    ladderize = TRUE,
    size = 0.2,
    color = "#454747"
  )
for (i in 1:nrow(dford)) {
  p <-
    p + geom_cladelabel(
      node = dford[i, "node"],
      label = NA, 
      #label = dford[i, "order"], # uncomment to plot order names
      hjust = 0.4,
      align = T,
      offset.text = 10,
      barsize = 1,
      fontsize = 2.5,
      fill = "grey20",
      color = "grey50"
    )
}
p

# Add order silhouettes
# (takes a while to display)

# Image directory
imgdir <- "DATA/phylopics/"

p1 = p +
  geom_image(
    x = 170,
    y = 148,
    image = paste0(
      imgdir,
      "Artiodactyla_PhyloPic.8b567be8.An-Ignorant-Atheist.Antilopinae.png"
    ),
    size = 0.04
  ) +
  geom_image(
    x = 170,
    y = 116,
    image = paste0(
      imgdir,
      "Perissodactlyla_PhyloPic.071ee517.David-Orr.Equus-ferus-caballus.png"
    ),
    size = 0.06
  ) +
  geom_image(
    x = 170,
    y = 100,
    image = paste0(
      imgdir,
      "Carnivora_PhyloPic.34e482b4.An-Ignorant-Atheist.Panthera.png"
    ),
    size = 0.06
  ) +
  geom_image(
    x = 170,
    y = 80,
    image = paste0(
      imgdir,
      "Chrioptera_PhyloPic.e7da460a.Margot-Michaud.Chiroptera_Eptesicus_Eptesicus-fuscus_Vespertilio-Noctilio_Vespertilionidae_Vespertilioniformes_Vespertilioninae_Vespertilionoidea.png"
    ),
    size = 0.06
  ) +
  geom_image(
    x = 170,
    y = 48,
    image = paste0(
      imgdir,
      "Primates_PhyloPic.071db0d0.Margot-Michaud.Papio_Papio-anubis.png"
    ),
    size = 0.06
  ) +
  geom_image(
    x = 170,
    y = 19,
    image = paste0(
      imgdir,
      "Rodentia_PhyloPic.570c7d9e.Alexandra-van-der-Geer.Rattus_Rattus-exulans.png"
    ),
    size = 0.07
  ) +
  geom_image(
    x = 170,
    y = 7,
    image = paste0(
      imgdir,
      "PhyloPic.b62bab6e.An-Ignorant-Atheist.Macropus-Macropus.png"
    ),
    size = 0.06
  )

p1

# Add barplot with life expectancy differences

# reorder dataframe to match the tree
d = fortify(tree_subset)
dd = subset(d, isTip)
ordered <- dd$label[order(dd$y, decreasing=TRUE)]
dif <- dif[match(ordered, row.names(dif)),]

# reformat data to long format
d1 <- dif %>% select(vertlife.species, starts_with("dif_"))
d1$order <- seq(1:nrow(d1))
d1 <- d1 %>% pivot_longer(cols = c(dif_surg_female, dif_surg_male,
                                   dif_horm_female, dif_horm_male))
d1 <- d1[which(!(is.na(d1$value))), ]
d1$sex <- ifelse(d1$name == "dif_surg_female" | d1$name == "dif_horm_female", "Female", "Male")
d1$method <- ifelse(d1$name == "dif_surg_female" | d1$name == "dif_surg_male", "Surgical", "Hormonal")
d1$bias <- ifelse(d1$value <= 0, "Contraception", "No contraception")

p2 <- ggplot(d1, aes(x=0, xend=value, y=reorder(vertlife.species, -order), yend= reorder(vertlife.species, -order), color = bias,
                     text = paste("species:", reorder(vertlife.species, -order)))) +
  geom_segment(linewidth = 1) +
  facet_wrap(method~sex, ncol = 4) +
  scale_color_manual(values = c("#336699", "#FF9933", NA), name = "") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.ticks  = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        panel.grid.minor.x = element_blank()) +
  xlab(bquote(~italic((l[n] - l[c]) / l[n]))) +
  ylab("")

# Plot both trees (takes a while)
p1 + p2
#ggsave("RESULTS/PLOTS/PlotRect.pdf", width = 25, height = 18, units = "cm")
dev.off()


#########
##########

# get class data
spp.class <- dplyr::distinct(dat, Phy, .keep_all = TRUE) %>% dplyr::select(Class, Phy)
names(spp.class)[2] <- "Species"

spp.class2 <- left_join(tip.label,spp.class,by = "Species")

# Actinopterygii and Osteichthyes all belong to fish, although there are differnt. In our dataset, only Actinopterygii only has one estimate. To make the figure easy, we change Actinopterygii to Osteichthyes 
spp.class2$Class[spp.class2$Class == "Actinopterygii"] <- "Osteichthyes"

# version 2
p1 <- ggtree(my.tr, branch.length = "branch.length") 
# https://yulab-smu.top/treedata-book/chapter7.html  
p2 <- p1 %<+% spp.class2 + geom_tiplab(aes(color  = Class),size=3,fontface = "italic") + geom_tippoint(aes(color=Class)) + xlim_expand(c(0,3), panel = "Tree") + guides(color="none")  # xlim issue https://stackoverflow.com/questions/48928250/ggtreefacet-plot-second-panel-uses-xlim-parameter-from-first-one
p3 <- facet_plot(p2, panel = 'N observations', data = re.spp.lnRR2, 
                 geom = geom_barh, 
                 mapping = aes(x = N_obs,fill=Class,color=Class), alpha = 0.4, stat='identity') +  guides(fill="none",color="none") + # #0072B2 https://yulab-smu.top/treedata-book/chapter12.html,
  geom_facet(panel = "Estimates of mean ratio (lnRR)", data = re.spp.lnRR2,
             geom = ggstance::geom_pointrangeh, # geom_pointrange
             mapping = aes(x = Mean, xmin = Lower_bound, xmax = Upper_bound, color = Class)) + # #CC79A7
  # add overall effect MLMA_lnRR$b
  geom_facet(panel = "Estimates of mean ratio (lnRR)", data = re.spp.lnRR2,
             geom = geom_vline, 
             mapping = aes(xintercept = -0.45),linetype="dashed",color="grey") +  # add multiple panels or multiple layers on the same panels (Figure 13.1) # https://yulab-smu.top/treedata-book/chapter7.html - at Chapter 13
  geom_facet(panel = "Estimates of CV ratio (lnCVR)", data = re.spp.lnCVR2,
             geom = ggstance::geom_pointrangeh, 
             mapping = aes(x = Mean, xmin = Lower_bound, xmax = Upper_bound, color = Class)) + 
  # add overall effect MLMA_lnCVR$b
  geom_facet(panel = "Estimates of CV ratio (lnCVR)", data = re.spp.lnRR2,
             geom = geom_vline, 
             mapping = aes(xintercept = 0.21), linetype="dashed",color="grey") + 
  theme_tree2() + theme(strip.background = element_rect(fill = "white"))

# lbs <- c(Tree = "Phylogenetic tree") # change facet name # https://guangchuangyu.github.io/cn/2018/09/facet_labeller/
# p3 <- ggtree::facet_labeller(p3,lbs)
p4 <- facet_widths(p3, c(Tree=0.7,`N observations`=0.2,`Estimates of mean ratio (lnRR)`=0.4,`Estimates of CV ratio (lnCVR)`=0.4)) # set the widths of specific panels (it also supports using a name vector)

png(filename = "Figure S7 (Phy est).png", width = 8, height = 6, units = "in", type = "windows", res = 400)
p4
dev.off()


# Figure S7
# getting blups (best linear predictors from the model) for lnRR
blups <- ranef(MLMA_lnRR) 
t.spp <- blups$Spp # this needs to be added
t.phy <- blups$Phy
t.study <- blups$StudyID # study
t.spp <- rownames_to_column(t.spp, var = "Species")
t.phy <- rownames_to_column(t.phy, var = "Species")
t.spp$Species <- t.phy$Species # uppercase
t.study <- rownames_to_column(t.study, var = "Study") # study

# average species
colnames(t.spp) <- c("Species", "Deviation", "SE", "Lower_bound", "Upper_bound")
colnames(t.phy) <- c("Species", "Deviation", "SE", "Lower_bound", "Upper_bound")
colnames(t.study) <- c("Study", "Deviation", "SE", "Lower_bound", "Upper_bound")

# sorting match between study and species
dat %>% group_by(StudyID) %>% summarise(Species = unique(Phy)) -> dat_match # warning information: https://stackoverflow.com/questions/62140483/how-to-interpret-dplyr-message-summarise-regrouping-output-by-x-override

index <- match(dat_match$Species, t.spp$Species)

spp.mean <- rep(MLMA_lnRR$b, dim(t.spp)[1]) + t.spp$Deviation + t.phy$Deviation
spp.se <- sqrt(MLMA_lnRR$se^2 + t.spp$SE^2 +  t.phy$SE^2) 
spp.lb <- spp.mean - spp.se * qnorm(0.975) 
spp.ub <- spp.mean + spp.se * qnorm(0.975)
re.spp <- tibble(Species = t.spp$Species, Mean = spp.mean, SE = spp.se, Lower_bound = spp.lb, Upper_bound = spp.ub) %>% arrange(Species)
N_obs <- dat %>% group_by(Phy) %>% summarise(N_obs=n()) # calculate sample size
names(N_obs)[1] <- "Species"
re.spp.lnRR <- left_join(re.spp,N_obs,by = "Species")
tip.label <- data.frame(Species = my.tr$tip.label) # extract tip label
re.spp.lnRR2 <- left_join(tip.label,re.spp.lnRR,by = "Species") 


# getting blups (best linear predictors from the model) for lnCVR
blups <- ranef(MLMA_lnCVR) 
t.spp <- blups$Spp # this needs to be added
t.phy <- blups$Phy
t.study <- blups$StudyID # study
t.spp <- rownames_to_column(t.spp, var = "Species")
t.phy <- rownames_to_column(t.phy, var = "Species")
t.spp$Species <- t.phy$Species # uppercase
t.study <- rownames_to_column(t.study, var = "Study") # study

# average species
colnames(t.spp) <- c("Species", "Deviation", "SE", "Lower_bound", "Upper_bound")
colnames(t.phy) <- c("Species", "Deviation", "SE", "Lower_bound", "Upper_bound")
colnames(t.study) <- c("Study", "Deviation", "SE", "Lower_bound", "Upper_bound")

# sorting match between study and species
dat %>% group_by(StudyID) %>% summarise(Species = unique(Phy)) -> dat_match # warning information: https://stackoverflow.com/questions/62140483/how-to-interpret-dplyr-message-summarise-regrouping-output-by-x-override

index <- match(dat_match$Species, t.spp$Species)

spp.mean <- rep(MLMA_lnCVR$b, dim(t.spp)[1]) + t.spp$Deviation + t.phy$Deviation
spp.se <- sqrt(MLMA_lnCVR$se^2 + t.spp$SE^2 +  t.phy$SE^2) 
spp.lb <- spp.mean - spp.se * qnorm(0.975) 
spp.ub <- spp.mean + spp.se * qnorm(0.975)
re.spp <- tibble(Species = t.spp$Species, Mean = spp.mean, SE = spp.se, Lower_bound = spp.lb, Upper_bound = spp.ub) %>% arrange(Species)
N_obs <- dat %>% group_by(Phy) %>% summarise(N_obs=n()) # calculate sample size
names(N_obs)[1] <- "Species"
re.spp.lnCVR <- left_join(re.spp,N_obs,by = "Species")
tip.label <- data.frame(Species = my.tr$tip.label) # extract tip label
re.spp.lnCVR2 <- left_join(tip.label,re.spp.lnCVR,by = "Species") 



# get class data
spp.class <- dplyr::distinct(dat, Phy, .keep_all = TRUE) %>% dplyr::select(Class, Phy)
names(spp.class)[2] <- "Species"

spp.class2 <- left_join(tip.label,spp.class,by = "Species")

# Actinopterygii and Osteichthyes all belong to fish, although there are differnt. In our dataset, only Actinopterygii only has one estimate. To make the figure easy, we change Actinopterygii to Osteichthyes 
spp.class2$Class[spp.class2$Class == "Actinopterygii"] <- "Osteichthyes"

# version 2
p1 <- ggtree(my.tr, branch.length = "branch.length") 
# https://yulab-smu.top/treedata-book/chapter7.html  
p2 <- p1 %<+% spp.class2 + geom_tiplab(aes(color  = Class),size=3,fontface = "italic") + geom_tippoint(aes(color=Class)) + xlim_expand(c(0,3), panel = "Tree") + guides(color="none")  # xlim issue https://stackoverflow.com/questions/48928250/ggtreefacet-plot-second-panel-uses-xlim-parameter-from-first-one
p3 <- facet_plot(p2, panel = 'N observations', data = re.spp.lnRR2, 
                 geom = geom_barh, 
                 mapping = aes(x = N_obs,fill=Class,color=Class), alpha = 0.4, stat='identity') +  guides(fill="none",color="none") + # #0072B2 https://yulab-smu.top/treedata-book/chapter12.html,
  geom_facet(panel = "Estimates of mean ratio (lnRR)", data = re.spp.lnRR2,
             geom = ggstance::geom_pointrangeh, # geom_pointrange
             mapping = aes(x = Mean, xmin = Lower_bound, xmax = Upper_bound, color = Class)) + # #CC79A7
  # add overall effect MLMA_lnRR$b
  geom_facet(panel = "Estimates of mean ratio (lnRR)", data = re.spp.lnRR2,
             geom = geom_vline, 
             mapping = aes(xintercept = -0.45),linetype="dashed",color="grey") +  # add multiple panels or multiple layers on the same panels (Figure 13.1) # https://yulab-smu.top/treedata-book/chapter7.html - at Chapter 13
  geom_facet(panel = "Estimates of CV ratio (lnCVR)", data = re.spp.lnCVR2,
             geom = ggstance::geom_pointrangeh, 
             mapping = aes(x = Mean, xmin = Lower_bound, xmax = Upper_bound, color = Class)) + 
  # add overall effect MLMA_lnCVR$b
  geom_facet(panel = "Estimates of CV ratio (lnCVR)", data = re.spp.lnRR2,
             geom = geom_vline, 
             mapping = aes(xintercept = 0.21), linetype="dashed",color="grey") + 
  theme_tree2() + theme(strip.background = element_rect(fill = "white"))

# lbs <- c(Tree = "Phylogenetic tree") # change facet name # https://guangchuangyu.github.io/cn/2018/09/facet_labeller/
# p3 <- ggtree::facet_labeller(p3,lbs)
p4 <- facet_widths(p3, c(Tree=0.7,`N observations`=0.2,`Estimates of mean ratio (lnRR)`=0.4,`Estimates of CV ratio (lnCVR)`=0.4)) # set the widths of specific panels (it also supports using a name vector)

png(filename = "Figure S7 (Phy est).png", width = 8, height = 6, units = "in", type = "windows", res = 400)
p4
dev.off()





