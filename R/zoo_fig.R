# figs for zoo data

rm(list = ls())
# pacakges
#remotes::install_github("GuangchuangYu/tidytree")
#if(!require("tidytree", quietly=TRUE))
#devtools::install_version("tidytree", version = "0.4.4")

library(tidyverse)

# Phylogenetic tree:

library(ape)
library(ggtree)
library(ggimage)
library(ggstance)
library(tidytree)
library(here)

# Plot layout:
library(patchwork)
library(cowplot)

dat_full <- read_csv(here("data", "zoo.csv"), na = c("", "NA"))

# trimming data

dim(dat_full)
# [1] "Chrysocyon_brachyurus" - maned wolf
# "Crocuta_crocuta" - spotted hyena
# "Panthera_uncia" - snow lepard
# "Neofelis_nebulosa" - clouded leopard
# "Pseudocheirus_peregrinus" - posumm


# dat$species == "zims.species"


# phylogeny == "vertlife.species"

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

dat$vertlife.species[which(dat$species == "Cervus canadensis")] <-"Cervus canadensis"
dat$phylogeny[which(dat$species == "Cervus canadensis")] <-"Cervus_canadensis"


dim(dat)




# getting relevant data
lifespan <- readRDS(here("Rdata", "lifespan_all.RDS"))



#setdiff(dat_full$species, lifespan$species)
#########

#lifespan_m <- lifespan %>% filter(sex == "male")
#lifespan_f <- lifespan %>% filter(sex == "female")

#dat$Phylogeny <- gsub("Perca_fluviatilis", "Lamperta_fluviatilis", dat$Phylogeny) #replace with the original name
# tree <- read.tree(here("data/tree_zoo4.tre"))
# 
# fortify(tree)
# # triming tree
# 
#to_drop <-
#  tree$tip.label[which(!(tree$tip.label %in% unique(lifespan$phylogeny)))]

#tree <- drop.tip(tree, to_drop)
# 
# # this works
# fortify(tree)

#write.tree(tree, here("data", "tree_zoo_all.tre"))


#tree <- read.tree( here("data", "tree_zoo_all.tre"))
tree <- read.tree( here("data", "tree_zoo4.tre"))

to_drop <-
  tree$tip.label[which(!(tree$tip.label %in% unique(lifespan$phylogeny)))]

tree <- drop.tip(tree, to_drop)

fortify(tree)

#tree <- compute.brlen(tree)
cor_tree <- vcv(tree, corr = TRUE)


# testing matching

match(dat$species, lifespan$species)
match(tree$tip.label, lifespan$phylogeny)
match(tree$tip.label, dat$phylogeny)
setdiff(tree$tip.label, dat$phylogeny)


# manipulating tree

tdat <- fortify(tree)
tdat2 <- subset(tdat, isTip)

# this ordering matches the tree
ordered_spp <- tdat2$label[order(tdat2$y, decreasing=TRUE)]



#############################
#  creating data to plot 
#############################

pdat <- dat %>% select(species, vertlife.species, order, 
                       Male_None_Mean, Female_None_Mean,
                       Male_Surgical_Mean,Female_Surgical_Mean, 
                       Male_Hormonal_Mean, Female_Hormonal_Mean, 
                       Male_Immunological_Mean, Female_Immunological_Mean
                       )

pdat %>% mutate(lrr_surgical_m = (Male_Surgical_Mean/Male_None_Mean - 1)*100,
                lrr_surgical_f = (Female_Surgical_Mean/Female_None_Mean - 1)*100,
                lrr_hormonal_m = (Male_Hormonal_Mean/Male_None_Mean - 1)*100,
                lrr_hormonal_f = (Female_Hormonal_Mean/Female_None_Mean - 1)*100,
                lrr_immunological_m = (Male_Immunological_Mean/Male_None_Mean - 1)*100,
                lrr_immunological_f = (Female_Immunological_Mean/Female_None_Mean - 1)*100,
                #) -> pdat
                phylogeny = gsub(" ","_", vertlife.species)) -> pdat


# matching the tree ordering and data ordering
pos_order <- match(ordered_spp, pdat$phylogeny)
pdat <- pdat[pos_order, ]
pdat$ordering <- seq(1:nrow(pdat))
#rownames(pdat) <- gsub(" ", "_", pdat$vertlife.species)

# reformat data to long format
pdat <- pdat %>% select(species, vertlife.species, order, phylogeny, ordering,
                        lrr_surgical_m, lrr_surgical_f,
                        lrr_hormonal_m, lrr_hormonal_f,
                        lrr_immunological_m, lrr_immunological_f)

pdat_long <- pdat %>% pivot_longer(cols = c(lrr_surgical_m, lrr_surgical_f,
                                            lrr_hormonal_m, lrr_hormonal_f,
                                            lrr_immunological_m, lrr_immunological_f))

pdat_long %>% mutate(sex = rep(c("male", "female"), nrow(pdat_long)/2),
                      type = rep(c("surgical", "surgical", 
                                   "hormonal", "hormonal", 
                                   "immunological", "immunological"), nrow(pdat_long)/6)) -> pdat_long

pdat_long %>% mutate(min_value = ifelse(is.na(value) == TRUE, NA, ifelse (value < 0, value, 0)), 
                     max_value = ifelse(is.na(value) == TRUE, NA, ifelse (value > 0, value, 0))) -> pdat_long

pdat_long$type <- factor(pdat_long$type, levels = c("surgical", "hormonal", "immunological"),
                         labels = c("Surgical", "Hormonal", "Immunological"))



# # use getMRCA() to obtain common ancestor nodes to position the order silhouttes
tree.tibble <- tidytree::as_tibble(tree)
ord <- unique(tax$order)

dford <- data.frame(order = ord, node = NA)

ldat <- vector("list", 10)

for (i in ord) {

  tip <-  as.vector(as.data.frame(pdat[which(pdat$order == i), "phylogeny"]))
  #tip <-  gsub(" ", "_", tip)
  if (length(tip) > 1) {
    tnode <- numeric(getMRCA(tree, tip = tip$phylogeny))

  } else 
    if (length(tip) == 1) {
   tnode <-
     match(tip$phylogeny, tree.tibble$label)
  }
  
  dford[[2]][dford[[1]] == i] <- tnode

}


tree <- as.ultrametric(tree)

p <-
  ggtree(
    tree,
    layout = "rectangular",
    ladderize = TRUE,
    size = 0.3,
    color = "#454747"
  )

pdat %>% select(phylogeny, order) -> odat

p0 <- p %<+% odat + 
  geom_tippoint(aes(color= order), size = 0.5) + xlim_expand(c(0,170), panel = "Tree") + 
  guides(color="none") 

# for (i in 1:nrow(dford)) {
#   p <-
#     p + geom_cladelabel2(
#       node = dford[i, "node"],
#       label = "", 
#       #label = dford[i, "order"], # uncomment to plot order names
#       hjust = 0.4,
#       align = T,
#       offset.text = 10,
#       barsize = 1,
#       fontsize = 2.5,
#       fill = "grey20",
#       color = "grey50",
#       angle = 0,
#       horizontal = TRUE,
#     )
# }
# p


# Add order silhouettes
# (takes a while to display)

# Image directory
imgdir <- "phylopics/"

p1 = p0 +
  geom_image(
    x = 170,
    y = 148,
    image = paste0(
      imgdir,
      "Artiodactyla_PhyloPic.8b567be8.An-Ignorant-Atheist.Antilopinae.png"
    ),
    size = 0.045
  ) +
  geom_image(
    x = 170,
    y = 116,
    image = paste0(
      imgdir,
      "Perissodactlyla_PhyloPic.071ee517.David-Orr.Equus-ferus-caballus.png"
    ),
    size = 0.05
  ) +
  geom_image(
    x = 170,
    y = 100,
    image = paste0(
      imgdir,
      "Carnivora_PhyloPic.34e482b4.An-Ignorant-Atheist.Panthera.png"
    ),
    size = 0.05
  ) +
  geom_image(
    x = 170,
    y = 80,
    image = paste0(
      imgdir,
      "Chrioptera_PhyloPic.e7da460a.Margot-Michaud.Chiroptera_Eptesicus_Eptesicus-fuscus_Vespertilio-Noctilio_Vespertilionidae_Vespertilioniformes_Vespertilioninae_Vespertilionoidea.png"
    ),
    size = 0.05
  ) +
  geom_image(
    x = 170,
    y = 48,
    image = paste0(
      imgdir,
      "Primates_PhyloPic.071db0d0.Margot-Michaud.Papio_Papio-anubis.png"
    ),
    size = 0.05
  ) +
  geom_image(
    x = 170,
    y = 19,
    image = paste0(
      imgdir,
      "Rodentia_PhyloPic.570c7d9e.Alexandra-van-der-Geer.Rattus_Rattus-exulans.png"
    ),
    size = 0.06
  ) +
  geom_image(
    x = 170,
    y = 7,
    image = paste0(
      imgdir,
      "PhyloPic.b62bab6e.An-Ignorant-Atheist.Macropus-Macropus.png"
    ),
    size = 0.05
  )

#p1



# putting into 
p2 <- ggplot(pdat_long,  
             aes(x = value, reorder(vertlife.species, -ordering))) +
  ggplot2::geom_errorbarh(aes(xmin = min_value, xmax = max_value, colour = sex), 
                 height = 0, show.legend = TRUE, size = 1, 
                 alpha = 0.8, position =position_dodge(width = 0.75)) +
  #geom_segment(linewidth = 1, alpha = 0.8, position =position_dodge(width = 0.95)) +
  facet_wrap(~type, ncol = 3, #scales = "free_x"
             ) +
  #facet_wrap(sex~type, ncol = 6) +
  #scale_color_manual(values = c("#FF9933", "#336699", NA), name = "") +
  scale_color_manual(values = c("#CC6677", "#88CCEE")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.ticks  = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor.x = element_blank()) +
  #guides(colour = colour) +
  theme(legend.position =  c(0.88, 0.1)) +
  #labs(x = "lnRR (effect size)", y = "", ) +
  xlab("Ratio: sterilization/control (%)") +
  ylab("") 

#p2
#library(ggpubr)
#p_phylo <- ggarrange(p1, p2)

p_phylo <-p1 + p2

##############