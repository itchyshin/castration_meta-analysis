# figures


library(tidyverse)
library(png)
library(metafor)
library(grid)

# blup plot

saveRDS(mod, here(Rdta, mod.rds))

mod <- readRDS(here(Rdta, mod.rds))

# getting blups (best linear predictors from the model)
blups <- ranef(mod) 
t.spp <- blups$Species_Latin # this needs to be added
t.phy <- blups$Phylogeny 
t.study <- blups$Study # study
t.spp <- rownames_to_column(t.spp, var = "Species")
t.phy <- rownames_to_column(t.phy, var = "Species")
t.study <-  rownames_to_column(t.study, var = "Study") # study

# average cultivar
colnames(t.spp) <- c("Species", "Deviation", "SE", "Lower_bound", "Upper_bound")
colnames(t.phy) <- c("Species", "Deviation", "SE", "Lower_bound", "Upper_bound")
colnames(t.study) <- c("Study", "Deviation", "SE", "Lower_bound", "Upper_bound")

# sorting match between study and species
dat %>% group_by(Study) %>% summarise(Species_Latin = unique(Species_Latin)) -> dat_match

dat_match$Species_Latin 
index <- match(dat_match$Species_Latin , t.spp$Species)


# knitr::kable(t.cultivar) col.names = c(’Cultivar’,
# ’Deviation’, ’SE’, ’Lower bound’, ’Upper bound’))
spp.mean <- rep(mod$b, dim(t.spp)[1]) + t.spp$Deviation + t.phy$Deviation
spp.se <- sqrt(mod$se^2 + t.spp$SE^2 +  t.phy$SE^2) 
spp.lb <- spp.mean - spp.se * qnorm(0.975) 
spp.ub <- spp.mean + spp.se * qnorm(0.975)
t.spp2 <- tibble(Species = t.spp$Species, Mean = spp.mean, SE = spp.se, Lower_bound = spp.lb, Upper_bound = spp.ub) %>% arrange(Species)

# study wise 

study.mean <- rep(mod$b, dim(t.study)[1]) + t.spp$Deviation[index] + t.phy$Deviation[index] + t.study$Deviation
study.se <- sqrt(mod$se^2 + t.spp$SE[index]^2 +  t.phy$SE[index]^2 + t.study$SE^2) 
study.lb <- study.mean - study.se * qnorm(0.975) 
study.ub <- study.mean + study.se * qnorm(0.975)
t.study2 <- tibble(Study = t.study$Study, Species = t.spp$Species[index], 
                   Study_Spp = as.factor(paste(Study, Species, sep = "_")), Mean = study.mean, 
                   SE = study.se, Lower_bound = study.lb, Upper_bound = study.ub) %>% arrange(Species)



# plotting

# getting all photos 
# reading in photos
filenames <- list.files("icons", pattern=".png", full.names=TRUE)
ldf <- lapply(filenames, readPNG)
names(ldf) <- substr(filenames, 7, 60)
#name <- substr(filenames, 7, 60)


spp.plot <- ggplot(data = t.spp2, aes(x = Mean, y = Species)) +
  geom_errorbarh(aes(xmin = Lower_bound, xmax = Upper_bound),  
                 height = 0, show.legend = FALSE, size = 0.5, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.5) +
  geom_vline(xintercept = mod$b, linetype = 1, colour = "red") +
  geom_point(aes(fill = Species), size = 3, shape = 21) + 
  xlim(-0, 0.25) +
  theme_bw() +
  labs(x = "lnRR (effect size)", y = "") + 
  theme(legend.position= c(0.99, 0.01), legend.justification = c(1, 0)) +
  theme(legend.title = element_text(size = 9)) +
  theme(legend.direction="horizontal") +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(size = 10, colour ="black",
                                   hjust = 0.5)) +
  guides(fill = "none")  + 
  #annotation_custom(rasterGrob(ldf), xmin = rep(0.01, 14), xmax = rep(0.05, 14), ymin = seq(0.5, 13.5, 1), ymax = seq(1.5, 14.5, 1)) 
  annotation_custom(rasterGrob(ldf[[1]]), xmin = 0.01, xmax = 0.05, ymin = 0.5, ymax = 1.5) +
  annotation_custom(rasterGrob(ldf[[2]]), xmin = 0.01, xmax = 0.05, ymin = 1.5, ymax = 2.5) +
  annotation_custom(rasterGrob(ldf[[3]]), xmin = 0.01, xmax = 0.05, ymin = 2.5, ymax = 3.5) +
  annotation_custom(rasterGrob(ldf[[4]]), xmin = 0.01, xmax = 0.05, ymin = 3.5, ymax = 4.5) +
  annotation_custom(rasterGrob(ldf[[5]]), xmin = 0.01, xmax = 0.05, ymin = 4.5, ymax = 5.5) +
  annotation_custom(rasterGrob(ldf[[6]]), xmin = 0.01, xmax = 0.05, ymin = 5.5, ymax = 6.5) +
  annotation_custom(rasterGrob(ldf[[7]]), xmin = 0.01, xmax = 0.05, ymin = 6.5, ymax = 7.5) +
  annotation_custom(rasterGrob(ldf[[8]]), xmin = 0.01, xmax = 0.05, ymin = 7.5, ymax = 8.5) +
  annotation_custom(rasterGrob(ldf[[9]]), xmin = 0.01, xmax = 0.05, ymin = 8.5, ymax = 9.5) +
  annotation_custom(rasterGrob(ldf[[10]]), xmin = 0.01, xmax = 0.05, ymin = 9.5, ymax = 10.5) +
  annotation_custom(rasterGrob(ldf[[11]]), xmin = 0.01, xmax = 0.05, ymin = 10.5, ymax = 11.5) +
  annotation_custom(rasterGrob(ldf[[12]]), xmin = 0.01, xmax = 0.05, ymin = 11.5, ymax = 12.5) +
  annotation_custom(rasterGrob(ldf[[13]]), xmin = 0.01, xmax = 0.05, ymin = 12.5, ymax = 13.5) +
  annotation_custom(rasterGrob(ldf[[14]]), xmin = 0.01, xmax = 0.05, ymin = 13.5, ymax = 14.5)

spp.plot
###################
# study specific effect
####################


study.plot <- ggplot(data = t.study2, aes(x = Mean, y = fct_reorder(Study, Species))) +
  geom_errorbarh(aes(xmin = Lower_bound, xmax = Upper_bound),  
                 height = 0, show.legend = FALSE, size = 0.5, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.5) +
  geom_vline(xintercept = mod$b, linetype = 1, colour = "red") +
  geom_point(aes(fill = Species), size = 3, shape = 21) + 
  xlim(-0.5, 0.8) +
  theme_bw() +
  labs(x = "lnRR (effect size)", y = "") + 
  theme(legend.position= c(0.99, 0.01), legend.justification = c(1, 0)) +
  theme(legend.title = element_text(size = 9)) +
  theme(legend.direction="horizontal") +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(size = 10, colour ="black",
                                   hjust = 0.5)) +
  guides(fill = "none")  + 
  #annotation_custom(rasterGrob(ldf), xmin = rep(0.01, 14), xmax = rep(0.05, 14), ymin = seq(0.5, 13.5, 1), ymax = seq(1.5, 14.5, 1)) 
  annotation_custom(rasterGrob(ldf[[1]]), xmin = -0.5, xmax = -0.4, ymin = 1.5, ymax = 3.5) +
  annotation_custom(rasterGrob(ldf[[2]]), xmin = -0.45, xmax = -0.35, ymin = 3+ 4, ymax = 5 + 4) +
  annotation_custom(rasterGrob(ldf[[3]]), xmin = -0.5, xmax = -0.4, ymin = 0 + 12, ymax = 2 + 12) +
  annotation_custom(rasterGrob(ldf[[4]]), xmin = -0.45, xmax = -0.35, ymin = 2 + 14, ymax = 4 + 14) +
  annotation_custom(rasterGrob(ldf[[5]]), xmin = -0.5, xmax = -0.4, ymin = 2.5 + 18, ymax = 4.5 + 18) +
  annotation_custom(rasterGrob(ldf[[6]]), xmin = -0.45, xmax = -0.35, ymin = 0.5 + 23, ymax = 2.5 + 23) +
  annotation_custom(rasterGrob(ldf[[7]]), xmin = -0.5, xmax = -0.4, ymin = 4 + 25, ymax = 6 + 25) +
  annotation_custom(rasterGrob(ldf[[8]]), xmin = -0.45, xmax = -0.35, ymin = 0 + 35, ymax = 2 + 35) +
  annotation_custom(rasterGrob(ldf[[9]]), xmin = -0.5, xmax = -0.4, ymin = 0.5 + 36, ymax = 2.5 + 36) +
  annotation_custom(rasterGrob(ldf[[10]]), xmin = -0.45, xmax = -0.35, ymin = 0 + 38, ymax = 2 + 38) +
  annotation_custom(rasterGrob(ldf[[11]]), xmin = -0.5, xmax = -0.4, ymin = 0 + 39, ymax = 2+ 39) +
  annotation_custom(rasterGrob(ldf[[12]]), xmin = -0.45, xmax = -0.35, ymin = 0.5 + 40, ymax = 2.5 + 40) +
  annotation_custom(rasterGrob(ldf[[13]]), xmin = -0.5, xmax = -0.4, ymin = 3.5 + 42, ymax = 5.5 + 42) +
  annotation_custom(rasterGrob(ldf[[14]]), xmin = -0.45, xmax = -0.35, ymin = 0 + 49, ymax = 2 + 49)

study.plot

