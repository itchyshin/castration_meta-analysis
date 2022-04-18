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
               here
)

########################
# data: life expectancy
########################

dat <- read_csv(here("data", "zoo.csv"), na = c("", "NA"))


glimpse(dat)
names(dat)

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
         species = species)

  
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

mod_m_surg <- rma.mv(yi, V = vi, 
                     random = list(~1|species),
                     data = dat_m_surg)
summary(mod_m_surg)
i2_ml(mod_m_surg)

orchard_plot(mod_m_surg,xlab = "lnRR (male surgical)", group = "species", data = dat_m_surg, g = FALSE)

# there is an outliner
dat_m_surg[which(dat_m_surg$yi == max(dat_m_surg$yi)), "species"]
#[1] "Pseudocheirus peregrinus"

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
         species = species)


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

mod_f_horm <- rma.mv(yi, V = vi, 
                     random = list(~1|species),
                     data = dat_f_horm)
summary(mod_f_horm)
i2_ml(mod_f_horm)

orchard_plot(mod_f_horm,xlab = "lnRR (female hormonal)", group = "species", data = dat_f_horm, g = FALSE)

######################
# female surgical data
######################
data_f_surg<- dat %>% filter(is.na(Female_Surgical_Mean) == FALSE, is.na(Female_None_Mean) == FALSE) %>% 
  mutate(F_control_m = Female_None_Mean,
         F_control_sd = sqrt(Female.None)*Female_None_SE,
         F_control_n = Female.None,
         F_surgical_m = Female_Surgical_Mean,
         F_surgical_sd = sqrt(Female.Surgical)*Female_Surgical_SE,
         F_surgical_n = Female.Surgical,
         species = species)


# getting effect size

data_f_surg <- escalc("ROM", 
                     m1i = F_surgical_m,
                     m2i = F_control_m,
                     sd1i = F_surgical_sd,
                     sd2i = F_control_sd,
                     n1i = F_surgical_n,
                     n2i = F_control_n,
                     data = data_f_surg,
)

mod_f_surg <- rma.mv(yi, V = vi, 
                     random = list(~1|species),
                     data = data_f_surg)
summary(mod_f_surg)
i2_ml(mod_f_surg)

orchard_plot(mod_f_surg, xlab = "lnRR (female surgical)", group = "species", data = data_f_surg, g = FALSE)

##############################
# data: gamma - mortality risk
##############################

dat2 <- read.csv(here("data", "gammaBaSTA.csv"))

#dat2 <- tibble(dat2)
glimpse(dat2)
names(dat2)

# filtering and getting SD

################
# male surgical
################

dat2_m_surg <- dat2 %>% filter(is.na(Male_Surgical_Mean) == FALSE) %>% 
  mutate(yi = Male_Surgical_Mean ,
         vi = Male_Surgical_SE^2,
         species = species)

hist(dat2_m_surg$yi)

mod2_m_surg <- rma.mv(yi, V = vi, 
                     random = list(~1|species),
                     data = dat2_m_surg)
summary(mod2_m_surg)
i2_ml(mod2_m_surg)

orchard_plot(mod2_m_surg,xlab = "lnHR (male surgical)", group = "species", data = dat2_m_surg, g = FALSE)

# there is an outliner
dat2_m_surg[which(dat2_m_surg$yi == min(dat2_m_surg$yi)), "species"]
#[1] "Pseudocheirus peregrinus"

#######################
# female hormonal data 
########################

dat2_f_horm <- dat2 %>% filter(is.na(Female_Hormonal_Mean) == FALSE) %>% 
  mutate(yi = Female_Hormonal_Mean ,
         vi = Female_Hormonal_SE^2,
         species = species)

hist(dat2_f_horm$yi)


mod2_f_horm <- rma.mv(yi, V = vi, 
                     random = list(~1|species),
                     data = dat2_f_horm)
summary(mod2_f_horm)
i2_ml(mod2_f_horm)

orchard_plot(mod2_f_horm,xlab = "lnHR (female hormonal)", group = "species", data = dat2_f_horm, g = FALSE)

######################
# female surgical data
######################
dat2_f_surg<- dat2 %>% filter(is.na(Female_Surgical_Mean) == FALSE) %>% 
  mutate(yi = Female_Surgical_Mean ,
         vi = Female_Surgical_SE^2,
         species = species)

hist(dat2_f_surg$yi)


mod2_f_surg <- rma.mv(yi, V = vi, 
                     random = list(~1|species),
                     data = dat2_f_surg)
summary(mod2_f_surg)
i2_ml(mod2_f_surg)

orchard_plot(mod2_f_surg, xlab = "lnHR (female surgical)", group = "species", data = dat2_f_surg, g = FALSE)
