# custom functions 

# function for getting lnRR for proportional data

lnrrp <- function(m1, m2, n1, n2) {
  # arcsine transforamtion
  asin_trans <- function(p) { asin(sqrt(p)) }
  # SD for arcsine distribution (see Wiki - https://en.wikipedia.org/wiki/Arcsine_distribution)
  var1 <- 1/8
  var2 <- 1/8
  # lnRR - with 2nd order correction
  lnrr <- log(asin_trans(m1)/asin_trans(m2)) + 
    0.5 * ((var1 / (n1 * asin_trans(m1)^2)) - (var2 / (n2 * asin_trans(m2)^2)))	
  
  var <- var1 / (n1 * asin_trans(m1)^2) + var1^2 / (2 * n1^2 * asin_trans(m1)^4)  + 
    var2 / (n2 * asin_trans(m2)^2) + var2^2 / (2 * n2^2 * asin_trans(m2)^4) 
  
  invisible(data.frame(yi = lnrr , vi = var))
}

# function for getting lnRR for mean data

lnrrm <- function(m1, m2, n1, n2, sd1, sd2) {
  # lnRR - with 2nd order correction
  lnrr <- log(m1/m2) + 
    0.5 * (((sd1^2) / (n1 * m1^2)) - ((sd2^2) / (n2 * m2^2)))	
  
  var <- sd1^2 / (n1 * m1^2) + sd1^4 / (2 * n1^2 * m2^4)  + 
    sd2^2 / (n2 * m2^2) + sd2^4 / (2 * n2^2 * m2^4) 
  
  invisible(data.frame(yi = lnrr , vi = var))
}


# function to get to 

lnrrm2 <- function(m1, m2, n1, n2, cv21, cv22) {
  # lnRR - with 2nd order correction
  lnrr <- log(m1/m2) + 
    0.5 * ((cv21 /n1) - (cv22 / n2))	
  
  var <- (cv21 / n1) + ((cv21^2) / (2 * n1^2))  + 
    (cv22/ n2) + ((cv22^2) / (2 * n2^2) )
  
  invisible(data.frame(yi = lnrr , vi = var))
}

# get 

# function for getting lnVR for mean data

lnvrm <- function(n1, n2, sd1, sd2){
  lnvr<-log(sd1 / sd2) + (1 / (2 * (n1 - 1))) - (1 / (2 * (n2 - 1))) 
  var<- 1 / (2 * (n1 - 1)) + 1 / (2 * (n1 - 1)^2) + 
    1 / (2 * (n2 - 1)) + 1 / (2 * (n2 - 1)^2)
  
  invisible(data.frame(yi = lnvr , vi = var))
}
# function for getting lnCVR for mean data

lncvrm <- function(m1, m2, n1, n2, sd1, sd2){
  yi<-log((sd1 / m1) / (sd2 / m2)) + (1 / (2 * (n1 - 1))) - (1 / (2 * (n2 - 1))) 
  lncvr<-yi + 0.5 * (((sd2^2) / (n2 * m2^2)) - ((sd1^2) / (n1 * m1^2)))	
  var<-sd1^2 / (n1 * m1^2) + sd1^4 / (2 * n1^2 * m1^4) +
    1 / (2 * (n1 - 1)) + 1 / (2 * (n1 - 1)^2) + 
    sd2^2 / (n2 * m2^2) + m2^4 / (2 * n2^2 * m2^4) +
    1 / (2 * (n2 - 1)) + 1 / (2 * (n2 - 1)^2)
  
  invisible(data.frame(yi = lncvr , vi = var))
}

# function for getting d for proportional data

smdp <- function(m1, m2, n1, n2) {
  mean_diff <- (car::logit(m1) - car::logit(m2))
  j_cor <- 1 - (3 /(4*(n1 + n2) -9))
  smd <- (mean_diff/(pi/sqrt(3)))*j_cor
  var <- ((n1 + n2)/(n1*n2) + (smd^2)/(2*(n1 + n2)))*(j_cor^2)
  invisible(data.frame(yi = smd , vi = var))
}

d_to_lnor <- function(smd){
  lnor <- smd * (sqrt(3)/pi)
  return(lnor)
} 

# funciont for getting d for traditional data

smdm <- function(m1, m2, n1, n2, sd1, sd2) {
  mean_diff <- m1 - m2
  sd_pool <- sqrt( ((n1 - 1)*sd1^2 + (n2 - 1)*sd2) / (n1 + n2 - 2) )
  j_cor <- 1 - (3 /(4*(n1 + n2) -9))
  smd <- (mean_diff/sd_pool)*j_cor
  var <- ((n1 + n2)/(n1*n2) + (smd^2)/(2*(n1 + n2)))*(j_cor^2)
  invisible(data.frame(yi = smd , vi = var))
}

d_to_lnor <- function(smd){
  lnor <- smd * (sqrt(3)/pi)
  return(lnor)
} 

lnor_to_d <- function(lnor){
  smd <- lnor*(pi/sqrt(3))
}
# invisible

# function to get lnrr for missing SD

# TODO - modify
# Below is the custom function to calculate the lnRR 
lnRR_func <- function(Mc, Nc, Me, Ne, aCV2c, aCV2e, rho = 0.5){
  lnRR <- log(Me/Mc) + 
    0.5 * ((aCV2e/Ne) - (aCV2c/Nc))	
  
  var_lnRR <- (aCV2c/Nc) + (aCV2e/Ne) - 
    2*rho*(sqrt(aCV2c)*sqrt(aCV2e)/(Nc)) 
  
  data.frame(lnRR,var_lnRR)
}


# getting CV2 

aCV2 <- dat %>% 
  group_by(Study_ID) %>%  # Group by study 
  summarise(CV2c = mean((SDc/Mc)^2, na.rm = T),  # Calculate the squared coefficient of variation for control and experimental groups
            CV2e = mean((SDe/Me)^2, na.rm = T)) %>% 
  ungroup() %>% # ungroup 
  summarise(aCV2c = mean(CV2c, na.rm = T), # Mean CV^2 for exp and control groups across studies
            aCV2e = mean(CV2e, na.rm = T)) 

effect <- lnRR_func(Mc = dat$Mc, 
                    Nc = dat$Nc, 
                    Me = dat$Me, 
                    Ne = dat$Ne, 
                    aCV2c = aCV2[[1]], 
                    aCV2e = aCV2[[2]],
                    rho = 0.8)  # Calculate effect sizes

# additional functions for something else

lnRR_arcsin <- function(m1, m2, n1, n2) { # m1 and m2 should be between 0 and 1
  # arcsine transformation
  asin_trans <- function(p) { asin(sqrt(p)) }
  # use the delta method - but keep observed
  #var1 <- sd1^2/(4*m1*(1-m1))
  #var2 <- sd2^2/(4*m2*(1-m2))
  
  # assuming theoritcal SD
  var1 <- 1/(8)
  var2 <- 1/(8)
  
  # lnRR - with 2nd order correction
  lnrr <- log(asin_trans(m1)/asin_trans(m2)) + 
    0.5 * ((var1 / (n1 * asin_trans(m1)^2)) - (var2 / (n2 * asin_trans(m2)^2)))	
  
  var <- var1 / (n1 * asin_trans(m1)^2) + var1^2 / (2 * n1^2 * asin_trans(m1)^4)  + 
    var2 / (n2 * asin_trans(m2)^2) + var2^2 / (2 * n2^2 * asin_trans(m2)^4) 
  
  invisible(data.frame(yi = lnrr , vi = var))
}


