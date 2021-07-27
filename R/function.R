# custom functions 

# function for getting lnRR for proportional data

lnrrp <- function(m1, m2, n1, n2) {
  # arcsine transforamtion
  asin_trans <- function(p) { asin(sqrt(p)) }
  # SD for arcsine distribution (see Wiki - https://en.wikipedia.org/wiki/Arcsine_distribution)
  sd1 <- sqrt(1/8)
  sd2 <- sqrt(1/8)
  # lnRR - with 2nd order correction
  lnrr <- log(asin_trans(m1)/asin_trans(m2)) + 
    0.5 * (((sd1^2) / (n1 * m1^2)) - ((sd2^2) / (n2 * m1^2)))	
  
  var <- sd1^2 / (n1 * m1^2) + sd1^4 / (2 * n1^2 * m2^4)  + 
    sd2^2 / (n2 * m2^2) + sd2^4 / (2 * n2^2 * m2^4) 
  
  invisible(data.frame(yi = lnrr , vi = var))
}

# function for getting lnRR for mean data

lnrrm <- function(m1, m2, n1, n2, sd1, sd2) {
  # lnRR - with 2nd order correction
  lnrr <- log(m1/m2) + 
    0.5 * (((sd1^2) / (n1 * m1^2)) - ((sd2^2) / (n2 * m1^2)))	
  
  var <- sd1^2 / (n1 * m1^2) + sd1^4 / (2 * n1^2 * m2^4)  + 
    sd2^2 / (n2 * m2^2) + sd2^4 / (2 * n2^2 * m2^4) 
  
  invisible(data.frame(yi = lnrr , vi = var))
}


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
