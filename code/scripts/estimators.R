# source("code/R functions/mid_p-randomized_pvalues.R")

# R function to compute the estimators 
####################################### non discrete #######################################
# return m0_hat and not pi0_hat

# PC (modified version of ZZD)
PC_ZZD <- function(pvalues, C, s) {
  return(C * min(length(pvalues), max(s, 2 * sum(pvalues))))
}

# PC (our modified PC) 
PC_new_DM <- function(pvalues) {
  return(2 * (1 + sum(pvalues)))
}

# Storey 
Storey <- function(pvalues, lambda) {
  sum_ <- sum(pvalues > lambda)
  cor <- 1 / (1 - lambda)
  return(cor * (1 +  sum_))
}

# hybrid PC-Storey

PC_Storey <- function(pvalues, lambda) {
  return((2 / (1 - lambda**2)) * (1 + sum(pvalues * (pvalues > lambda))))
  
}

# hybrid Quad-Storey

Quad_Storey <- function(pvalues, lambda) {
  return((3 / (1 - lambda**3)) * (1 + sum(pvalues**2 * (pvalues > lambda))))
}

quad_PC <- function(pvalues) {
  return(3 * (1 + sum(pvalues **2)))
}


# Polynomial Storey-PC estimator

poly_PC_STOREY <- function(pvalues, lambda, r) {
  return(( (r + 1) / (1 - (lambda**(r + 1)))) * (1 + sum((pvalues**r) * (pvalues > lambda))))
}

new_poly <- function(pvalues, lambda, r) {
 return( (3/(2*(1 - lambda))) * (1 + sum(as.numeric(((1- ((1 - pvalues) / (1 - lambda))**2) * (pvalues > lambda))))))
}

# Poi <- function(pvalues, lambda) {
#   return((1 + sum(qpois(pvalues, lambda))) / length(pvalues))
# }


####################################### generic discrete adjustements functions ###########################################

# Rescaling function that can be used with natural or mid p-values

resc_estim <- function(pvalues, pCDFlist, g_func, lambda = NULL, r = NULL, pvaluesupportlist = NULL) {
  
  nu_vec <- c()
  
  for (i in 1:length(pvalues)) {
    
    pvaluesproba <- diff(c(0, unique(pCDFlist[[i]])))
    if(is.null(pvaluesupportlist)) nu <- p_nullexpectation_gfunc(unique(pCDFlist[[i]]), pvaluesproba, g_func, lambda=lambda, r=r)
    else nu <- p_nullexpectation_gfunc(pvaluesupportlist[[i]], pvaluesproba, g_func, lambda=lambda, r=r)
    
    nu_vec <- c(nu_vec, nu)
  }
  # print(which(nu_vec == 0))
  
  if (length(formals(g_func)) == 3) return(sum(g_func(pvalues, lambda, r) * (1/nu_vec)) + 1/min(nu_vec))
  else return(sum(g_func(pvalues) * (1/nu_vec)) + 1/min(nu_vec))
  
}





# # Storey rescaled estimator 
# Storey_res <- function(pvalues, lambda, CDFs) {
#   stepf <- lapply(CDFs, function(x) stepfun(x, c(0, x)))
#   m0_hat = 0
#   for (i in 1:length(pvalues)) {
#     m0_hat = m0_hat + (1 / (1 - stepf[[i]](lambda))) * (pvalues[i] > lambda)
#   }
#   m0_hat = m0_hat + (1 / (1 - lambda))
#   
#   pi0_hat = m0_hat / length(pvalues)
#   return(pi0_hat)
# }
# 
# 
# # PC on mid pvalues
# PC_midp <- function(pvalues, CDF){
#   mid_pvalues <- midp.values(CDF, pvalues)
#   return((2 + 2 * sum(mid_pvalues)) / length(pvalues))
# }

# randomized estimator
  ## PC randomized 
PC_randomized <- function(pvalues, CDFs, n.sim=100){
  results <- replicate(n.sim, 1 / PC_new_DM(randomizedp.values(CDFs, pvalues)))
  return(1 / (mean(results)))
}

## Storey randomized
Storey_randomized <- function(pvalues, CDFs, lambda, n.sim=100){
  results <- replicate(n.sim, 1 / Storey(randomizedp.values(CDFs, pvalues), lambda)) 
  return(1 / (mean(results)))
}

## Polynomial randomized
Poly_randomized <- function(pvalues, CDFs, lambda, r, n.sim = 100){
  results <- replicate(n.sim, 1 / poly_PC_STOREY(randomizedp.values(CDFs, pvalues), lambda, r))
  return(1 / (mean(results)))
}
