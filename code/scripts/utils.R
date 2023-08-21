
generate_pvalues <- function(m, mu1, pi0, seed=NULL) {
  #-----------------------------------------------------------------------------
  if (!is.null(seed)){
    set.seed(seed)
  }
  #-----------------------------------------------------------------------------
  
  true_indices <- rbinom(m, 1, 1-pi0)  # c(rep(1, m * (1 - pi0)), rep(0, m * pi0))
  # mu_0 <- 0
  # max.mu.1 = 3
  # runif(m, 0, max.mu.1)
  observations <- rnorm(m) + mu1 * true_indices
  raw.pvalues <- pnorm(observations, lower.tail = FALSE)
  
  data <- list(true_indices = true_indices, observations = observations, raw = raw.pvalues)
  return(data)
}


BH_proc <- function(pvalues, alpha, adaptive = FALSE, hat_pi0 = NULL) {
  
  m = length(pvalues)
  if(!adaptive) hat_pi0 <- 1 
  
  if (length(which(sort(pvalues) <= ((1:m) * alpha)/ (m * hat_pi0))) > 0) {
    k_hat <- max(which(sort(pvalues) <= ((1:m) * alpha)/ (m * hat_pi0)))
    rejections_index <- order(pvalues)[1:k_hat]
  }
  else {
    k_hat <- 0
    rejections_index <- c()
  }
  
  
  output <- list(k_hat = k_hat, rej = rejections_index)
  return(output)
}


proc_power <- function(true_indices, rejection_set) {
  
  alternative_indices <- which(true_indices == 1)
  nb_true_discoveries <- sum(rejection_set %in% alternative_indices)
  power <- nb_true_discoveries / length(alternative_indices)
  
  return(power)
}


# bias variance and mse for the oracle case where we would know the 
# expectation and variance, of the generic g function, 
# under the null and the alternative 

bias <- function(m, m0, g_alt_exp, nu) {
  return((1/nu) * (1 + (m - m0) * g_alt_exp))
}

variance <- function(m, m0, g_null_var, g_alt_var, nu) {
  return((1/nu**2) * (m0 * g_null_var + (m - m0) * g_alt_var))
}

mse <- function(m, m0, g_alt_exp, g_null_var, g_alt_var, nu) {
  bias_ <- bias(m, m0, g_alt_exp, nu) / m
  var_ <- variance(m, m0, g_null_var, g_alt_var, nu) / m**2
  return(bias_**2 + var_)
}

g_storey <- function(u, lambda) {
  return((u > lambda)*1)
}

g_pc <- function(u) {
  return(u)
}

g_sto_pc <- function(u , lambda) {
  return((u > lambda) * u)
}

g_quad <- function(u) {
  return(u**2)
}

g_sto_quad <- function(u, lambda) {
  return((u > lambda) * u**2)
}

g_poly <- function(u, lambda, r) {
  return((u > lambda) * u**r)
}

g_new_poly <- function(u, lambda, r = 1) {
  return(as.numeric((( 1 - ((1 - u) / (1 - lambda))**2) * (u >lambda))))
  
}

g_alt_gaussian_density_2 <- function(x, mu) {
  return(dnorm(qnorm(1 - x) - mu) * (1 / dnorm(qnorm(1 - x))))
}

g_alt_gaussian_density <- function(x, mu) {
  return(exp(-mu * qnorm(x) - (mu**2) / 2))
}

g_unif_density <- function(x) {
  if((x <= 1) & (x >=0)) return(1)
  else return(0)
  
}

integrand_null_exp <- function(x, g_func, lambda = NULL, r = NULL) {
  if (length(formals(g_func)) == 2) return(g_unif_density(x) * g_func(x, lambda))
  else if (length(formals(g_func)) == 3) return(g_unif_density(x) * g_func(x, lambda, r))
  else return(g_unif_density(x) * g_func(x))
}


integrand_null_var <- function(x, g_func, lambda = NULL, r = NULL) {
  expectation <- integrate(integrand_null_exp, 0, 1, g_func, lambda, r)$value
  if (length(formals(g_func)) == 2) return(g_unif_density(x) * (g_func(x, lambda) - expectation)**2)
  else if (length(formals(g_func)) == 3) return(g_unif_density(x) * (g_func(x, lambda, r) - expectation)**2)
  else return(g_unif_density(x) * (g_func(x) - expectation)**2)
}

integrand_exp <- function(x, mu, g_func, lambda = NULL, r = NULL) {
  if (length(formals(g_func)) == 2) return(g_alt_gaussian_density(x, mu) * g_func(x, lambda))
  else if (length(formals(g_func)) == 3) return(g_alt_gaussian_density(x, mu) * g_func(x, lambda, r))
  else return(g_alt_gaussian_density(x, mu) * g_func(x))
}

integrand_var <- function(x, mu, g_func, lambda = NULL, r = NULL) {
  expectation <- integrate(integrand_exp, 0, 1, mu, g_func, lambda, r)$value
  if (length(formals(g_func)) == 2) return(g_alt_gaussian_density(x, mu) * (g_func(x, lambda) - expectation)**2)
  else if (length(formals(g_func)) == 3) return(g_alt_gaussian_density(x, mu) * (g_func(x, lambda, r) - expectation)**2)
  else return(g_alt_gaussian_density(x, mu) * (g_func(x) - expectation)**2)
}

p_nullexpectation_gfunc <- function(pvaluesupport, pvalueproba, g_func, lambda = NULL, r = NULL) {
  
  if (length(formals(g_func)) == 3) return(sum(g_func(pvaluesupport, lambda, r) * pvalueproba))
  else return(sum(g_func(pvaluesupport) * pvalueproba))
  
  
}





# g_alt_gaussian_exp <- function(g_func, mu, lambda) {
#         
#   return(integrate(integrand_exp, 0, 1, mu, g_func, lambda)$value)
# }
# 
# g_alt_gaussian_var <- function(g_func, mu, lambda) {
#   return(integrate(integrand_var, 0, 1, mu, g_func, lambda)$value)
# }
