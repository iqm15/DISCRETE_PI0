## remote server paths 
setwd("/users/home/meah/phd_projects/pi0_discrete/code/scripts")
source("mid_p-randomized_pvalues.R")
source("estimators.R")
source("utils.R")

## local mac paths 
# setwd("/Users/iqraa/Documents/phd_projects/discrete_pi0/code/scripts/")
# source("mid_p-randomized_pvalues.R")
# source("estimators.R")
# source("utils.R")

library(lubridate)
library(svMisc)
library(rjson)

#------------------ load parameters from json file -----------------------------
param.list <- fromJSON(file = "../config_files/xp_compare_cont.json") 
#-------------------------------------------------------------------------------

nb_methods = 6 ##

alpha <- param.list$alpha

# parameter of interest
nulls_prop <- seq(param.list$nulls_proportion$begin, param.list$nulls_proportion$end,
                  by = param.list$nulls_proportion$by)  ##

mu1_range <- seq(param.list$mu1$begin, param.list$mu$end, by = param.list$mu1$by)

#------------------ Run xp -----------------------------------------------------
for (mu1 in mu1_range){
  
  estimates <- matrix(data = NA, nrow = nb_methods * param.list$nb_run,  ##
                      ncol = length(nulls_prop), byrow = FALSE)  ##
  
  rejections <- matrix(data = NA, nrow = (nb_methods + 2) * param.list$nb_run,  ##
                       ncol = length(nulls_prop), byrow = FALSE)  ##
  
  power <- matrix(data = NA, nrow = (nb_methods + 2) * param.list$nb_run,  ##
                       ncol = length(nulls_prop), byrow = FALSE)  ##  
  
  for (i in 1:length(nulls_prop)) {  
    
    
    vec_1 <- c()
    vec_2 <- c()
    vec_3 <- c()
    
    for (j in 1:param.list$nb_run) {
      
      
      # generate data
      data <- generate_pvalues(param.list$m, mu1, nulls_prop[i])
      pvalues <- data$raw
      
      # compute estimators for the different methods 
      Storey_estim <- Storey(pvalues, param.list$lambda) / param.list$m
      PC_estim <- PC_new_DM(pvalues) / param.list$m
      PC_STOREY_estim <- PC_Storey(pvalues, param.list$lambda) / param.list$m
      QUAD_estim <- quad_PC(pvalues) / param.list$m
      PC_ZZD_estim <- PC_ZZD(pvalues, param.list$C, param.list$s) / param.list$m
      QUAD_STOREY_estim <- Quad_Storey(pvalues, param.list$lambda) / param.list$m
      
      estimates_ <- c(Storey_estim, PC_estim, PC_STOREY_estim, QUAD_estim, PC_ZZD_estim, QUAD_STOREY_estim)
      
      # Run adaptive BH procedure with each estimator and also raw BH and oracle adaptive BH
      ABH_STOREY <- BH_proc(pvalues, alpha, TRUE, Storey_estim)
      ABH_PC <- BH_proc(pvalues, alpha, TRUE, PC_estim)
      ABH_PC_STOREY <- BH_proc(pvalues, alpha, TRUE, PC_STOREY_estim)
      ABH_QUAD <- BH_proc(pvalues, alpha, TRUE, QUAD_estim)
      ABH_PC_ZZD <- BH_proc(pvalues, alpha, TRUE, PC_ZZD_estim)
      ABH_QUAD_STOREY <- BH_proc(pvalues, alpha, TRUE, QUAD_STOREY_estim)
      BH <- BH_proc(pvalues, alpha)
      oracle_BH <- BH_proc(pvalues, alpha, TRUE, nulls_prop[i])
   
      rejections_ <- c(ABH_STOREY$k_hat, ABH_PC$k_hat, ABH_PC_STOREY$k_hat, ABH_QUAD$k_hat, ABH_PC_ZZD$k_hat, ABH_QUAD_STOREY$k_hat, BH$k_hat, oracle_BH$k_hat)
      
      power_ <- c(
        proc_power(data$true_indices, ABH_STOREY$rej),
        proc_power(data$true_indices, ABH_PC$rej),
        proc_power(data$true_indices, ABH_PC_STOREY$rej),
        proc_power(data$true_indices, ABH_QUAD$rej),
        proc_power(data$true_indices, ABH_PC_ZZD$rej),
        proc_power(data$true_indices, ABH_QUAD_STOREY$rej),
        proc_power(data$true_indices, BH$rej),
        proc_power(data$true_indices, oracle_BH$rej)
      )
      
      vec_1 <- c(vec_1, estimates_)
      vec_2 <- c(vec_2, rejections_)
      vec_3 <- c(vec_3, power_)
    }
    
    estimates[, i] <- vec_1
    rejections[, i] <- vec_2
    power[, i] <- vec_3
 }
    
  #-------------------- save data frame for estimators  ----------------------------
  method = rep(c("STOREY", "PC", "PC_STOREY", "QUAD", "PC_ZZD", "QUAD_STOREY"), param.list$nb_run) ##

  df_estimates <- data.frame(cbind(method, estimates))

  file_name = gsub(" " , "", paste("../../xp_data/gaussian/estimators/", gsub(" ", "_", paste("data_estimates_mu1", as.character(mu1), now(), sep="_")), ".csv"))

  write.csv(df_estimates, file_name, row.names = FALSE)
  
  #-------------------- save data frame for rejection nb ABH procedures  ----------------------------  

  method = rep(c("STOREY", "PC", "PC_STOREY", "QUAD", "PC_ZZD", "QUAD_STOREY", "BH", "oracle_ABH"), param.list$nb_run) ##

  df_estimates <- data.frame(cbind(method, rejections))

  file_name = gsub(" " , "", paste("../../xp_data/gaussian/rejections/", gsub(" ", "_", paste("data_rejections_mu1", as.character(mu1), now(), sep="_")), ".csv"))

  write.csv(df_estimates, file_name, row.names = FALSE)
  
  #-------------------- save data frame for power of ABH procedure  ---------------------------- 
  
  method = rep(c("STOREY", "PC", "PC_STOREY", "QUAD", "PC_ZZD", "QUAD_STOREY", "BH", "oracle_ABH"), param.list$nb_run) ##
  
  df_estimates <- data.frame(cbind(method, power))
  
  file_name = gsub(" " , "", paste("../../xp_data/gaussian/power/", gsub(" ", "_", paste("data_power_mu1", as.character(mu1), now(), sep="_")), ".csv"))
  
  write.csv(df_estimates, file_name, row.names = FALSE)
  
}