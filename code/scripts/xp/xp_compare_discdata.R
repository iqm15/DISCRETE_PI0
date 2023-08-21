## remote server paths 
setwd("/users/home/meah/phd_projects/pi0_discrete/code/scripts")
source("/users/home/meah/phd_projects/pi0_discrete/code/scripts/mid_p-randomized_pvalues.R")
source("/users/home/meah/phd_projects/pi0_discrete/code/scripts/estimators.R")
source("/users/home/meah/phd_projects/pi0_discrete/code/scripts/utils.R")

## local mac paths 
# setwd("/Users/iqraa/Documents/phd_projects/discrete_pi0/code/scripts/")
# source("/Users/iqraa/Documents/phd_projects/discrete_pi0/code/scripts/mid_p-randomized_pvalues.R")
# source("/Users/iqraa/Documents/phd_projects/discrete_pi0/code/scripts/estimators.R")
# source("/Users/iqraa/Documents/phd_projects/discrete_pi0/code/scripts/utils.R")

library(OnlineSuperUnif)
library(lubridate)
library(svMisc)
library(rjson)

#------------------ load parameters from json file -----------------------------
param.list <- fromJSON(file = "../config_files/xp_compare_disc.json") 
#-------------------------------------------------------------------------------

nb_methods = 12 ##

# parameter of interest
non_nulls_prop <- seq(param.list$non_nulls_proportion$begin, param.list$non_nulls_proportion$end,
                      by = param.list$non_nulls_proportion$by)  ##

signal_strength <- seq(param.list$p3$begin, param.list$p3$end, by = param.list$p3$by)

estimates <- matrix(data = NA, nrow = nb_methods * param.list$nb_run,  ##
                               ncol = length(non_nulls_prop), byrow = FALSE)  ##

#------------------ Run xp -----------------------------------------------------
for (p3 in signal_strength) {  
  
  
  for (i in 1:length(non_nulls_prop)){
    vec <- c()
    for (j in 1:param.list$nb_run) {
      
        data <- data_simulation(param.list$N, 
                                param.list$m, 
                                non_nulls_prop[i], 
                                p3, 
                                param.list$cluster_option)
        
        p_value_data <- pvalues_simulation(data$data)
        p_values <- p_value_data$raw
        CDF <- p_value_data$support
        
        
        midp_values <- midp.values(CDF, p_values)
        
        estimates_ <- c(
                        Storey(p_values, param.list$lambda) / length(p_values), 
                        PC_new_DM(p_values) / length(p_values),
                        poly_PC_STOREY(p_values,  param.list$lambda, param.list$r) / length(p_values),
                        
                        resc_estim(p_values, CDF, g_poly, lambda=param.list$lambda, r=0) / length(p_values),
                        resc_estim(p_values, CDF, g_pc) / length(p_values),
                        resc_estim(p_values, CDF, g_poly, lambda=param.list$lambda, r=param.list$r) / length(p_values),
                        
                        resc_estim(midp_values$obs_midp, CDF, g_poly, lambda=param.list$lambda_midp, r=0, pvaluesupportlist=midp_values$midpsupport) / length(p_values),
                        resc_estim(midp_values$obs_midp, CDF, g_pc, pvaluesupportlist=midp_values$midpsupport) / length(p_values),
                        resc_estim(midp_values$obs_midp, CDF, g_poly, lambda=param.list$lambda_midp, r=param.list$r, pvaluesupportlist=midp_values$midpsupport) / length(p_values),
                        
                        Storey_randomized(p_values, CDF, lambda=param.list$lambda) / length(p_values),
                        PC_randomized(p_values, CDF) / length(p_values),
                        Poly_randomized(p_values, CDF, lambda=param.list$lambda, r=param.list$r) / length(p_values)
                        )
    
        
        vec <- c(vec, estimates_)
      }
  
      estimates[, i] <- vec 
    }
    
    #-------------------- save data frame  ----------------------------
    method = rep(c("STOREY", "PC", "POLY", 
                   "RESC_STOREY", "RESC_PC", "RESC_POLY", 
                   "RESC_STOREY_midp", "RESC_PC_midp", "RESC_POLY_midp",
                   "STOREY_RANDOM", "PC_RANDOM", "POLY_RANDOM"), param.list$nb_run)
    
    df_estimates <- data.frame(cbind(method, estimates))
    
    file_name = gsub(" " , "", paste("../../xp_data/FET/estimators/", gsub(" ", "_", paste("data_estimates_signal_strength", as.character(p3), now(), sep="_")), ".csv"))
    
    write.csv(df_estimates, file_name, row.names = FALSE)


}