## remote server paths 
setwd("/users/home/meah/phd_projects/pi0_discrete/code/scripts")
source("/users/home/meah/phd_projects/pi0_discrete/code/scripts/mid_p-randomized_pvalues.R")
source("/users/home/meah/phd_projects/pi0_discrete/code/scripts/estimators.R")

## local mac paths 
# setwd("/Users/iqraa/Documents/phd_projects/discrete_pi0/code/scripts/")
# source("/Users/iqraa/Documents/phd_projects/discrete_pi0/code/scripts/mid_p-randomized_pvalues.R")
# source("/Users/iqraa/Documents/phd_projects/discrete_pi0/code/scripts/estimators.R")
library(OnlineSuperUnif)
library(lubridate)
library(svMisc)
library(rjson)

#------------------ load parameters from json file -----------------------------
param.list <- fromJSON(file = "../config_files/xp_compare_DM_vs_ZZD.json") 
#-------------------------------------------------------------------------------

nb_methods = 2 ##

# parameter of interest
nulls_prop <- seq(param.list$nulls_proportion$begin, param.list$nulls_proportion$end,
                      by = param.list$nulls_proportion$by)  ##


estimates <- matrix(data = NA, nrow = nb_methods * param.list$nb_run,  ##
                    ncol = length(nulls_prop), byrow = FALSE)  ##
# print(estimates)

#------------------ Run xp -----------------------------------------------------
for (i in nulls_prop) {  
  vec <- c()
  
  for (j in 1:param.list$nb_run) {
    
    pvalues <- generate_pvalues(param.list$m, param.list$mu1, i)$raw
    
    estimates_ <- c( PC_new_DM(pvalues) / param.list$m, 
                     PC_ZZD(pvalues, param.list$C, param.list$s) / param.list$m)

    vec <- c(vec, estimates_)
  }
  
  estimates[, i] <- vec 
}

#-------------------- save data frame  ----------------------------
method = rep(c("PC_DM", "PC_ZZD"), length.out = nb_methods * param.list$nb_run) 

df_estimates <- data.frame(cbind(method, estimates))

file_name = gsub(" " , "", paste("../../xp_data/", gsub(" ", "_", paste("compare_DM_vs_ZZD", param.list$param_interest, "study", now(), sep="_")), ".csv"))

write.csv(df_estimates, file_name, row.names = FALSE)


