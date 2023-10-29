# getting the MRM results

setwd("/Users/miriam/Desktop/revised_kelp_code/Code") # needs to be changed
rm(list=ls())
# location for storing the results
resloc <- "../Results/"
if (!dir.exists(resloc)){
  dir.create(resloc, recursive=TRUE)
}
# loading the data
datloc <- "../Data/"

load(file=paste0(resloc, "new_data_processing.RData"))

library(ncdf4)
library(R.matlab)
library(geosphere)
library(plyr)
library(wsyn)
library(ncf)
library(hexbin)
library(ggplot2)
library(gtools)
library(ecodist)
library(PBSmapping)
library(rjson)
library(raptr)
library(ggsn)



diag(sm) <- NA
logit_sm <- logit(sm, min=-1, max=1)
diag(sm) <- 1

source("altered_mrm_function.R")

# transport, linear, logit
t_lin_logit <- my_mrm(as.dist(logit_sm) ~ as.dist(symProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# transport, linear, not logit
t_lin_nlogit <- my_mrm(as.dist(sm) ~ as.dist(symProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# transport, log, logit
t_log_logit <- my_mrm(as.dist(logit_sm) ~ as.dist(logSymProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# transport, log, not logit
t_log_nlogit <- my_mrm(as.dist(sm) ~ as.dist(logSymProbMat) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# connectivity, linear, logit
c_lin_logit <- my_mrm(as.dist(logit_sm) ~ as.dist(symProbMatF) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# connectivity, linear, not logit
c_lin_nlogit <- my_mrm(as.dist(sm) ~ as.dist(symProbMatF) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# connectivity, log, logit
c_log_logit <- my_mrm(as.dist(logit_sm) ~ as.dist(logSym) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))
# connectivity, log, not logit
c_log_nlogit <- my_mrm(as.dist(sm) ~ as.dist(logSym) + as.dist(newDistMat) + as.dist(smWaves) + as.dist(smNO3))

predictors <- c("dispersal", "distance", "waves", "nitrate")

all_sr_0.1 <- data.frame(predictors,
                          t_lin_logit$coef[7:10],
                          t_lin_nlogit$coef[7:10],
                          t_log_logit$coef[7:10],
                          t_log_nlogit$coef[7:10],
                          c_lin_logit$coef[7:10],
                          c_lin_nlogit$coef[7:10],
                          c_log_logit$coef[7:10],
                          c_log_nlogit$coef[7:10])
saveRDS(all_sr_0.1, file="all_sr_0.1.Rds")






# MRM for only islands (sites 50-117)

sm_islands <- sm[-(1:50), -(1:50)]
symProbMat_islands <- symProbMat[-(1:50), -(1:50)]
logSymProbMat_islands <- logSymProbMat[-(1:50), -(1:50)]
symProbMatF_islands <- symProbMatF[-(1:50), -(1:50)]
logSym_islands <- logSym[-(1:50), -(1:50)]
newDistMat_islands <- newDistMat[-(1:50), -(1:50)]
smWaves_islands <- smWaves[-(1:50), -(1:50)]
smNO3_islands <- smNO3[-(1:50), -(1:50)]


diag(sm_islands) <- NA
logit_sm_islands <- logit(sm_islands, min=-1, max=1)
diag(sm_islands) <- 1


# transport, linear, logit
t_lin_logit <- my_mrm(as.dist(logit_sm_islands) ~ as.dist(symProbMat_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# transport, linear, not logit
t_lin_nlogit <- my_mrm(as.dist(sm_islands) ~ as.dist(symProbMat_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# transport, log, logit
t_log_logit <- my_mrm(as.dist(logit_sm_islands) ~ as.dist(logSymProbMat_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# transport, log, not logit
t_log_nlogit <- my_mrm(as.dist(sm_islands) ~ as.dist(logSymProbMat_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# connectivity, linear, logit
c_lin_logit <- my_mrm(as.dist(logit_sm_islands) ~ as.dist(symProbMatF_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# connectivity, linear, not logit
c_lin_nlogit <- my_mrm(as.dist(sm_islands) ~ as.dist(symProbMatF_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# connectivity, log, logit
c_log_logit <- my_mrm(as.dist(logit_sm_islands) ~ as.dist(logSym_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))
# connectivity, log, not logit
c_log_nlogit <- my_mrm(as.dist(sm_islands) ~ as.dist(logSym_islands) + as.dist(newDistMat_islands) + as.dist(smWaves_islands) + as.dist(smNO3_islands))

predictors <- c("dispersal", "distance", "waves", "nitrate")

islands_sr_0.1 <- data.frame(predictors,
                             t_lin_logit$coef[7:10],
                             t_lin_nlogit$coef[7:10],
                             t_log_logit$coef[7:10],
                             t_log_nlogit$coef[7:10],
                             c_lin_logit$coef[7:10],
                             c_lin_nlogit$coef[7:10],
                             c_log_logit$coef[7:10],
                             c_log_nlogit$coef[7:10])
saveRDS(islands_sr_0.1, file="islands_sr_0.1.Rds")




# MRM for only mainland (sites 1-49)

sm_main <- sm[-(1:50), -(50:117)]
symProbMat_main <- symProbMat[-(50:117), -(50:117)]
logSymProbMat_main <- logSymProbMat[-(50:117), -(50:117)]
symProbMatF_main <- symProbMatF[-(50:117), -(50:117)]
logSym_main <- logSym[-(50:117), -(50:117)]
newDistMat_main <- newDistMat[-(50:117), -(50:117)]
smWaves_main <- smWaves[-(50:117), -(50:117)]
smNO3_main <- smNO3[-(50:117), -(50:117)]


diag(sm_main) <- NA
logit_sm_main <- logit(sm_main, min=-1, max=1)
diag(sm_main) <- 1

# transport, linear, logit
t_lin_logit <- my_mrm(as.dist(logit_sm_main) ~ as.dist(symProbMat_main) + as.dist(newDistMat_main) + as.dist(smWaves_main) + as.dist(smNO3_main))
# transport, linear, not logit
t_lin_nlogit <- my_mrm(as.dist(sm_main) ~ as.dist(symProbMat_main) + as.dist(newDistMat_main) + as.dist(smWaves_main) + as.dist(smNO3_main))
# transport, log, logit
t_log_logit <- my_mrm(as.dist(logit_sm_main) ~ as.dist(logSymProbMat_main) + as.dist(newDistMat_main) + as.dist(smWaves_main) + as.dist(smNO3_main))
# transport, log, not logit
t_log_nlogit <- my_mrm(as.dist(sm_main) ~ as.dist(logSymProbMat_main) + as.dist(newDistMat_main) + as.dist(smWaves_main) + as.dist(smNO3_main))
# connectivity, linear, logit
c_lin_logit <- my_mrm(as.dist(logit_sm_main) ~ as.dist(symProbMatF_main) + as.dist(newDistMat_main) + as.dist(smWaves_main) + as.dist(smNO3_main))
# connectivity, linear, not logit
c_lin_nlogit <- my_mrm(as.dist(sm_main) ~ as.dist(symProbMatF_main) + as.dist(newDistMat_main) + as.dist(smWaves_main) + as.dist(smNO3_main))
# connectivity, log, logit
c_log_logit <- my_mrm(as.dist(logit_sm_main) ~ as.dist(logSym_main) + as.dist(newDistMat_main) + as.dist(smWaves_main) + as.dist(smNO3_main))
# connectivity, log, not logit
c_log_nlogit <- my_mrm(as.dist(sm_main) ~ as.dist(logSym_main) + as.dist(newDistMat_main) + as.dist(smWaves_main) + as.dist(smNO3_main))

predictors <- c("dispersal", "distance", "waves", "nitrate")

main_sr_0.1 <- data.frame(predictors,
                             t_lin_logit$coef[7:10],
                             t_lin_nlogit$coef[7:10],
                             t_log_logit$coef[7:10],
                             t_log_nlogit$coef[7:10],
                             c_lin_logit$coef[7:10],
                             c_lin_nlogit$coef[7:10],
                             c_log_logit$coef[7:10],
                             c_log_nlogit$coef[7:10])
saveRDS(main_sr_0.1, file="main_sr_0.1.Rds")

