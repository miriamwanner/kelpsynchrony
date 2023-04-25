library(ggplot2)

setwd("/Users/miriam/Documents/Github/kelpsynchrony/Code/") # needs to be changed


col1 <- c(("region"), rep("all", 8), 
          rep("island", 8), rep("main", 8))
col2 <- c(("dispersal metric"), rep(c(rep("transport", 4), 
                                      rep("connectivity", 4)), 3))
col3 <- c(("scale of connectivity"), rep(c(rep("linear", 2), 
                                           rep("log", 2)), 6))
col4 <- c(("scale of synchrony"), rep(c("logit", "not logit"), 12))

all_0.1 <- readRDS(file = "all_sr_0.1.Rds")
all_0.02 <- readRDS(file = "all_sr_0.02.Rds")
all_0.5 <- readRDS(file = "all_sr_0.5.Rds")

islands_0.1 <- readRDS(file = "islands_sr_0.1.Rds")
islands_0.02 <- readRDS(file = "islands_sr_0.02.Rds")
islands_0.5 <- readRDS(file = "islands_sr_0.5.Rds")

main_0.1 <- readRDS(file = "main_sr_0.1.Rds")
main_0.02 <- readRDS(file = "main_sr_0.02.Rds")
main_0.5 <- readRDS(file = "main_sr_0.5.Rds")

col5_results_0.1_dispersal <- c("dispersal", 
                                all_0.1[1,2], all_0.1[1,3],
                                all_0.1[1,4], all_0.1[1,5],
                                all_0.1[1,6], all_0.1[1,7],
                                all_0.1[1,8], all_0.1[1,9],
                                islands_0.1[1,2], islands_0.1[1,3],
                                islands_0.1[1,4], islands_0.1[1,5],
                                islands_0.1[1,6], islands_0.1[1,7],
                                islands_0.1[1,8], islands_0.1[1,9],
                                main_0.1[1,2], main_0.1[1,3],
                                main_0.1[1,4], main_0.1[1,5],
                                main_0.1[1,6], main_0.1[1,7],
                                main_0.1[1,8], main_0.1[1,9])
col6_results_0.1_distance <- c("distance",
                               all_0.1[2,2], all_0.1[2,3],
                               all_0.1[2,4], all_0.1[2,5],
                               all_0.1[2,6], all_0.1[2,7],
                               all_0.1[2,8], all_0.1[2,9],
                               islands_0.1[2,2], islands_0.1[2,3],
                               islands_0.1[2,4], islands_0.1[2,5],
                               islands_0.1[2,6], islands_0.1[2,7],
                               islands_0.1[2,8], islands_0.1[2,9],
                               main_0.1[2,2], main_0.1[2,3],
                               main_0.1[2,4], main_0.1[2,5],
                               main_0.1[2,6], main_0.1[2,7],
                               main_0.1[2,8], main_0.1[2,9])
col7_results_0.1_waves <- c("waves",
                            all_0.1[3,2], all_0.1[3,3],
                            all_0.1[3,4], all_0.1[3,5],
                            all_0.1[3,6], all_0.1[3,7],
                            all_0.1[3,8], all_0.1[3,9],
                            islands_0.1[3,2], islands_0.1[3,3],
                            islands_0.1[3,4], islands_0.1[3,5],
                            islands_0.1[3,6], islands_0.1[3,7],
                            islands_0.1[3,8], islands_0.1[3,9],
                            main_0.1[3,2], main_0.1[3,3],
                            main_0.1[3,4], main_0.1[3,5],
                            main_0.1[3,6], main_0.1[3,7],
                            main_0.1[3,8], main_0.1[3,9])
col8_results_0.1_nitrate <- c("nitrate",
                              all_0.1[4,2], all_0.1[4,3],
                              all_0.1[4,4], all_0.1[4,5],
                              all_0.1[4,6], all_0.1[4,7],
                              all_0.1[4,8], all_0.1[4,9],
                              islands_0.1[4,2], islands_0.1[4,3],
                              islands_0.1[4,4], islands_0.1[4,5],
                              islands_0.1[4,6], islands_0.1[4,7],
                              islands_0.1[4,8], islands_0.1[4,9],
                              main_0.1[4,2], main_0.1[4,3],
                              main_0.1[4,4], main_0.1[4,5],
                              main_0.1[4,6], main_0.1[4,7],
                              main_0.1[4,8], main_0.1[4,9])
results_0.1 <- data.frame(col1, col2, col3, col4, 
                          col5_results_0.1_dispersal,
                          col6_results_0.1_distance,
                          col7_results_0.1_waves,
                          col8_results_0.1_nitrate)

col5_results_0.02_dispersal <- c("dispersal", 
                                 all_0.02[1,2], all_0.02[1,3],
                                 all_0.02[1,4], all_0.02[1,5],
                                 all_0.02[1,6], all_0.02[1,7],
                                 all_0.02[1,8], all_0.02[1,9],
                                 islands_0.02[1,2], islands_0.02[1,3],
                                 islands_0.02[1,4], islands_0.02[1,5],
                                 islands_0.02[1,6], islands_0.02[1,7],
                                 islands_0.02[1,8], islands_0.02[1,9],
                                 main_0.02[1,2], main_0.02[1,3],
                                 main_0.02[1,4], main_0.02[1,5],
                                 main_0.02[1,6], main_0.02[1,7],
                                 main_0.02[1,8], main_0.02[1,9])
col6_results_0.02_distance <- c("distance",
                                all_0.02[2,2], all_0.02[2,3],
                                all_0.02[2,4], all_0.02[2,5],
                                all_0.02[2,6], all_0.02[2,7],
                                all_0.02[2,8], all_0.02[2,9],
                                islands_0.02[2,2], islands_0.02[2,3],
                                islands_0.02[2,4], islands_0.02[2,5],
                                islands_0.02[2,6], islands_0.02[2,7],
                                islands_0.02[2,8], islands_0.02[2,9],
                                main_0.02[2,2], main_0.02[2,3],
                                main_0.02[2,4], main_0.02[2,5],
                                main_0.02[2,6], main_0.02[2,7],
                                main_0.02[2,8], main_0.02[2,9])
col7_results_0.02_waves <- c("waves",
                             all_0.02[3,2], all_0.02[3,3],
                             all_0.02[3,4], all_0.02[3,5],
                             all_0.02[3,6], all_0.02[3,7],
                             all_0.02[3,8], all_0.02[3,9],
                             islands_0.02[3,2], islands_0.02[3,3],
                             islands_0.02[3,4], islands_0.02[3,5],
                             islands_0.02[3,6], islands_0.02[3,7],
                             islands_0.02[3,8], islands_0.02[3,9],
                             main_0.02[3,2], main_0.02[3,3],
                             main_0.02[3,4], main_0.02[3,5],
                             main_0.02[3,6], main_0.02[3,7],
                             main_0.02[3,8], main_0.02[3,9])
col8_results_0.02_nitrate <- c("nitrate",
                               all_0.02[4,2], all_0.02[4,3],
                               all_0.02[4,4], all_0.02[4,5],
                               all_0.02[4,6], all_0.02[4,7],
                               all_0.02[4,8], all_0.02[4,9],
                               islands_0.02[4,2], islands_0.02[4,3],
                               islands_0.02[4,4], islands_0.02[4,5],
                               islands_0.02[4,6], islands_0.02[4,7],
                               islands_0.02[4,8], islands_0.02[4,9],
                               main_0.02[4,2], main_0.02[4,3],
                               main_0.02[4,4], main_0.02[4,5],
                               main_0.02[4,6], main_0.02[4,7],
                               main_0.02[4,8], main_0.02[4,9])
results_0.02 <- data.frame(col1, col2, col3, col4, 
                          col5_results_0.02_dispersal,
                          col6_results_0.02_distance,
                          col7_results_0.02_waves,
                          col8_results_0.02_nitrate)

col5_results_0.5_dispersal <- c("dispersal", 
                                all_0.5[1,2], all_0.5[1,3],
                                all_0.5[1,4], all_0.5[1,5],
                                all_0.5[1,6], all_0.5[1,7],
                                all_0.5[1,8], all_0.5[1,9],
                                islands_0.5[1,2], islands_0.5[1,3],
                                islands_0.5[1,4], islands_0.5[1,5],
                                islands_0.5[1,6], islands_0.5[1,7],
                                islands_0.5[1,8], islands_0.5[1,9],
                                main_0.5[1,2], main_0.5[1,3],
                                main_0.5[1,4], main_0.5[1,5],
                                main_0.5[1,6], main_0.5[1,7],
                                main_0.5[1,8], main_0.5[1,9])
col6_results_0.5_distance <- c("distance",
                               all_0.5[2,2], all_0.5[2,3],
                               all_0.5[2,4], all_0.5[2,5],
                               all_0.5[2,6], all_0.5[2,7],
                               all_0.5[2,8], all_0.5[2,9],
                               islands_0.5[2,2], islands_0.5[2,3],
                               islands_0.5[2,4], islands_0.5[2,5],
                               islands_0.5[2,6], islands_0.5[2,7],
                               islands_0.5[2,8], islands_0.5[2,9],
                               main_0.5[2,2], main_0.5[2,3],
                               main_0.5[2,4], main_0.5[2,5],
                               main_0.5[2,6], main_0.5[2,7],
                               main_0.5[2,8], main_0.5[2,9])
col7_results_0.5_waves <- c("waves",
                            all_0.5[3,2], all_0.5[3,3],
                            all_0.5[3,4], all_0.5[3,5],
                            all_0.5[3,6], all_0.5[3,7],
                            all_0.5[3,8], all_0.5[3,9],
                            islands_0.5[3,2], islands_0.5[3,3],
                            islands_0.5[3,4], islands_0.5[3,5],
                            islands_0.5[3,6], islands_0.5[3,7],
                            islands_0.5[3,8], islands_0.5[3,9],
                            main_0.5[3,2], main_0.5[3,3],
                            main_0.5[3,4], main_0.5[3,5],
                            main_0.5[3,6], main_0.5[3,7],
                            main_0.5[3,8], main_0.5[3,9])
col8_results_0.5_nitrate <- c("nitrate",
                              all_0.5[4,2], all_0.5[4,3],
                              all_0.5[4,4], all_0.5[4,5],
                              all_0.5[4,6], all_0.5[4,7],
                              all_0.5[4,8], all_0.5[4,9],
                              islands_0.5[4,2], islands_0.5[4,3],
                              islands_0.5[4,4], islands_0.5[4,5],
                              islands_0.5[4,6], islands_0.5[4,7],
                              islands_0.5[4,8], islands_0.5[4,9],
                              main_0.5[4,2], main_0.5[4,3],
                              main_0.5[4,4], main_0.5[4,5],
                              main_0.5[4,6], main_0.5[4,7],
                              main_0.5[4,8], main_0.5[4,9])
results_0.5 <- data.frame(col1, col2, col3, col4, 
                          col5_results_0.5_dispersal,
                          col6_results_0.5_distance,
                          col7_results_0.5_waves,
                          col8_results_0.5_nitrate)






# plotting

ggplot(results_0.1, aes(x=col1[2:25])) +
  geom_line(aes(y=col5_results_0.1_dispersal[2:25]), color="red")




  