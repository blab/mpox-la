library(doParallel)
library(foreach)
library(tidyverse)
library(viridis)
#library(metR)
library(RColorBrewer)

source('./scripts/utils_inference.R')
source('./scripts/utils_cluster_alloc.R')
source('./scripts/utils_distrib.R')


## Input files 

cluster_alloc <-read_csv("../mpox_la/scripts/cluster_distribution_2023onwards.csv")
#cluster_alloc <-read_csv("../mpox_la/scripts/cluster_distribution_2024only.csv")
#cluster_alloc <-read_csv("../mpox_la/scripts/cluster_distribution_2023only.csv")

#cluster_alloc <- readRDS(file_cluster_alloc)

## Probability that transmission occurs before mutation
# p_trans_before_mut <- readRDS(file_proba_trans_before_mut) %>% 
#   ungroup() %>% 
#   filter(pathogen == 'MERS') %>% 
#   select(p_trans_before_mut) %>%
#   as.numeric()

p_trans_before_mut <- 1

## Scenario for the proportion of infection sequenced
n_sequences <- sum(cluster_alloc$cluster_size * cluster_alloc$count)
n_suspected_cases <- 328
prop_cases_sequenced <- n_sequences/n_suspected_cases
vec_p_detect_cases <- c(0.05, 0.1, 0.25, 0.5) # Proportion of infections detected (as cases)

## Setting up the boundary for the optimization
R_min <- 0.01
R_max <- 10.0
k_min <- 0.001
k_max <- 10.0

## Setting up the inference
max_cluster_size_inference <- 10000
n_cores <- 4

## Running the grid search
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)

t0 <- Sys.time()
df_inference <- Reduce('bind_rows',
                       foreach(i_p_detect = 1:length(vec_p_detect_cases),
                               .packages = c('tidyverse')) %dopar% {
                                 
                                 ## Define the function to minimize
                                 minus_log_lik <- function(param){
                                   R <- param[1]
                                   k <- param[2]
                                   
                                   ll <- get_loglik_obs(vec_cluster_size = cluster_alloc$cluster_size,
                                                        vec_n_clusters = cluster_alloc$count,
                                                        R = R, k = k,
                                                        p = p_trans_before_mut,
                                                        p_detect = vec_p_detect_cases[i_p_detect] * prop_cases_sequenced, 
                                                        max_cluster_size = max_cluster_size_inference)
                                   
                                   
                                   return(-ll)
                                 }
                                 
                                 mle_estim <- optim(par = c(1.0, 0.1),
                                                    fn = minus_log_lik,
                                                    method = 'L-BFGS-B',
                                                    lower = c(R_min, k_min), upper = c(R_max, k_max))
                                 
                                 CI_R <- get_CI(1, minus_log_lik, mle_estim,
                                                param_min = R_min, param_max = R_max,
                                                increment_param = 0.01)
                                 CI_k <- get_CI(2, minus_log_lik, mle_estim,
                                                param_min = k_min, param_max = k_max,
                                                increment_param = 0.001)
                                 
                                 bind_cols(tibble(mle_estim = mle_estim$par,
                                                  param = c('R', 'k'),
                                                  p_detect = rep(vec_p_detect_cases[i_p_detect], 2)),
                                           bind_rows(CI_R, CI_k))
                                 
                               })
t1 <- Sys.time()
print(t1 - t0)
stopCluster(cl)
write.csv(df_inference, "cluster_inference_2023onwards_lessdetect.csv", row.names = FALSE)
print(df_inference)


# 
# p <- 0.7
# vec_R <- seq(0.5, 2.0, 0.5)
# vec_k <- c(0.01, 0.1, 1.)
# max_cluster_size <- 10000
# 
# # Assuming your dataframe is named `df`
# df_transformed <- df_inference %>%
#   pivot_wider(
#     names_from = param,      # Use the 'param' column for column names
#     values_from = everything(),  # Spread all other columns
#     id_cols = p_detect       # Use 'p_detect' as the index
#   )
# 
# 
# ### Generating the distribution
# df_theoretical_distrib <- Reduce('bind_rows', lapply(vec_R, FUN = function(R){
#   Reduce('bind_rows', lapply(vec_k, FUN = function(k){
#     theoretical_distrib_clusters(R, k, p, max_cluster_size = max_cluster_size) %>% 
#       mutate(R = R, k = k, p = p)
#   })) 
# }))
# 
# 
# # Create a list of tuples (R, k)
# vec_Rk <- expand.grid(R = , k = vec_k) %>% split(1:nrow(.))
# 
# # Apply function using a single lapply() over tuples
# df_theoretical_distrib <- Reduce('bind_rows', lapply(vec_Rk, FUN = function(Rk){
#   theoretical_distrib_clusters(Rk$R, Rk$k, p, max_cluster_size = max_cluster_size) %>%
#     mutate(R = Rk$R, k = Rk$k, p = p)
# }))
# 
# 
# 
# ### Plotting the distribution
# plt_theoretical_distrib_cluster_size <- df_theoretical_distrib %>%
#   mutate(k = paste0('k = ', k)) %>%
#   ggplot(aes(x = cluster_size,
#              y = prob,
#              group = interaction(k, R),
#              colour = as.factor(R))) +
#   geom_step() +
#   facet_wrap(.~ k, nrow = 1) +
#   scale_x_continuous(name = 'Size of cluster of identical sequences',
#                      breaks = c(1, seq(10, 50, 10)),
#                      expand = expansion(mult = c(0., 0.05))) +
#   scale_y_continuous(name = 'Cumulative distribution\nfunction') +
#   scale_colour_manual(name = 'R',
#                       breaks = seq(0.5, 2.0, 0.5),
#                       labels = c('0.5', '1.0', '1.5', '2.0'),
#                       values = RColorBrewer::brewer.pal(5, 'BuPu')[-1]) +
#   theme_bw() +
#   coord_cartesian(xlim = c(1, 50), ylim = c(0.0, 1.0)) +
#   theme(panel.grid.minor.x = element_blank(),
#         strip.background = element_rect(fill = 'gray22'),
#         strip.text = element_text(colour = 'white'),
#         panel.spacing = unit(1, 'lines'))
# 
# plot(plt_theoretical_distrib_cluster_size)
