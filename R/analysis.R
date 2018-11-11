library(nichefillr)
library(dplyr)
library(ggplot2)
library(drake)

#' First some quick tests that the nichefillr package is working as expected

data(example_parms)

nichefillr:::plot_K_contour(example_parms$K_parms, c(-5, 5), c(-5, 5))

example_parms$macro_parms$tot_time <- 100000
test_sim <- sim_radiation(example_parms)
plot(test_sim)

test_sim$sim_object$traits[ , test_sim$sim_object$extant] %>%
  t %>%
  as.data.frame() %>%
  ggplot(aes(V1, V2)) +
  geom_hex()

sim_animation(test_sim, file_name = "test_sim", view = FALSE)

#' Now to test the new compiled fitness function

nichefillr:::temp_fitness_interface(x = c(5, 0), d = example_parms$macro_parms$d, 
                                    m = 2, 
                                    u = length(example_parms$K_parms$hz), 
                                    a = example_parms$K_parms$a,
                                    h0 = example_parms$K_parms$h0,
                                    h_z = example_parms$K_parms$hz,
                                    P0 = example_parms$K_parms$P0i[1],
                                    sigma0_i = example_parms$K_parms$sig0i,
                                    P_z = example_parms$K_parms$Piz[1, , drop = TRUE],
                                    D0 = example_parms$a_parms$D_i[1],
                                    b_iz = example_parms$K_parms$biz,
                                    state = c(10, 10, 10, 10, 0, 0), V_gi = example_parms$macro_parms$V_gi,
                                    sigma_iz = example_parms$K_parms$sigiz, 
                                    gamma_i = example_parms$a_parms$gamma_i,
                                    C = example_parms$a_parms$C)

#' First we create our "spoke and hub" fitness landscape in two dimensions and
#' test if it worked as expected.

hub_and_spoke_K <- list(h0 = 0.5,
                        hz = 1*c(1, 1, 1, 1, 1, 0.75, 0.75, 0.75, 0.75),
                        biz = cbind(c(0, 0), c(0, 1), c(1, 0), c(0, -1), c(-1, 0),
                                    c(0, 0.55), c(0.55, 0), c(0, -0.55), c(-0.55, 0)),
                        sigiz = cbind(c(0.15, 0.15), replicate(4, c(0.1, 0.1)), 
                                      c(0.02, 0.2), c(0.2, 0.02), c(0.02, 0.2), c(0.2, 0.02)),
                        Piz = rbind(c(1.1, 1.1, 1.1, 1.1, 1.1, 4, 4, 4, 4),
                                    c(1, 1, 1, 1, 1, 4, 4, 4, 4)),
                        sig0i = c(5, 5),
                        P0i = c(6, 6),
                        a = 0.15)

nichefillr:::plot_K_contour(hub_and_spoke_K, c(-2, 2), c(-2, 2), bins = 10)
full_params <- list(K_parms = c(hub_and_spoke_K, c_var = 0.02),
                    a_parms = list(gamma_i = c(0.1, 0.1),
                                   D_i = c(1, 1),
                                   C = 1),
                    macro_parms = list(b_rate = 0.005,
                                       init_traits = c(0.1, 0.1),
                                       e_var = c(0.1, 0.1),
                                       init_Ns = c(0.001, 0.001),
                                       init_br = 2,
                                       check_extinct = 0.005,
                                       tot_time = 100000,
                                       V_gi = c(0.02, 0.02),
                                       m = 2L, d = 2L,
                                       save_tree = TRUE,
                                       save_tree_interval = 100,
                                       progress = TRUE,
                                       mult = 1))

sim_hub_and_spoke_1 <- sim_radiation(full_params)
plot(sim_hub_and_spoke_1)
