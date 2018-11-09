library(nichefillr)
library(dplyr)
library(ggplot2)
library(drake)

#' First some quick tests that the nichefillr package is working as expected

data(example_parms)

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

