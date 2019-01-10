library(diffeqr)
diffeq_setup()
library(nichefillr)
library(readr)

test_params <- generate_landscape2(10000, 10, num_peaks = 50,
                                  h_to_sig_ratio = 2,
                                  P_min_max = c(0.75, 1.5),
                                  a = 0.01,
                                  dirichlet_param = 10)

nichefillr:::plot_K_contour(test_params, x_lims = c(-60, 60), y_lims = c(-60, 60),
                            contour_res = 200, bins = 12)

test_params2 <- generate_landscape2(10000, 10, num_peaks = 50,
                                   h_to_sig_ratio = 2,
                                   P_min_max = c(0.75, 1.5),
                                   a = 0.01,
                                   dirichlet_param = 10,
                                   d = 3)

nichefillr:::plot_K_contour3D(test_params, x_lims = c(-55, 55), y_lims = c(-55, 55),
                              z_lims = c(-50, 50), 
                            contour_res = 200, bins = 12)



full_params_low <- list(K_parms = c(test_params, c_var = 0),
                    a_parms = list(gamma_i = c(0.25, 0.25),
                                   D_i = c(1, 1),
                                   C = 1),
                    macro_parms = list(b_rate = 0.0005,
                                       init_traits = c(0.01, 0.01),
                                       e_var = c(0.1, 0.1),
                                       init_Ns = c(0.05, 0.05),
                                       init_br = 10,
                                       check_extinct = 0.01,
                                       tot_time = 50000,
                                       V_gi = c(0.01, 0.01),
                                       m = 2L, d = 2L,
                                       save_tree = TRUE,
                                       save_tree_interval = 100,
                                       progress = TRUE,
                                       mult = 2))

test_sim_low <- sim_radiation(full_params_low, trait_hist_prop = 0.25)
write_rds(test_sim_low, "anim/test_sim_low.rds")
plot(test_sim_low, type = "trace")


full_params_med1 <- list(K_parms = c(test_params, c_var = 0),
                    a_parms = list(gamma_i = c(1.5, 1.5),
                                   D_i = c(1, 1),
                                   C = 1),
                    macro_parms = list(b_rate = 0.001,
                                       init_traits = c(0.01, 0.01),
                                       e_var = c(0.1, 0.1),
                                       init_Ns = c(0.05, 0.05),
                                       init_br = 10,
                                       check_extinct = 0.01,
                                       tot_time = 100000,
                                       V_gi = c(0.01, 0.01),
                                       m = 2L, d = 2L,
                                       save_tree = TRUE,
                                       save_tree_interval = 100,
                                       progress = TRUE,
                                       mult = 2))

test_sim_med1 <- sim_radiation(full_params_med1, trait_hist_prop = 0.25)
write_rds(test_sim_med1, "anim/test_sim_med1.rds")
plot(test_sim_med1, type = "trace")

full_params_med2 <- list(K_parms = c(test_params, c_var = 0),
                         a_parms = list(gamma_i = c(2.5, 2.5),
                                        D_i = c(1, 1),
                                        C = 1),
                         macro_parms = list(b_rate = 0.001,
                                            init_traits = c(0.01, 0.01),
                                            e_var = c(0.1, 0.1),
                                            init_Ns = c(0.05, 0.05),
                                            init_br = 10,
                                            check_extinct = 0.01,
                                            tot_time = 100000,
                                            V_gi = c(0.01, 0.01),
                                            m = 2L, d = 2L,
                                            save_tree = TRUE,
                                            save_tree_interval = 100,
                                            progress = TRUE,
                                            mult = 2))

test_sim_med2 <- sim_radiation(full_params_med2, trait_hist_prop = 0.25)
write_rds(test_sim_med2, "anim/test_sim_med2.rds")
plot(test_sim_med2, type = "trace")

full_params_high <- list(K_parms = c(test_params, c_var = 0),
                        a_parms = list(gamma_i = c(5, 5),
                                       D_i = c(1, 1),
                                       C = 1),
                        macro_parms = list(b_rate = 0.001,
                                           init_traits = c(0.01, 0.01),
                                           e_var = c(0.1, 0.1),
                                           init_Ns = c(0.05, 0.05),
                                           init_br = 10,
                                           check_extinct = 0.01,
                                           tot_time = 100000,
                                           V_gi = c(0.01, 0.01),
                                           m = 2L, d = 2L,
                                           save_tree = TRUE,
                                           save_tree_interval = 100,
                                           progress = TRUE,
                                           mult = 2))

test_sim_high <- sim_radiation(full_params_high, trait_hist_prop = 0.25)
write_rds(test_sim_high, "anim/test_sim_high.rds")
plot(test_sim_high, type = "trace")

sim_animation(test_sim_low, file_name = "anim/test_sim_low.gif",
              contour_res = 200, x_lims = c(-60, 60), y_lims = c(-60, 60),
              num_frames = 1000)

sim_animation(test_sim_med1, file_name = "anim/test_sim_med1.gif",
              contour_res = 200, x_lims = c(-60, 60), y_lims = c(-60, 60),
              num_frames = 1000)

sim_animation(test_sim_med2, file_name = "anim/test_sim_med2.gif",
              contour_res = 200, x_lims = c(-60, 60), y_lims = c(-60, 60),
              num_frames = 1000)

sim_animation(test_sim_high, file_name = "anim/test_sim_high.gif",
              contour_res = 200, x_lims = c(-60, 60), y_lims = c(-60, 60),
              num_frames = 1000)



hub_and_spoke_K <- list(h0 = 0.5,
                        hz = 10*c(1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 0.1),
                        biz = 20*cbind(c(0, 0), c(0, 1), c(1, 0), c(0, -1), c(-1, 0),
                                    c(0, 0.55), c(0.55, 0), c(0, -0.55), c(-0.55, 0)),
                        sigiz = 20*cbind(c(0.15, 0.15), replicate(4, c(0.1, 0.1)), 
                                      c(0.04, 0.2), c(0.2, 0.04), c(0.02, 0.4), c(0.4, 0.02)),
                        Piz = rbind(c(1.1, 1.1, 1.1, 1.1, 1.1, 4, 4, 4, 4),
                                    c(1, 1, 1, 1, 1, 4, 4, 4, 4)),
                        sig0i = c(50, 50),
                        P0i = c(6, 6),
                        a = 0.002)

nichefillr:::plot_K_contour(hub_and_spoke_K, c(-30, 30), c(-30, 30), bins = 12)
full_params <- list(K_parms = c(hub_and_spoke_K, c_var = 0.02),
                    a_parms = list(gamma_i = c(1.5, 1.5),
                                   D_i = c(1, 1),
                                   C = 1),
                    macro_parms = list(b_rate = 0.005,
                                       init_traits = c(0.01, 0.01),
                                       e_var = c(0.1, 0.1),
                                       init_Ns = c(0.05, 0.05),
                                       init_br = 10,
                                       check_extinct = 0.005,
                                       tot_time = 100000,
                                       V_gi = c(0.01, 0.01),
                                       m = 2L, d = 2L,
                                       save_tree = TRUE,
                                       save_tree_interval = 100,
                                       progress = TRUE,
                                       mult = 2))

sim_hub_and_spoke_1 <- sim_radiation(full_params, trait_hist_prop = 0.25)
plot(sim_hub_and_spoke_1, type = "trace")

write_rds(sim_hub_and_spoke_1, "anim/sim_hub_and_spoke_1.rds")

sim_hub_and_spoke_1 <- read_rds("anim/sim_hub_and_spoke_1.rds")

sim_animation(sim_hub_and_spoke_1, file_name = "anim/hub_and_spoke_example2.gif",
              contour_res = 200, x_lims = c(-30, 30), y_lims = c(-30, 30),
              num_frames = 1000, include_phylogeny = TRUE, expand_factor = 0.05,
              max_size = 24, width = 8, height = 10)

########## set up and run "thank you" simulation ##########

long <- 0.6
short <- 0.1
steep <- c(1.5, 1.5)
thanks <- list(h0 = 1,
                        hz = 20*c(rep(1, 23), 0.5),
                        biz = 20*cbind(c(-2.5, 1), c(-2.5, 0), c(-2, 0), c(-1.5, 0), c(-1, -0.5),
                                       c(-0.25, 0), c(0.5, 0), c(1, -0.25),
                                       c(1.5, 0), c(2, -0.5), c(2.5, 0), c(3, 0),
                                       c(3.25, 0.25), c(3.25, -0.25),
                                       c(3.3, 0.5), c(3.3, -0.5),
                                       c(3.8, 0.25), c(4.8, -0.25),
                                       c(4.3, 0), c(4.3, 0.5), c(4.3, -0.5),
                                       c(5.3, 0.25), c(5.3, -0.75),
                                       c(0.75, -0.75)) - c(20, 0),
                        sigiz = 20*cbind(c(long, short), c(short, long), c(short, long), c(long/2, short), 
                                         c(short, long/2),
                                         c(long*3/4, long*3/4), c(short, long),
                                         c(short, long*3/4), c(long/2, short),
                                         c(short, long/2), c(short, long),
                                         c(short*6/4, short),
                                         c(short, short), c(short, short),
                                         c(short, short*3/2), c(short, short*3/2),
                                         c(short*3/2, short*3/2), c(short*3/2, short*3/2),
                                         c(long*2/4, short), c(long*3/4, short), c(long*3/4, short),
                                         c(short, long*3/4), c(short*3/2, short*3/2),
                                         c(short, short)),
                        Piz = cbind(steep, steep, steep, steep, steep,
                                    c(1.2, 1.2), steep, steep, steep, steep, c(1.2, 1.2),
                                    c(1, 1), c(1, 1),
                                    steep, steep, c(1.5, 1.5), c(1.5, 1.5),
                                    steep, steep,
                                    steep, steep, steep, 
                                    c(1.5, 1.5),
                                    c(1, 1)),
                        sig0i = c(80, 20),
                        P0i = c(6, 6),
                        a = 5)

nichefillr:::plot_K_contour(thanks, c(-200, 300), c(-60, 60), bins = 12, contour_res = 200)

full_params_thanks <- list(K_parms = c(thanks, c_var = 0),
                    a_parms = list(gamma_i = c(4, 4),
                                   D_i = c(1, 1),
                                   C = 1),
                    macro_parms = list(b_rate = 0.01,
                                       init_traits = c(10, 0),
                                       e_var = c(0.1, 0.1),
                                       init_Ns = c(0.05, 0.05),
                                       init_br = 5,
                                       check_extinct = 0.005,
                                       tot_time = 50000,
                                       V_gi = c(0.01, 0.01),
                                       m = 2L, d = 2L,
                                       save_tree = TRUE,
                                       save_tree_interval = 100,
                                       progress = TRUE,
                                       mult = 2))

thanks_sim <- sim_radiation(full_params_thanks, trait_hist_prop = 0.25)

plot(thanks_sim, type = "trace")

test <- sim_animation(thanks_sim, file_name = "anim/thanks_animation.gif",
              contour_res = 200, x_lims = c(-100, 100), y_lims = c(-40, 40),
              num_frames = 1000, include_phylogeny = FALSE, expand_factor = 0.05,
              max_size = 6, width = 10, height = 4, test_frame = 500, end_col = 0.95)

rstudioapi::viewer(test)

test <- sim_animation(thanks_sim, file_name = "anim/thanks_animation.gif",
                      contour_res = 200, x_lims = c(-100, 100), y_lims = c(-40, 40),
                      num_frames = 1000, include_phylogeny = FALSE, expand_factor = 0.05,
                      max_size = 6, width = 10, height = 4, end_col = 0.95)
