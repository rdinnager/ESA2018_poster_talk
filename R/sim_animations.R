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