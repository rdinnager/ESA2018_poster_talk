library(diffeqr)
library(nichefillr)
library(dplyr)
library(ggplot2)
library(drake)

diffeq_setup()

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

### Now let's try a random peaks landscape

test_params <- generate_landscape2(1000, 5, num_peaks = 50,
                                  h_to_sig_ratio = 1,
                                  P_min_max = c(0.75, 1.5),
                                  a = 0.01,
                                  dirichlet_param = 10)

nichefillr:::plot_K_contour(test_params, x_lims = c(-20, 20), y_lims = c(-20, 20),
                            contour_res = 200, bins = 12)

vis_K_landscape(test_params, test_params$niche_rad)

test_params <- generate_landscape2(1000, 5, num_peaks = 50,
                                   h_to_sig_ratio = 3,
                                   P_min_max = c(0.75, 1.5),
                                   a = 0.01,
                                   dirichlet_param = 10,
                                   d = 3)

nichefillr:::plot_K_contour3D(test_params, x_lims = c(-10, 10), y_lims = c(-10, 10),
                              z_lims = c(-6, 6),
                            contour_res = 200, bins = 12)


full_params <- list(K_parms = c(test_params, c_var = 0.5),
                    a_parms = list(gamma_i = c(1, 1),
                                   D_i = c(1, 1),
                                   C = 1),
                    macro_parms = list(b_rate = 0.001,
                                       init_traits = c(0.1, 0.1),
                                       e_var = c(0.1, 0.1),
                                       init_Ns = c(0.1, 0.1),
                                       init_br = 10,
                                       check_extinct = 0.005,
                                       tot_time = 100000,
                                       V_gi = c(0.01, 0.01),
                                       m = 2L, d = 2L,
                                       save_tree = TRUE,
                                       save_tree_interval = 100,
                                       progress = TRUE,
                                       mult = 2))

test_sim <- sim_radiation(full_params)
plot(test_sim)
sim_animation(test_sim, file_name = "test_sim_big", view = FALSE)

######### generate params ##############

n_runs <- 4000

gen_params <- function(d = 2) {
  #comp_vol <- rexp(1, 1/40)
  
  comp_1d <- rexp(1, 1/1.25)
  gamma <- comp_1d * (d ^ (-0.5))
  test_params <- generate_landscape2(1000, 10, num_peaks = 50,
                                     h_to_sig_ratio = 2,
                                     P_min_max = c(0.75, 1.5),
                                     a = 0.01,
                                     dirichlet_param = 10,
                                     d = d)
  
  full_params <- list(K_parms = c(test_params, c_var = 0.5),
                      a_parms = list(gamma_i = rep(gamma, d),
                                     D_i = rep(1, d),
                                     C = 1),
                      macro_parms = list(b_rate = 0.001,
                                         init_traits = rep(0.01, d),
                                         e_var = rep(0.1, d),
                                         init_Ns = rep(0.1, 2),
                                         init_br = 10,
                                         check_extinct = 0.005,
                                         tot_time = 200000,
                                         V_gi = rep(0.01, d),
                                         m = 2L, d = as.integer(d),
                                         save_tree = TRUE,
                                         save_tree_interval = 100,
                                         progress = FALSE,
                                         mult = 2))
  
  full_params
}

library(pbapply)
library(parallel)

ds <- sample(4, n_runs, replace = TRUE) + 1

#test <- gen_params(ds[1])

cl <- makeCluster(8L)
clusterEvalQ(cl, {library(nichefillr)})

param_sets <- pblapply(ds, gen_params, cl = cl)

library(readr)
write_rds(param_sets, "results/params4/param_sets.rds")

stopCluster(cl)

######## run models ##########
library(diffeqr)
diffeq_setup()
library(nichefillr)
library(dplyr)
library(ggplot2)
library(readr)
param_sets <- read_rds("results/params3/param_sets.rds")

## fix error from above

#param_sets <- lapply(param_sets, function(x) {x$macro_parms$init_Ns <- c(0.1, 0.1); x})

# sims_test <- sim_radiation_multi(param_sets[1:9], save_folder = "results/sims",
#                                  ncpus = 9, save_prefix = "sim_set_1_")

sims_test <- sim_radiation(param_sets[[1]])

sims_null <- sim_radiation_multi(param_sets, save_folder = "results/sims3",
                                 ncpus = 9, save_prefix = "sim_set_3_")

######## Analyse models #########
library(nichefillr)
library(dplyr)
library(ggplot2)
library(readr)
library(ape)
library(FD)

mod_files <- list.files("results/sims3", full.names = TRUE)

# mod_file <- mod_files[2]
extract_div_data <- function(mod_file) {
  sim <- read_rds(mod_file)
  if(class(sim) != "try-error") {
    SR <- sum(sim$sim_object$extant)
    extant_list <- sim$sim_object$extant_list$as.list()
    times <- sim$sim_object$tree_times$as.list()
    SR_trend <- sapply(extant_list, length)
    SR_trend <- data_frame(Time = times, SR = SR_trend)
    traits <- t(sim$sim_object$traits)[sim$sim_object$extant, , drop = FALSE]
    rownames(traits) <- sim$sim_object$phylo$tip.label[sim$sim_object$extant]
    if(SR > 1) {
      capt <- capture.output(FD <- dbFD(traits, stand.x = FALSE))
      FD <- as.data.frame(FD[1:8])
    } else {
      FD <- NA
    }
    params <- c(d = sim$sim_params$macro_parms$d, gamma = sim$sim_params$a_parms$gamma_i[1])
    #print(mod_file)
    rm(sim)
    gc()
    return(list(SR = SR, SR_trend = SR_trend, FD = FD, params = params))
  } else {
    return(NULL)
  }
}

library(pbapply)

sim_dat <- pblapply(mod_files, extract_div_data)
write_rds(sim_dat, "results/sim_dat3.rds")


sim_dat <- read_rds("results/sim_dat2.rds")
sim_dat <- sim_dat[sapply(sim_dat, function(x) !is.null(x))]


dat <- sim_dat[[1]]
make_df <- function(dat) {
  if(is.na(dat$FD)) {
    dat$FD <- data_frame(nbsp = NA, sing.sp = NA,
                         FRic = NA, qual.FRic = NA,
                         FEve = NA, FDiv = NA,
                         FDis = NA, RaoQ = NA)
  }
  dat$FD %>%
    mutate(SR = dat$SR, gamma = dat$params["gamma"], d = dat$params["d"])
}

dat_df <- lapply(sim_dat, make_df) %>%
  bind_rows

calc_vol <- function(sigg, d) {
  ((2*pi)^(d/2))*(sigg^d)
}

dat_df <- dat_df %>%
  rowwise %>%
  mutate(comp_vol = calc_vol(gamma, d))

ggplot(dat_df, aes(comp_vol, SR)) +
  geom_point(aes(colour = as.factor(d))) +
  geom_smooth(aes(colour = as.factor(d), fill = as.factor(d))) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 50, 100, 200, 400, 600)) +
  scale_color_brewer(name = "Dimensionality", palette = "Dark2") +
  scale_fill_brewer(name = "Dimensionality", palette = "Dark2") +
  xlab(expression(paste("Competition Strength (", gamma, ")"))) +
  ylab("Equilibrium Species Richness") +
  theme_minimal() +
  theme(legend.position = c(0.75, 0.75),
        panel.grid = element_blank())

ggplot(dat_df, aes(gamma, SR)) +
  geom_point(aes(colour = as.factor(d))) +
  geom_smooth(aes(colour = as.factor(d), fill = as.factor(d))) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 0.35, 0.5, 1, 2, 4, 6)) +
  scale_color_brewer(name = "Dimensionality", palette = "Dark2") +
  scale_fill_brewer(name = "Dimensionality", palette = "Dark2") +
  xlab(expression(paste("Competition Strength (", gamma, ")"))) +
  ylab("Equilibrium Species Richness") +
  theme_minimal() +
  theme(legend.position = c(0.75, 0.75),
        panel.grid = element_blank())

dat_df2 <- dat_df %>%
  dplyr::filter(d == 2)


library(scico)
library(scales)
ggplot(dat_df2, aes(SR, RaoQ)) +
  geom_density2d(colour = "grey30", bins = 20) +
  geom_point(aes(colour = gamma), alpha = 0.8) +
  scale_color_scico(name = "Competition\nStrength", palette = "tofino", trans = "log1p", breaks = c(0.25, 1.0, 2.5, 5.0, 7.5)) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  ylab("Functional Diversity (Rao's Quadratic Entropy)") +
  xlab("Species Richness")

sim_dat <- read_rds("results/sim_dat.rds")

grab_summ <- function(dat) {
  if(!is.null(dat)) {
    data_frame(SR = dat$SR, gamma = dat$params["gamma"], d = dat$params["d"])
  } else {
    data_frame(SR = NA, gamma = NA, d = NA)
  }
}


sim_summ <- lapply(sim_dat, grab_summ) %>%
  bind_rows

low_comp <- sim_summ[which((sim_summ$SR < 10) & (sim_summ$gamma < 0.25) & (sim_summ$d == 2)), ]

low_comp_num <- which((sim_summ$SR < 10) & (sim_summ$gamma < 0.25) & (sim_summ$d == 2))[4]
low_comp_params <- read_rds(mod_files[low_comp_num])$sim_params
low_comp_sim <- read_rds(mod_files[low_comp_num])

library(nichefillr)
library(readr)
library(ape)
library(ggtree)

plot(low_comp_sim, type = "trace", expand_factor = 8)
plot(low_comp_sim$sim_object$phylo)

ggtree(low_comp_sim$sim_object$phylo)
low_comp_alive <- drop.tip(low_comp_sim$sim_object$phylo, which(!low_comp_sim$sim_object$extant))
ggtree(low_comp_alive)

med_comp <- sim_summ[which((sim_summ$SR > 57) & (sim_summ$SR < 63) & (sim_summ$gamma > 0.95) & (sim_summ$gamma < 1.1) & (sim_summ$d == 2)), ]

med_comp_num <- which((sim_summ$SR > 57) & (sim_summ$SR < 63) & (sim_summ$gamma > 0.95) & (sim_summ$gamma < 1.1) & (sim_summ$d == 2))[1]
med_comp_params <- read_rds(mod_files[med_comp_num])$sim_params
med_comp_sim <- read_rds(mod_files[med_comp_num])

plot(med_comp_sim, type = "trace", expand_factor = 1)

ggtree(med_comp_sim$sim_object$phylo)
med_comp_alive <- drop.tip(med_comp_sim$sim_object$phylo, which(!med_comp_sim$sim_object$extant))
ggtree(med_comp_alive)


high_comp <- sim_summ[which((sim_summ$SR < 17) & (sim_summ$gamma > 2.5) & (sim_summ$d == 2)), ]

high_comp_num <- which((sim_summ$SR < 17) & (sim_summ$gamma > 2.5) & (sim_summ$d == 2))[29]
high_comp_params <- read_rds(mod_files[high_comp_num])$sim_params
high_comp_sim <- read_rds(mod_files[high_comp_num])

plot(high_comp_sim, type = "trace", expand_factor = 0.8)

ggtree(high_comp_sim$sim_object$phylo)
high_comp_alive <- drop.tip(high_comp_sim$sim_object$phylo, which(!high_comp_sim$sim_object$extant))
ggtree(high_comp_alive)


ode_init <- c(as.vector(high_comp_sim$sim_object$traits[ , high_comp_sim$sim_object$extant][,1:6]), high_comp_sim$sim_object$Ns[high_comp_sim$sim_object$extant][1:6])

# get_fit <- function(x) {
#   nichefillr:::temp_fitness_interface(x = x, 
#                                       d = as.integer(high_comp_sim$sim_params$macro_parms$d), 
#                                       m = as.integer(length(ode_init)/3), 
#                                       u = as.integer(length(high_comp_sim$sim_params$K_parms$hz)), 
#                                       a = high_comp_sim$sim_params$K_parms$a,
#                                       h0 = high_comp_sim$sim_params$K_parms$h0, 
#                                       h_z = high_comp_sim$sim_params$K_parms$hz,
#                                       P0 = high_comp_sim$sim_params$K_parms$P0i[1],
#                                       sigma0_i = high_comp_sim$sim_params$K_parms$sig0i,
#                                       P_z = high_comp_sim$sim_params$K_parms$Piz[1, , drop = TRUE],
#                                       D0 = high_comp_sim$sim_params$a_parms$D_i[1],
#                                       b_iz = high_comp_sim$sim_params$K_parms$biz,
#                                       state = ode_init, 
#                                       V_gi = high_comp_sim$sim_params$macro_parms$V_gi,
#                                       sigma_iz = high_comp_sim$sim_params$K_parms$sigiz, 
#                                       gamma_i = high_comp_sim$sim_params$a_parms$gamma_i,
#                                       C = 1)
#   
# }

# get_fit(c(0, 0))

library(tidyr)
contour_df <- crossing(Niche_Axis_1 = seq(-30, 30, length.out = 100),
                       Niche_Axis_2 = seq(-30, 30, length.out = 100))

spec_traits <- matrix(ode_init[1:12], nrow = 2)
spec_Ns <- ode_init[13:18] / 1

K_parms <- high_comp_sim$sim_params$K_parms
a_parms <- high_comp_sim$sim_params$a_parms
x <- c(5, -3)
xmat <- apply(spec_traits, 2, function(y) (((y - x)^2) / (2 * a_parms$gamma_i))^a_parms$D_i) 
comp <- exp(-apply(xmat, 2, sum)) * spec_Ns
sum(comp)


z <- apply(contour_df %>% as.matrix, 1, 
           function(y) nichefillr:::fit_func(y, high_comp_sim$sim_params$K_parms, 
                                            high_comp_sim$sim_params$a_parms, 
                                            spec_traits, spec_Ns, res_x = c(0, 0), d = 2))

contour_df <- contour_df %>%
  mutate(K = z) %>%
  transform(K = ifelse(K < -1, NA, K))

library(scico)
pp <- ggplot(contour_df, aes(Niche_Axis_1, Niche_Axis_2)) +
  geom_raster(aes(fill = K)) +
  geom_contour(aes(z = K), data = contour_df, colour = "grey20", bins = 10) +
  scale_fill_scico(palette = "bilbao") +
  coord_equal() +
  theme_minimal()
pp

library(diffeqr)
diffeq_setup()
library(nichefillr)
low_comp_sim <- sim_radiation(low_comp_params, trait_hist_prop = 0.25)

test1 <- read_rds("results/sims/sim_set_1_1008.rds")


######## make plot of species richness traces ########
library(pbapply)
library(readr)
library(dplyr)
sim_dat <- read_rds("results/sim_dat.rds")
sim_dat <- sim_dat[sapply(sim_dat, function(x) !is.null(x))]

#num <- 100
get_SR_trend <- function(num) {
  x <- sim_dat[[num]]
  if(!is.null(x)) {
    df <- x$SR_trend %>% 
      mutate(sim_num = num) %>%
      mutate(d = x$params["d"], gamma = x$params["gamma"]) %>%
      transform(Time = unlist(Time))
  } else {
    df <- data_frame(Time = NULL, SR = NULL, sim_num = NULL, d = NULL, gamma = NULL)
  }
  df
}

SR_trends <- pblapply(seq_along(sim_dat), get_SR_trend) %>%
  bind_rows %>%
  mutate(gamma_interval = cut_interval(log(gamma + 1), 9))

ggplot(SR_trends %>% filter(d == 2), aes(Time, SR)) +
  geom_line(aes(colour = gamma, group = sim_num), alpha = 0.05) +
  facet_wrap(~gamma_interval, nrow = 3, scales = "free") +
  theme_minimal()

ggplot(SR_trends %>% filter(d == 2) %>% group_by(sim_num) %>% mutate(SR = SR / max(SR)), aes(Time, SR)) +
  geom_line(aes(group = sim_num), alpha = 0.05) +
  ylab("Maximum Species Richness") +
  xlab("Time Step") +
  theme_minimal() +
  theme(panel.grid = element_blank())

sample_sim_num <- SR_trends %>%
  filter(d == 2) %>%
  group_by(sim_num) %>%
  summarise(gamma_interval = gamma_interval[1], gamma = gamma[1]) %>%
  ungroup() %>%
  group_by(gamma_interval) %>%
  sample_frac(0.075)

library(scico)
ggplot(SR_trends %>% filter(sim_num %in% sample_sim_num$sim_num), 
       aes(Time, SR)) +
  geom_line(aes(colour = gamma, group = sim_num), alpha = 0.5) +
  scale_color_scico(name = "Competition\nStrength", palette = "tofino", trans = "log1p", breaks = c(0.25, 1.0, 2.5, 5.0)) +
  ylab("Species Richness") +
  xlab("Time Step") +
  theme_minimal() +
  theme(panel.grid = element_blank())

####### run RPANDAS diversification models on simulated phylogenies ##############
library(RPANDA)

library(nichefillr)
library(dplyr)
library(ggplot2)
library(readr)
library(ape)
library(FD)

mod_files <- list.files("results/sims", full.names = TRUE)

# mod_file <- mod_files[2]
run_div_models <- function(mod_file) {
  sim <- read_rds(mod_file)
  if(class(sim) != "try-error") {
    if(sum(sim$sim_object$extant) > 2) {
      phyl <- drop.tip(sim$sim_object$phylo, which(!sim$sim_object$extant))
      phyl$edge.length <- phyl$edge.length / max(phyl$edge.length)
      model_1 <- try(fit_coal_cst(phyl, cst.rate = TRUE))
      model_2 <- try(fit_coal_cst(phyl, cst.rate = FALSE))
      model_5 <- try(fit_coal_var(phyl, cst.lamb = TRUE, mu.0 = TRUE))
      model_6 <- try(fit_coal_var(phyl, cst.lamb = FALSE, mu.0 = TRUE))
      return(list(model_1 = model_1, model_2 = model_2, model_5 = model_5, model_6 = model_6,
                  d = sim$sim_params$macro_parms$d, gamma = sim$sim_params$a_parms$gamma_i[1]))
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

library(pbapply)

sim_mod_res <- pblapply(mod_files, run_div_models)
write_rds(sim_mod_res, "results/sim_mod_res.rds")

sim_mod_res <- read_rds("results/sim_mod_res.rds")

#mod_num <- 1
get_data_from_mods <- function(mod_num) {
  mod <- sim_mod_res[[mod_num]]
  if(!is.null(mod)) {
    
    if(class(mod$model_1) != "try-error") {
      df1 <- data_frame(model = "model_1", aicc = mod$model_1$aicc, tau0 = mod$model_1$tau0)
    } else {
      df1 <- data_frame(model = "model_1", aicc = NA, tau0 = NA)
    }
    if(class(mod$model_2) != "try-error") {
      df2 <- data_frame(model = "model_2", aicc = mod$model_2$aicc, tau0 = mod$model_2$tau0, gamma = mod$model_2$gamma)
    } else {
      df2 <- data_frame(model = "model_2", aicc = NA, tau0 = NA, gamma = NA)
    }
    if(class(mod$model_5) != "try-error") {
      df5 <- data_frame(model = "model_5", aicc = mod$model_5$aicc, lamb0 = mod$model_5$lamb0)
    } else {
      df5 <- data_frame(model = "model_5", aicc = NA, lamb0 = NA)
    }
    if(class(mod$model_6) != "try-error") {
      df6 <- data_frame(model = "model_6", aicc = mod$model_6$aicc, lamb0 = mod$model_6$lamb0, alpha = mod$model_6$alpha)
    } else {
      df6 <- data_frame(model = "model_6", aicc = NA, lamb0 = NA, alpha = NA)
    }
    return(bind_rows(df1, df2, df5, df6) %>%
             mutate(sim_num = mod_num) %>%
             mutate(d = mod$d, gamma2 = mod$gamma))
  } else {
    data_frame(model = NA, aicc = NA)
  }
}

model_df <- pblapply(seq_along(sim_mod_res), get_data_from_mods) %>%
  bind_rows

ggplot(model_df %>% filter(model == "model_1"), aes(gamma2, tau0)) +
  geom_point()

ggplot(model_df %>% filter(model == "model_2"), aes(gamma2, tau0)) +
  geom_point()

ggplot(model_df %>% filter(model == "model_2"), aes(gamma2, gamma)) +
  geom_point() + ylim(c(0, 40))

ggplot(model_df %>% filter(model == "model_5"), aes(gamma2, lamb0)) +
  geom_point()

ggplot(model_df %>% filter(model == "model_6"), aes(gamma2, lamb0)) +
  geom_point()

ggplot(model_df %>% filter(model == "model_6"), aes(gamma2, alpha)) +
  geom_point() + ylim(c(-10, 40))

best_mods <- model_df %>%
  filter(is.finite(aicc)) %>%
  group_by(sim_num) %>%
  summarise(best_mod = model[which.min(aicc)], gamma2 = gamma2[1]) %>%
  mutate(model_1 = ifelse(best_mod == "model_1", 1, 0),
         model_2 = ifelse(best_mod == "model_2", 1, 0),
         model_5 = ifelse(best_mod == "model_5", 1, 0),
         model_6 = ifelse(best_mod == "model_6", 1, 0))

ggplot(best_mods, aes(gamma2, model_1)) +
  geom_point() +
  geom_smooth()

ggplot(best_mods, aes(gamma2, model_2)) +
  geom_point() +
  geom_smooth()

ggplot(best_mods, aes(gamma2, model_5)) +
  geom_point() +
  geom_smooth()

ggplot(best_mods, aes(gamma2, model_6)) +
  geom_point() +
  geom_smooth()