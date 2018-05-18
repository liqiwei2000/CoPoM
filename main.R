# Data: X with n observations and p variables
# observation_supervised == FALSE: unsupervised; 
#  observation_supervised == TRUE: supervised with known cluster allocation.
# variable_supervised == FALSE: unsupervised; 
#  variable_supervised == TRUE: supervised with known variable identification.
# s_constraint == FALSE: gamma prior with parameters a_s = b_s = 1; 
#  s_constraint == TRUE: Dirichlet process prior;  
# g_constraint == FALSE: Gamma prior with parameters a_g = b_g = 1;
#  g_constraint == TRUE: Dirichlet process prior; 

# Load libraries and functions ------------------------------------------------------------------
source("Data_Files/functions.R");
source("Data_Files/mcmc.R");
library(mcclust);
library(lattice);

# Input data ------------------------------------------------------------------------------------
load("Data_Files/simulated_data.RData");
# load("Data_Files/real_data.RData");
n <- dim(X)[1];
p <- dim(X)[2];

# Load hyperparameters --------------------------------------------------------------------------
a <- 1;
b <- 1;
alpha <- 0.1;
omega <- 0.01;
d_0 <- 1;
M <- floor(n/2);
L <- floor(p/2);
p_0 <- 0.5;
hyperparameters <- c(a, b, alpha, omega, d_0, M, L, p_0);

# Load environment variables -------------------------------------------------------------------
observation_supervised <- FALSE;
variable_supervised <- FALSE;
s_constraint <- TRUE;
g_constraint <- TRUE;
alpha_supervised <- FALSE;
omega_supervised <- FALSE;
p_0_supervised <- FALSE;
ZIP <- TRUE;
environment <- c(observation_supervised, variable_supervised, s_constraint, g_constraint, alpha_supervised, omega_supervised, ZIP, p_0_supervised);

# Load algorithm settings -----------------------------------------------------------------------
iter <- 10000;
burn <- 1000;
start_K <- 2;
start_gamma_sum <- 2;

# Store the samples -----------------------------------------------------------------------------
z_store <- matrix(0L, iter, n);
K_store <- rep(0L, iter);
alpha_store <- rep(0L, iter);
ari_store <- rep(0L, iter);
gamma_sum_store <- rep(NA, iter);
gamma_ppi_store <- rep(0L, p);
omega_store <- rep(0L, iter);
s_store <- matrix(0, iter, n);
g_store <- matrix(0, iter, p);
map_z_store <- rep(0, iter);
map_gamma_store <- rep(0, iter);
d_mean <- matrix(0, n, p);

# Generate starting parameters ------------------------------------------------------------------
parameters <- starting(start_K, start_gamma_sum, hyperparameters, environment, X);

# Start MCMC algorithms -------------------------------------------------------------------------
wall_time <- 0;
system_time <- 0;
count <- 0;
for (i in 1:iter) {
  parameters <- core(parameters, hyperparameters, environment, X);
  
  # Monitor the process
  if (count == floor(i/iter*100)) {
    print(paste0(count, "% has been done"));
    count <- count + 1;
  }
  
  # Output the parameters into log
  wall_time <- wall_time + parameters$time[1];
  system_time <- system_time + parameters$time[2];
  z_store[i,] <- parameters$z;
  K_store[i] <- length(unique(parameters$z));
  gamma_sum_store[i] <- sum(parameters$gamma);
  s_store[i,] <- parameters$s;
  g_store[i,] <- parameters$g;
  alpha_store[i] <- parameters$alpha;
  omega_store[i] <- parameters$omega;
  map_z_store[i] <- parameters$map_z;
  map_gamma_store[i] <- parameters$map_gamma;
  
  # Obtain the MAP estimate
  if (i == burn) {
    z_map <- parameters$z;
    gamma_map <- parameters$gamma;
    d_map <- parameters$d;
    map_z_temp <- parameters$map_z;
    map_gamma_temp <- parameters$map_gamma;
  } else if (i > burn) {
    gamma_ppi_store <- gamma_ppi_store + parameters$gamma;
    if (parameters$map_z > map_z_temp) {
      z_map <- parameters$z;
      map_z_temp <- parameters$map_z;
    }
    if (parameters$map_gamma > map_gamma_temp) {
      gamma_map <- parameters$gamma;
      map_gamma_temp <- parameters$map_gamma;
    }
  }
}
print(paste0("Running time = ", round(wall_time/60, 3), " minutes"));

# Obtain the MP estimate
mpv <- gamma_ppi_store/(iter - burn);
gamma_ppi <- (mpv >= 0.5);
ppm <- matrix(0, nrow = n, ncol = n);
for (ii in burn:iter) {
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
       if (z_store[ii, i] == z_store[ii, j]) {
         ppm[i, j] <- ppm[i, j] + 1;
         ppm[j, i] <- ppm[j, i] + 1;
       }
    }
  }
}
ppm <- ppm/(iter - burn + 1);
diag(ppm) <- rep(1, n);
z_ppm <- minbinder(ppm, method = "comp")$cl;
