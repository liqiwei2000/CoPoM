# **********************************************************************************************
# **********************************************************************************************
# Name: Constrained ZIP Mixtures (CoZIPm)
# Version: 3.0x
# Created: 04/04/2014
# Last update: 03/29/2016
# Author: Qiwei Li, Ph.D. student, Department of Statistics, Rice University
# **********************************************************************************************
# **********************************************************************************************

starting = function(start_K, start_gamma_sum, hyperparameters, environment, X) {
  # Load data, hyperparameters, and environment information ------------------------------------
  n <- dim(X)[1];
  p <- dim(X)[2];
  
  a <- hyperparameters[1];
  b <- hyperparameters[2];
  alpha <- hyperparameters[3];
  omega <- hyperparameters[4];
  d_0 <- hyperparameters[5];
  M <- hyperparameters[6];
  L <- hyperparameters[7];
  p_0 <- hyperparameters[8];
  
  observation_supervised <- environment[1];
  variable_supervised <- environment[2];
  s_constraint <- environment[3];
  g_constraint <- environment[4];
  alpha_supervised <- environment[5];
  omega_supervised <- environment[6];
  ZIP <- environment[7];
  p_0_supervised <- environment[8];
  
  # Initialize Gamma ---------------------------------------------------------------------------
  gamma <- rep(0, p);
  if (variable_supervised) {
    gamma <- gammat;
  } else {
    gamma[sample.int(p, start_gamma_sum, replace = FALSE)] <- 1;
    # gamma[1:start_gamma_sum] <- 1;
  }
  if (!omega_supervised) {
    a_omega <- 0.2;
    b_omega <- 1.8;
    omega <- rbeta(1, a_omega + sum(gamma == 1), b_omega + sum(gamma == 0));
  }
  
  # Initialize S -------------------------------------------------------------------------------
  if (s_constraint) {
    tau_eta <- 1;
    sigma_s <- 0.1;
    tau_s <- 1;
    c_s <- 0;
    t <- runif(M, 0, 1);
    eta <- rnorm(M, 0, tau_eta);
    phi <- rep(1/M, M);
    xis <- sample.int(M, n, replace = TRUE, prob = phi);
    psis <- sample.int(2, n, replace = TRUE) - 1;
    s <- exp(rnorm(n, psis*eta[xis] + (1 - psis)*(c_s - t[xis]*eta[xis])/(1 - t[xis]), sigma_s));
  } else {
    a_s <- 1;
    b_s <- 1;
    s <- rgamma(n, shape = a_s, rate = b_s);
  }
  s <- rep(1, n);
  
  # Initialize G -------------------------------------------------------------------------------
  if (g_constraint) {
    tau_mu <- 1;
    sigma_g <- 0.1;
    tau_g <- 1;
    c_g <- 0;
    q <- runif(L, 0, 1);
    mu <- rnorm(L, 0, tau_mu);
    pi <- rep(1/L, L);
    xi <- sample.int(L, p, replace = TRUE, prob = pi);
    psi <- sample.int(2, p, replace = TRUE) - 1;
    g <- exp(rnorm(p, psi*mu[xi] + (1 - psi)*(c_g - q[xi]*mu[xi])/(1 - q[xi]), sigma_g));
  } else {
    a_g <- 1;
    b_g <- 1;
    g <- rgamma(p, shape = a_g, rate = b_g);
  }
  g <- rep(1, p);
  
  # Initialize Z -------------------------------------------------------------------------------
  if (observation_supervised) {
    z <- zt;
  } else {
    z <- organize_label(sample.int(start_K, n, replace = TRUE));
  }
  K <- length(unique(z));
  if (!alpha_supervised) {
    etaeta <- rbeta(1, alpha + 1, n);
    a_alpha <- 1;
    b_alpha <- 10;
    prob <- c(a_alpha + K - 1, n*(b_alpha - log(etaeta)));
    prob <- prob/sum(prob);
    if (runif(1) < prob[1]) {
      flag <- 0;
      alpha <- rgamma(1, shape = a_alpha + K, rate = b_alpha - log(etaeta));
    } else {
      flag <- 1;
      alpha <- rgamma(1, shape = a_alpha + K - 1, rate = b_alpha - log(etaeta));
    }
  }
  
  # Initialize d_0 ----------------------------------------------------------------------------
  # d_0 <- rgamma(1, shape = a + sum(X[, which(gamma == 0)]), rate = b + sum(g[which(gamma == 0)])*sum(s));
  # d_0 <- 1;
  
  # Initialize theta ----------------------------------------------------------------------------
  theta <- matrix(0, nrow = n, ncol = p);
  if (ZIP) {
    # d <- matrix(d_0, nrow = n, ncol = p);
    # for (k in 1:K) {
    #   for (j in which(gamma == 1)) {
    #     d[z == k, j] <- rgamma(1, shape = a + sum(X[z == k, j]), rate = b + g*sum(s[z == k]))
    #   }
    # }
    zero_list <- which(X == 0, arr.ind = TRUE);
    theta[cbind(zero_list[, 1], zero_list[, 2])] <- rbinom(dim(zero_list)[1], 1, p_0);
    if (!p_0_supervised) {
      a_p_0 <- 1;
      b_p_0 <- 1;
      p_0 <- rbeta(1, a_p_0 + sum(theta == 1), b_p_0 + sum(theta == 0));
    }
  }
  
  # Return the results -------------------------------------------------------------------------
  if (s_constraint & g_constraint) {
    return(list(theta = theta, p_0 = p_0, z = z, alpha = alpha, gamma = gamma, omega = omega, d_0 = d_0, s = s, g = g, t = t, eta = eta, phi = phi, xis = xis, psis = psis, q = q, mu = mu, pi = pi, xi = xi, psi = psi));
  } else if (s_constraint & !g_constraint) {
    return(list(theta = theta, p_0 = p_0, z = z, alpha = alpha, gamma = gamma, omega = omega, d_0 = d_0, s = s, g = g, t = t, eta = eta, phi = phi, xis = xis, psis = psis));
  } else if (!s_constraint & g_constraint) {
    return(list(theta = theta, p_0 = p_0, z = z, alpha = alpha, gamma = gamma, omega = omega, d_0 = d_0, s = s, g = g, q = q, mu = mu, pi = pi, xi = xi, psi = psi));
  } else {
    return(list(theta = theta, p_0 = p_0, z = z, alpha = alpha, gamma = gamma, omega = omega, d_0 = d_0, s = s, g = g));
  } 
}

# **********************************************************************************************
# **********************************************************************************************

core = function(parameters, hyperparameters, environment, X) {
  # Load data, parameters, hyperparameters, and environment information ------------------------
  n <- dim(X)[1];
  p <- dim(X)[2];
  
  observation_supervised <- environment[1];
  variable_supervised <- environment[2];
  s_constraint <- environment[3];
  g_constraint <- environment[4];
  alpha_supervised <- environment[5];
  omega_supervised <- environment[6];
  ZIP <- environment[7];
  p_0_supervised <- environment[8];
  
  a <- hyperparameters[1];
  b <- hyperparameters[2];
  M <- hyperparameters[6];
  L <- hyperparameters[7];
  if (s_constraint) {
    tau_eta <- 1;
    sigma_s <- 0.1;
    tau_s <- 1;
    c_s <- 0;
  } else {  
    a_s <- 1;
    b_s <- 1;
  }
  if (g_constraint) {
    tau_mu <- 1;
    sigma_g <- 0.1;
    tau_g <- 1;
    c_g <- 0;
  } else {
    a_g <- 1;
    b_g <- 1;
  }
  if (!omega_supervised) {
    a_omega <- 0.2;
    b_omega <- 1.8;
  }
  if (ZIP) {
    zero_list <- which(X == 0, arr.ind = TRUE);
    if (!p_0_supervised) {
      a_p_0 <- 1;
      b_p_0 <- 1;
    }
  }
  
  theta <- parameters$theta;
  p_0 <- parameters$p_0;
  z <- parameters$z;
  alpha <- parameters$alpha;
  K <- length(unique(z));
  gamma <- parameters$gamma;
  omega <- parameters$omega;
  d <- parameters$d;
  d_0 <- parameters$d_0;
  s <- parameters$s;
  g <- parameters$g;
  if (s_constraint) {
    t <- parameters$t;
    eta <- parameters$eta;
    phi <- parameters$phi;
    xis <- parameters$xis;
    psis <- parameters$psis;
  } 
  if (g_constraint) {
    q <- parameters$q;
    mu <- parameters$mu;
    pi <- parameters$pi;
    xi <- parameters$xi;
    psi <- parameters$psi;
  }
  
  # Start MCMC sampling -----------------------------------------------------------------------
  start_time <- proc.time();
  
  # Sampling for updating S -------------------------------------------------------------------
  if (!s_constraint) {
    # Metropolis-Hastings sampling for updating S
    for (i in 1:n) {
      s_temp <- rgamma(1, shape = s[i], rate = 1);
      if (s_temp != 0) {
        hastings <- loglikelihood_n2s(X[i, which(theta[i,] == 0)], gamma[which(theta[i,] == 0)], d_0, s_temp, g[which(theta[i,] == 0)], a, b) + dgamma(s_temp, shape = a_s, rate = b_s, log = TRUE) - loglikelihood_n2s(X[i, which(theta[i,] == 0)], gamma[which(theta[i,] == 0)], d_0, s[i], g[which(theta[i,] == 0)], a, b) - dgamma(s[i], shape = a_s, rate = b_s, log = TRUE);  
        if (!is.na(hastings)) {
          if (hastings >= log(runif(1))) {
            s[i] <- s_temp;
          }
        }
      }
    }
  } else {
    # Metropolis-Hastings sampling for updating S
    for (i in 1:n) {
      s_temp <- exp(rnorm(1, mean = log(s[i]), sd = tau_s));
      hastings <- loglikelihood_n2s(X[i, which(theta[i,] == 0)], gamma[which(theta[i,] == 0)], d_0, s_temp, g[which(theta[i,] == 0)], a, b) + dnorm(log(s_temp), psis[i]*eta[xis[i]] + (1 - psis[i])*(c_s - t[xis[i]]*eta[xis[i]])/(1 - t[xis[i]]), sigma_s, log = TRUE) - loglikelihood_n2s(X[i, which(theta[i,] == 0)], gamma[which(theta[i,] == 0)], d_0, s[i], g[which(theta[i,] == 0)], a, b) - dnorm(log(s[i]), psis[i]*eta[xis[i]] + (1 - psis[i])*(c_s - t[xis[i]]*eta[xis[i]])/(1 - t[xis[i]]), sigma_s, log = TRUE);
      if (hastings >= log(runif(1))) {
        s[i] <- s_temp;
      }
    }
    # Gibbs sampling for updating eta
    for (m in 1:M) {
      a_temp <- (sum(xis == m & psis == 1) + sum(xis == m & psis == 0)*(t[m]/(1 - t[m]))^2)/sigma_s^2 + 1/tau_eta^2;
      b_temp <- (sum(log(s[which(xis == m & psis == 1)])) - (t[m]/(1 - t[m]))*sum(log(s[which(xis == m & psis == 0)]) - c_s/(1 - t[m])))/sigma_s^2;
      eta[m] <- rnorm(1, mean = b_temp/a_temp, sd = sqrt(1/a_temp));
    }
    # Gibbs sampling for updating Xi_s
    for (i in 1:n) {
      prob <- dnorm(log(s[i]), (psis[i]*eta + (1 - psis[i])*(c_s - t*eta)/(1 - t)), sigma_s)*phi;
      if (sum(prob) == 0) {
        prob <- rep(1/M, M);
      } else {
        prob <- prob/sum(prob);
      }
      xis[i] <- sample.int(M, 1, replace = TRUE, prob);
    }
    # Gibbs sampling for updating Phi
    v <- rep(0, M);
    temp <- 1;
    for (m in 1:M) {
      v[m] <- rbeta(1, 1 + sum(xis == m), 1 + sum(xis > m));
      if (m == 1) {
        phi[m] <- v[1];
      } else {
        temp <- temp*(1 - v[m - 1]);
        phi[m] <- temp*v[m];
      }
    }
    # Gibbs sampling for updating Psi_s
    for (i in 1:n){
      prob <- rep(0, 2);
      prob[2] <- dnorm(log(s[i]), mean = eta[xis[i]], sd = sigma_s)*t[xis[i]];
      prob[1] <- dnorm(log(s[i]), mean = (c_s - t[xis[i]]*eta[xis[i]])/(1 - t[xis[i]]), sd = sigma_s)*(1 - t[xis[i]]);
      if (sum(prob) == 0) {
        prob <- c(0.5, 0.5);
      } else {
        prob <- prob/sum(prob);
      }
      psis[i] <- sample.int(2, 1, replace = TRUE, prob) - 1
    }
    # Gibbs sampling for updating T
    for (m in 1:M) {
      t[m] <- rbeta(1, 1 + sum(xis == m & psis == 1), 1 + sum(xis == m & psis == 0));
    }
  }
  
  # Sampling for updating G -------------------------------------------------------------------
  if (!g_constraint) {
    # Gibbs sampling for updating G
    g[gamma == 0] <- rgamma(p - sum(gamma), shape = a_g + colSums(X)[gamma == 0], rate = b_g + d_0*sum(s));
    # Metropolis-Hastings sampling for updating G
    for (j in which(gamma == 1)) {
      g_temp <- rgamma(1, shape = g[j], rate = 1);
      if (g_temp != 0) {
        hastings <- loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], gamma[j], d_0, s[which(theta[, j] == 0)], g_temp, a, b) + dgamma(g_temp, shape = a_g, rate = b_g, log = TRUE) - loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], gamma[j], d_0, s[which(theta[, j] == 0)], g[j], a, b) - dgamma(g[j], shape = a_g, rate = b_g, log = TRUE);  
        if (!is.na(hastings)) {
          if (hastings >= log(runif(1))) {
            g[j] <- g_temp;
          }
        }
      }
    }
  } else {
    # Metropolis-Hastings sampling for updating G
    for (j in 1:p) {
      g_temp <- exp(rnorm(1, mean = log(g[j]), sd = tau_g));
      hastings <- loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], gamma[j], d_0, s[which(theta[, j] == 0)], g_temp, a, b) + dnorm(log(g_temp), (psi[j]*mu[xi[j]] + (1 - psi[j])*(c_g - q[xi[j]]*mu[xi[j]])/(1 - q[xi[j]])), sigma_g, log = TRUE) - loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], gamma[j], d_0, s[which(theta[, j] == 0)], g[j], a, b) - dnorm(log(g[j]), (psi[j]*mu[xi[j]] + (1 - psi[j])*(c_g - q[xi[j]]*mu[xi[j]])/(1 - q[xi[j]])), sigma_g, log = TRUE);  
      if (hastings >= log(runif(1))) {
        g[j] <- g_temp;
      }
    }
    # Gibbs sampling for updating M
    for (l in 1:L) {
      a_temp <- (sum(xi == l & psi == 1) + sum(xi == l & psi == 0)*(q[l]/(1 - q[l]))^2)/sigma_g^2 + 1/tau_mu^2;
      b_temp <- (sum(log(g[which(xi == l & psi == 1)])) - (q[l]/(1 - q[l]))*sum(log(g[which(xi == l & psi == 0)]) - c_g/(1 - q[l])))/sigma_g^2;
      mu[l] <- rnorm(1, mean = b_temp/a_temp, sd = sqrt(1/a_temp));
    } 
    # Gibbs sampling for updating Xi_g
    for (j in 1:p) {
      prob <- dnorm(log(g[j]), psi[j]*mu + (1 - psi[j])*(c_g - q*mu)/(1 - q), sigma_g)*pi;
      if (sum(prob) == 0) {
        prob <- rep(1/L, L);
      } else {
        prob <- prob/sum(prob);
      }
      xi[j] <- sample.int(L, 1, replace = TRUE, prob);
    }
    # Gibbs sampling for updating Pi
    v <- rep(0, L);
    temp <- 1;
    for (l in 1:L) {
      v[l] <- rbeta(1, 1 + sum(xi == l), 1 + sum(xi > l));
      if (l == 1) {
        pi[l] <- v[1];
      } else {
        temp <- temp*(1 - v[l - 1]);
        pi[l] <- temp*v[l];
      }
    }
    # Gibbs sampling for updating Psi_g
    for (j in 1:p){
      prob <- rep(0, 2);
      prob[2] <- dnorm(log(g[j]), mean = mu[xi[j]], sd = sigma_g)*q[xi[j]];
      prob[1] <- dnorm(log(g[j]), mean = (c_g - q[xi[j]]*mu[xi[j]])/(1 - q[xi[j]]), sd = sigma_g)*(1 - q[xi[j]]);
      if (sum(prob) == 0) {
        prob <- c(0.5, 0.5);
      } else {
        prob <- prob/sum(prob);
      }
      psi[j] <- sample.int(2, 1, replace = TRUE, prob) - 1
    }
    # Gibbs sampling for updating Q
    for (l in 1:L) {
      q[l] <- rbeta(1, 1 + sum(xi == l & psi == 1), 1 + sum(xi == l & psi == 0));
    }
  }
  # g <- colSums(X)
  
  # Metropolis-Hastings sampling for updating Gamma -------------------------------------------
  map_gamma <- 0;
  if (!variable_supervised) {
    for (e in 1:20) {
      temp <- sample.int(3, 1);
      if (temp == 1) {
        # Add
        if (sum(gamma) != length(gamma)) {
          j <- sample(which(gamma == 0 & which(colSums(theta) != n)), 1);
          # hastings <- loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], 1, d_0, s[which(theta[, j] == 0)], g[j], a, b) - loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], 0, d_0, s[which(theta[, j] == 0)], g[j], a, b) + log(p - sum(gamma)) - log(sum(gamma) + 1);
          hastings <- loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], 1, d_0, s[which(theta[, j] == 0)], g[j], a, b) - loglikelihood_p2s(X[which(theta[, j] == 0), j], rep(1, length(which(theta[, j] == 0))), 1, d_0, s[which(theta[, j] == 0)], g[j], a, b) + log(p - sum(gamma)) - log(sum(gamma) + 1);
          if (omega_supervised) {
            hastings <- hastings + log(omega) - log(1 - omega); 
          } else {
            hastings <- hastings + lgamma(a_omega + sum(gamma) + 1) + lgamma(b_omega + p - sum(gamma) - 1) - lgamma(a_omega + sum(gamma)) - lgamma(b_omega + p - sum(gamma));
          }
          if (hastings >= log(runif(1))) {
            gamma[j] <- 1;
            # acceptance[1] <- acceptance[1] + 1;
          }
        }
      } else if (temp == 2) {
        # Delete
        if (sum(gamma) > 1) {
          j <- sample(which(gamma == 1 & which(colSums(theta) != n)), 1);
          # hastings <- loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], 0, d_0, s[which(theta[, j] == 0)], g[j], a, b) - loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], 1, d_0, s[which(theta[, j] == 0)], g[j], a, b) + log(sum(gamma)) - log(p - sum(gamma) + 1);
           hastings <- loglikelihood_p2s(X[which(theta[, j] == 0), j], rep(1, length(which(theta[, j] == 0))), 1, d_0, s[which(theta[, j] == 0)], g[j], a, b) - loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], 1, d_0, s[which(theta[, j] == 0)], g[j], a, b) + log(sum(gamma)) - log(p - sum(gamma) + 1);
          if (omega_supervised) {
            hastings <- hastings + log(1 - omega) - log(omega); 
          } else {
            hastings <- hastings + lgamma(a_omega + sum(gamma) - 1) + lgamma(b_omega + p - sum(gamma) + 1) - lgamma(a_omega + sum(gamma)) - lgamma(b_omega + p - sum(gamma));
          }
          if (hastings >= log(runif(1))) {
            gamma[j] <- 0;
            # acceptance[2] <- acceptance[2] + 1;
          }
        }
      } else if (temp == 3) {
        # Swap
        if (sum(gamma) != 0 && sum(gamma) != length(gamma)) {
          j_0 <- sample(which(gamma == 0 & which(colSums(theta) != n)), 1);
          j_1 <- sample(which(gamma == 1 & which(colSums(theta) != n)), 1);
          # hastings <- loglikelihood_p2s(X[which(theta[, j_0] == 0), j_0], z[which(theta[, j_0] == 0)], 1, d_0, s[which(theta[, j_0] == 0)], g[j_0], a, b) + loglikelihood_p2s(X[which(theta[, j_1] == 0), j_1], z[which(theta[, j_1] == 0)], 0, d_0, s[which(theta[, j_1] == 0)], g[j_1], a, b) - loglikelihood_p2s(X[which(theta[, j_0] == 0), j_0], z[which(theta[, j_0] == 0)], 0, d_0, s[which(theta[, j_0] == 0)], g[j_0], a, b) - loglikelihood_p2s(X[which(theta[, j_1] == 0), j_1], z[which(theta[, j_1] == 0)], 1, d_0, s[which(theta[, j_1] == 0)], g[j_1], a, b);
          hastings <- loglikelihood_p2s(X[which(theta[, j_0] == 0), j_0], z[which(theta[, j_0] == 0)], 1, d_0, s[which(theta[, j_0] == 0)], g[j_0], a, b) + loglikelihood_p2s(X[which(theta[, j_1] == 0), j_1], rep(1, length(which(theta[, j_1] == 0))), 1, d_0, s[which(theta[, j_1] == 0)], g[j_1], a, b) - loglikelihood_p2s(X[which(theta[, j_0] == 0), j_0], rep(1, length(which(theta[, j_0] == 0))), 1, d_0, s[which(theta[, j_0] == 0)], g[j_0], a, b) - loglikelihood_p2s(X[which(theta[, j_1] == 0), j_1], z[which(theta[, j_1] == 0)], 1, d_0, s[which(theta[, j_1] == 0)], g[j_1], a, b);
          if (hastings >= log(runif(1))) {
            gamma[j_0] <- 1;
            gamma[j_1] <- 0;
            # acceptance[3] <- acceptance[3] + 1;
          }
        }
      }
    }
    
    if (!omega_supervised) {
      map_gamma <- map_gamma + dbeta(omega, a_omega, b_omega, log = TRUE);
      omega <- rbeta(1, a_omega + sum(gamma == 1), b_omega + sum(gamma == 0));
    }
    
    # map_gamma <- loglikelihood_p2s(X, z, gamma, d_0, s, g, a, b) + dbinom(sum(gamma), p, omega, log = TRUE);
    map_gamma <- dbinom(sum(gamma), p, omega, log = TRUE);
    for (j in 1:p) {
      map_gamma <- map_gamma + loglikelihood_p2s(X[which(theta[, j] == 0), j], z[which(theta[, j] == 0)], 1, d_0, s[which(theta[, j] == 0)], g[j], a, b);
    }
  }

  # Gibbs sampling for updating d_0 -----------------------------------------------------------
  # d_0 <- rgamma(1, shape = a + sum(X[, which(gamma == 0)]), rate = b + sum(g[which(gamma == 0)])*sum(s));
  
  # Gibbs sampling for updating Z ------------------------------------------------------------
  map_z <- 0;
  if (!observation_supervised) {
    for (i in 1:n) {
      valid_gamma_index <- intersect(which(theta[i,] == 0), which(gamma == 1));
      if (length(valid_gamma_index) != 0) {
        prob <- rep(0, K + 1);
        for (k in 1:K) {
          z_temp <- z;
          z_temp[i] <- k;
          z_temp2 <- z;
          z_temp2[i] <- 0;
          if (length(valid_gamma_index) == 1) {
            if (sum(z_temp2 == k) == 0) {
              prob[k] <- (a)*log(b) - lgamma(a) + lgamma(a + (X[z_temp == k, valid_gamma_index])) - (a + (X[z_temp == k, valid_gamma_index]))*log(b + g[valid_gamma_index]*sum(s[z_temp == k])) + log(sum(z_temp2 == k)) - log(n + alpha - 1); 
            } else if (sum(z_temp2 == k) == 1) {
              prob[k] <- (a + (X[z_temp2 == k, valid_gamma_index]))*log(b + g[valid_gamma_index]*sum(s[z_temp2 == k])) - lgamma(a + (X[z_temp2 == k, valid_gamma_index])) + lgamma(a + sum(X[z_temp == k, valid_gamma_index])) - (a + sum(X[z_temp == k, valid_gamma_index]))*log(b + g[valid_gamma_index]*sum(s[z_temp == k])) + log(sum(z_temp2 == k)) - log(n + alpha - 1); 
            } else {
              prob[k] <- (a + sum(X[z_temp2 == k, valid_gamma_index]))*log(b + g[valid_gamma_index]*sum(s[z_temp2 == k])) - lgamma(a + sum(X[z_temp2 == k, valid_gamma_index])) + lgamma(a + sum(X[z_temp == k, valid_gamma_index])) - (a + sum(X[z_temp == k, valid_gamma_index]))*log(b + g[valid_gamma_index]*sum(s[z_temp == k])) + log(sum(z_temp2 == k)) - log(n + alpha - 1);
            }
          } else {
            if (sum(z_temp2 == k) == 0) {
              prob[k] <- sum((a)*log(b) - lgamma(a) + lgamma(a + (X[z_temp == k, valid_gamma_index])) - (a + (X[z_temp == k, valid_gamma_index]))*log(b + g[valid_gamma_index]*sum(s[z_temp == k]))) + log(sum(z_temp2 == k)) - log(n + alpha - 1); 
            } else if (sum(z_temp2 == k) == 1) {
              prob[k] <- sum((a + (X[z_temp2 == k, valid_gamma_index]))*log(b + g[valid_gamma_index]*sum(s[z_temp2 == k])) - lgamma(a + (X[z_temp2 == k, valid_gamma_index])) + lgamma(a + colSums(X[z_temp == k, valid_gamma_index])) - (a + colSums(X[z_temp == k, valid_gamma_index]))*log(b + g[valid_gamma_index]*sum(s[z_temp == k]))) + log(sum(z_temp2 == k)) - log(n + alpha - 1); 
            } else {
              prob[k] <- sum((a + colSums(X[z_temp2 == k, valid_gamma_index]))*log(b + g[valid_gamma_index]*sum(s[z_temp2 == k])) - lgamma(a + colSums(X[z_temp2 == k, valid_gamma_index])) + lgamma(a + colSums(X[z_temp == k, valid_gamma_index])) - (a + colSums(X[z_temp == k, valid_gamma_index]))*log(b + g[valid_gamma_index]*sum(s[z_temp == k]))) + log(sum(z_temp2 == k)) - log(n + alpha - 1); 
            }
          }
        }
        if (sum(z_temp2 == k) != 0) {
          prob[K + 1] <- sum(a*log(b) - lgamma(a) + lgamma(a + X[i, valid_gamma_index]) - (a + X[i, valid_gamma_index])*log(b + s[i]*g[valid_gamma_index])) + log(alpha) - log(n + alpha - 1);
        } else {
          prob <- prob[1:K];
        }
        prob <- exp(prob - max(prob));
        prob <- prob/sum(prob);
        z[i] <- sample.int(length(prob), 1, prob = prob);
        z <- organize_label(z);
        K <- length(unique(z));
      }
    }
    for (i in 1:n) {
      z_temp <- z;
      z_temp[i] <- 0;
      a_temp <- rep(NA, p);
      b_temp <- rep(NA, p);
      for (j in which(gamma == 1)) {
        a_temp[j] <- a + sum(X[z_temp == z[i], j]);
        b_temp[j] <- b + sum(s[z_temp == z[i]])*g[j];
      }
      if (sum(z_temp == z[i]) == 0) {
        map_z <- map_z + loglikelihood_n2s(X[i,], gamma, d_0, s[i], g, a, b) + log(alpha) - log(n + alpha - 1);
      } else {
        map_z <- map_z + loglikelihood_n2(X[i,], gamma, d_0, s[i], g, a_temp, b_temp) + log(sum(z_temp == z[i])) - log(n + alpha - 1);
      }
    }
    if (!alpha_supervised) {
      etaeta <- rbeta(1, alpha + 1, n);
      a_alpha <- 1;
      b_alpha <- 10;
      map_z <- map_z + dgamma(alpha, shape = a_alpha, rate = b_alpha, log = TRUE);
      prob <- c(a_alpha + K - 1, n*(b_alpha - log(etaeta)));
      prob <- prob/sum(prob);
      if (runif(1) <= prob[1]) {
        alpha <- rgamma(1, shape = a_alpha + K, rate = b_alpha - log(etaeta));
      } else {
        alpha <- rgamma(1, shape = a_alpha + K - 1, rate = b_alpha - log(etaeta));
      }
    }
  }
  
  # Gibbs sampling for updating Theta -------------------------------------------
  if (ZIP) {
    d <- matrix(d_0, nrow = n, ncol = p);
    for (k in 1:K) {
      for (j in which(gamma == 1)) {
        d[z == k, j] <- rgamma(1, shape = a + sum(X[intersect(which(theta[, j] == 0), which(z == k)), j]), rate = b + g*sum(s[intersect(which(theta[, j] == 0), which(z == k))]));
        # d[z == k, j] <- rgamma(1, shape = a + sum(X[z == k, j]), rate = b + g*sum(s[z == k]));
      }
    }
    # zero_list <- which(X == 0, arr.ind = TRUE);
    theta[cbind(zero_list[, 1], zero_list[, 2])] <- rbinom(dim(zero_list)[1], 1, p_0/(p_0 + (1 - p_0)*exp(-s[zero_list[, 1]]*g[zero_list[, 2]]*d[cbind(zero_list[, 1], zero_list[, 2])])));
    if (!p_0_supervised) {
      # a_p_0 <- 1;
      # b_p_0 <- 1;
      p_0 <- rbeta(1, a_p_0 + sum(theta == 1), b_p_0 + sum(theta == 0));
    }
  }
  
  end_time <- proc.time();
  
  # Return the results ------------------------------------------------------------------------
  time <- end_time - start_time;
  if (s_constraint & g_constraint) {
    return(list(time = time, theta = theta, p_0 = p_0, z = z, d_0 = d_0, gamma = gamma, alpha = alpha, omega = omega, d = d, s = s, g = g, map_gamma = map_gamma, map_z = map_z, t = t, eta = eta, phi = phi, xis = xis, psis = psis, q = q, mu = mu, pi = pi, xi = xi, psi = psi));
  } else if (s_constraint & !g_constraint) {
    return(list(time = time, theta = theta, p_0 = p_0, z = z, d_0 = d_0, gamma = gamma, alpha = alpha, omega = omega, d = d, s = s, g = g, map_gamma = map_gamma, map_z = map_z, t = t, eta = eta, phi = phi, xis = xis, psis = psis));
  } else if (!s_constraint & g_constraint) {
    return(list(time = time, theta = theta, p_0 = p_0, z = z, d_0 = d_0, gamma = gamma, alpha = alpha, omega = omega, d = d, s = s, g = g, map_gamma = map_gamma, map_z = map_z, q = q, mu = mu, pi = pi, xi = xi, psi = psi));
  } else {
    return(list(time = time, theta = theta, p_0 = p_0, z = z, d_0 = d_0, gamma = gamma, alpha = alpha, omega = omega, d = d, s = s, g = g, map_gamma = map_gamma, map_z = map_z));
  } 
}