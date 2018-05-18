# Generate simulated data
simulated_data_generator = function(z, gamma, s, g, mu, sigma) {
  n <- length(z);
  p <- length(gamma);
  K <- length(unique(z));
  d <- matrix(NA, nrow = n, ncol = p);
  d[, gamma == 0] <- 1;
  for (j in which(gamma == 1)) {
    temp <- exp(rnorm(K, mu, sigma));
    for (k in 1:K) {
      d[z == k, j] <- temp[k];
    }
  }
  X <- matrix(rpois(n*p, d*s%*%t(g)), nrow = n, ncol = p);
  return (list(X = X, d = d));
}

# Calculate the full data log-likelihood
loglikelihood = function(X, z, gamma, d, d_0, s, g) {
  return (sum(dpois(X, d*s%*%t(g), log = TRUE)));
}

# Calculate the collapsed data log-likelihood by observation
loglikelihood_n2 = function(X, gamma, d_0, s, g, a, b) {
  if (is.matrix(X)) {
    n <- dim(X)[1];
    loglklh <- 0;
    for (i in 1:n) {
      loglklh <- loglklh + loglikelihood_n2(X[i,], gamma, d_0, s[i], g, a, b);
    }
  } else {
    loglklh <- sum(X)*log(s) + sum(X*log(g), na.rm = TRUE) - sum(lfactorial(X));
    loglklh <- loglklh - s*d_0*sum(g[gamma == 0]) + sum(X[gamma == 0])*log(d_0);
    if (length(a) == 1) {
      for (j in which(gamma == 1)) {
        loglklh <- loglklh + a*log(b) - lgamma(a) + lgamma(a + X[j]) - (a + X[j])*log(b + s*g[j]);
      }
    } else {
      for (j in which(gamma == 1)) {
        loglklh <- loglklh + a[j]*log(b[j]) - lgamma(a[j]) + lgamma(a[j] + sum(X[j])) - (a[j] + sum(X[j]))*log(b[j] + s*sum(g[j]));
      }
    }
  }
  return (loglklh);
}

loglikelihood_n2s = function(X, gamma, d_0, s, g, a, b) {
  if (is.matrix(X)) {
    loglklh <- sum(X*log(s%*%t(g))) - sum(lfactorial(X)) + sum(X[, gamma == 0])*log(d_0) - sum((s%*%t(g[gamma == 0]))*d_0) + length(z)*sum(gamma == 1)*(a*log(b) - lgamma(a)) + sum(lgamma(a + X[, gamma == 1])) - sum((a + X[, gamma == 1])*log(b + s%*%t(g[gamma == 1])));
  } else {
    loglklh <- sum(X*log(s*g)) - sum(lfactorial(X)) + sum(X[gamma == 0])*log(d_0) - sum((s*g[gamma == 0])*d_0) + sum(gamma == 1)*(a*log(b) - lgamma(a)) + sum(lgamma(a + X[gamma == 1])) - sum((a + X[gamma == 1])*log(b + s*g[gamma == 1]));
  }
  return (loglklh);
}

# Calculate the collapsed data log-likelihood by variable
loglikelihood_p2 = function(X, z, gamma, d_0, s, g, a, b) {
  if (is.matrix(X)) {
    p <- dim(X)[2];
    loglklh <- 0;
    for (j in 1:p) {
      loglklh <- loglklh + loglikelihood_p2(X[, j], z, gamma[j], d_0, s, g[j], a, b);
    }
  } else {
    K <- length(unique(z));
    loglklh <- sum(X)*log(g) + sum(X*log(s)) - sum(lfactorial(X));
    if (gamma == 0) {
      loglklh <- loglklh - sum(s)*d_0*g + sum(X)*log(d_0);
    } else {
      # for (i in 1:n) {
      #   if (length(a) == 1) {
      #     loglklh <- loglklh + a*log(b) - lgamma(a) + lgamma(a + sum(X[i])) - (a + sum(X[i]))*log(b + g*sum(s[i]));
      #   } else {
      #     loglklh <- loglklh + a[i]*log(b[i]) - lgamma(a[i]) + lgamma(a[i] + sum(X[i])) - (a[i] + sum(X[i]))*log(b[i] + g*sum(s[i]));
      #   }
      # }
      for (k in 1:K) {
        loglklh <- loglklh + a*log(b) - lgamma(a) + lgamma(a + sum(X[z == k])) - (a + sum(X[z == k]))*log(b + g*sum(s[z == k]));
      }
    }
  }
  return (loglklh);
}

loglikelihood_p2s = function(X, z, gamma, d_0, s, g, a, b) {
  if (is.matrix(X)) {
    loglklh <- sum(X*log(s%*%t(g))) - sum(lfactorial(X)) + sum(X[, gamma == 0])*log(d_0) - sum((s%*%t(g[gamma == 0]))*d_0) + length(unique(z))*sum(gamma == 1)*(a*log(b) - lgamma(a));
    for (k in unique(z)) {
      X_temp <- X[z == k, gamma == 1];
      if (is.matrix(X_temp)) {
        loglklh <- loglklh + sum(lgamma(a + colSums(X_temp))) - sum((a + colSums(X_temp))*log(b + g[gamma == 1]*sum(s[z == k])));
      } else {
        loglklh <- loglklh + sum(lgamma(a + X_temp)) - sum((a + X_temp)*log(b + g[gamma == 1]*sum(s[z == k])));
      }
    }
  } else {
    loglklh <- sum(X*log(s*g) - lfactorial(X));
    if (gamma == 0) {
      loglklh <- loglklh + sum(X)*log(d_0) - sum(s*g*d_0);
    } else {
      loglklh <- loglklh + length(unique(z))*(a*log(b) - lgamma(a));
      for (k in unique(z)) {
        loglklh <- loglklh + lgamma(a + sum(X[z == k])) - (a + sum(X[z == k]))*log(b + g*sum(s[z == k]));
      }
    }
  }
  return (loglklh);
}

# Evaluate results using clustering error rate (CER)
cer = function(zt, z) {
  cer <- 0
  n <- length(z);
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (zt[i] == zt[j] && z[i] != z[j]) {
        cer = cer + 1;
      } else if (zt[i] != zt[j] && z[i] == z[j]) {
        cer = cer + 1;
      }
    }
  }
  return (cer/choose(n, 2));
}

# Evaluate results using adjusted Rand index (ARI)
ari <- function(zt, z) {
  library(mclust);
  return (adjustedRandIndex(zt, z));
}

# Visualize the simulated data
data_visualizer = function(X, z, gamma, s, g, mean = FALSE) {
  n <- dim(X)[1];
  p <- dim(X)[2];
  if (mean == FALSE) {
    X <- X/(s%*%t(g));
  } else {
    d <- rep(1, length(unique(z)));
    K <- length(d);
    for (k in 1:K) {
      d[k] <- sum(X[z == k, gamma == 1])/sum(s[z == k])/sum(g[gamma == 1]);
      X[z == k, gamma == 1] <- d[k];
    }
    d_0 <- sum(X[, gamma == 0])/sum(s)/sum(g[gamma == 0]);
    X[, gamma == 0] <- d_0;
  }
  X <- X[, sort.int(gamma, decreasing = TRUE, index.return = TRUE)$ix];
  X <- X[sort.int(z, decreasing = TRUE, index.return = TRUE)$ix,];
  min <- min(X, na.rm = TRUE);
  max <- max(X, na.rm = TRUE);
  ColorRamp <- rgb(seq(0.95, 0.99, length = 50), seq(0.95, 0.05, length = 50), seq(0.95, 0.05, length = 50));
  ColorLevels <- seq(min, max, length = length(ColorRamp));
  layout(matrix(data = c(1, 2), nrow = 1, ncol = 2), widths = c(4, 1), heights = c(1, 1));
  par(mar = c(5,5, 2.5, 1), font = 2);
  image(1:ncol(X), 1:nrow(X), t(X), col = ColorRamp, xlab = "Variable", ylab = "Observation", axes = FALSE, zlim = c(min, max), main = NA);
  box();
  axis(side = 2, at = seq(1, n, 1), labels = sort(z, decreasing = TRUE), cex.axis = 1.0);
  axis(side = 1, at = seq(1,p,1), labels = sort(gamma, decreasing = TRUE), las = 1, cex.axis = 1);
}

# Switching label
switch_label = function(X, z, gamma) {
  K <- length(unique(z));
  z_new <- rep(NA, 1, length(z));
  mu <- rep(NA, 1, K);
  for(k in 1:K) {
    mu[k] <- mean(X[z == k, gamma == 1]);
  }
  temp <- order(mu, decreasing = TRUE);
  for(k in 1:K) {
    z_new[z == temp[k]] <- k;
  }
  return (z_new);
}

# Organizing label
organize_label = function(z) {
  z_new <- rep(NA, 1, length(z));
  count <- 1;
  for(k in unique(z)) {
    z_new[z == k] <- count;
    count <- count + 1;
  }
  return (z_new);
}

#
tabulate_error = function(gamma_true, gamma) {
  table = matrix(0L, 2, 2);
  p <- length(gamma_true);
  for (i in 1:p) {
    table[gamma[i] + 1, gamma_true[i] + 1] <- table[gamma[i] + 1, gamma_true[i] + 1] + 1;
  }
  return (table);
}

#
logprob = function(prob) {
  n <- length(prob);
  temp <- rep(0, n);
  for (i in 1:n) {
    temp[i] <- 1/sum(exp(prob - prob[i]));
  }
  return (temp);
}

#
pz = function(z, alpha) {
  K <- length(unique(z));
  temp <- K*log(alpha) + lgamma(alpha);
  for (k in 1:K) {
    temp <- temp + lgamma(sum(z == k));
  }
  temp <- temp - lgamma(alpha + length(z));
  return (temp);
}