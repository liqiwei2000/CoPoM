### Set working directory --------------------------------------------------------------
setwd("~/Dropbox/Qiwei");

### Load function ----------------------------------------------------------------------
source("./code/functions.R");

### Input settings ---------------------------------------------------------------------
# Features
p <- 1000; # Number of features
p_gamma <- 50; # Number of discriminatory features
gammat <- rep(0, p); 
gammat[1:p_gamma] <- 1;
gammat <- gammat[sample(p, size = p, replace = FALSE)]; # True gamma

# Observations
zt <- c(rep(1, 6), rep(2, 15), rep(3, 9)); # True z
n <- length(zt); # Number of observations
K <- length(unique(zt)); # Number of groups

# Scaling factors
st <- runif(n, 0.5, 1.5); # True s
gt <- rexp(p, 1/3); # True g

# Group means
mu <- 0;
sigma <- 1; # Between-group standard deviation

### Generate simulated Poisson data -----------------------------------------------------
temp <- simulated_data_generator(zt, gammat, st, gt, mu, sigma);
X <- temp$X;
dt <- temp$d; # True d
# Remove those features with all zero counts
if (length(which(colSums(X) == 0)) > 0) {
  gt <- gt[-which(colSums(X) == 0)];
  gammat <- gammat[-which(colSums(X) == 0)];
  X <- X[, -which(colSums(X) == 0)];
}

### Generate simulated zero-inflated Poisson data --------------------------------------
# Percentage of zeros
p_0 <- 0.25;
thetat <- matrix(sample(c(0, 1), dim(X)[1]*dim(X)[2], replace = TRUE, prob = c(1 - p_0, p_0)), nrow = dim(X)[1], ncol = dim(X)[2]);
X[which(thetat == 1, arr.ind = TRUE)] <- 0;
# Remove those features with all zero counts
if (length(which(colSums(X) == 0)) > 0) {
  gt <- gt[-which(colSums(X) == 0)];
  gammat <- gammat[-which(colSums(X) == 0)];
  thetat <- thetat[, -which(colSums(X) == 0)];
  X <- X[, -which(colSums(X) == 0)];
}

### Save simulated data ----------------------------------------------------------------
save(X, zt, gammat, st, gt, file = "Data_Files/simulated_data.RData");
