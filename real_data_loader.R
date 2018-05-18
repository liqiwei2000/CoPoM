### Set working directory --------------------------------------------------------------
setwd("~/Dropbox/Qiwei");

### Load data --------------------------------------------------------------------------
load("Data_Files/real_data.RData");

### Print the summary of real data ------------------------------------------------------
print(paste0("Number of observations: ", dim(X)[1]));
print(paste0("Number of features: ", dim(X)[2]));
print(paste0("Percentage of zeros: ", round(sum(X == 0)/dim(X)[2]/dim(X)[1]*100, 3)));
