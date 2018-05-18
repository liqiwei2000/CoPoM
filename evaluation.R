# Visualize the results -----------------------------------------------------------------------
par(mfrow = c(2, 2));
plot(c(start_K, K_store[1:iter]), type = "l", xlab = "Iterations", ylab = "Number of clusters", ylim = c(0, max(start_K, K_store)), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5);
hist(K_store[burn:iter], xlab = "Number of clusters", main = "", breaks = c(0:max(K_store)), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5);
plot(gamma_sum_store[1:iter], type = "l", xlab = "Iterations", ylab = "Number of inc. variables", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5);
plot(mpv, ylim = c(0, 1), type = "h", xlab = "Variable index", ylab = "Frequency", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5);
par(mfrow = c(1, 1));

levelplot(ppm, col.regions = gray(0:100/100)[100:1], xlab = "", ylab = "", , cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5);

par(mfrow = c(1, 2));
par(pty = "s");
plot(st, s_store[iter,], xlab = "True", ylab = "Estimate", main = "S", pch = 19, cex = 0.5, xlim = c(0, max(c(st, s_store[iter,]))), ylim = c(0, max(c(st, s_store[iter,]))));
abline(a = 0, b = 1, col = 2);
plot(gt, g_store[iter,], xlab = "True", ylab = "Estimate", main = "G", pch = 19, cex = 0.5, xlim = c(0, max(c(gt, g_store[iter,]))), ylim = c(0, max(c(gt, g_store[iter,]))));
abline(a = 0, b = 1, col = 2);
par(pty = "m");
par(mfrow = c(1, 1));

par(mfrow = c(1, 2));
plot(omega_store[1:iter], type = "l", xlab = "Iterations", ylab = "omega");
plot(alpha_store[1:iter], type = "l", xlab = "Iterations", ylab = "alpha");
par(mfrow = c(1, 1));

tabulate_error(gammat, gamma_ppi);

tabulate_error(gammat, gamma_map);

ari(zt, z_ppm);

ari(zt, z_map);
