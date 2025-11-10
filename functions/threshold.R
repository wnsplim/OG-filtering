function(og, maxthres, diff, outgroup = NULL) {
  
  if (maxthres <= 0 || maxthres != floor(maxthres)) {
    stop("Error: maxthres must be a positive integer.")
  }
  
  mean <- numeric()
  sd <- numeric()
  cv <- numeric()
  mincov <- numeric()
  
  pb <- txtProgressBar(min = diff, max = max(seq(diff, maxthres, by = diff)), style = 3)
  
  for (i in seq(diff, maxthres, by = diff)) {
    temp <- og_filter(og = og, threshold = i, outgroup = outgroup, quiet = TRUE, statistics = TRUE)
    mean <- c(mean, temp[[2]])
    sd <- c(sd, temp[[3]])
    cv <- c(cv, temp[[4]])
    mincov <- c(mincov, temp[[5]])
    setTxtProgressBar(pb, i)
  }
  cat("\n")
  
  
  min_ind_sd <- sort(union(which(cv %in% cv[which(diff(sign(diff(sd))) == 2) + 1]), which(sd == min(sd))))
  min_ind_cv <- sort(union(which(cv %in% cv[which(diff(sign(diff(cv))) == 2) + 1]), which(cv == min(cv))))
  max_ind_mincov <- sort(union(which(cv %in% cv[which(diff(sign(diff(-mincov))) == 2) + 1]), which(mincov == max(mincov))))
  
  plot(mean, type = "l", col = "black", xlab = "Threshold", xaxt = "n", ylab = "", ylim = c(min(c(min(mean),min(mincov)))*0.95, max(c(max(mean),max(mincov)))*1.05))
  x_labels <- seq(from = 0, to = max(seq(diff, maxthres, by = diff)), by = floor(max(seq(diff, maxthres, by = diff)))/10)
  axis(1, at = x_labels/diff, labels = x_labels)
  
  lines(mincov, col = "blue")
  points(max_ind_mincov, mincov[max_ind_mincov], col = "blue", pch = 15)
  
  par(new = TRUE)
  
  plot(cv, type = "l", col = "orange", axes = FALSE, xlab = "", ylab = "", xaxt = "n", ylim = c(min(c(min(cv),min(sd)))*0.95, max(c(max(cv),max(sd)))*1.05))
  lines(sd, col = "red")
  
  axis(4)
  
  points(min_ind_sd, sd[min_ind_sd], col = "red", pch = 15)
  points(min_ind_cv, cv[min_ind_cv], col = "orange", pch = 15)
  
  legend("topright", legend = c("mean", "mincov", "sd", "CV"), col = c("black", "blue", "red", "orange"), pch = 15, bty = "n")
  
  abline(v = x_labels/diff, col = "gray", lty = 2)
  
  plot <- recordPlot()
  cat("Local minima (sd):", paste((min_ind_sd) * diff, collapse = ", "), "\n")
  cat("Local minima (CV):", paste((min_ind_cv) * diff, collapse = ", "), "\n")
  cat("Local Maxima (minimum coverage):", paste((max_ind_mincov) * diff, collapse = ", "), "\n")
  
  return(list(plot = plot))
}
