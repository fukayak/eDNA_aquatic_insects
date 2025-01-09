plot_fig1 <- function(phi, theta, psi) {

  if (!require("RColorBrewer")) install.packages("RColorBrewer")

  phi   <- t(phi)
  theta <- t(theta)
  psi   <- t(psi)

  list_order <- c("Ephemeroptera", "Plecoptera",
                  "Trichoptera", "Odonata",
                  "Diptera", "Coleoptera", "others")

  ord <- set_order(psi[, 1], as.character(data_raw$order), list_order)

  layout(
    matrix(
      c(rep(1, 4), rep(2, 1), rep(3, 4), rep(4, 1), rep(5, 4), rep(6, 1)),
      5, 3)
  )
  mar1 <- c(0, 2.5, 4, 1)
  mar2 <- c(5, 2.5, 1, 1)

  colsx <- RColorBrewer::brewer.pal(8, "Set1")[c(2, 5, 1, 3, 4, 8, 7)]
  cols <- rep(colsx[7], nrow(psi))
  cols[data_raw$order == "Ephemeroptera"] <- colsx[1]
  cols[data_raw$order == "Plecoptera"]    <- colsx[2]
  cols[data_raw$order == "Trichoptera"]   <- colsx[3]
  cols[data_raw$order == "Odonata"]       <- colsx[4]
  cols[data_raw$order == "Diptera"]       <- colsx[5]
  cols[data_raw$order == "Coleoptera"]    <- colsx[6]

  # psi
  par(mar = mar1)
  plot(0, type = "n", xlim = c(0, 1), ylim = c(1, nrow(psi)),
       xaxt = "n", yaxt = "n", ylab = "", xlab = "",
       main = "(a) Site occupancy probability")
  arrows(y0 = 1:nrow(psi), x0 = psi[ord, 2], x1 = psi[ord, 3],
         length = 0, col = cols[ord])
  points(y = 1:nrow(psi), x = psi[ord, 1], pch = 21,
         bg = cols[ord], cex = 1.5)
  axis(side = 2, at = seq(0, nrow(psi), 10), labels = NA, las = 2)
  axis(side = 1, at = seq(0, 1, by = 0.2), labels = NA)

  par(mar = mar2)
  hist(psi[, 1], col = rgb(0.7, 0.7, 0.7, 0.5), xlim = c(0, 1),
       breaks = seq(0, 1, 0.1), yaxt = "n", xlab = "", ylab = "",
       main = "", cex.axis = 1.2)
  axis(side = 2, at = seq(0, 120, by = 20), las = 2, cex.axis = 1.2)
  box()

  # theta
  par(mar = mar1)
  plot(0, type = "n", xlim = c(0, 1), ylim = c(1, nrow(theta)),
       xaxt = "n", yaxt = "n", ylab = "", xlab = "",
       main = "(b) Sequence capture probability")
  arrows(y0 = 1:nrow(theta), x0 = theta[ord, 2], x1 = theta[ord, 3],
         length = 0, col = cols[ord])
  points(y = 1:nrow(theta), x = theta[ord, 1], pch = 21,
         bg = cols[ord], cex = 1.5)
  axis(side = 2, at = seq(0, nrow(theta), 10), labels = NA, las = 2)
  axis(side = 1, at = seq(0, 1, by = 0.2), labels = NA)

  par(mar = mar2)
  hist(theta[, 1], col = rgb(0.7, 0.7, 0.7, 0.5), xlim = c(0, 1),
       breaks = seq(0, 1, 0.1),
       yaxt = "n", xlab = "", ylab = "", main = "", cex.axis = 1.2)
  axis(side = 2, at = seq(0, 120, by = 20), las = 2, cex.axis = 1.2)
  box()

  # phi
  par(mar = mar1)
  plot(0, type = "n", xlim = range(c(0, phi)), ylim = c(1, nrow(phi)),
       xaxt = "n", yaxt = "n", ylab = "", xlab = "",
       main = "(c) Sequence relative dominance")
  arrows(y0 = 1:nrow(phi), x0 = phi[ord, 2], x1 = phi[ord, 3],
         length = 0, col = cols[ord])
  points(y = 1:nrow(phi), x = phi[ord, 1], pch = 21,
         bg = cols[ord], cex = 1.5)
  axis(side = 2, at = seq(0, nrow(phi), 10), labels = NA, las = 2)
  axis(side = 1, at = seq(0, 2.5, by = 0.5), labels = NA)
  legend("bottomright", "topright", legend = list_order, title = "Order",
         bty = "n", pch = 19, col = colsx)

  par(mar = mar2)
  hist(phi[, 1], col = rgb(0.7, 0.7, 0.7, 0.5), xlim = range(c(0, phi)),
       breaks = seq(0, 2.1, 0.1),
       yaxt = "n", xlab = "", ylab = "", main = "", cex.axis = 1.2)
  axis(side = 2, at = seq(0, 120, by = 20), las = 2, cex.axis = 1.2)
  box()
}

set_order <- function(param, order, list_order) {

  out <- vector(length = length(param))
  N_last <- 0

  for (i in seq_along(list_order)) {
    is_selected <- order == list_order[i]
    ord <- order(param[is_selected], decreasing = TRUE)
    N_first <- N_last + 1
    N_last  <- N_last + sum(is_selected)
    out[N_first:N_last] <- which(is_selected)[ord]
  }

  rev(out)
}

get_cor_mat <- function(phi, theta, psi) {

  log_phi     <- log(phi)
  logit_theta <- qlogis(theta)
  logit_psi   <- qlogis(psi)

  post_cor_mat <- array(dim = c(nrow(phi), 3, 3))
  for (i in 1:nrow(phi)) {
    post_cor_mat[i, 2, 1] <- cor(logit_theta[i, ], log_phi[i, ])
    post_cor_mat[i, 3, 1] <- cor(logit_psi[i, ],   log_phi[i, ])
    post_cor_mat[i, 3, 2] <- cor(logit_psi[i, ],   logit_theta[i, ])
  }

  out <- array(dim = c(3, 3, 3))
  dimnames(out) <- list(c("log(phi)", "logit(theta)", "logit(psi)"),
                        c("log(phi)", "logit(theta)", "logit(psi)"),
                        c("Posterior median", "95% CI Lower", "95% CI Upper"))

  out[, , 1] <- apply(post_cor_mat, c(2, 3), median)
  out[, , 2] <- apply(post_cor_mat, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
  out[, , 3] <- apply(post_cor_mat, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)

  out
}

plot_fig2 <- function(result1000, result800, result600, result400) {

  if (!require("ggplot2")) install.packages("ggplot2")
  library(ggplot2)

  result1000$vol <- 1000
  result800$vol  <- 800
  result600$vol  <- 600
  result400$vol  <- 400

  cols <- c("#0080FF", "#3399FF", "#66B3FF", "#FFB366", "#FF9933", "#FF8000")

  d <- rbind(result1000, result800, result600, result400)
  d$K <- as.factor(d$K); d$vol <- as.factor(d$vol)
  colnames(d) <- c("replicate", "reads", "detection", "vol")
  g <- ggplot(d, aes(x = reads, y = detection, colour = replicate))
  g + geom_line() + geom_point() + theme_light() +
      scale_colour_manual(values = cols,
                          name = "Number of replicates") +
      facet_grid(. ~ vol) + ylim(0, max(d$detection)) +
      ylab("Expected number of species detected") +
      xlab("Sequencing depth")
}

plot_fig3 <- function(result1000, result800, result600, result400) {

  if (!require("RColorBrewer")) install.packages("RColorBrewer")
  if (!require("ggplot2")) install.packages("ggplot2")
  library(RColorBrewer)
  library(ggplot2)

  result1000$vol <- 1000
  result800$vol  <- 800
  result600$vol  <- 600
  result400$vol  <- 400

  d <- rbind(result1000, result800, result600, result400)
  d$vol <- factor(d$vol, levels = c("1000", "800", "600", "400"))
  colnames(d) <- c("B", "lambda1", "lambda2", "replicate", "reads",
                   "detection", "vol")
  g <- ggplot(d, aes(x = replicate, y = detection, colour = vol))
  g + geom_line() + geom_point() + theme_light() +
      scale_colour_manual(values = brewer.pal(8, "Set2")[c(1:4)],
                          name = "Filtered volume (ml)") +
      ylim(0, max(d$detection)) +
      ylab("Expected number of species detected") +
      xlab("Number of replicates")
}

