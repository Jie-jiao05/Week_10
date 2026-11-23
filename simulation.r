## set simulation parameters
n <- 500
nsim <- 1000
significance.level <- 0.05
pz <- 0.2
alpha0 <- 0
alpha1.values <- c(0, 1, 2)
beta0 <- -3
beta1 <- 0
beta2 <- 2

set.seed(2025)

## store Type I error rates
results <- data.frame(
  alpha1 = alpha1.values,
  adjusted = NA,
  unadjusted = NA
)

## repeat for each alpha1 scenario
for (j in 1:length(alpha1.values)) {
  
  alpha1 <- alpha1.values[j]
  
  p.adj <- numeric(nsim)
  p.unadj <- numeric(nsim)
  
  ## repeat simulation 1000 times
  for (i in 1:nsim) {
    
    z <- rbinom(n, size = 1, prob = pz)
    px <- exp(alpha0 + alpha1 * z) / (1 + exp(alpha0 + alpha1 * z))
    x <- rbinom(n, size = 1, prob = px)
    py <- exp(beta0 + beta1 * x + beta2 * z) /(1 + exp(beta0 + beta1 * x + beta2 * z))
    y <- rbinom(n, size = 1, prob = py)

    ## combine three random variables into a data frame 
    dat <- data.frame(lung = y, coffee = x, smoke = z)

    ## fit unadjusted logistic regression mode
    unadj.mod <- glm(lung ~ coffee, data = dat, family = "binomial")
    p.unadj[i] <- summary(unadj.mod)$coef["coffee", "Pr(>|z|)"]

    ## fit adjusted logistic regression model
    adj.mod <- glm(lung ~ coffee + smoke, data = dat, family = "binomial")
    p.adj[i] <- summary(adj.mod)$coef["coffee", "Pr(>|z|)"]
  }
  
  results$unadjusted[j] <- mean(p.unadj < significance.level)
  results$adjusted[j] <- mean(p.adj < significance.level)
}

## create and save the results table
table_display <- rbind(
  c("", "Adjusted", "Unadjusted"),
  cbind(results$alpha1,
        sprintf("%.1f%%", 100 * results$adjusted),
        sprintf("%.1f%%", 100 * results$unadjusted))
)
png("simulation_results_table.png", width = 1200, height = 800, res = 150)
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.07, 0.92, "Table 1. Estimated Type I Error Rates",
     adj = 0, font = 2, cex = 1.8)
text(0.07, 0.85, sprintf("Sample size = %d  |  Simulations = %d  |  Significance level = %.2f",
                         n, nsim, significance.level),
     adj = 0, cex = 0.95, col = "gray30")
x <- seq(0.15, 0.82, length.out = ncol(table_display))
y <- seq(0.70, 0.24, length.out = nrow(table_display))
for (i in 1:nrow(table_display))
  text(x, y[i], table_display[i, ], font = ifelse(i == 1, 2, 1), cex = 1.1)
text(x[1], y[1], expression(alpha[1]), font = 2, cex = 1.1)
row_gap <- abs(y[2] - y[1])
line_y <- c(y[1] + row_gap / 2, y - row_gap / 2)
segments(0.07, line_y, 0.93, line_y, lwd = 1.4, col = "gray25")
text(0.07, min(line_y) - 0.08,
     "Note: The adjusted model controls for smoking; the unadjusted model does not.",
     adj = 0, cex = 0.85, col = "gray30")
dev.off()
