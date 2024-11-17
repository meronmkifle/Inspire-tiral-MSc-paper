# Article title 
# The prognostic utility of soluble fms-like tyrosine kinase-1 (sFlt-1) and placental growth factor (PlGF)
# biomarkers for predicting preeclampsia: A secondary analysis of data from the INSPIRE trial data   

# Libraries
# install.libraries('pacman')
library(pacman)
pacman::p_load(readxl, tidyverse, rms, reshape2, mfp, plyr, Hmisc, MASS, pROC, DescTools, ggrepel, gridExtra, boot, janitor)

# Load dataset
Inspire <- read_excel("Inspire.xls")
attach(Inspire)

# Descriptive analysis
Inspire %>% group_by(PE = pet_after_visit1_within7days) %>% 
  summarise(across(c(age_at_recruitment, gestweeks_recruitment, bmi, logsflt1vi, logpigfv1, log10_sFlt1_PIGF_atvisit1), 
                   list(mean = mean, sd = sd, median = median, quart25 = ~quantile(., 0.25), quart75 = ~quantile(., 0.75))))

# PE outcome by biomarker values (log transformed)
par(mfrow = c(1, 3))
boxplot(logsflt1vi ~ pet_after_visit1_within7days, xlab = "Preeclampsia", ylab = "sFlt values (pg/mL)")
boxplot(logpigfv1 ~ pet_after_visit1_within7days, xlab = "Preeclampsia", ylab = "PIGF values (pg/mL)")
boxplot(log10_sFlt1_PIGF_atvisit1 ~ pet_after_visit1_within7days, xlab = "Preeclampsia", ylab = "sflt/PIGF ratio values (pg/mL)")

# sflt-1/PIGF cut off proportion by PE within 7 days
Inspire %>% janitor::tabyl(ratio_class_atvisit1, pet_after_visit1_within7days)

# Function to fit, validate and plot logistic models
fit_validate_plot <- function(mod, pred_var, title) {
  summary(mod)
  
  prob <- predict(mod, type = "response")
  lin_pred <- predict(mod, type = "link")
  
  plot(prob ~ get(pred_var), xlab = pred_var, ylab = "Predicted probabilities of PE",
       main = paste("Observed", pred_var, "values Vs predicted probabilities"), pch = 20)
  
  print(PseudoR2(mod, "Nagelkerke"))
  print(BIC(mod))
  
  mod_log_1 <- glm(pet_after_visit1_within7days2 ~ offset(lin_pred), family = "binomial")
  print(mod_log_1$coefficients)
  
  mod_log_2 <- glm(pet_after_visit1_within7days2 ~ lin_pred, family = "binomial", x = TRUE, y = TRUE)
  print(mod_log_2$coef)
  
  c1 <- roc(pet_after_visit1_within7days2 ~ prob, ci = TRUE, plot = TRUE, legacy.axes = TRUE, print.auc = TRUE)
  print(c1)
  
  groups <- cut(prob, breaks = quantile(prob, prob = seq(0, 1, 0.1)), labels = 1:10, include.lowest = TRUE)
  gpdata <- cbind(Inspire, groups, prob)
  obs <- ddply(gpdata, ~groups, summarise, mean = mean(as.numeric(pet_after_visit1_within7days2)))[, 2]
  exp <- ddply(gpdata, ~groups, summarise, mean = mean(prob))
  obsn <- table(pet_after_visit1_within7days2, groups)[1,]
  
  lci <- pmax(0, (obs - (1.96 * (((obs * (1 - obs)) / obsn) ^ 0.5))))
  uci <- pmin(1, (obs + (1.96 * (((obs * (1 - obs)) / obsn) ^ 0.5))))
  
  par(pty = "s")
  plot(obs ~ exp[, 2], xlim = c(0, 1), ylim = c(0, 1), col = "red", ylab = "Observed", xlab = "Predicted")
  lines(c(0, 1), c(0, 1), lty = 2)
  for (i in 1:10) { lines(c(exp[i, 2], exp[i, 2]), c(lci[i], uci[i]), col = "green") }
  
  h <- hist(prob, breaks = 50, plot = FALSE)
  for (i in seq_along(h$mids)) {
    lines(c(h$mids[i], h$mids[i]), c(rep(1, length(h$mids))[i], 1 - ((h$counts[i] / max(h$counts)) / 10)))
  }
  
  obs_all <- predict(loess(pet_after_visit1_within7days2 ~ prob, span = 1))
  lines_data <- data.frame(prob, obs_all)
  lines_data2 <- lines_data[order(prob), ]
  lines(lines_data2[, 1], lines_data2[, 2], col = "blue")
  
  legend(0, 0.999, c("Risk groups", "Reference line", "95% CI", "Loess"),
         col = c("red", "black", "green", "blue"), lty = c(0, 2, 1, 1), pch = c(1, NA, NA, NA), bty = "n")
}

# Models
sflt_mod <- glm(pet_after_visit1_within7days2 ~ factor(trial_arm) + logsflt1vi, family = "binomial")
fit_validate_plot(sflt_mod, "logsflt1vi", "sFlt-1")

PIGF_mod <- glm(pet_after_visit1_within7days2 ~ factor(trial_arm) + logpigfv1, family = "binomial")
fit_validate_plot(PIGF_mod, "logpigfv1", "PIGF")

sflt_PIGF_mod <- glm(pet_after_visit1_within7days2 ~ factor(trial_arm) + log10_sFlt1_PIGF_atvisit1, family = "binomial")
fit_validate_plot(sflt_PIGF_mod, "log10_sFlt1_PIGF_atvisit1", "sflt-1/PIGF ratio")

cutoff_mod <- glm(pet_after_visit1_within7days2 ~ factor(trial_arm) + ratio_class_atvisit1, family = "binomial")
fit_validate_plot(cutoff_mod, "ratio_class_atvisit1", "sflt-1/PIGF cutoff")

# Delong test for statistical difference of models AUC curves
roc_sflt <- roc(pet_after_visit1_within7days2 ~ predict(sflt_mod, type = "response"))
roc_plgf <- roc(pet_after_visit1_within7days2 ~ predict(PIGF_mod, type = "response"))
roc_sflt_plgf_ratio <- roc(pet_after_visit1_within7days2 ~ predict(sflt_PIGF_mod, type = "response"))
roc_sflt_plgf_cutoff <- roc(pet_after_visit1_within7days2 ~ predict(cutoff_mod, type = "response"))

par(mfrow = c(2, 2))
roc_tests <- list(
  SFLT_PLGF = roc.test(roc_sflt, roc_plgf, method = "delong"),
  SFLT_RATIO = roc.test(roc_sflt, roc_sflt_plgf_ratio, method = "delong"),
  SFLT_CUTOFF = roc.test(roc_sflt, roc_sflt_plgf_cutoff, method = "delong"),
  PLGF_RATIO = roc.test(roc_plgf, roc_sflt_plgf_ratio, method = "delong"),
  PLGF_CUTOFF = roc.test(roc_plgf, roc_sflt_plgf_cutoff, method = "delong"),
  RATIO_CUTOFF = roc.test(roc_sflt_plgf_ratio, roc_sflt_plgf_cutoff, method = "delong")
)

roc_tests
