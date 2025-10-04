# ==============================================================================
# Prognostic Utility of sFlt-1 and PlGF Biomarkers for Predicting Preeclampsia
# Secondary Analysis of INSPIRE Trial Data
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup and Data Loading
# ------------------------------------------------------------------------------

# Install and load required packages
# install.packages('pacman')  # Run this line if pacman is not installed
library(pacman)
pacman::p_load(
  readxl,      # For reading Excel files
  tidyverse,   # For data manipulation and visualization
  rms,         # Regression modeling strategies
  reshape2,    # Data reshaping
  mfp,         # Multivariable fractional polynomials
  plyr,        # Data manipulation
  Hmisc,       # Harrell miscellaneous functions
  MASS,        # Modern applied statistics
  pROC,        # ROC curve analysis
  DescTools,   # Descriptive statistics tools
  ggrepel,     # Better label positioning
  gridExtra,   # Arranging multiple plots
  boot,        # Bootstrap methods
  janitor      # Data cleaning
)

# Load the dataset
inspire_data <- read_excel("Inspire.xls")
attach(inspire_data)

# ------------------------------------------------------------------------------
# 2. Descriptive Statistics
# ------------------------------------------------------------------------------

# Calculate summary statistics by preeclampsia outcome (PE within 7 days)
descriptive_stats <- inspire_data %>% 
  group_by(PE = pet_after_visit1_within7days) %>% 
  summarise(
    across(
      c(age_at_recruitment, gestweeks_recruitment, bmi, 
        logsflt1vi, logpigfv1, log10_sFlt1_PIGF_atvisit1),
      list(
        mean = mean,
        sd = sd,
        median = median,
        quart25 = ~quantile(., 0.25),
        quart75 = ~quantile(., 0.75)
      )
    )
  )

print(descriptive_stats)

# ------------------------------------------------------------------------------
# 3. Visual Exploration of Biomarkers by PE Outcome
# ------------------------------------------------------------------------------

# Create three side-by-side boxplots comparing biomarker distributions
par(mfrow = c(1, 3))

# sFlt-1 levels by PE outcome
boxplot(logsflt1vi ~ pet_after_visit1_within7days, 
        xlab = "Preeclampsia", 
        ylab = "sFlt-1 values (log pg/mL)",
        main = "sFlt-1 by PE Status",
        col = c("lightblue", "coral"))

# PlGF levels by PE outcome
boxplot(logpigfv1 ~ pet_after_visit1_within7days, 
        xlab = "Preeclampsia", 
        ylab = "PlGF values (log pg/mL)",
        main = "PlGF by PE Status",
        col = c("lightgreen", "coral"))

# sFlt-1/PlGF ratio by PE outcome
boxplot(log10_sFlt1_PIGF_atvisit1 ~ pet_after_visit1_within7days, 
        xlab = "Preeclampsia", 
        ylab = "sFlt-1/PlGF ratio (log10)",
        main = "sFlt-1/PlGF Ratio by PE Status",
        col = c("lightyellow", "coral"))

# Reset plot layout
par(mfrow = c(1, 1))

# ------------------------------------------------------------------------------
# 4. Cross-tabulation: sFlt-1/PlGF Cut-off vs PE Outcome
# ------------------------------------------------------------------------------

# Examine how many patients fall into each ratio category by PE status
ratio_table <- inspire_data %>% 
  janitor::tabyl(ratio_class_atvisit1, pet_after_visit1_within7days)

print(ratio_table)

# ------------------------------------------------------------------------------
# 5. Custom Function for Model Fitting, Validation, and Visualization
# ------------------------------------------------------------------------------

fit_validate_plot <- function(model, predictor_var, model_title) {
  
  # Print model summary
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Model:", model_title, "\n")
  cat(rep("=", 70), "\n\n", sep = "")
  print(summary(model))
  
  # Generate predictions
  prob <- predict(model, type = "response")    # Predicted probabilities
  lin_pred <- predict(model, type = "link")    # Linear predictor
  
  # Plot observed values vs predicted probabilities
  plot(prob ~ get(predictor_var), 
       xlab = predictor_var, 
       ylab = "Predicted Probability of PE",
       main = paste("Predicted Probabilities vs", predictor_var),
       pch = 20, col = rgb(0, 0, 1, 0.3))
  
  # Model fit statistics
  cat("\nNagelkerke R-squared:", PseudoR2(model, "Nagelkerke"), "\n")
  cat("BIC:", BIC(model), "\n\n")
  
  # Calibration check: intercept-only model with offset
  calibration_intercept <- glm(pet_after_visit1_within7days2 ~ offset(lin_pred), 
                               family = "binomial")
  cat("Calibration intercept (should be ~0):\n")
  print(calibration_intercept$coefficients)
  
  # Calibration slope check (should be ~1 for perfect calibration)
  calibration_slope <- glm(pet_after_visit1_within7days2 ~ lin_pred, 
                          family = "binomial", x = TRUE, y = TRUE)
  cat("\nCalibration slope (should be ~1):\n")
  print(calibration_slope$coef)
  
  # ROC curve analysis
  cat("\nROC Analysis:\n")
  roc_curve <- roc(pet_after_visit1_within7days2 ~ prob, 
                   ci = TRUE, 
                   plot = TRUE, 
                   legacy.axes = TRUE, 
                   print.auc = TRUE,
                   main = paste("ROC Curve -", model_title))
  print(roc_curve)
  
  # Calibration plot with 10 risk groups
  groups <- cut(prob, 
                breaks = quantile(prob, prob = seq(0, 1, 0.1)), 
                labels = 1:10, 
                include.lowest = TRUE)
  
  gpdata <- cbind(inspire_data, groups, prob)
  
  # Calculate observed and expected proportions by group
  obs <- ddply(gpdata, ~groups, 
               summarise, 
               mean = mean(as.numeric(pet_after_visit1_within7days2)))[, 2]
  exp <- ddply(gpdata, ~groups, 
               summarise, 
               mean = mean(prob))
  obsn <- table(pet_after_visit1_within7days2, groups)[1, ]
  
  # Calculate 95% confidence intervals
  lci <- pmax(0, obs - (1.96 * sqrt((obs * (1 - obs)) / obsn)))
  uci <- pmin(1, obs + (1.96 * sqrt((obs * (1 - obs)) / obsn)))
  
  # Create calibration plot
  par(pty = "s")  # Square plot
  plot(obs ~ exp[, 2], 
       xlim = c(0, 1), ylim = c(0, 1), 
       col = "red", 
       pch = 16,
       ylab = "Observed Probability", 
       xlab = "Predicted Probability",
       main = paste("Calibration Plot -", model_title))
  
  # Add reference line (perfect calibration)
  lines(c(0, 1), c(0, 1), lty = 2, lwd = 2)
  
  # Add confidence intervals for each risk group
  for (i in 1:10) {
    lines(c(exp[i, 2], exp[i, 2]), c(lci[i], uci[i]), col = "green", lwd = 2)
  }
  
  # Add histogram showing distribution of predicted probabilities
  h <- hist(prob, breaks = 50, plot = FALSE)
  for (i in seq_along(h$mids)) {
    lines(c(h$mids[i], h$mids[i]), 
          c(1, 1 - ((h$counts[i] / max(h$counts)) / 10)),
          col = "gray", lwd = 0.5)
  }
  
  # Add smoothed calibration curve using loess
  obs_all <- predict(loess(pet_after_visit1_within7days2 ~ prob, span = 1))
  lines_data <- data.frame(prob, obs_all)
  lines_data_sorted <- lines_data[order(prob), ]
  lines(lines_data_sorted[, 1], lines_data_sorted[, 2], col = "blue", lwd = 2)
  
  # Add legend
  legend("topleft", 
         c("Risk groups", "Perfect calibration", "95% CI", "Loess smooth"),
         col = c("red", "black", "green", "blue"), 
         lty = c(0, 2, 1, 1), 
         pch = c(16, NA, NA, NA), 
         bty = "n")
  
  cat("\n", rep("=", 70), "\n\n", sep = "")
}

# ------------------------------------------------------------------------------
# 6. Fit and Validate Logistic Regression Models
# ------------------------------------------------------------------------------

# Model 1: sFlt-1 as predictor
cat("\n### Model 1: sFlt-1 ###\n")
sflt_model <- glm(pet_after_visit1_within7days2 ~ factor(trial_arm) + logsflt1vi, 
                  family = "binomial")
fit_validate_plot(sflt_model, "logsflt1vi", "sFlt-1")

# Model 2: PlGF as predictor
cat("\n### Model 2: PlGF ###\n")
plgf_model <- glm(pet_after_visit1_within7days2 ~ factor(trial_arm) + logpigfv1, 
                  family = "binomial")
fit_validate_plot(plgf_model, "logpigfv1", "PlGF")

# Model 3: sFlt-1/PlGF ratio as predictor
cat("\n### Model 3: sFlt-1/PlGF Ratio ###\n")
ratio_model <- glm(pet_after_visit1_within7days2 ~ factor(trial_arm) + log10_sFlt1_PIGF_atvisit1, 
                   family = "binomial")
fit_validate_plot(ratio_model, "log10_sFlt1_PIGF_atvisit1", "sFlt-1/PlGF Ratio")

# Model 4: Categorical cutoff approach
cat("\n### Model 4: sFlt-1/PlGF Cutoff ###\n")
cutoff_model <- glm(pet_after_visit1_within7days2 ~ factor(trial_arm) + ratio_class_atvisit1, 
                    family = "binomial")
fit_validate_plot(cutoff_model, "ratio_class_atvisit1", "sFlt-1/PlGF Cutoff")

# ------------------------------------------------------------------------------
# 7. Statistical Comparison of Model Performance (DeLong Test)
# ------------------------------------------------------------------------------

# Generate ROC curves for all models
roc_sflt <- roc(pet_after_visit1_within7days2 ~ predict(sflt_model, type = "response"))
roc_plgf <- roc(pet_after_visit1_within7days2 ~ predict(plgf_model, type = "response"))
roc_ratio <- roc(pet_after_visit1_within7days2 ~ predict(ratio_model, type = "response"))
roc_cutoff <- roc(pet_after_visit1_within7days2 ~ predict(cutoff_model, type = "response"))

# Perform pairwise DeLong tests to compare AUCs
cat("\n", rep("=", 70), "\n", sep = "")
cat("Pairwise Comparisons of Model AUCs (DeLong Test)\n")
cat(rep("=", 70), "\n\n", sep = "")

roc_comparisons <- list(
  SFLT_vs_PLGF = roc.test(roc_sflt, roc_plgf, method = "delong"),
  SFLT_vs_RATIO = roc.test(roc_sflt, roc_ratio, method = "delong"),
  SFLT_vs_CUTOFF = roc.test(roc_sflt, roc_cutoff, method = "delong"),
  PLGF_vs_RATIO = roc.test(roc_plgf, roc_ratio, method = "delong"),
  PLGF_vs_CUTOFF = roc.test(roc_plgf, roc_cutoff, method = "delong"),
  RATIO_vs_CUTOFF = roc.test(roc_ratio, roc_cutoff, method = "delong")
)

# Print all comparison results
for (comparison_name in names(roc_comparisons)) {
  cat("\n", comparison_name, ":\n", sep = "")
  print(roc_comparisons[[comparison_name]])
  cat("\n")
}

# ------------------------------------------------------------------------------
# 8. Summary Visualization: Compare All ROC Curves
# ------------------------------------------------------------------------------

plot(roc_sflt, col = "red", lwd = 2, main = "Comparison of All Models")
plot(roc_plgf, col = "blue", add = TRUE, lwd = 2)
plot(roc_ratio, col = "green", add = TRUE, lwd = 2)
plot(roc_cutoff, col = "purple", add = TRUE, lwd = 2)

legend("bottomright", 
       legend = c(
         paste("sFlt-1 (AUC:", round(auc(roc_sflt), 3), ")"),
         paste("PlGF (AUC:", round(auc(roc_plgf), 3), ")"),
         paste("Ratio (AUC:", round(auc(roc_ratio), 3), ")"),
         paste("Cutoff (AUC:", round(auc(roc_cutoff), 3), ")")
       ),
       col = c("red", "blue", "green", "purple"),
       lwd = 2,
       bty = "n")

# ------------------------------------------------------------------------------
# End of Analysis
# ------------------------------------------------------------------------------
