# =============================================================================
# Cohort 11 — Kaplan-Meier Overall Survival by B-cell Threshold
# Reproducing Figure 4E: Chang et al. (2025), Annals of Oncology
# DOI: https://doi.org/10.1016/j.annonc.2024.09.014
#
# Stratification: B-cell abundance ≤ 5.5% vs > 5.5% of PBMCs
# Outcome:        Overall survival (OS)
# Statistics:     Log-rank test + Cox proportional hazards HR
# =============================================================================

library(readxl)
library(survival)
library(survminer)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

df <- read_excel(
  "Data_of_in-house_cohorts.xlsx",
  sheet = "Cohort11",
  skip  = 1        # skip merged top-header row
)

# Inspect column names — adjust if your read produces different strings
head(df)
colnames(df)

# -----------------------------------------------------------------------------
# 2. Compute B-cell abundance
#    Definition: B cells / (Lymphocytes + Monocytes) * 100
# -----------------------------------------------------------------------------

df$B_cell_abundance <- df[["Flow B Cells Count"]] /
  (df[["Flow Lymphocytes Count"]] + df[["Flow Monocytes Count"]]) * 100

# -----------------------------------------------------------------------------
# 3. Stratify by 5.5% threshold (external cutoff from paper, derived Cohort 7)
# -----------------------------------------------------------------------------

df$B_group <- ifelse(df$B_cell_abundance > 5.5, "B > 5.5%", "B \u2264 5.5%")
df$B_group <- factor(df$B_group, levels = c("B \u2264 5.5%", "B > 5.5%"))

# Verify group sizes — expected: 26 (≤5.5%) and 34 (>5.5%)
table(df$B_group)

# -----------------------------------------------------------------------------
# 4. Build survival object
# -----------------------------------------------------------------------------

# OS status: 1 = event (death), 0 = censored
surv_obj <- Surv(
  time  = df[["OS months"]],
  event = df[["OS status"]]
)

# -----------------------------------------------------------------------------
# 5. Log-rank test
# -----------------------------------------------------------------------------

logrank_fit <- survdiff(surv_obj ~ B_group, data = df)
p_val       <- 1 - pchisq(logrank_fit$chisq, df = 1)
cat("Log-rank p-value:", round(p_val, 3), "\n")
# Expected: ~0.030

# -----------------------------------------------------------------------------
# 6. Cox proportional hazards — HR with 95% CI
# -----------------------------------------------------------------------------

cox_fit <- coxph(surv_obj ~ B_group, data = df)
summary(cox_fit)

HR    <- exp(coef(cox_fit))
HR_CI <- exp(confint(cox_fit))
cat(sprintf("HR = %.2f (%.2f\u2013%.2f)\n", HR, HR_CI[1], HR_CI[2]))
# Expected: HR = 0.54 (0.31–0.94)

# -----------------------------------------------------------------------------
# 7. Kaplan-Meier fit
# -----------------------------------------------------------------------------

km_fit <- survfit(surv_obj ~ B_group, data = df)
summary(km_fit, times = c(0, 10, 20, 30, 40, 50))  # Number at risk check

# -----------------------------------------------------------------------------
# 8. Plot — reproducing Figure 4E
# -----------------------------------------------------------------------------

# Annotation string for HR and p-value
hr_label <- sprintf(
  "HR = %.2f (%.2f\u2013%.2f)\nP = %.3f",
  HR, HR_CI[1], HR_CI[2], p_val
)

km_plot <- ggsurvplot(
  km_fit,
  data          = df,
  
  # Lines
  palette       = c("#1C3E8A", "#C0392B"),  # blue = ≤5.5%, red = >5.5%
  size          = 0.8,
  censor.shape  = c("|", "+"),              # match paper symbols
  censor.size   = 4,
  
  # Axes
  xlab          = "Time (months)",
  ylab          = "OS probability",
  xlim          = c(0, 54),
  ylim          = c(0, 1),
  break.x.by    = 10,
  break.y.by    = 0.25,
  
  # Legend
  legend        = c(0.72, 0.88),
  legend.title  = "",
  legend.labs   = c("B \u2264 5.5%", "B > 5.5%"),
  
  # Number at risk table
  risk.table         = TRUE,
  risk.table.col     = "strata",
  risk.table.height  = 0.22,
  risk.table.y.text  = FALSE,    # show colored lines instead of text labels
  tables.theme       = theme_cleantable(),
  
  # Panel label
  title         = "",
  
  # Theme
  ggtheme       = theme_classic(base_size = 12) +
    theme(
      axis.line   = element_line(linewidth = 0.5),
      plot.margin = margin(10, 10, 5, 10)
    ),
  
  # Confidence intervals (set TRUE to show, FALSE to hide as in paper)
  conf.int      = FALSE,
  
  # p-value (use custom annotation below instead for exact formatting)
  pval          = FALSE
)

# Add HR and p-value annotation (matching paper position)
km_plot$plot <- km_plot$plot +
  annotate(
    "text",
    x     = 1.5,
    y     = 0.32,
    label = hr_label,
    hjust = 0,
    size  = 3.5
  ) +
  # Panel label "E"
  annotate(
    "text",
    x      = -1,
    y      = 1.04,
    label  = "E",
    fontface = "bold",
    size   = 5,
    hjust  = 1
  )

# -----------------------------------------------------------------------------
# 9. Save
# -----------------------------------------------------------------------------

# Display in RStudio viewer
print(km_plot)

# Save as PNG (adjust width/height/dpi as needed)
ggsave(
  filename = "cohort11_KM_OS.png",
  plot     = print(km_plot),
  width    = 5.5,
  height   = 5.5,
  dpi      = 300
)

# =============================================================================
# Session info
# =============================================================================
sessionInfo()
