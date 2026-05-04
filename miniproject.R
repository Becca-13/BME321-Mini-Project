# load libraries¨
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(pROC)
library(janitor)
library(stringr)

# load file
setwd("~/Desktop/BME321/miniproject")

data <- read_csv2("Cohort11.csv", skip = 1)

data <- clean_names(data)

head(data)
colnames(data)

# create responder column
data <- data %>%
  mutate(
    responder = case_when(
      best_response %in% c("CR", "PR") ~ 1,
      best_response %in% c("SD", "PD") ~ 0,
      TRUE ~ NA_real_
    )
  )

table(data$responder)

# create B-cell and T-cell abundance variables
data <- data %>%
  mutate(
    flow_lymphocytes_count = as.numeric(flow_lymphocytes_count),
    flow_b_cells_count = as.numeric(flow_b_cells_count),
    flow_t_cd4_count = as.numeric(flow_t_cd4_count),
    flow_t_cd8_count = as.numeric(flow_t_cd8_count),
    
    flow_t_cell_count = flow_t_cd4_count + flow_t_cd8_count,
    
    b_cell_percent =
      flow_b_cells_count / flow_lymphocytes_count * 100,
    
    t_cell_percent =
      flow_t_cell_count / flow_lymphocytes_count * 100
  )

# keep only rows for ROC and no NA
data_clean <- data %>%
  filter(
    !is.na(responder),
    !is.na(b_cell_percent),
    !is.na(t_cell_percent)
  )

dim(data_clean)
table(data_clean$responder)

# B-cell distribution by response
ggplot(data_clean,
       aes(x = factor(responder),
           y = b_cell_percent)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(
    x = "Responder",
    y = "B-cell percentage",
    title = "B-cell abundance and ICB response"
  )

# ROC curve for B cells
roc_b <- roc(
  data_clean$responder,
  data_clean$b_cell_percent
)

# ROC curve for T cells
roc_t <- roc(
  data_clean$responder,
  data_clean$t_cell_percent
)

auc(roc_b)
auc(roc_t)

# plot
par(mar = c(5,5,4,2))

plot(
  roc_b,
  legacy.axes = TRUE,
  col = "lightblue",
  lwd = 2,
  main = "ROC curves — ICB response prediction",
  xlim = c(0,1),
  ylim = c(0,1)
)

plot(
  roc_t,
  legacy.axes = TRUE,
  col = "red",
  lwd = 3,
  add = TRUE
)

abline(a = 0, b = 1, lty = 2, col = "gray")

legend(
  "bottomright",
  inset = 0.03, 
  legend = c(
    paste0("B cells (AUC = ", round(auc(roc_b), 3), ")"),
    paste0("T cells (AUC = ", round(auc(roc_t), 3), ")"),
    "Random (AUC = 0.5)"
  ),
  col = c("lightblue", "red", "black"),
  lwd = c(3, 3, 1),
  lty = c(1, 1, 2),
  bty = "n",
  cex = 0.8
)

# add TNM stage --> 1, 2+3, 4
data_clean <- data_clean %>%
  mutate(
    tnm_group = case_when(
      tnm_stage %in% c("I") ~ "1",
      
      tnm_stage %in% c(
        "II",
        "IIA",
        "III"
      ) ~ "2_3",
      
      tnm_stage %in% c(
        "IV",
        "IVA",
        "IVB",
        "IVC"
      ) ~ "4",
      
      TRUE ~ NA_character_
    )
  )

data_clean$tnm_group <- as.factor(data_clean$tnm_group)

# new dataset for tnm
data_tnm <- data_clean %>%
  filter(
    !is.na(responder),
    !is.na(b_cell_percent),
    !is.na(tnm_group)
  )

table(data_tnm$responder)
table(data_tnm$tnm_group)

# make a regression model 
model_b <- glm(
  responder ~ b_cell_percent,
  data = data_tnm,
  family = binomial
)

model_b_tnm <- glm(
  responder ~ b_cell_percent + tnm_group,
  data = data_tnm,
  family = binomial
)

# predictions
pred_b <- predict(model_b, type = "response")
pred_b_tnm <- predict(model_b_tnm, type = "response")

roc_b <- roc(data_tnm$responder, pred_b)
roc_b_tnm <- roc(data_tnm$responder, pred_b_tnm)


# roc
par(mar = c(5,5,4,2))

plot(
  roc_b,
  legacy.axes = TRUE,
  col = "lightblue",
  lwd = 3,
  main = "ROC curves — ICB response prediction (Cohort 11)",
  xlim = c(0,1),
  ylim = c(0,1)
)

plot(
  roc_b_tnm,
  legacy.axes = TRUE,
  col = "red",
  lwd = 3,
  add = TRUE
)

abline(a = 0, b = 1, lty = 2, col = "gray")

legend(
  "topleft",
  legend = c(
    paste0("B cells (AUC = ", round(auc(roc_b), 3), ")"),
    paste0("B cells + TNM (AUC = ", round(auc(roc_b_tnm), 3), ")"),
    "Random (AUC = 0.5)"
  ),
  col = c("lightblue", "red", "gray"),
  lwd = c(3, 3, 1),
  lty = c(1, 1, 2),
  bty = "n"
)

