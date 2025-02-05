# load libraries
library(R.utils)
library(plyr)
library(raster)
library(sf)
library(stars)
library(zoo)
library(ggplot2)
library(rasterVis)
library(ggpubr)
library(tidyverse)
library(writexl)
library(dplyr)
library(readxl)
library(beepr)
library(mgcv)
library(reshape)
library(dismo)
library(patchwork)
library(modEvA)
library(devtools)
options(digits=4)  # Fewer decimal places - works better for instructor
library(sperrorest)
library(raster)
library(RColorBrewer)
library(pROC)


## Data import and cutoff ----
# load file
Area1 <-  read.csv(file="Area1_Rmodel.csv", header = TRUE, sep = ";")
Area2 <-  read.csv(file="Area2_Rmodel.csv", header = TRUE, sep = ";")

Area1 <- na.omit(Area1)
Area2 <- na.omit(Area2)

# make presence/absence a factor
Area1$Crollo <- factor(Area1$Crollo)
Area1$Geologia <- factor(Area1$Geologia)

Area2$Crollo <- factor(Area2$Crollo)
Area2$Geologia <- factor(Area2$Geologia)


# cutoff per area1


Area1.trafo <- function(x) {
  
  x$Slope[x$Slope > 71] <- 71
  x$PlanarC[x$PlanarC < (-5)] <- (-5)
  x$PlanarC[x$PlanarC > (5)] <- 5
  x$ProfileC[x$ProfileC > (5)] <- 5
  x$ProfileC[x$ProfileC < (-5)] <- (-5)
  x$NGSS[x$NGSS < 55] <- 55
  x$NGSS[x$NGSS > 118] <- 118
  x$FT[x$FT < 17] <- 17
  x$FT[x$FT > 75] <- 75
  
  
  return(x)
}

# cutoff per area2
Area2.trafo <- function(x) {
  
  x$Slope[x$Slope > 71] <- 71
  x$Slope[x$Slope < 10] <- 10
  x$PlanarC[x$PlanarC < (-5)] <- (-5)
  x$PlanarC[x$PlanarC > (5)] <- 5
  x$ProfileC[x$ProfileC > (5)] <- 5
  x$ProfileC[x$ProfileC < (-5)] <- (-5)
  x$TWI[x$TWI > (4.5)] <- 4.5
  x$TWI[x$TWI < (0)] <-0
  x$WD[x$WD >(280)] <-(280)
  x$WD[x$WD <(130)] <-(130)
  x$NGSS[x$NGSS < 55] <- 55
  x$NGSS[x$NGSS > 118] <- 118
  x$FT[x$FT < 17] <- 17
  x$FT[x$FT > 75] <- 75
  x$EWI[x$EWI > 17] <- 17
  x$EWI[x$EWI < 14] <- 14
  
  return(x)
}

# aggiusta dati con le funzioni di cutoff appena create
Area1 <- Area1.trafo(Area1)
Area2 <- Area2.trafo(Area2)

# General formulas to fit presence1 
fo.1 <- Crollo ~ Cord_Y + Slope + 
  East+ PlanarC + ProfileC + WD + EWI + FT

sfo.1 <- Crollo ~ s(Cord_Y) + s(Slope) + 
  s(East) + s(PlanarC) + s(ProfileC) + s(WD) + s(EWI) + s(FT)


# ANALISI PRELIMINARE DELLE CSF
# definisci funzione generale
my.gam <- function(formula, data, family = binomial,smooth_formula,method="REML",select=TRUE) {
  fit <- gam(smooth_formula, data, family = family,method=method,select=TRUE)
  return(fit)
}

#fitta modello per presence1
fit.Area1.1 <- my.gam(fo.1,Area1,smooth_formula=sfo.1,method="REML",select=TRUE)
fit.Area2.1 <- my.gam(fo.1,Area2,smooth_formula=sfo.1,method="REML",select=TRUE)


## Area1 ----

### Non-spatial CV ----

keeps = c("Cord_Y", "Slope", "East", "PlanarC", "ProfileC", "WD", "EWI", "FT")

res <- sperrorest(
  formula = fo.1,
  data = Area1,
  model_fun = my.gam,
  model_args = list(smooth_formula = sfo.1, method = "REML", select = TRUE),
  pred_args = list(type = "response"),
  smp_fun = partition_cv,
  smp_args = list(
    nfold = 5,
    repetition = 1:10
  )
)


get_split_data <- function(res, repetition, fold, data) {
  indexes = res[["represampling"]][[repetition]][[fold]][[data]]
  return(Area1[indexes, keeps, drop = FALSE])
}

get_auroc_values <- function(res, repetition, fold, data) {
  auroc = res[[2]][[repetition]][[fold]][[data]][["auroc"]]
  return(auroc)
}

calculate_mess_cv <- function(res) {
  # Initialize an empty list to collect rows (better performance than growing a dataframe)
  results_list <- list()
  
  # Loop through repetitions
  for (repetition in seq_along(res$represampling)) {
    # Loop through folds
    for (fold in seq_along(res$represampling[[repetition]])) {
      # Get training and testing data
      train_dt <- get_split_data(res, repetition, fold, "train")
      test_dt <- get_split_data(res, repetition, fold, "test")
      
      # Calculate MESS
      mess_table <- MESS(P = test_dt, V = train_dt)
      
      # Summary statistics
      min_mess <- min(mess_table[,"TOTAL"], na.rm = TRUE)
      max_mess <- max(mess_table[,"TOTAL"], na.rm = TRUE)
      mean_mess <- mean(mess_table[,"TOTAL"], na.rm = TRUE)
      median_mess <- median(mess_table[,"TOTAL"], na.rm = TRUE)
      q10_mess <- quantile(mess_table[,"TOTAL"], probs = 0.10, na.rm = TRUE)
      q90_mess <- quantile(mess_table[,"TOTAL"], probs = 0.90, na.rm = TRUE)
      
      # Adding auroc
      auroc_train = get_auroc_values(res, repetition, fold, "train")
      auroc_test = get_auroc_values(res, repetition, fold, "test")
      diff_perf = auroc_train - auroc_test
      
      # Add results to the list
      results_list[[length(results_list) + 1]] <- data.frame(
        repetition = repetition,
        fold = fold,
        min_mess = min_mess,
        max_mess = max_mess,
        mean_mess = mean_mess,
        median_mess = median_mess,
        q10_mess = q10_mess,
        q90_mess = q90_mess,
        auroc_train = auroc_train,
        auroc_test = auroc_test,
        diff_perf = diff_perf
      )
    }
  }
  
  # Combine all rows into a dataframe
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  # Saving
  write.csv(results_df, "Area1/mess_cv.csv", row.names = FALSE)
  
  # Return the final dataframe
  return(results_df)
}

res_mess <- calculate_mess_cv(res)

# Create the individual plots
plot_mean <- ggplot(res_mess, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Mean_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_median <- ggplot(res_mess, aes(x = median_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Median_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_q10 <- ggplot(res_mess, aes(x = q10_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Q10_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_q90 <- ggplot(res_mess, aes(x = q90_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Q90_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine the plots using patchwork and add a title
combined_plot <- (plot_mean | plot_median) / (plot_q10 | plot_q90) +
  plot_annotation(
    title = "Random Cross-Validation in Area1: MESS vs diff_perf",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

print(combined_plot)
ggsave("plots/MESS_Area1/CV_mess_vs_diff_perf.jpg", plot = combined_plot, width = 16, height = 12, dpi = 300)

write.csv(res_mess, "Area1/mess_auroc_CV.csv", row.names = FALSE)
# View the combined dataframe
print(res_mess)

# Create the boxplot for train and test AUROC scores
boxplot_auroc <- ggplot(res_mess, aes(x = factor(1), y = auroc_train)) +
  geom_boxplot(aes(fill = "Train"), color = "goldenrod1", alpha = 0.7) +
  geom_boxplot(aes(x = factor(2), y = auroc_test, fill = "Test"), color = "firebrick", alpha = 0.7) +
  scale_fill_manual(values = c("Train" = "goldenrod1", "Test" = "firebrick")) +
  theme_minimal() +
  labs(
    title = "Train vs Test AUROC Scores",
    x = "Data Partition",
    y = "AUROC Score"
  ) +
  scale_x_discrete(labels = c("Train", "Test")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Save the boxplot as an image
ggsave("plots/MESS_Area1/boxplot_auroc_scores.jpg", plot = boxplot_auroc, width = 8, height = 6, dpi = 300)

# Display the plot
print(boxplot_auroc)

# Filter the data for different ranges of auroc_train
res_mess_above_0.9 <- subset(res_mess, auroc_train > 0.9)
res_mess_between_0.85_0.9 <- subset(res_mess, auroc_train > 0.85 & auroc_train <= 0.9)
res_mess_between_0.8_0.85 <- subset(res_mess, auroc_train > 0.8 & auroc_train <= 0.85)
res_mess_below_0.8 <- subset(res_mess, auroc_train <= 0.8)

# Create individual scatter plots
plot_above_0.9 <- ggplot(res_mess_above_0.9, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train > 0.9", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_between_0.85_0.9 <- ggplot(res_mess_between_0.85_0.9, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train 0.85 - 0.9", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_between_0.8_0.85 <- ggplot(res_mess_between_0.8_0.85, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train 0.8 - 0.85", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_below_0.8 <- ggplot(res_mess_below_0.8, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train < 0.8", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots into a 2x2 layout using patchwork
final_plot <- (plot_above_0.9 + plot_between_0.85_0.9) / 
  (plot_between_0.8_0.85 + plot_below_0.8)

# Combine plots into a 1x2 layout using patchwork
final_plot <- (plot_between_0.85_0.9 | plot_between_0.8_0.85)  +
  plot_annotation(
    title = "Cross-Validation: MESS vs diff_perf suddivisi per Train Auroc",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# Display the combined plot
print(final_plot)

ggsave("plots/MESS_Area1/CV_divided_train.jpg", plot = final_plot, width = 16, height = 12, dpi = 300)


### Spatial CV ----

keeps = c("Cord_Y", "Slope", "East", "PlanarC", "ProfileC", "WD", "EWI", "FT")

res <- sperrorest(
  coords = c("Cord_X", "Cord_Y"),
  formula = fo.1,
  data = Area1,
  model_fun = my.gam,
  model_args = list(smooth_formula = sfo.1, method = "REML", select = TRUE),
  pred_args = list(type = "response"),
  smp_fun = partition_kmeans,
  smp_args = list(
    nfold = 5,
    repetition = 1:10
    #balancing_steps = 1
  )
)

res_mess <- calculate_mess_cv(res)

# Create the individual plots
plot_mean <- ggplot(res_mess, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Mean_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_median <- ggplot(res_mess, aes(x = median_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Median_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_q10 <- ggplot(res_mess, aes(x = q10_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Q10_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_q90 <- ggplot(res_mess, aes(x = q90_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Q90_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine the plots using patchwork and add a title
combined_plot <- (plot_mean | plot_median) / (plot_q10 | plot_q90) +
  plot_annotation(
    title = "K-mean CV (in Area1): MESS vs diff_perf",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
ggsave("plots/MESS_Area1/kmean_mess_vs_diff_perf.jpg", plot = combined_plot, width = 16, height = 12, dpi = 300)

write.csv(res_mess, "Area1/mess_auroc_CV_kmean.csv", row.names = FALSE)
# View the combined dataframe
print(res_mess)

# Create the boxplot for train and test AUROC scores
boxplot_auroc <- ggplot(res_mess, aes(x = factor(1), y = auroc_train)) +
  geom_boxplot(aes(fill = "Train"), color = "goldenrod1", alpha = 0.7) +
  geom_boxplot(aes(x = factor(2), y = auroc_test, fill = "Test"), color = "firebrick", alpha = 0.7) +
  scale_fill_manual(values = c("Train" = "goldenrod1", "Test" = "firebrick")) +
  theme_minimal() +
  labs(
    title = "Kmean Train vs Test AUROC Scores",
    x = "Data Partition",
    y = "AUROC Score"
  ) +
  scale_x_discrete(labels = c("Train", "Test")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Save the boxplot as an image
ggsave("plots/MESS_Area1/boxplot_auroc_scores_kmean.jpg", plot = boxplot_auroc, width = 8, height = 6, dpi = 300)

# Display the plot
print(boxplot_auroc)

# Filter the data for different ranges of auroc_train
res_mess_above_0.9 <- subset(res_mess, auroc_train > 0.9)
res_mess_between_0.85_0.9 <- subset(res_mess, auroc_train > 0.85 & auroc_train <= 0.9)
res_mess_between_0.8_0.85 <- subset(res_mess, auroc_train > 0.8 & auroc_train <= 0.85)
res_mess_below_0.8 <- subset(res_mess, auroc_train <= 0.8)

# Create individual scatter plots
plot_above_0.9 <- ggplot(res_mess_above_0.9, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train > 0.9", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_between_0.85_0.9 <- ggplot(res_mess_between_0.85_0.9, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train 0.85 - 0.9", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_between_0.8_0.85 <- ggplot(res_mess_between_0.8_0.85, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train 0.8 - 0.85", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_below_0.8 <- ggplot(res_mess_below_0.8, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train < 0.8", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots into a 2x2 layout using patchwork
final_plot <- (plot_above_0.9 + plot_between_0.85_0.9) / 
  (plot_between_0.8_0.85 + plot_below_0.8) +
  plot_annotation(
    title = "K-mean (Area1): MESS vs diff_perf suddivisi per Train Auroc",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )


# Save the boxplot as an image
ggsave("plots/MESS_Area1/A1_kmean_divided_train.jpg", plot = final_plot, width = 8, height = 6, dpi = 300)


# Display the combined plot
print(final_plot)


## Area2 ----
### Non-spatial CV ----

keeps = c("Cord_Y", "Slope", "East", "PlanarC", "ProfileC", "WD", "EWI", "FT")

res <- sperrorest(
  formula = fo.1,
  data = Area2,
  model_fun = my.gam,
  model_args = list(smooth_formula = sfo.1, method = "REML", select = TRUE),
  pred_args = list(type = "response"),
  smp_fun = partition_cv,
  smp_args = list(
    nfold = 5,
    repetition = 1:10
  )
)


get_split_data <- function(res, repetition, fold, data) {
  indexes = res[["represampling"]][[repetition]][[fold]][[data]]
  return(Area2[indexes, keeps, drop = FALSE])
}

get_auroc_values <- function(res, repetition, fold, data) {
  auroc = res[[2]][[repetition]][[fold]][[data]][["auroc"]]
  return(auroc)
}

calculate_mess_cv <- function(res) {
  # Initialize an empty list to collect rows (better performance than growing a dataframe)
  results_list <- list()
  
  # Loop through repetitions
  for (repetition in seq_along(res$represampling)) {
    # Loop through folds
    for (fold in seq_along(res$represampling[[repetition]])) {
      # Get training and testing data
      train_dt <- get_split_data(res, repetition, fold, "train")
      test_dt <- get_split_data(res, repetition, fold, "test")
      
      # Calculate MESS
      mess_table <- MESS(P = test_dt, V = train_dt)
      
      # Summary statistics
      min_mess <- min(mess_table[,"TOTAL"], na.rm = TRUE)
      max_mess <- max(mess_table[,"TOTAL"], na.rm = TRUE)
      mean_mess <- mean(mess_table[,"TOTAL"], na.rm = TRUE)
      median_mess <- median(mess_table[,"TOTAL"], na.rm = TRUE)
      q10_mess <- quantile(mess_table[,"TOTAL"], probs = 0.10, na.rm = TRUE)
      q90_mess <- quantile(mess_table[,"TOTAL"], probs = 0.90, na.rm = TRUE)
      
      # Adding auroc
      auroc_train = get_auroc_values(res, repetition, fold, "train")
      auroc_test = get_auroc_values(res, repetition, fold, "test")
      diff_perf = auroc_train - auroc_test
      
      # Add results to the list
      results_list[[length(results_list) + 1]] <- data.frame(
        repetition = repetition,
        fold = fold,
        min_mess = min_mess,
        max_mess = max_mess,
        mean_mess = mean_mess,
        median_mess = median_mess,
        q10_mess = q10_mess,
        q90_mess = q90_mess,
        auroc_train = auroc_train,
        auroc_test = auroc_test,
        diff_perf = diff_perf
      )
    }
  }
  
  # Combine all rows into a dataframe
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  # Saving
  write.csv(results_df, "Area2/mess_cv.csv", row.names = FALSE)
  
  # Return the final dataframe
  return(results_df)
}

res_mess <- calculate_mess_cv(res)

# Create the individual plots
plot_mean <- ggplot(res_mess, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Mean_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_median <- ggplot(res_mess, aes(x = median_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Median_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_q10 <- ggplot(res_mess, aes(x = q10_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Q10_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_q90 <- ggplot(res_mess, aes(x = q90_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Q90_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine the plots using patchwork and add a title
combined_plot <- (plot_mean | plot_median) / (plot_q10 | plot_q90) +
  plot_annotation(
    title = "Random Cross-Validation in Area2: MESS vs diff_perf",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

print(combined_plot)
ggsave("plots/MESS_Area2/A2_CV_mess_vs_diff_perf.jpg", plot = combined_plot, width = 16, height = 12, dpi = 300)

write.csv(res_mess, "Area2/mess_auroc_CV.csv", row.names = FALSE)
# View the combined dataframe
print(res_mess)

# Create the boxplot for train and test AUROC scores
boxplot_auroc <- ggplot(res_mess, aes(x = factor(1), y = auroc_train)) +
  geom_boxplot(aes(fill = "Train"), color = "goldenrod1", alpha = 0.7) +
  geom_boxplot(aes(x = factor(2), y = auroc_test, fill = "Test"), color = "firebrick", alpha = 0.7) +
  scale_fill_manual(values = c("Train" = "goldenrod1", "Test" = "firebrick")) +
  theme_minimal() +
  labs(
    title = "Area2 - Train vs Test AUROC Scores",
    x = "Data Partition",
    y = "AUROC Score"
  ) +
  scale_x_discrete(labels = c("Train", "Test")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Save the boxplot as an image
ggsave("plots/MESS_Area2/boxplot_auroc_scores.jpg", plot = boxplot_auroc, width = 8, height = 6, dpi = 300)

# Display the plot
print(boxplot_auroc)

# Filter the data for different ranges of auroc_train
res_mess_above_0.9 <- subset(res_mess, auroc_train > 0.9)
res_mess_between_0.85_0.9 <- subset(res_mess, auroc_train > 0.85 & auroc_train <= 0.9)
res_mess_between_0.8_0.85 <- subset(res_mess, auroc_train > 0.8 & auroc_train <= 0.85)
res_mess_below_0.8 <- subset(res_mess, auroc_train <= 0.8)

# Create individual scatter plots
plot_above_0.9 <- ggplot(res_mess_above_0.9, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train > 0.9", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_between_0.85_0.9 <- ggplot(res_mess_between_0.85_0.9, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train 0.85 - 0.9", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_between_0.8_0.85 <- ggplot(res_mess_between_0.8_0.85, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train 0.8 - 0.85", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_below_0.8 <- ggplot(res_mess_below_0.8, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train < 0.8", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots into a 2x2 layout using patchwork
final_plot <- (plot_above_0.9 + plot_between_0.85_0.9) +
  plot_annotation(
    title = "Cross-Validation: MESS vs diff_perf suddivisi per Train Auroc",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

ggsave("plots/MESS_Area2/A2_CV_divided_train.jpg", plot = final_plot, width = 16, height = 12, dpi = 300)


# Display the combined plot
print(final_plot)

### Spatial CV ----

keeps = c("Cord_Y", "Slope", "East", "PlanarC", "ProfileC", "WD", "EWI", "FT")

res <- sperrorest(
  coords = c("Cord_X", "Cord_Y"),
  formula = fo.1,
  data = Area2,
  model_fun = my.gam,
  model_args = list(smooth_formula = sfo.1, method = "REML", select = TRUE),
  pred_args = list(type = "response"),
  smp_fun = partition_kmeans,
  smp_args = list(
    nfold = 5,
    repetition = 1:10
    #balancing_steps = 1
  )
)

res_mess <- calculate_mess_cv(res)

# Create the individual plots
plot_mean <- ggplot(res_mess, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Mean_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_median <- ggplot(res_mess, aes(x = median_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Median_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_q10 <- ggplot(res_mess, aes(x = q10_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Q10_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

plot_q90 <- ggplot(res_mess, aes(x = q90_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Q90_mess vs diff_perf") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine the plots using patchwork and add a title
combined_plot <- (plot_mean | plot_median) / (plot_q10 | plot_q90) +
  plot_annotation(
    title = "K-mean CV (in Area2): MESS vs diff_perf",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

print(combined_plot)
ggsave("plots/MESS_Area2/kmean_mess_vs_diff_perf.jpg", plot = combined_plot, width = 16, height = 12, dpi = 300)

write.csv(res_mess, "Area2/mess_auroc_CV_kmean.csv", row.names = FALSE)
# View the combined dataframe
print(res_mess)

# Create the boxplot for train and test AUROC scores
boxplot_auroc <- ggplot(res_mess, aes(x = factor(1), y = auroc_train)) +
  geom_boxplot(aes(fill = "Train"), color = "goldenrod1", alpha = 0.7) +
  geom_boxplot(aes(x = factor(2), y = auroc_test, fill = "Test"), color = "firebrick", alpha = 0.7) +
  scale_fill_manual(values = c("Train" = "goldenrod1", "Test" = "firebrick")) +
  theme_minimal() +
  labs(
    title = "Kmean Train vs Test AUROC Scores",
    x = "Data Partition",
    y = "AUROC Score"
  ) +
  scale_x_discrete(labels = c("Train", "Test")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Save the boxplot as an image
ggsave("plots/MESS_Area2/boxplot_auroc_scores_kmean.jpg", plot = boxplot_auroc, width = 8, height = 6, dpi = 300)

# Display the plot
print(boxplot_auroc)

# Filter the data for different ranges of auroc_train
res_mess_above_0.9 <- subset(res_mess, auroc_train > 0.9)
res_mess_between_0.85_0.9 <- subset(res_mess, auroc_train > 0.85 & auroc_train <= 0.9)
res_mess_between_0.8_0.85 <- subset(res_mess, auroc_train > 0.8 & auroc_train <= 0.85)
res_mess_below_0.8 <- subset(res_mess, auroc_train <= 0.8)

# Create individual scatter plots
plot_above_0.9 <- ggplot(res_mess_above_0.9, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train > 0.9", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_between_0.85_0.9 <- ggplot(res_mess_between_0.85_0.9, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train 0.85 - 0.9", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_between_0.8_0.85 <- ggplot(res_mess_between_0.8_0.85, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train 0.8 - 0.85", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

plot_below_0.8 <- ggplot(res_mess_below_0.8, aes(x = mean_mess, y = diff_perf)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "AUROC Train < 0.8", x = "Mean Mess", y = "Delta AUROC") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots into a 2x2 layout using patchwork
final_plot <- (plot_above_0.9 + plot_between_0.85_0.9) + 
  plot_annotation(
    title = "K-mean: MESS vs diff_perf suddivisi per Train Auroc",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# Display the combined plot
print(final_plot)

ggsave("plots/MESS_Area2/A2_kmean_divided_train.jpg", plot = final_plot, width = 8, height = 6, dpi = 300)



## Mess con raster Area 2 ####
keeps = c("Slope", "East", "PlanarC", "ProfileC", "WD", "EWI",  "FT")

train_data <- Area2[ , keeps, drop = FALSE]

## RASTER MESS
import_A2 <- function() {
  library(raster)
  
  # Load raster layers for Area2
  slope <- raster("Area2/Slope_A2.tif")
  east <- raster("Area2/East_A2.tif")
  planarC <- raster("Area2/PlanarC_A2.tif")
  profileC <- raster("Area2/ProfileC_A2.tif")
  ft <- raster("Area2/FT_RES_A2.tif")
  ewi <- raster("Area2/EWI_A2.tif")
  wd <- raster("Area2/WD_A2.tif")
  
  # Apply thresholding
  slope[slope > 90] <- 30.67
  slope[slope > 71] <- 71
  east[east > 4] <- 2.55
  planarC[planarC > 5] <- 5
  planarC[planarC < -5] <- -5
  profileC[profileC > 5] <- 5
  profileC[profileC < -5] <- -5
  ft[ft < 17] <- 17
  ft[ft > 75] <- 75
  ewi[ewi > 17] <- 17
  ewi[ewi < 14] <- 14
  wd[wd > 280] <- 280
  wd[wd < 130] <- 130
  
  # Align rasters to match the slope raster
  ft <- resample(ft, slope, method = "bilinear")
  ewi <- resample(ewi, slope, method = "bilinear")
  wd <- resample(wd, slope, method = "bilinear")
  
  # Stack all raster layers into a RasterStack
  s <- stack(slope, east, planarC, profileC, wd, ewi, ft)
  names(s) <- c("Slope", "East", "PlanarC", "ProfileC", "WD", "EWI", "FT")
  
  return(s)
}

# Example usage:
s2_agg <- import_A2()

# We try to reduce the size for increase speed and a better view
s2_agg <- aggregate(s2_agg, fact = 1, fun = mean)


mess <- mess(s2_agg, train_data, full = TRUE, filename = "plots/MESS_raster/mess_Area2.tiff", overwrite=TRUE)
names(mess) <- c("Slope", "East", "PlanarC", "ProfileC", "WD", "EWI", "FT", "TOTAL")

### Individual plot MESS ----
plot(mess)
plot(mess[["FT"]], main="FT")
plot(mess[["Slope"]], main="Slope")
plot(mess[["East"]], main="East")
plot(mess[["PlanarC"]], main="PlanarC")
plot(mess[["ProfileC"]], main="ProfileC")
plot(mess[["WD"]], main="WD")
plot(mess[["EWI"]], main="EWI")
#plot(mess[["NGSS"]], main="NGSS")
plot(mess[["FT"]], main="FT")
plot(mess[["TOTAL"]], main="TOTAL")

## Mess con raster Area 2 ####
keeps = c("Slope", "East", "PlanarC", "ProfileC", "WD", "EWI",  "FT")

train_data <- Area2[ , keeps, drop = FALSE]

## RASTER MESS
s1_agg

# Example usage:
s1_agg <- import_A1()

s1_agg <- dropLayer(s1_agg, which(names(s1_agg) == "Cord_Y"))

# We try to reduce the size for increase speed and a better view
s1_agg <- aggregate(s1_agg, fact = 1, fun = mean)


mess <- mess(s2_agg, train_data, full = TRUE, filename = "plots/MESS_raster/mess_Area1.tiff", overwrite=TRUE)
names(mess) <- c("Slope", "East", "PlanarC", "ProfileC", "WD", "EWI", "FT", "TOTAL")

### Individual plot MESS ----
plot(mess)
plot(mess[["FT"]], main="FT")
plot(mess[["Slope"]], main="Slope")
plot(mess[["East"]], main="East")
plot(mess[["PlanarC"]], main="PlanarC")
plot(mess[["ProfileC"]], main="ProfileC")
plot(mess[["WD"]], main="WD")
plot(mess[["EWI"]], main="EWI")
#plot(mess[["NGSS"]], main="NGSS")
plot(mess[["FT"]], main="FT")
plot(mess[["TOTAL"]], main="TOTAL")


## Apply model ####
# Re-defining the functions
# General formulas to fit presence1 
fo.1 <- Crollo ~ Cord_Y + Slope + 
  East+ PlanarC + ProfileC + WD + EWI + FT

sfo.1 <- Crollo ~ s(Cord_Y) + s(Slope) + 
  s(East) + s(PlanarC) + s(ProfileC) + s(WD) + s(EWI) + s(FT)

# Function to train the GAM model
my.gam <- function(formula, data, family = binomial, smooth_formula, method = "REML", select = TRUE) {
  fit <- gam(smooth_formula, data = data, family = family, method = method, select = select)
  return(fit)
}

#### Area1 per predire Area1 (train) ----
import_A1 <- function() {
  
  # Load raster layers
  slope <- raster("Area1/SLOPEA1.tif")
  east <- raster("Area1/EASTA1.tif")
  planarC <- raster("Area1/PLANCA1.tif")
  profileC <- raster("Area1/PROFCA1.tif")
  ft <- raster("Area1/FTA1.tif")
  ewi <- raster("Area1/EWIA1.tif")
  wd <- raster("Area1/WDA1.tif")
  
  # Apply thresholds and transformations
  slope[slope > 90] <- 27.43
  slope[slope > 71] <- 71
  east[east > 4] <- 0.002
  planarC[planarC > 5] <- 5
  planarC[planarC < -5] <- -5
  profileC[profileC > 5] <- 5
  profileC[profileC < -5] <- -5
  ft[ft < 17] <- 17
  ft[ft > 75] <- 75
  ewi[ewi > 17] <- 17
  ewi[ewi < 14] <- 14
  wd[wd > 280] <- 280
  wd[wd < 130] <- 130
  
  # Align raster resolutions using resampling
  ft <- resample(ft, slope, method = "bilinear")
  ewi <- resample(ewi, slope, method = "bilinear")
  wd <- resample(wd, slope, method = "bilinear")
  
  # Stack rasters into a multi-layer object
  s <- stack(slope, east, planarC, profileC, wd, ewi, ft)
  names(s) <- c("Slope", "East", "PlanarC", "ProfileC", "WD", "EWI", "FT")
  
  # Extract the y-coordinates
  y_coords <- coordinates(s)[, 2]
  y_raster <- raster(s[["Slope"]])
  values(y_raster) <- y_coords
  
  s <- stack(y_raster, s)
  # Rename raster layers to match the model variables
  names(s) <- c("Cord_Y", "Slope", "East", "PlanarC", "ProfileC", "WD", "EWI", "FT")
  
  return(s)
}

# Example usage:
s1 <- import_A1()


# We try to reduce the size for increase speed
# Change here to increase granularity -> Fact 1 = real size
s1_agg <- aggregate(s1, fact = 1, fun = mean)

# Train the GAM model (if not already trained)
fit.Area1.1 <- my.gam(formula = fo.1, data = Area1, family = binomial, smooth_formula = sfo.1)

# ROC curve train
predictions <- predict(fit.Area1.1, newdata = Area1, type = "response")
par(pty = "s")
roc_obj_train <- roc(Area1$Crollo, predictions, 
               percent=TRUE, 
               plot=TRUE, 
               legacy.axes = TRUE,
               xlab = "Percentuale di area comulata [%]",
               ylab = "Percentuale comulata dei fenomeni di frana [%]",
               print.auc = TRUE
               )

# Extract Landslide Points from Area1
landslide_points <- Area1[Area1$Crollo == 1, c("Cord_X", "Cord_Y")]

# Generate Susceptibility Map with Landslide Points
# Predict probabilities for raster
prob_raster <- predict(s1_agg, model = fit.Area1.1, type = "response")

# Save the full raster stack/brick as a TIFF
writeRaster(prob_raster, filename = "plots/Predictions/Area1_in_Area1_suscettibilita.tif", 
            format = "GTiff", overwrite = TRUE)

# Convert raster to dataframe
df <- as.data.frame(rasterToPoints(prob_raster))
colnames(df) <- c("x", "y", "probability")

# Create susceptibility map and overlay landslide points
raster_plot <- ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = probability)) +
  scale_fill_gradientn(colors = c("green", "yellow", "red"), name = "Susceptibility", limits = c(0, 1)) +
  geom_point(data = landslide_points, aes(x = Cord_X, y = Cord_Y), color = "black", size = 1.5, shape = 16, alpha = 0.7) +  # Overlay landslide points
  labs(title = "Mappa di Suscettibilità di Area1 con Frane") +
  theme_minimal() +
  coord_fixed() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

print(raster_plot)


#### Area1 per predire Area2 (test) ----
# Extract the y-coordinates
s2 <- import_A2()
y_coords <- coordinates(s2)[, 2]
y_raster <- raster(s2[["Slope"]])
values(y_raster) <- y_coords

s2 <- stack(y_raster, s2)
# Rename raster layers to match the model variables
names(s2) <- c("Cord_Y", "Slope", "East", "PlanarC", "ProfileC", "WD", "EWI", "FT")

s2_agg <- aggregate(s2, fact = 1, fun = mean)

fit.Area1.1 <- my.gam(formula = fo.1, data = Area1, family = binomial, smooth_formula = sfo.1)

predictions <- predict(fit.Area1.1, newdata = Area2, type = "response")
par(pty = "s")
roc_obj_train <- roc(Area2$Crollo, predictions, 
                     percent=TRUE, 
                     plot=TRUE, 
                     legacy.axes = TRUE,
                     xlab = "Percentuale di area comulata [%]",
                     ylab = "Percentuale comulata dei fenomeni di frana [%]",
                     print.auc = TRUE
)

# Extract Landslide Points from Area1
landslide_points <- Area2[Area2$Crollo == 1, c("Cord_X", "Cord_Y")]

# Generate Susceptibility Map with Landslide Points
# Predict probabilities for raster
prob_raster <- predict(s2_agg, model = fit.Area1.1, type = "response")

# Save the full raster stack/brick as a TIFF
writeRaster(prob_raster, filename = "plots/Predictions/Area1_in_Area2_suscettibilita.tif", 
            format = "GTiff", overwrite = TRUE)

# Convert raster to dataframe
df <- as.data.frame(rasterToPoints(prob_raster))
colnames(df) <- c("x", "y", "probability")

# Create susceptibility map and overlay landslide points
raster_plot <- ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = probability)) +
  scale_fill_gradientn(colors = c("green", "yellow", "red"), name = "Susceptibility", limits = c(0, 1)) +
  geom_point(data = landslide_points, aes(x = Cord_X, y = Cord_Y), color = "black", size = 1.5, shape = 16, alpha = 0.7) +  # Overlay landslide points
  labs(title = "Mappa di Suscettibilità di Area2 con Frane") +
  theme_minimal() +
  coord_fixed() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

print(raster_plot)

#### Area2 per predire Area2 (train) ----
fit.Area1.1 <- my.gam(formula = fo.1, data = Area2, family = binomial, smooth_formula = sfo.1)

predictions <- predict(fit.Area1.1, newdata = Area2, type = "response")
par(pty = "s")
roc_obj_train <- roc(Area2$Crollo, predictions, 
                     percent=TRUE, 
                     plot=TRUE, 
                     legacy.axes = TRUE,
                     xlab = "Percentuale di area comulata [%]",
                     ylab = "Percentuale comulata dei fenomeni di frana [%]",
                     print.auc = TRUE
)

# Extract Landslide Points from Area1
landslide_points <- Area2[Area2$Crollo == 1, c("Cord_X", "Cord_Y")]

# Generate Susceptibility Map with Landslide Points
# Predict probabilities for raster
prob_raster <- predict(s2_agg, model = fit.Area1.1, type = "response")

# Save the full raster stack/brick as a TIFF
writeRaster(prob_raster, filename = "plots/Predictions/Area2_in_Area2_suscettibilita.tif", 
            format = "GTiff", overwrite = TRUE)

# Convert raster to dataframe
df <- as.data.frame(rasterToPoints(prob_raster))
colnames(df) <- c("x", "y", "probability")

# Create susceptibility map and overlay landslide points
raster_plot <- ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = probability)) +
  scale_fill_gradientn(colors = c("green", "yellow", "red"), name = "Susceptibility", limits = c(0, 1)) +
  geom_point(data = landslide_points, aes(x = Cord_X, y = Cord_Y), color = "black", size = 1.5, shape = 16, alpha = 0.7) +  # Overlay landslide points
  labs(title = "Mappa di Suscettibilità di Area2 con Frane") +
  theme_minimal() +
  coord_fixed() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

print(raster_plot)

#### Area2 per predire Area1 (test) ----
fit.Area2.1 <- my.gam(formula = fo.1, data = Area2, family = binomial, smooth_formula = sfo.1)

predictions <- predict(fit.Area2.1, newdata = Area1, type = "response")
par(pty = "s")
roc_obj_train <- roc(Area1$Crollo, predictions, 
                     percent=TRUE, 
                     plot=TRUE, 
                     legacy.axes = TRUE,
                     xlab = "Percentuale di area comulata [%]",
                     ylab = "Percentuale comulata dei fenomeni di frana [%]",
                     print.auc = TRUE
)

# Extract Landslide Points from Area1
landslide_points <- Area1[Area1$Crollo == 1, c("Cord_X", "Cord_Y")]

# Generate Susceptibility Map with Landslide Points
# Predict probabilities for raster
prob_raster <- predict(s1_agg, model = fit.Area1.1, type = "response")

# Save the full raster stack/brick as a TIFF
writeRaster(prob_raster, filename = "plots/Predictions/Area2_in_Area1_suscettibilita.tif", 
            format = "GTiff", overwrite = TRUE)

# Convert raster to dataframe
df <- as.data.frame(rasterToPoints(prob_raster))
colnames(df) <- c("x", "y", "probability")

# Create susceptibility map and overlay landslide points
raster_plot <- ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = probability)) +
  scale_fill_gradientn(colors = c("green", "yellow", "red"), name = "Susceptibility", limits = c(0, 1)) +
  geom_point(data = landslide_points, aes(x = Cord_X, y = Cord_Y), color = "black", size = 1.5, shape = 16, alpha = 0.7) +  # Overlay landslide points
  labs(title = "Mappa di Suscettibilità di Area1 con Frane") +
  theme_minimal() +
  coord_fixed() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

print(raster_plot)