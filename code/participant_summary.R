# Participant summary ----

setwd(".")
set.seed(0509)

library(tidyverse)
library(rstatix)
library(GGally)
library(patchwork)

participants <- read_tsv("data/metadata/patient_data.tsv") %>%
  filter(SubjectID != "NW046") # Missing metadata

## Analysis ----

summary_stats <- participants %>%
  group_by(Location) %>%
  summarise(
    n = n(),
    age_median = median(Age, na.rm = TRUE),
    age_q1 = quantile(Age, 0.25, na.rm = TRUE),
    age_q3 = quantile(Age, 0.75, na.rm = TRUE),
    n_female = sum(Sex == "Female"),
    p_female = n_female / n * 100,
    n_male = sum(Sex == "Male"),
    p_male = n_male / n * 100,
    res_time_median = median(TimeNW, na.rm = TRUE),
    res_time_q1 = quantile(TimeNW, 0.25, na.rm = TRUE),
    res_time_q3 = quantile(TimeNW, 0.75, na.rm = TRUE),
    fi_score_median = median(FIScore, na.rm = TRUE),
    fi_score_q1 = quantile(FIScore, 0.25, na.rm = TRUE),
    fi_score_q3 = quantile(FIScore, 0.75, na.rm = TRUE),
    n_0.3_0.4 = sum(FIRange == "0.3-0.4"),
    p_0.3_0.4 = n_0.3_0.4 / n * 100,
    n_0.4_0.5 = sum(FIRange == "0.4-0.5"),
    p_0.4_0.5 = n_0.4_0.5 / n * 100,
    n_0.5_0.6 = sum(FIRange == "0.5-0.6"),
    p_0.5_0.6 = n_0.5_0.6 / n * 100,
    n_0.6_0.7 = sum(FIRange == "0.6-0.7"),
    p_0.6_0.7 = n_0.6_0.7 / n * 100,
    n_0.7_0.8 = sum(FIRange == "0.7-0.8"),
    p_0.7_0.8 = n_0.7_0.8 / n * 100
  )

write_tsv(summary_stats, "results/participant_summary_stats.tsv")

# Manor vs Center
participants %>%
  wilcox_test(Age ~ Location)
participants %>%
  wilcox_test(TimeNW ~ Location)
participants %>%
  wilcox_test(FIScore ~ Location)
chisq_test(table(participants$Location, participants$Sex))

# Males vs Females
participants %>%
  wilcox_test(Age ~ Sex)
participants %>%
  wilcox_test(TimeNW ~ Sex)
participants %>%
  wilcox_test(FIScore ~ Sex)

# Figure S1: Participant Baseline Characteristics & Sex Differences ----

# Custom correlation function to show r and p
cor_func <- function(data, mapping, size = 4, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  test <- cor.test(x, y)

  r_val <- round(test$estimate, 2)
  p_val <- test$p.value
  p_text <- if (p_val < 0.001) "p < 0.001" else paste0("p = ", format.pval(p_val, digits = 3))

  display_text <- paste0("r = ", r_val, "\n", p_text)

  ggplot(data = data.frame()) +
    annotate("text", x = 0.5, y = 0.5, label = display_text, size = size) +
    theme_void()
}

# Custom function for lower triangle to show points and smooth
lower_func <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.4, size = 1, color = "grey30") +
    geom_smooth(method = "lm", color = "black", fill = "black", alpha = 0.2, linewidth = 0.5)
}

# Prepare data with lowercase names for plotting
plot_data <- participants %>%
  rename(age = Age, sex = Sex, res_time = TimeNW, fi_score = FIScore)

# Panel A: Correlation Matrix
cor_theme <- theme_classic() +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 11),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

panel_a <- ggpairs(
  plot_data %>% select(age, res_time, fi_score),
  columnLabels = c("Age", "Residence Time", "Frailty (FI)"),
  upper = list(continuous = wrap(cor_func, size = 4)),
  lower = list(continuous = wrap(lower_func)),
  diag = list(continuous = wrap("densityDiag", fill = "gray50", alpha = 0.4))
) + cor_theme

# Panel B: Box plots by Sex
create_sex_boxplot <- function(data, var, ylabel) {
  stat_test <- data %>%
    wilcox_test(as.formula(paste(var, "~ sex"))) %>%
    add_significance()

  ggplot(data, aes(x = sex, y = !!sym(var), fill = sex)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
    labs(x = NULL, y = ylabel) +
    annotate("text",
      x = 1.5, y = max(data[[var]], na.rm = TRUE) * 1.05,
      label = paste0("p = ", format.pval(stat_test$p, digits = 3)),
      size = 4, fontface = "italic"
    ) +
    scale_fill_manual(values = c("Female" = "#BC8F8F", "Male" = "#2E8B57")) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 10),
      axis.text.x = element_text(size = 10, face = "bold")
    )
}

p_age <- create_sex_boxplot(plot_data, "age", "Age (years)")
p_res <- create_sex_boxplot(plot_data, "res_time", "Residence time (years)")
p_fi <- create_sex_boxplot(plot_data, "fi_score", "FI")

panel_b <- (p_age / p_res / p_fi)

# Combine and Save
fig_1_composite <- wrap_elements(ggmatrix_gtable(panel_a)) + panel_b +
  plot_layout(widths = c(1.8, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30")
    )
  )

ggsave("results/Supp1_participants_baseline_vars.pdf", fig_1_composite, width = 12, height = 8)


# Frailty associations (FIScore only)
run_frailty_glm <- function(data) {
  covariates <- names(data)[which(names(data) == "anemia"):which(names(data) == "vitd")]
  formula_str <- paste("FIScore ~", paste(covariates, collapse = " + "))

  model <- glm(as.formula(formula_str), data = data, family = quasibinomial)

  # Extract coefficients
  summary_res <- as.data.frame(summary(model)$coefficients)
  colnames(summary_res) <- c("estimate", "std_error", "statistic", "p_value")
  summary_res$term <- rownames(summary_res)

  return(summary_res)
}

frailty_results <- run_frailty_glm(participants) %>%
  filter(term != "(Intercept)")

write_tsv(frailty_results, "results/frailty_glm.tsv")

## Figures ----

# Bar chart of associations
diseases <- c(
  "anemia", "arthritis", "atherosclerosis", "copd", "dementia",
  "depression", "diabetes", "dyslipidemia", "hypertension",
  "hypothyroidism", "osteoporosis", "renal_disease"
)

frailty_associations_figure <- frailty_results %>%
  mutate(term_clean = case_when(
    term == "copd" ~ "COPD",
    term == "insulin" ~ "Insulin",
    term == "iron_supplements" ~ "Iron supplements",
    term == "renal_disease" ~ "Renal disease",
    term == "dementia" ~ "Dementia",
    term == "metformin" ~ "Metformin",
    term == "atherosclerosis" ~ "Atherosclerosis",
    term == "arthritis" ~ "Arthritis",
    term == "ccbs" ~ "Calcium channel blockers",
    term == "acetaminophen" ~ "Acetaminophen",
    term == "nitrates" ~ "Nitrates",
    term == "benzodiazepines" ~ "Benzodiazepines",
    term == "statins" ~ "Statins",
    term == "diabetes" ~ "Diabetes",
    term == "diuretics" ~ "Diuretics",
    term == "ppis" ~ "Proton pump inhibitors",
    term == "laxatives" ~ "Laxatives",
    term == "aaps" ~ "Atypical antipsychotics",
    term == "anemia" ~ "Anemia",
    term == "hypertension" ~ "Hypertension",
    term == "levothyroxine" ~ "Levothyroxine",
    term == "dyslipidemia" ~ "Dyslipidemia",
    term == "depression" ~ "Depression",
    term == "tcas" ~ "Tricyclic antidepressants",
    term == "hypothyroidism" ~ "Hypothyroidism",
    term == "antiplatelets" ~ "Antiplatelets",
    term == "osteoporosis" ~ "Osteoporosis",
    term == "beta_blockers" ~ "Beta blockers",
    term == "vitb12" ~ "Vitamin B12",
    term == "vitd" ~ "Vitamin D",
    term == "ace_inhibitors" ~ "ACE inhibitors",
    term == "ssris" ~ "SSRIs",
    term == "opioids" ~ "Opioids",
    TRUE ~ term
  )) %>%
  mutate(type = ifelse(term %in% diseases, "Comorbidity", "Medication")) %>%
  mutate(p_signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )) %>%
  # Order by estimate descending
  mutate(term_clean = fct_reorder(term_clean, estimate, .desc = TRUE)) %>%
  # Position labels slightly above or below the bars
  mutate(
    label_pos = ifelse(estimate >= 0, estimate + 0.02, estimate - 0.02),
    v_just = ifelse(estimate >= 0, 0.5, 1)
  ) %>%
  ggplot(aes(x = term_clean, y = estimate, fill = type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  scale_fill_manual(values = c("Comorbidity" = "#FFF68F", "Medication" = "#BC8F8F")) +
  geom_text(aes(y = label_pos, label = p_signif, vjust = v_just), size = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.direction = "horizontal",
    legend.key.size = unit(0.3, "cm"),
    panel.grid = element_blank()
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Estimate (FI)",
    fill = NULL
  )

ggsave(
  "results/Main1_Frailty_Associations.pdf",
  frailty_associations_figure,
  width = 6, height = 3
)
