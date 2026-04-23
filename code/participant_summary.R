# Participant summary ----

setwd(".")
set.seed(0509)

library(tidyverse)
library(rstatix)

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
diseases <- c("anemia", "arthritis", "atherosclerosis", "copd", "dementia",
              "depression", "diabetes", "dyslipidemia", "hypertension",
              "hypothyroidism", "osteoporosis", "renal_disease")

frailty_associations_figure <- frailty_results %>%
  mutate(type = ifelse(term %in% diseases, "Comorbidity", "Medication")) %>%
  mutate(p_signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )) %>%
  # Order by estimate descending
  mutate(term = fct_reorder(term, estimate, .desc = TRUE)) %>%
  # Position labels slightly above or below the bars
  mutate(
    label_pos = ifelse(estimate >= 0, estimate + 0.02, estimate - 0.02),
    v_just = ifelse(estimate >= 0, 0.5, 1)
  ) %>%

  ggplot(aes(x = term, y = estimate, fill = type)) +
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
    y = "Effect size",
    fill = NULL
  )

ggsave(
  "results/frailty_associations.pdf",
  frailty_associations_figure,
  width = 6, height = 3)
