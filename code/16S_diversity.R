# Alpha and beta diversity analysis ----

setwd(".")
set.seed(0509)

library(tidyverse)
library(rstatix)
library(vegan)
library(cowplot)
library(scales)
library(interactions)
library(lme4)
library(lmerTest)

## Load data ----

participants <- read_tsv("data/metadata/patient_data.tsv") %>%
  filter(SubjectID != "NW046")

alpha_diversity <-
  map(
    list.files("data/16S/alpha_diversity", pattern = "\\.tsv$", full.names = TRUE),
    function(file) {
      df <- read_tsv(file, show_col_types = FALSE)
      colnames(df)[1] <- "sample"
      df
    }
  ) %>%
  reduce(full_join, by = "sample") %>%
  pivot_longer(cols = -sample, names_to = "metric", values_to = "value") %>%
  mutate(SubjectID = str_extract(sample, "^[^-]+")) %>%
  filter(SubjectID != "NW046")

distance_matrices_files <- list.files("data/16S/beta_diversity", pattern = "\\_distance.tsv$", full.names = TRUE)
distance_matrices <- map(distance_matrices_files, function(file) {
  dm <- read_tsv(file, show_col_types = FALSE)
  colnames(dm)[1] <- "sample"
  dm
})
names(distance_matrices) <- basename(distance_matrices_files) %>% str_remove("\\_distance.tsv$")

## Calculate microbiome stability based on beta diversity distances ----
beta_stability <- map(names(distance_matrices), function(metric_name) {
  distance_matrices[[metric_name]] %>%
    pivot_longer(cols = -sample, names_to = "sample_b", values_to = "distance") %>%
    rename(sample_a = sample) %>%
    mutate(
      SubjectID = str_extract(sample_a, "^[^-]+"),
      subject_b = str_extract(sample_b, "^[^-]+")
    ) %>%
    filter(SubjectID == subject_b, sample_a < sample_b, SubjectID != "NW046") %>%
    select(sample_a, sample_b, SubjectID, !!sym(metric_name) := distance)
}) %>%
  reduce(full_join, by = c("sample_a", "sample_b", "SubjectID"))

## Generate Figure 2A -----

alpha_range <- alpha_diversity %>%
  filter(metric == "faith_pd") %>%
  pull(value) %>%
  range(na.rm = TRUE)
beta_range <- beta_stability %>%
  pull(rpca) %>%
  range(na.rm = TRUE)
fi_range <- c(0, 1)

p_a_data <- bind_rows(
  # Alpha diversity facet
  alpha_diversity %>%
    filter(metric == "faith_pd") %>%
    select(SubjectID, value) %>%
    mutate(metric_type = "Alpha", facet = "Alpha diversity"),
  participants %>%
    select(SubjectID, value = FIScore) %>%
    mutate(
      value = rescale(value, to = alpha_range, from = fi_range),
      metric_type = "Frailty", facet = "Alpha diversity"
    ),
  # Beta diversity facet
  beta_stability %>%
    select(SubjectID, value = rpca) %>%
    mutate(metric_type = "Beta", facet = "Beta diversity"),
  participants %>%
    select(SubjectID, value = FIScore) %>%
    mutate(
      value = rescale(value, to = beta_range, from = fi_range),
      metric_type = "Frailty", facet = "Beta diversity"
    )
) %>%
  mutate(
    SubjectID = factor(SubjectID, levels = participants %>%
      arrange(FIScore) %>%
      pull(SubjectID)),
    facet = factor(facet, levels = c("Alpha diversity", "Beta diversity"), labels = c("Phylogenetic diversity (alpha)", "Phylogenetic RPCA (beta)"))
  )

p_a <- ggplot(p_a_data, aes(x = value, y = SubjectID, fill = metric_type)) +
  geom_point(
    data = subset(p_a_data, metric_type == "Frailty"),
    shape = 23, size = 1.2, color = "black", stroke = 0.4, alpha = 0.8
  ) +
  geom_point(
    data = subset(p_a_data, metric_type != "Frailty"),
    shape = 21, size = 1, alpha = 0.8
  ) +
  facet_wrap(~facet, scales = "free_x") +
  # Custom "legend" in each facet (box + diamond + text)
  geom_rect(
    data = data.frame(
      facet = factor(c("Phylogenetic diversity (alpha)", "Phylogenetic RPCA (beta)"), levels = c("Phylogenetic diversity (alpha)", "Phylogenetic RPCA (beta)")),
      xmin = c(max(alpha_range) - 0.2 * diff(alpha_range), max(beta_range) - 0.2 * diff(beta_range)),
      xmax = c(max(alpha_range) + 0.02 * diff(alpha_range), max(beta_range) + 0.02 * diff(beta_range)),
      ymin = 3,
      ymax = 5
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white", color = "black", linewidth = 0.2, inherit.aes = FALSE
  ) +
  geom_point(
    data = data.frame(
      facet = factor(c("Phylogenetic diversity (alpha)", "Phylogenetic RPCA (beta)"), levels = c("Phylogenetic diversity (alpha)", "Phylogenetic RPCA (beta)")),
      x = c(max(alpha_range) - 0.05 * diff(alpha_range), max(beta_range) - 0.05 * diff(beta_range)),
      y = 4
    ),
    aes(x = x, y = y), shape = 23, fill = "green4", size = 2.2, stroke = 0.4, inherit.aes = FALSE
  ) +
  geom_text(
    data = data.frame(
      facet = factor(c("Phylogenetic diversity (alpha)", "Phylogenetic RPCA (beta)"), levels = c("Phylogenetic diversity (alpha)", "Phylogenetic RPCA (beta)")),
      x = c(max(alpha_range) - 0.08 * diff(alpha_range), max(beta_range) - 0.08 * diff(beta_range)),
      y = 4
    ),
    aes(x = x, y = y, label = "FI"), size = 2.8, hjust = 1.5, inherit.aes = FALSE
  ) +
  scale_x_continuous(name = NULL) +
  scale_fill_manual(
    values = c("Alpha" = "#FFA07A", "Beta" = "#7AC5CD", "Frailty" = "green4")
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 5),
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(size = 7),
    panel.grid.major = element_line(linewidth = 0.1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.2),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL)

## Diversity and stability correlations ----

# Calculate means
mean_alpha_diversity <- alpha_diversity %>%
  group_by(SubjectID, metric) %>%
  summarise(across(value, mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = metric, values_from = value, names_prefix = "alpha_diversity_")

mean_beta_stability <- beta_stability %>%
  group_by(SubjectID) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE), .names = "beta_stability_{.col}"))

participants <- participants %>%
  left_join(mean_alpha_diversity) %>%
  left_join(mean_beta_stability)

# Find extremes for mean Faith's PD and mean RPCA stability
participants %>%
  select(SubjectID, alpha_diversity_faith_pd, beta_stability_rpca) %>%
  pivot_longer(cols = -SubjectID, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  filter(value == max(value, na.rm = TRUE) | value == min(value, na.rm = TRUE)) %>%
  arrange(metric, desc(value)) %>%
  print()

# Diversity associations with metadata ----
metrics <- c(
  grep("^alpha_diversity_", names(participants), value = TRUE),
  grep("^beta_stability_", names(participants), value = TRUE)
)

# Binary variables (Wilcoxon test)
map_dfr(c("Sex", "Location"), function(var_name) {
  map_dfr(metrics, function(metric_name) {
    wilcox_test(participants, as.formula(paste(metric_name, "~", var_name))) %>%
      mutate(metric = metric_name, variable = var_name)
  })
})

# Continuous variables (Spearman correlation)
map_dfr(c("Age", "FIScore", "TimeNW"), function(var_name) {
  map_dfr(metrics, function(metric_name) {
    cor.test(participants[[metric_name]], participants[[var_name]],
      method = "spearman", exact = FALSE
    ) %>%
      broom::tidy() %>%
      mutate(metric = metric_name, variable = var_name)
  })
})

cor_results <-
  expand.grid(
    alpha = names(participants) %>% grep("^alpha_diversity_", ., value = TRUE),
    beta = names(participants) %>% grep("^beta_stability_", ., value = TRUE),
    stringsAsFactors = FALSE
  ) %>%
  rowwise() %>%
  mutate(
    cor_test_res = list(
      cor.test(participants[[alpha]], participants[[beta]],
        method = "spearman", exact = FALSE
      ) %>% broom::tidy()
    )
  ) %>%
  unnest(cor_test_res)

write_tsv(cor_results, "results/diversity_stability_correlations.tsv")

## Beta diversity analysis for Panel C and results ----

# Load PCoA coordinates
pcoa_file <- "data/16S/beta_diversity_grouped/rpca_pcoa.tsv"
pcoa_lines <- readLines(pcoa_file)
site_start <- grep("^Site\\t", pcoa_lines)[1]
site_n <- as.numeric(str_split(pcoa_lines[site_start], "\\t")[[1]][2])
pcoa_coords <- read_tsv(pcoa_file, skip = site_start, n_max = site_n, col_names = c("SubjectID", "PC1", "PC2", "PC3"), show_col_types = FALSE)

eigvals_idx <- grep("^Eigvals\\t", pcoa_lines)[1]
eigvals <- as.numeric(str_split(pcoa_lines[eigvals_idx + 1], "\\t")[[1]])
eigvals <- eigvals[!is.na(eigvals)]
prop_idx <- grep("^Proportion explained\\t", pcoa_lines)[1]
prop_explained <- as.numeric(str_split(pcoa_lines[prop_idx + 1], "\\t")[[1]])
prop_explained <- prop_explained[!is.na(prop_explained)]

pc1_eig <- eigvals[1]
pc1_prop <- prop_explained[1] * 100
pc2_eig <- eigvals[2]
pc2_prop <- prop_explained[2] * 100

# Run PERMANOVA
grouped_dms <- list.files("data/16S/beta_diversity_grouped", pattern = "\\_distance.tsv$", full.names = TRUE) %>%
  set_names(nm = basename(.) %>% str_remove("\\_distance.tsv$")) %>%
  map(~ read_tsv(.x, show_col_types = FALSE) %>% rename(sample = 1))

variables_to_test <- c("Age", "Sex", "FIScore", "TimeNW", "Location")
target_samples <- participants$SubjectID

permanova_results <- map_dfr(names(grouped_dms), function(metric_name) {
  dm_df <- grouped_dms[[metric_name]]
  common_samples <- intersect(target_samples, dm_df$sample)

  dm_matrix <- dm_df %>%
    filter(sample %in% common_samples) %>%
    arrange(match(sample, common_samples)) %>%
    select(all_of(common_samples)) %>%
    as.matrix()

  metadata_subset <- participants %>%
    filter(SubjectID %in% common_samples) %>%
    arrange(match(SubjectID, common_samples))

  set.seed(0509)
  form <- as.formula(paste("as.dist(dm_matrix) ~", paste(variables_to_test, collapse = " + ")))
  adonis2(form, data = metadata_subset, permutations = 999, by = "margin") %>%
    broom::tidy() %>%
    filter(!term %in% c("Residual", "Total")) %>%
    mutate(metric = metric_name)
})

write_tsv(permanova_results, "results/beta_permanova_results.tsv")

# Panel B: Faith's PD vs Microbiome Instability (Weighted UniFrac)
p_b_data <- participants %>% filter(!is.na(alpha_diversity_faith_pd), !is.na(beta_stability_weighted_unifrac))
cor_test <- cor.test(p_b_data$alpha_diversity_faith_pd, p_b_data$beta_stability_weighted_unifrac, method = "spearman")

p_b <- ggplot(p_b_data, aes(x = alpha_diversity_faith_pd, y = beta_stability_weighted_unifrac, color = FIScore)) +
  geom_point(alpha = 0.8, size = 1) +
  geom_smooth(method = "lm", color = "black", fill = "grey100", linewidth = 0.4) +
  scale_color_viridis_c(
    option = "plasma", name = "FI",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "left",
      label.position = "bottom",
      barheight = unit(0.12, "cm"),
      barwidth = unit(2.5, "cm"),
      title.vjust = 1,
      title.hjust = 1
    )
  ) +
  annotate("text",
    x = Inf, y = Inf, label = sprintf("rho = %.2f\np = %.3f", cor_test$estimate, cor_test$p.value),
    hjust = 1.1, vjust = 1.35, fontface = "bold", size = 2.5
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1,
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  ) +
  labs(x = "Mean PD", y = "Mean distance (WU)")

# Interaction models (Results only)
interaction_results <- expand.grid(
  alpha = names(participants) %>% grep("^alpha_diversity_", ., value = TRUE),
  beta = names(participants) %>% grep("^beta_stability_", ., value = TRUE),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    model_res = list({
      form <- as.formula(paste0("FIScore ~ ", alpha, " * ", beta))
      m <- glm(form, data = participants, family = quasibinomial)
      broom::tidy(m)
    })
  ) %>%
  unnest(model_res)

write_tsv(interaction_results, "results/frailty_interaction_models.tsv")

# Panel C: RPCA PCoA colored by Location
p_c_permanova <- permanova_results %>%
  filter(metric == "rpca", term == "Location") %>%
  slice(1)

p_c_data <- pcoa_coords %>%
  inner_join(participants %>% select(SubjectID, Location), by = "SubjectID")

p_c <- ggplot(p_c_data, aes(x = PC1, y = PC2, color = Location)) +
  stat_ellipse(linewidth = 0.3, alpha = 0.6) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1,
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 0, 0, 0),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6)
  ) +
  labs(
    x = sprintf("PC1 (%.2f%%)", pc1_prop),
    y = sprintf("PC2 (%.2f%%)", pc2_prop)
  ) +
  annotate("text",
    x = Inf, y = Inf,
    label = sprintf("R2 = %.3f\np = %.3f", p_c_permanova$R2, p_c_permanova$p.value),
    hjust = 1.1, vjust = 1.35, size = 2.5, fontface = "bold"
  )


# Final composite plot
bottom_row <- plot_grid(p_b, p_c, ncol = 2, labels = c("B", "C"), label_size = 10, align = "h", axis = "bt")
composite_plot <- plot_grid(p_a, bottom_row, ncol = 1, labels = c("A", ""), rel_heights = c(2, 1), label_size = 10)

ggsave("results/Main2_Alpha_Beta_Diversity.pdf", composite_plot, width = 5, height = 7)

# PERMANOVA for diseases and medications
covariates <- c(disease_cols, med_cols)
formula_str <- paste("as.dist(dm_matrix) ~", paste(covariates, collapse = " + "))

beta_covariate_results <- map(names(grouped_dms), function(metric_name) {
  dm_df <- grouped_dms[[metric_name]]

  dm_matrix <- dm_df %>%
    filter(sample %in% target_samples) %>%
    arrange(match(sample, target_samples)) %>%
    select(all_of(target_samples)) %>%
    as.matrix()

  # Use by = "margin" to test each covariate while controlling for others
  res <- adonis2(as.formula(formula_str), data = participants, permutations = 999, by = "margin")
  broom::tidy(res) %>%
    mutate(metric = metric_name) %>%
    filter(!term %in% c("Residual", "Total"))
}) %>% bind_rows()

write_tsv(beta_covariate_results, "results/beta_covariates_permanova_results.tsv")

## Group comparisons: Diseases and Medications ----

disease_cols <- c(
  "anemia", "arthritis", "atherosclerosis", "copd", "dementia", "depression",
  "diabetes", "dyslipidemia", "hypertension", "hypothyroidism", "osteoporosis", "renal_disease"
)

med_cols <- c(
  "ace_inhibitors", "acetaminophen", "aaps", "benzodiazepines", "beta_blockers",
  "antiplatelets", "ccbs", "diuretics", "insulin", "iron_supplements", "laxatives",
  "levothyroxine", "metformin", "nitrates", "opioids", "ppis", "ssris", "statins",
  "tcas", "vitb12", "vitd"
)

metrics <- c(
  grep("^alpha_diversity_", names(participants), value = TRUE)
)

run_diversity_glm <- function(data, alpha_metric) {
  covariates <- names(data)[which(names(data) == "anemia"):which(names(data) == "vitd")]
  formula_str <- paste(alpha_metric, "~", paste(covariates, collapse = " + "))

  model <- glm(as.formula(formula_str), data = data)

  # Extract coefficients
  summary_res <- as.data.frame(summary(model)$coefficients)
  colnames(summary_res) <- c("estimate", "std_error", "statistic", "p_value")
  summary_res$term <- rownames(summary_res)
  summary_res$alpha_metric <- alpha_metric # Add the alpha metric to the results

  return(summary_res)
}

alpha_metrics_for_glm <- grep("^alpha_diversity_", names(participants), value = TRUE)

diversity_results <- map_dfr(alpha_metrics_for_glm, function(metric_name) {
  run_diversity_glm(participants, metric_name)
}) %>%
  filter(term != "(Intercept)")

write_tsv(diversity_results, "results/diversity_glm_results.tsv")

# Bar chart of associations
diseases <- c(
  "anemia", "arthritis", "atherosclerosis", "copd", "dementia",
  "depression", "diabetes", "dyslipidemia", "hypertension",
  "hypothyroidism", "osteoporosis", "renal_disease"
)

diversity_associations_figure <- diversity_results %>%
  filter(alpha_metric == "alpha_diversity_faith_pd") %>%
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
    TRUE ~ term # Default to original term if no match
  )) %>%
  mutate(type = ifelse(term %in% diseases, "Comorbidity", "Medication")) %>%
  mutate(p_signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    # p_value < 0.1 ~ "•",
    TRUE ~ ""
  )) %>%
  # Order by estimate descending
  mutate(term_clean = fct_reorder(term_clean, estimate, .desc = TRUE)) %>%
  # Position labels slightly above or below the bars
  mutate(
    label_pos = ifelse(estimate >= 0, estimate + 0.25, estimate - 0.25),
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
    y = "Estimate (PD)",
    fill = NULL
  )

ggsave(
  "results/Main2_Alpha_Associations.pdf",
  diversity_associations_figure,
  width = 6, height = 3
)
