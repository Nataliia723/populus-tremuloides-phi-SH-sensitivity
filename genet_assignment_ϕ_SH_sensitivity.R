# ============================================================
# GENET ASSIGNMENT — ϕ AND SH SENSITIVITY ANALYSES
# Populus tremuloides
#
# Purpose:
# 1. Evaluate sensitivity of genet assignments to alternative
#    pairwise relatedness (ϕ) and shared heterozygosity (SH)
#    thresholds.
# 2. Compare within-genet and between-genet distributions.
# 3. Verify that samples not assigned to the 25 multi-sample
#    genets are genetically distinct from those genets and
#    from one another.
#
# Final dataset:
#   421 samples
#   25 multi-sample genets = 413 samples
#   8 singleton genets
#   33 genets in total
# ============================================================


# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------

library(openxlsx)
library(dplyr)
library(tidyr)


# ------------------------------------------------------------
# 2. FILE PATHS
# ------------------------------------------------------------

# Change only this path if the shared folder is moved
base_dir <- "C:/Users/Nataliia/Desktop/MAÎTRISE NK/ANALYSES_PEUPLIER/ANALYSES_ROOS"

clone_file <- file.path(
  base_dir,
  "clonality-20250630T150532Z-1-001",
  "clonal_groups.xlsx"
)

sh_file <- file.path(
  base_dir,
  "clonality-20250630T150532Z-1-001",
  "SH_index_table.xlsx"
)

phi_file <- file.path(
  base_dir,
  "relatedness-20250624T012216Z-1-001",
  "relatedness_output.xlsx"
)


# ------------------------------------------------------------
# 3. IMPORT DATA
# ------------------------------------------------------------

clone_groups <- read.xlsx(clone_file)

sh_table <- read.xlsx(sh_file)

relatedness <- read.xlsx(phi_file)


# ------------------------------------------------------------
# 4. BASIC CHECKS
# ------------------------------------------------------------

# Number of samples assigned to multi-sample genets
n_distinct(clone_groups$Sample)

# Number of multi-sample genets
n_distinct(clone_groups$Clone_Group)

# Expected:
# 413 samples
# 25 multi-sample genets


# ------------------------------------------------------------
# 5. PREPARE PHI DATA
# ------------------------------------------------------------

# Remove self-comparisons
rel_pairs <- relatedness %>%
  filter(INDV1 != INDV2)

# Create an unordered pair ID so reciprocal comparisons
# (A-B and B-A) are counted only once
rel_unique <- rel_pairs %>%
  mutate(
    Pair_ID = ifelse(
      INDV1 < INDV2,
      paste(INDV1, INDV2, sep = "__"),
      paste(INDV2, INDV1, sep = "__")
    )
  ) %>%
  distinct(Pair_ID, .keep_all = TRUE)


# ------------------------------------------------------------
# 6. ADD GENET IDENTITIES TO PHI PAIRS
# ------------------------------------------------------------

clone_lookup_1 <- clone_groups %>%
  select(Sample, Clone_Group) %>%
  rename(
    INDV1 = Sample,
    Genet1 = Clone_Group
  )

clone_lookup_2 <- clone_groups %>%
  select(Sample, Clone_Group) %>%
  rename(
    INDV2 = Sample,
    Genet2 = Clone_Group
  )

phi_assigned <- rel_unique %>%
  inner_join(clone_lookup_1, by = "INDV1") %>%
  inner_join(clone_lookup_2, by = "INDV2") %>%
  mutate(
    Comparison = ifelse(
      Genet1 == Genet2,
      "Within_genet",
      "Between_genets"
    )
  )


# Check number of unique comparisons
table(phi_assigned$Comparison)

# Expected:
# Within_genet    = 9,109
# Between_genets  = 75,969


# ------------------------------------------------------------
# 7. PHI DISTRIBUTIONS
# ------------------------------------------------------------

phi_distribution <- phi_assigned %>%
  group_by(Comparison) %>%
  summarise(
    N = n(),
    Minimum = min(RELATEDNESS_PHI, na.rm = TRUE),
    Q01 = quantile(RELATEDNESS_PHI, 0.01, na.rm = TRUE),
    Q05 = quantile(RELATEDNESS_PHI, 0.05, na.rm = TRUE),
    Median = median(RELATEDNESS_PHI, na.rm = TRUE),
    Q95 = quantile(RELATEDNESS_PHI, 0.95, na.rm = TRUE),
    Q99 = quantile(RELATEDNESS_PHI, 0.99, na.rm = TRUE),
    Maximum = max(RELATEDNESS_PHI, na.rm = TRUE),
    .groups = "drop"
  )

phi_distribution

# Expected ranges:
# Within genets:   0.378 to 0.498
# Between genets: -0.250 to -0.007


# ------------------------------------------------------------
# 8. PHI THRESHOLD SENSITIVITY
# ------------------------------------------------------------

phi_thresholds <- c(
  0.350,
  0.375,
  0.400,
  0.425,
  0.450,
  0.475,
  0.490
)

phi_sensitivity <- bind_rows(
  lapply(
    phi_thresholds,
    function(threshold) {
      
      within <- phi_assigned %>%
        filter(Comparison == "Within_genet")
      
      between <- phi_assigned %>%
        filter(Comparison == "Between_genets")
      
      data.frame(
        Metric = "Pairwise relatedness (phi)",
        Threshold = threshold,
        
        Within_total = nrow(within),
        
        Within_above = sum(
          within$RELATEDNESS_PHI >= threshold,
          na.rm = TRUE
        ),
        
        Within_percent = 100 * mean(
          within$RELATEDNESS_PHI >= threshold,
          na.rm = TRUE
        ),
        
        Between_total = nrow(between),
        
        Between_above = sum(
          between$RELATEDNESS_PHI >= threshold,
          na.rm = TRUE
        ),
        
        Between_percent = 100 * mean(
          between$RELATEDNESS_PHI >= threshold,
          na.rm = TRUE
        )
      )
    }
  )
)

phi_sensitivity


# ------------------------------------------------------------
# 9. PREPARE SH DATA
# ------------------------------------------------------------

# SH_index_table already contains each unordered pair once.

sh_lookup_1 <- clone_groups %>%
  select(Sample, Clone_Group) %>%
  rename(
    X1 = Sample,
    Genet1 = Clone_Group
  )

sh_lookup_2 <- clone_groups %>%
  select(Sample, Clone_Group) %>%
  rename(
    X2 = Sample,
    Genet2 = Clone_Group
  )

sh_assigned <- sh_table %>%
  inner_join(sh_lookup_1, by = "X1") %>%
  inner_join(sh_lookup_2, by = "X2") %>%
  mutate(
    Comparison = ifelse(
      Genet1 == Genet2,
      "Within_genet",
      "Between_genets"
    )
  )


# Check number of comparisons
table(sh_assigned$Comparison)

# Expected:
# Within_genet    = 9,109
# Between_genets  = 75,969


# ------------------------------------------------------------
# 10. SH DISTRIBUTIONS
# ------------------------------------------------------------

sh_distribution <- sh_assigned %>%
  group_by(Comparison) %>%
  summarise(
    N = n(),
    Minimum = min(SH, na.rm = TRUE),
    Q01 = quantile(SH, 0.01, na.rm = TRUE),
    Q05 = quantile(SH, 0.05, na.rm = TRUE),
    Median = median(SH, na.rm = TRUE),
    Q95 = quantile(SH, 0.95, na.rm = TRUE),
    Q99 = quantile(SH, 0.99, na.rm = TRUE),
    Maximum = max(SH, na.rm = TRUE),
    .groups = "drop"
  )

sh_distribution

# Expected ranges:
# Within genets:  0.881 to 0.998
# Between genets: 0.163 to 0.293


# ------------------------------------------------------------
# 11. SH THRESHOLD SENSITIVITY
# ------------------------------------------------------------

sh_thresholds <- c(
  0.800,
  0.850,
  0.900,
  0.925,
  0.950,
  0.975,
  0.980
)

sh_sensitivity <- bind_rows(
  lapply(
    sh_thresholds,
    function(threshold) {
      
      within <- sh_assigned %>%
        filter(Comparison == "Within_genet")
      
      between <- sh_assigned %>%
        filter(Comparison == "Between_genets")
      
      data.frame(
        Metric = "Shared heterozygosity (SH)",
        Threshold = threshold,
        
        Within_total = nrow(within),
        
        Within_above = sum(
          within$SH >= threshold,
          na.rm = TRUE
        ),
        
        Within_percent = 100 * mean(
          within$SH >= threshold,
          na.rm = TRUE
        ),
        
        Between_total = nrow(between),
        
        Between_above = sum(
          between$SH >= threshold,
          na.rm = TRUE
        ),
        
        Between_percent = 100 * mean(
          between$SH >= threshold,
          na.rm = TRUE
        )
      )
    }
  )
)

sh_sensitivity


# ------------------------------------------------------------
# 12. FINAL SENSITIVITY TABLE
# ------------------------------------------------------------

sensitivity_table <- bind_rows(
  phi_sensitivity,
  sh_sensitivity
) %>%
  select(
    Metric,
    Threshold,
    Within_above,
    Within_percent,
    Between_above,
    Between_percent
  ) %>%
  mutate(
    Within_percent = round(Within_percent, 2),
    Between_percent = round(Between_percent, 2)
  )

sensitivity_table


# ------------------------------------------------------------
# 13. IDENTIFY SAMPLES NOT ASSIGNED TO THE 25 MULTI-SAMPLE
#     GENETS
# ------------------------------------------------------------

# Obtain all 421 samples from SH table
all_samples <- sort(
  unique(
    c(
      sh_table$X1,
      sh_table$X2
    )
  )
)

assigned_samples <- clone_groups$Sample

unassigned_samples <- setdiff(
  all_samples,
  assigned_samples
)

unassigned_samples
length(unassigned_samples)

# Expected:
# ptrem_014
# ptrem_080
# ptrem_174
# ptrem_293
# ptrem_349
# ptrem_410
# ptrem_539
# ptrem_540
#
# n = 8


# ------------------------------------------------------------
# 14. COMPARE EACH SINGLETON WITH THE 25 MULTI-SAMPLE GENETS:
#     PHI
# ------------------------------------------------------------

clone_lookup <- clone_groups %>%
  select(Sample, Clone_Group)

phi_singletons <- rel_unique %>%
  mutate(
    Singleton = case_when(
      INDV1 %in% unassigned_samples &
        INDV2 %in% assigned_samples ~ INDV1,
      
      INDV2 %in% unassigned_samples &
        INDV1 %in% assigned_samples ~ INDV2,
      
      TRUE ~ NA_character_
    ),
    
    Assigned_sample = case_when(
      INDV1 %in% unassigned_samples &
        INDV2 %in% assigned_samples ~ INDV2,
      
      INDV2 %in% unassigned_samples &
        INDV1 %in% assigned_samples ~ INDV1,
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Singleton)) %>%
  left_join(
    clone_lookup,
    by = c("Assigned_sample" = "Sample")
  )

best_phi_singleton <- phi_singletons %>%
  group_by(Singleton) %>%
  slice_max(
    RELATEDNESS_PHI,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    Singleton,
    Best_PHI_sample = Assigned_sample,
    Best_PHI_genet = Clone_Group,
    Max_PHI = RELATEDNESS_PHI
  )


# ------------------------------------------------------------
# 15. COMPARE EACH SINGLETON WITH THE 25 MULTI-SAMPLE GENETS:
#     SH
# ------------------------------------------------------------

sh_singletons <- sh_table %>%
  mutate(
    Singleton = case_when(
      X1 %in% unassigned_samples &
        X2 %in% assigned_samples ~ X1,
      
      X2 %in% unassigned_samples &
        X1 %in% assigned_samples ~ X2,
      
      TRUE ~ NA_character_
    ),
    
    Assigned_sample = case_when(
      X1 %in% unassigned_samples &
        X2 %in% assigned_samples ~ X2,
      
      X2 %in% unassigned_samples &
        X1 %in% assigned_samples ~ X1,
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Singleton)) %>%
  left_join(
    clone_lookup,
    by = c("Assigned_sample" = "Sample")
  )

best_sh_singleton <- sh_singletons %>%
  group_by(Singleton) %>%
  slice_max(
    SH,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    Singleton,
    Best_SH_sample = Assigned_sample,
    Best_SH_genet = Clone_Group,
    Max_SH = SH
  )


# ------------------------------------------------------------
# 16. SINGLETON SUMMARY
# ------------------------------------------------------------

singleton_summary <- best_phi_singleton %>%
  full_join(
    best_sh_singleton,
    by = "Singleton"
  ) %>%
  mutate(
    PHI_ge_045 = Max_PHI >= 0.45,
    SH_ge_090 = Max_SH >= 0.90
  ) %>%
  arrange(Singleton)

singleton_summary

# Expected:
# No singleton should have:
# PHI >= 0.45
# or
# SH >= 0.90


# ------------------------------------------------------------
# 17. COMPARE THE 8 SINGLETONS WITH ONE ANOTHER
# ------------------------------------------------------------

phi_singleton_pairs <- rel_unique %>%
  filter(
    INDV1 %in% unassigned_samples,
    INDV2 %in% unassigned_samples
  ) %>%
  select(
    INDV1,
    INDV2,
    RELATEDNESS_PHI
  ) %>%
  arrange(desc(RELATEDNESS_PHI))

sh_singleton_pairs <- sh_table %>%
  filter(
    X1 %in% unassigned_samples,
    X2 %in% unassigned_samples
  ) %>%
  select(
    X1,
    X2,
    SH
  ) %>%
  arrange(desc(SH))


# Check: 8 choose 2 = 28 comparisons
stopifnot(nrow(phi_singleton_pairs) == 28)
stopifnot(nrow(sh_singleton_pairs) == 28)

phi_singleton_pairs
sh_singleton_pairs


# Check whether any singleton pair reaches the thresholds
any(phi_singleton_pairs$RELATEDNESS_PHI >= 0.45)
any(sh_singleton_pairs$SH >= 0.90)

# Expected:
# FALSE
# FALSE


# Highest similarity among singleton pairs
phi_singleton_pairs %>%
  slice_max(
    RELATEDNESS_PHI,
    n = 1,
    with_ties = FALSE
  )

sh_singleton_pairs %>%
  slice_max(
    SH,
    n = 1,
    with_ties = FALSE
  )

# Expected highest pair:
# ptrem_539 / ptrem_540
#
# phi = 0.364288
# SH  = 0.7895
#
# These values remain below the high-similarity thresholds.


# ------------------------------------------------------------
# 18. FINAL CONSISTENCY CHECKS
# ------------------------------------------------------------

n_multi_sample_genets <- n_distinct(
  clone_groups$Clone_Group
)

n_singleton_genets <- length(
  unassigned_samples
)

n_total_genets <- n_multi_sample_genets +
  n_singleton_genets

cat(
  "\nMulti-sample genets:",
  n_multi_sample_genets,
  "\nSingleton genets:",
  n_singleton_genets,
  "\nTOTAL GENETS:",
  n_total_genets,
  "\n"
)

stopifnot(n_multi_sample_genets == 25)
stopifnot(n_singleton_genets == 8)
stopifnot(n_total_genets == 33)


# ------------------------------------------------------------
# 19. EXPORT RESULTS
# ------------------------------------------------------------

output_file <- file.path(
  base_dir,
  "genet_assignment_phi_SH_sensitivity.xlsx"
)

write.xlsx(
  list(
    "Sensitivity_table" = sensitivity_table,
    "PHI_distribution" = phi_distribution,
    "SH_distribution" = sh_distribution,
    "Singleton_summary" = singleton_summary,
    "Singleton_pairs_PHI" = phi_singleton_pairs,
    "Singleton_pairs_SH" = sh_singleton_pairs
  ),
  file = output_file,
  overwrite = TRUE
)

cat(
  "\nResults exported to:\n",
  output_file,
  "\n"
)


# ------------------------------------------------------------
# 20. SESSION INFORMATION
# ------------------------------------------------------------

sessionInfo()