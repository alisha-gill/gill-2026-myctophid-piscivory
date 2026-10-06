# 
# Gill et al. (in review) — A global meta-analysis of piscivory among
# lanternfish (Myctophidae)
# 
# Description:
#   Fits Bayesian mixed-effects logistic regression models to lanternfish
#   diet data and produces all manuscript figures, tables, and supplementary
#   materials reported in the paper.
#
# Input:
#   data/Myctophid_diet.csv          — main diet dataset
#   data/Myctophid_traits.csv        — species-level functional traits
#   data/dist_to_shore_cache.csv     — pre-computed coastline distances
#                                      (generated on first run; reload to skip)
#
# Output:
#   models/   — fitted brms model objects and cached LOO results (.rds)
#   figures/  — all manuscript and supplementary figures (.pdf and .tiff)
#   tables/   — all manuscript and supplementary tables (.docx)
#
# Analysis structure:
#   Section 2  — data preparation and collinearity checks
#   Section 3A — candidate models; model selection in two stages:
#                (1) fixed effects on a shared complete-case subset (n = 190)
#                (2) random-effects structure among the top-ranked genus x size
#                    models on their common data (n = 226);
#                best model reported from its fit on complete data
#   Section 3B — trait-based models; traits vs. genus x size on a shared
#                trait subset; best trait model reported on complete data
#   Section 3C — sensitivity analyses of the best-supported model
#   Section 4  — tables and figures
#
# NOTE on spatial data:
#   Distance to shore is computed by querying an external coastline service
#   (rnaturalearth). This query can be slow and server-dependent. After the
#   first successful run, the result is cached as a .csv file in data/ so
#   that re-running the script does not require re-querying. If you have
#   received this cache file with the repository, the download block will be
#   skipped automatically. 
#
# NOTE on cached models:
#   fit_or_load() and loo_or_load() load an existing .rds file whenever one
#   exists, regardless of the formula or data passed. If a model's formula or
#   data is changed, delete its .rds file (and its loo_*.rds file) so it is
#   refit.

# 
# SECTION 1: SETUP##############################################################
# 

## Clear environment -----------------------------------------------------------
rm(list = ls())

## Libraries -------------------------------------------------------------------

# Core data wrangling & visualisation
library(readr)
library(tidyverse)
library(broom.mixed)
library(patchwork)

# Spatial
library(maps)
library(sf)
library(rnaturalearth)

# Bayesian modelling & diagnostics
library(brms)
library(tidybayes)
library(bayesplot)
library(loo)
library(car)

# Tables & output
library(flextable)
library(officer)

## Output directories ----------------------------------------------------------
# Create all output directories up front 
dir.create("models",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)
dir.create("data",    showWarnings = FALSE)

# 
# SECTION 2: DATA###############################################################
# 

## 2a: Diet data ---------------------------------------------------------------
diet <- read_csv("Myctophid_diet_expanded.csv")

# restrict to adult specimens with explicitly-reported fish presence/absence
diet <- diet %>%
  filter(explicit == "yes")

# Capitalise genus names
diet <- diet %>%
  mutate(genusconfirmed = str_to_sentence(genusconfirmed))

# Set reference level for genus (Diaphus: most-sampled, intermediate piscivory)
diet$genusconfirmed <- relevel(factor(diet$genusconfirmed), ref = "Diaphus")

# Ensure ocean and environment are factors
diet$ocean       <- factor(diet$ocean)
diet$environment <- factor(diet$environment)

# Identification method, collapsed to two categories for use as a covariate.
diet <- diet %>%
  mutate(
    method_clean = case_when(
      method == "gut" ~ "visual_only",
      method %in% c("dna", "sia", "gut and dna", "gut and sia") ~ "molecular_isotope_inclusive",
      TRUE ~ NA_character_
    ),
    method_clean = factor(method_clean, levels = c("visual_only", "molecular_isotope_inclusive"))
  )

message(sprintf("method_clean coverage: %d / %d (%.0f%%)",
                sum(!is.na(diet$method_clean)), nrow(diet),
                100 * mean(!is.na(diet$method_clean))))

## 2b: Traits data -------------------------------------------------------------
traits <- read_csv("Myctophid_traits.csv")

# Clean and derive trait variables from the trait table. Matched at the
# species level only (spconfirmed <-> species) 

traits_clean <- traits %>%
  mutate(
    sp_lower = str_squish(str_to_lower(species)),
    genus_lower = str_squish(str_to_lower(genus)),
    jaw_class_sp = factor(jaw_class_sp, levels = c("short", "moderate", "long")),
    raptorial_sp = as.integer(raptorial_sp),
    # Centrobranchus lacks a comparable gill-raker structure entirely (lath-like
    # rakers are replaced by patches of small teeth) -- this is a different
    # morphology, not "zero rakers" on the same continuum as every other
    # genus. Set to NA here so it's
    # excluded from the continuous gill-raker model.
    raker_structure = if_else(genus_lower == "centrobranchus",
                              "tooth-patch", "typical-raker"),
    gillrakers_sp = if_else(genus_lower == "centrobranchus",
                            NA_real_, gillrakers_sp)
  )

message(sprintf("Centrobranchus: %d species flagged with tooth-patch raker structure (excluded from continuous gill-raker model, retained as raker_structure)",
                sum(traits_clean$raker_structure == "tooth-patch")))

# Match diet observations to trait data by species (spconfirmed <-> species) 

diet <- diet %>%
  mutate(sp_fixed = str_squish(str_to_lower(spconfirmed)))

diet <- diet %>%
  left_join(
    traits_clean %>%
      select(sp_lower, gillrakers_sp, jaw_class_sp, raptorial_sp, raker_structure),
    by = c("sp_fixed" = "sp_lower")
  ) %>%
  rename(
    gillrakers      = gillrakers_sp,
    jaw_class       = jaw_class_sp,
    raptorial_teeth = raptorial_sp
  ) %>%
  mutate(
    trait_source    = if_else(!is.na(jaw_class) | !is.na(raptorial_teeth),
                              "species-matched", "unmatched"),
    gillrakers_z    = as.numeric(scale(gillrakers)),
    raptorial_teeth = factor(raptorial_teeth, levels = c(0, 1),
                             labels = c("villiform-only", "raptorial"))
  )

## 2c: Map data (for Figure 1) -------------------------------------------------
world <- map_data("world")

## 2d: Distance to shore -------------------------------------------------------
# Pre-computed distances are cached after the first run to avoid re-querying

dist_cache <- "data/dist_to_shore_cache.csv"

if (file.exists(dist_cache)) {
  
  message("Loading cached distance-to-shore data from: ", dist_cache)
  dist_cached <- read_csv(dist_cache)
  diet <- diet %>%
    left_join(dist_cached, by = c("ref", "species", "yearstart", "locality"))
  
} else {
  
  message("Computing distance to shore — this may take 1-2 minutes...")
  
  diet_with_coords <- diet %>%
    filter(!is.na(longmid) & !is.na(latmid))
  
  fish_points <- st_as_sf(diet_with_coords,
                          coords = c("longmid", "latmid"),
                          crs    = 4326)
  
  coastline     <- ne_coastline(scale = "medium", returnclass = "sf")
  distances     <- st_distance(fish_points, coastline)
  min_distances <- apply(distances, 1, min)
  
  diet_with_coords$dist_to_shore <- min_distances
  
  # Cache result
  dist_cached <- diet_with_coords %>%
    select(ref, species, yearstart, locality, dist_to_shore)
  write_csv(dist_cached, dist_cache)
  message("Distance-to-shore cached to: ", dist_cache)
  
  diet <- diet %>%
    left_join(dist_cached, by = c("ref", "species", "yearstart", "locality"))
}

## 2e: Standardise predictors --------------------------------------------------
# Body size (sizez) is standardised globally across all genera using maximum
# standard length, placing all genera on a common scale. 
# All other continuous predictors are also standardised globally.

diet <- diet %>%
  mutate(
    sizez                = as.numeric(scale(maxsize)),
    depthz               = as.numeric(scale(maxdepth)),
    dist_to_shore_km     = as.numeric(dist_to_shore) / 1000,
    dist_to_shore_scaled = as.numeric(scale(dist_to_shore_km))
  )

# Pre-compute size back-transform constants (used in Figures 4 & 5 axis labels)
size_mean <- mean(diet$maxsize, na.rm = TRUE)
size_sd   <- sd(diet$maxsize,   na.rm = TRUE)

# Stomach-count weighting columns (used in the weighting sensitivity, 3C-ii)
diet <- diet %>%
  mutate(
    nstomachs_weight    = ifelse(is.na(nstomachs), 1, nstomachs),
    nstomachs_logweight = ifelse(is.na(nstomachs), log1p(1), log1p(nstomachs))
  )
stopifnot(all(c("nstomachs_weight", "nstomachs_logweight") %in% names(diet)))

## 2f: Collinearity checks -----------------------------------------------------
# Traits only add information if they're not just re-deriving `sizez` or
# each other. Checked pairwise below, plus a formal VIF check on the
# combined fixed-effects specification.

diet_traits_complete <- diet %>%
  filter(!is.na(gillrakers_z), !is.na(jaw_class), !is.na(raptorial_teeth), !is.na(sizez))

message(sprintf("Complete cases for collinearity checks: %d / %d",
                nrow(diet_traits_complete), nrow(diet)))

# (a) Jaw class vs. body size
message("\n--- (a) Jaw class vs. body size (maxsize, mm) ---")
diet_traits_complete %>%
  group_by(jaw_class) %>%
  summarize(n = n(), mean_maxsize = mean(maxsize, na.rm = TRUE),
            sd_maxsize = sd(maxsize, na.rm = TRUE), .groups = "drop") %>%
  print()

aov_jaw_size <- aov(maxsize ~ jaw_class, data = diet_traits_complete)
print(summary(aov_jaw_size))
cat("Eta-squared:", summary.lm(aov_jaw_size)$r.squared, "\n")

# (b) Gill rakers vs. body size
message("\n--- (b) Gill rakers vs. body size ---")
print(cor.test(diet_traits_complete$gillrakers, diet_traits_complete$maxsize))

# (c) Gill rakers vs. jaw class
message("\n--- (c) Gill rakers vs. jaw class ---")
diet_traits_complete %>%
  group_by(jaw_class) %>%
  summarize(n = n(), mean_gr = mean(gillrakers, na.rm = TRUE),
            sd_gr = sd(gillrakers, na.rm = TRUE), .groups = "drop") %>%
  print()
print(summary(aov(gillrakers ~ jaw_class, data = diet_traits_complete)))

# (d) Raptorial teeth vs. jaw class
message("\n--- (d) Raptorial teeth vs. jaw class ---")
tab_teeth_jaw <- table(diet_traits_complete$jaw_class, diet_traits_complete$raptorial_teeth)
print(tab_teeth_jaw)
print(chisq.test(tab_teeth_jaw))

# (e) Raptorial teeth vs. body size
message("\n--- (e) Raptorial teeth vs. body size ---")
print(t.test(maxsize ~ raptorial_teeth, data = diet_traits_complete))

# (f) GVIF on the trait specification (traits + size). No model includes genus
# together with traits: traits are almost entirely nested within genus (see
# nesting check below), so a GVIF including genus is not informative.
vif_check <- glm(fish == "yes" ~ sizez + gillrakers_z + jaw_class + raptorial_teeth,
                 data = diet_traits_complete, family = binomial)
message("\n--- (f) GVIF check on trait specification ---")
print(vif(vif_check))

# How far traits are nested within genus
diet_traits_complete %>%
  group_by(genusconfirmed) %>%
  summarise(jaw_levels = n_distinct(jaw_class),
            tooth_levels = n_distinct(raptorial_teeth), .groups = "drop") %>%
  summarise(genera = n(),
            single_jaw_class = sum(jaw_levels == 1),
            single_tooth_type = sum(tooth_levels == 1)) %>%
  print()

# (g) Genus vs. body size -- can each genus inform its own size slope?
message("\n--- (g) Genus vs. body size ---")
size_by_genus <- diet %>%
  filter(!is.na(sizez)) %>%
  group_by(genusconfirmed) %>%
  summarise(n = n(), n_sizes = n_distinct(maxsize),
            sd_sizez = round(sd(sizez), 2), .groups = "drop") %>%
  arrange(n_sizes, n)
print(size_by_genus, n = Inf)
cat("Eta-squared (size ~ genus):",
    round(summary(lm(sizez ~ genusconfirmed, data = diet))$r.squared, 3), "\n")

# (h) Genus vs. latitude / ocean -- confounding of genus with region
message("\n--- (h) Genus vs. latitude and ocean ---")
cat("Eta-squared (latitude ~ genus):",
    round(summary(lm(latmid ~ genusconfirmed, data = diet))$r.squared, 3), "\n")
tab_genus_ocean <- table(diet$genusconfirmed, diet$ocean)
print(tab_genus_ocean)
chi <- suppressWarnings(chisq.test(tab_genus_ocean))   # sparse cells: use descriptively
cat("Cramer's V (genus x ocean):",
    round(sqrt(chi$statistic / (sum(tab_genus_ocean) * (min(dim(tab_genus_ocean)) - 1))), 3), "\n")

# (i) Environmental covariates against each other (MultiModelF specification)
message("\n--- (i) Environmental covariates ---")
env_complete <- diet %>%
  filter(!is.na(latmid), !is.na(depthz), !is.na(dist_to_shore_scaled),
         !is.na(season), !is.na(sizez))
print(round(cor(env_complete[, c("sizez", "latmid", "depthz", "dist_to_shore_scaled")]), 2))
print(vif(glm(fish == "yes" ~ sizez + latmid + season + depthz + dist_to_shore_scaled,
              data = env_complete, family = binomial)))

# 
# SECTION 3A: CANDIDATE MODELS AND MODEL SELECTION##############################
# 

## 3A-i: Set up ----------------------------------------------------------------

# Prior specifications 
# priorA = moderate (main analysis)
# priorB = diffuse  (sensitivity analysis)
# priorC = tight    (sensitivity analysis)

priorA <- c(prior(normal(0, 1),   class = "b"),
            prior(cauchy(0, 1),   class = "sd"))

priorB <- c(prior(normal(0, 5),   class = "b"),
            prior(cauchy(0, 2),   class = "sd"))

priorC <- c(prior(normal(0, 0.5), class = "b"),
            prior(cauchy(0, 0.5), class = "sd"))

# Set brm() defaults 
brm_defaults <- list(
  family    = bernoulli(link = "logit"),
  data      = diet,
  chains    = 4,
  iter      = 4000,
  warmup    = 2000,
  control   = list(adapt_delta = 0.95, max_treedepth = 12),
  cores     = 4,
  save_pars = save_pars(all = TRUE),
  seed      = 1234
)

# fit_or_load() 
# Fits a brms model and saves it as an .rds file, or loads the saved file if it
# already exists. This means the models only need to be fitted once — re-running
# the script after the initial fit takes seconds rather than hours.

fit_or_load <- function(file, formula, data = brm_defaults$data, prior = priorA) {
  if (file.exists(file)) {
    readRDS(file)
  } else {
    m <- do.call(brm, c(list(formula = formula, data = data, prior = prior),
                        brm_defaults[setdiff(names(brm_defaults), "data")]))
    saveRDS(m, file)
    m
  }
}

# loo_or_load()
# moment_match = TRUE is only applied where the unstabilised loo() actually
# flags influential points (Pareto k > 0.7). Results are cached like the model
# fits themselves, so this only runs once rather than on every script
# execution. cores = 4 parallelises the pointwise computation.

loo_or_load <- function(file, model, cores = 4) {
  if (file.exists(file)) return(readRDS(file))
  l <- loo(model, cores = cores)
  if (any(l$diagnostics$pareto_k > 0.7)) {
    message("Refitting with moment matching: ", file)
    l <- loo(model, moment_match = TRUE, cores = cores)
  }
  saveRDS(l, file)
  l
}

## 3A-ii: Univariable models ---------------------------------------------------
GeneraModel   <- fit_or_load("models/GeneraModel.rds",
                             fish ~ genusconfirmed + (1 | ref))
SizeModel      <- fit_or_load("models/SizeModel.rds",
                              fish ~ sizez + (1 | ref))
SeasonModel    <- fit_or_load("models/SeasonModel.rds",
                              fish ~ season + (1 | ref))
LocationModel  <- fit_or_load("models/LocationModel.rds",
                              fish ~ latmid + (1 | ref))
DepthModel     <- fit_or_load("models/DepthModel.rds",
                              fish ~ depthz + (1 | ref))
DistShoreModel <- fit_or_load("models/DistShoreModel.rds",
                              fish ~ dist_to_shore_scaled + (1 | ref))
MethodModel    <- fit_or_load("models/MethodModel.rds",
                              fish ~ method_clean + (1 | ref))

## 3A-iii: Multivariable models ------------------------------------------------
MultiModelA <- fit_or_load("models/MultiModelA.rds",
                           fish ~ sizez + genusconfirmed + (1 | ref))
MultiModelB <- fit_or_load("models/MultiModelB.rds",
                           fish ~ sizez + genusconfirmed + latmid + (1 | ref))
MultiModelC <- fit_or_load("models/MultiModelC.rds",
                           fish ~ sizez + genusconfirmed + season + (1 | ref))
MultiModelD <- fit_or_load("models/MultiModelD.rds",
                           fish ~ sizez + genusconfirmed + dist_to_shore_scaled + (1 | ref))
MultiModelE <- fit_or_load("models/MultiModelE.rds",
                           fish ~ sizez + genusconfirmed + depthz + (1 | ref))
MultiModelF <- fit_or_load("models/MultiModelF.rds",
                           fish ~ sizez + genusconfirmed + latmid + season +
                             depthz + dist_to_shore_scaled + (1 | ref))
MultiModelG <- fit_or_load("models/MultiModelG.rds",
                           fish ~ sizez + genusconfirmed + method_clean + (1 | ref))

## 3A-iv: Interaction and random-slope models ----------------------------------
MultiModelA_interaction <- fit_or_load("models/MultiModelA_interaction.rds",
                                       fish ~ sizez * genusconfirmed + (1 | ref))

MultiModelA_interaction_method <- fit_or_load("models/MultiModelA_interaction_method.rds",
                                              fish ~ sizez * genusconfirmed + method_clean + (1 | ref))

MultiModelA_rs <- fit_or_load("models/MultiModelA_rs.rds",
                              fish ~ sizez * genusconfirmed + (1 + sizez | ref))

MultiModelA_rsnc <- fit_or_load("models/MultiModelA_rsnc.rds",
                                fish ~ sizez * genusconfirmed + (1 + sizez || ref))

MultiModelA_rs_method <- fit_or_load("models/MultiModelA_rs_method.rds",
                                     fish ~ sizez * genusconfirmed + method_clean + (1 + sizez | ref))

## 3A-v: LOO on each model's own complete data ---------------------------------
# These full-data LOO results are reused for the best model's diagnostics
# (Section 4) and the sensitivity analyses (Section 3C). They are NOT used
# for model selection, since the models are fit on different n (see 3A-vi).

model_list_all <- list(
  GeneraModel = GeneraModel, SizeModel = SizeModel, SeasonModel = SeasonModel,
  LocationModel = LocationModel, DepthModel = DepthModel, DistShoreModel = DistShoreModel,
  MethodModel = MethodModel,
  MultiModelA = MultiModelA, MultiModelB = MultiModelB, MultiModelC = MultiModelC,
  MultiModelD = MultiModelD, MultiModelE = MultiModelE, MultiModelF = MultiModelF,
  MultiModelG = MultiModelG,
  MultiModelA_interaction = MultiModelA_interaction,
  MultiModelA_interaction_method = MultiModelA_interaction_method,
  MultiModelA_rs = MultiModelA_rs, MultiModelA_rsnc = MultiModelA_rsnc,
  MultiModelA_rs_method = MultiModelA_rs_method
)

loo_list <- purrr::imap(model_list_all, function(m, name) {
  loo_or_load(paste0("models/loo_", name, ".rds"), m)
})

## 3A-vi: Model selection on the shared complete-case subset -------------------
# All 19 candidate models refit on the rows complete for every predictor used
# across the candidate set (n = 221), so LOO is directly comparable across all
# of them. The selected model is then reported from its fit on its own
# complete data (n = 271; see 3A-vii).

shared_vars <- c("fish", "ref", "genusconfirmed", "sizez", "latmid", "season",
                 "depthz", "dist_to_shore_scaled", "method_clean")

diet_shared <- diet %>%
  filter(if_all(all_of(shared_vars), ~ !is.na(.x))) %>%
  mutate(genusconfirmed = droplevels(genusconfirmed))   # Diaphus stays reference

shared_models <- imap(model_list_all, function(m, name) {
  fit_or_load(paste0("models/", name, "_shared.rds"),
              formula(m), data = diet_shared)
})

loo_list_shared <- imap(shared_models, function(m, name) {
  loo_or_load(paste0("models/loo_", name, "_shared.rds"), m)
})

## 3A-vii: Random-effects structure and selected model -------------------------
# Stage 1 (3A-vi, n = 221) selects the fixed effects: the genus x size
# interaction is retained; environmental covariates are not. Stage 2 compares
# random-effects structures among the five top-ranked genus x size models on
# their common data (n = 271, all fit on identical rows). Random size slopes
# improved predictive performance; the intercept-slope correlation was not
# identified, so the uncorrelated structure is selected. Identification method
# is constant within study and is absorbed by the study random effects (models
# with and without method were tied), so it is not included.

re_candidates <- c("MultiModelA_rsnc", "MultiModelA_rs", "MultiModelA_rs_method",
                   "MultiModelA_interaction", "MultiModelA_interaction_method")

stopifnot(n_distinct(map_int(model_list_all[re_candidates], nobs)) == 1)

best_name  <- "MultiModelA_rsnc"
best_model <- model_list_all[[best_name]]   # fish ~ sizez * genusconfirmed + (1 + sizez || ref)

# 
# SECTION 3B: TRAIT-BASED MODELS################################################
# 
# Tests whether species-level functional traits explain piscivory better than,
# or in addition to, genus and body size. Framed as testing *why* genus
# differs (a mechanistic follow-up)
#
# NOTE ON SAMPLE SIZE: the single-trait models in 3B-i use per-model complete
# cases (each on its own maximal data). gillrakers_z is NA for all
# Centrobranchus observations (tooth-patch morphology, no comparable raker
# count -- see Section 2b/2f), so any model containing it drops that genus,
# while models built only on jaw_class/raptorial_teeth retain it. LOO is
# therefore NOT comparable across 3B-i models; trait comparisons are made on
# a shared subset in 3B-ii.

## 3B-i: Trait models on their own complete data -------------------------------

GillRakerModel_cc <- fit_or_load(
  "models/GillRakerModel_cc.rds",
  fish ~ gillrakers_z + (1 | ref),
  data = diet %>% filter(!is.na(gillrakers_z))
)

JawClassModel_cc <- fit_or_load(
  "models/JawClassModel_cc.rds",
  fish ~ jaw_class + (1 | ref),
  data = diet %>% filter(!is.na(jaw_class))
)

ToothTypeModel_cc <- fit_or_load(
  "models/ToothTypeModel_cc.rds",
  fish ~ raptorial_teeth + (1 | ref),
  data = diet %>% filter(!is.na(raptorial_teeth))
)

trait_cc_data <- diet %>%
  filter(!is.na(gillrakers_z), !is.na(jaw_class), !is.na(raptorial_teeth))

message(sprintf("Combined trait model complete cases: %d / %d",
                nrow(trait_cc_data), nrow(diet)))

CombinedTraitModel <- fit_or_load(
  "models/CombinedTraitModel.rds",
  fish ~ gillrakers_z + jaw_class + raptorial_teeth + (1 | ref),
  data = trait_cc_data
)

JawToothModel <- fit_or_load(
  "models/MouthTraitModel.rds",
  fish ~ jaw_class + raptorial_teeth + (1 | ref),
  data = trait_cc_data
)

## 3B-ii: Traits vs. genus x size on a shared trait subset ----------------------
# All trait models plus the best-supported genus x size model, refit on rows
# complete for every predictor among them (n = 223), so LOO is directly
# comparable. Centrobranchus drops out here (no gill-raker count; Section 2b).
# The best trait model is then reported from its fit on its own complete data
# (ToothTypeModel_cc, n = 333; Section 4).

trait_vars <- c("fish", "ref", "genusconfirmed", "sizez", "method_clean",
                "gillrakers_z", "jaw_class", "raptorial_teeth")

diet_traits <- diet %>%
  filter(if_all(all_of(trait_vars), ~ !is.na(.x))) %>%
  mutate(across(c(genusconfirmed, jaw_class, raptorial_teeth), droplevels))

message(sprintf("Trait comparison subset: n = %d rows, %d studies, %d genera",
                nrow(diet_traits), n_distinct(diet_traits$ref),
                n_distinct(diet_traits$genusconfirmed)))

trait_formulas <- list(
  GillRakerModel                 = fish ~ gillrakers_z + (1 | ref),
  JawClassModel                  = fish ~ jaw_class + (1 | ref),
  ToothTypeModel                 = fish ~ raptorial_teeth + (1 | ref),
  JawToothModel                  = fish ~ jaw_class + raptorial_teeth + (1 | ref),
  CombinedTraitModel             = fish ~ gillrakers_z + jaw_class + raptorial_teeth + (1 | ref)
)
trait_formulas[[best_name]] <- formula(best_model)   # selected genus x size model (3A-vii)

trait_models_shared <- imap(trait_formulas, function(f, name) {
  fit_or_load(paste0("models/", name, "_traitset.rds"), f, data = diet_traits)
})

stopifnot(all(map_int(trait_models_shared, nobs) == nrow(diet_traits)))

loo_list_traits <- imap(trait_models_shared, function(m, name) {
  loo_or_load(paste0("models/loo_", name, "_traitset.rds"), m)
})

## 3B-iii: Full-data trait LOO by comparability group (not used in tables) -----
# Retained from the earlier per-model approach. Superseded by the shared-subset
# comparison in 3B-ii (Table: model comparison), which is the one reported.

trait_model_list <- list(
  GeneraModel = GeneraModel, SizeModel = SizeModel,
  GillRakerModel_cc = GillRakerModel_cc, JawClassModel_cc = JawClassModel_cc,
  ToothTypeModel_cc = ToothTypeModel_cc, CombinedTraitModel = CombinedTraitModel
)

trait_loo_list <- list(
  GeneraModel = loo_list$GeneraModel,   # reused from Section 3A
  SizeModel   = loo_list$SizeModel,     # reused from Section 3A
  GillRakerModel_cc = loo_or_load("models/loo_GillRakerModel_cc.rds", 
                                  GillRakerModel_cc),
  JawClassModel_cc  = loo_or_load("models/loo_JawClassModel_cc.rds", JawClassModel_cc),
  ToothTypeModel_cc = loo_or_load("models/loo_ToothTypeModel_cc.rds", ToothTypeModel_cc),
  CombinedTraitModel = loo_or_load("models/loo_CombinedTraitModel.rds", CombinedTraitModel)
)

trait_loo_df <- imap_dfr(trait_loo_list, function(l, name) {
  m <- trait_model_list[[name]]
  tibble(
    Model = name,
    n = nobs(m),
    LOOIC = round(l$estimates["looic", "Estimate"], 1),
    SE = round(l$estimates["looic", "SE"], 1),
    elpd_loo = round(l$estimates["elpd_loo", "Estimate"], 1),
    elpd_loo_per_n = round(l$estimates["elpd_loo", "Estimate"] / nobs(m), 4),
    Pareto_k_high = sum(l$diagnostics$pareto_k > 0.7)
  )
}) %>% arrange(n, LOOIC)

# 
# SECTION 3C: SENSITIVITY ANALYSES OF THE BEST-SUPPORTED MODEL##################
#

## 3C-i: Prior sensitivity ------------------------------------------------------
# All sensitivity refits take their formula from best_model (3A-vii), and their
# cache file names include best_name, so changing the selected model refits
# them rather than silently loading fits of an earlier model.

best_formula <- formula(best_model)
best_rhs     <- best_model$formula$formula[[3]]   # right-hand side, for weighted refits

best_priorB <- fit_or_load(paste0("models/", best_name, "_priorB.rds"),
                           best_formula, prior = priorB)

best_priorC <- fit_or_load(paste0("models/", best_name, "_priorC.rds"),
                           best_formula, prior = priorC)

extract_key_params <- function(model, prior_label) {
  fixef(model, probs = c(0.025, 0.975)) %>%
    as.data.frame() %>%
    rownames_to_column("Parameter") %>%
    filter(Parameter %in% c("Intercept", "sizez",
                            "genusconfirmedGymnoscopelus",
                            "sizez:genusconfirmedGymnoscopelus")) %>%
    mutate(
      Prior    = prior_label,
      OR       = round(exp(Estimate), 2),
      OR_lower = round(exp(Q2.5),    2),
      OR_upper = round(exp(Q97.5),   2),
      OR_fmt   = sprintf("%.2f [%.2f, %.2f]", OR, OR_lower, OR_upper),
      Parameter_clean = case_when(
        Parameter == "Intercept"                         ~ "Intercept (Diaphus, mean size)",
        Parameter == "sizez"                              ~ "Size (Diaphus)",
        Parameter == "genusconfirmedGymnoscopelus"        ~ "Gymnoscopelus",
        Parameter == "sizez:genusconfirmedGymnoscopelus"  ~ "Size × Gymnoscopelus"
      )
    ) %>%
    select(Prior, Parameter_clean, OR_fmt)
}

sensitivity_df <- bind_rows(
  extract_key_params(best_model,  "Moderate: normal(0,1), Cauchy(0,1)"),
  extract_key_params(best_priorB, "Diffuse:  normal(0,5), Cauchy(0,2)"),
  extract_key_params(best_priorC, "Tight:    normal(0,0.5), Cauchy(0,0.5)")
)

loo_priorB   <- loo_or_load(paste0("models/loo_", best_name, "_priorB.rds"), best_priorB)
loo_priorC   <- loo_or_load(paste0("models/loo_", best_name, "_priorC.rds"), best_priorC)
loo_moderate <- loo_list[[best_name]]   # computed in 3A-v

loo_sensitivity <- tibble(
  Prior = c("Moderate: normal(0,1), Cauchy(0,1)",
            "Diffuse:  normal(0,5), Cauchy(0,2)",
            "Tight:    normal(0,0.5), Cauchy(0,0.5)"),
  LOOIC = c(round(loo_moderate$estimates["looic", "Estimate"], 1),
            round(loo_priorB$estimates["looic",   "Estimate"], 1),
            round(loo_priorC$estimates["looic",   "Estimate"], 1)),
  Pareto_k_high = c(sum(loo_moderate$diagnostics$pareto_k > 0.7),
                    sum(loo_priorB$diagnostics$pareto_k   > 0.7),
                    sum(loo_priorC$diagnostics$pareto_k   > 0.7))
)

## 3C-ii: Random-effects structure ----------------------------------------------
# Key estimates under each of the five candidate random-effects structures
# compared in 3A-vii (predictive comparison is reported in the model
# comparison table; this shows whether the choice changes the conclusions).

re_labels <- c(
  MultiModelA_rsnc               = "Intercept + size slope, uncorrelated",
  MultiModelA_rs                 = "Intercept + size slope, correlated",
  MultiModelA_rs_method          = "Intercept + size slope, correlated; + method",
  MultiModelA_interaction        = "Intercept only",
  MultiModelA_interaction_method = "Intercept only; + method"
)
re_labels[best_name] <- paste(re_labels[best_name], "(selected)")
stopifnot(setequal(names(re_labels), re_candidates))

re_structure_df <- imap_dfr(model_list_all[re_candidates],
                            ~ extract_key_params(.x, re_labels[[.y]]))

## 3C-iii: Stomach-count weighting (attempted; not adopted) ---------------------
# not viable here: nstomachsfish populated for <2% of observations). 
# Two weighting schemes were tested against the primary model structure; 
# both produced numerically unstable fits relative to the
# unweighted primary model, so neither is adopted. 

weighted_raw <- fit_or_load(
  paste0("models/", best_name, "_weighted_raw.rds"),
  as.formula(bquote(fish | weights(nstomachs_weight) ~ .(best_rhs)))
)

weighted_log <- fit_or_load(
  paste0("models/", best_name, "_weighted_log.rds"),
  as.formula(bquote(fish | weights(nstomachs_logweight) ~ .(best_rhs)))
)

loo_weighted_raw <- loo_or_load(paste0("models/loo_", best_name, "_weighted_raw.rds"), weighted_raw)
loo_weighted_log <- loo_or_load(paste0("models/loo_", best_name, "_weighted_log.rds"), weighted_log)

re_sd <- function(m, term) VarCorr(m)$ref$sd[term, "Estimate"]

weighting_diagnostics <- tibble(
  Model = c("Unweighted (primary model)", "Weighted: raw nstomachs", "Weighted: log1p(nstomachs)"),
  SD_intercept_study = c(re_sd(best_model, "Intercept"),
                         re_sd(weighted_raw, "Intercept"),
                         re_sd(weighted_log, "Intercept")),
  SD_slope_study     = c(re_sd(best_model, "sizez"),
                         re_sd(weighted_raw, "sizez"),
                         re_sd(weighted_log, "sizez")),
  Pareto_k_high = c(sum(loo_moderate$diagnostics$pareto_k > 0.7),
                    sum(loo_weighted_raw$diagnostics$pareto_k > 0.7),
                    sum(loo_weighted_log$diagnostics$pareto_k > 0.7))
)

## 3C-iv: Small-n genera sensitivity --------------------------------------------
# Refits the primary
# model structure excluding genera with fewer than 5 observations, to
# confirm results are not driven by sparsely-sampled taxa.

genus_n <- diet %>% count(genusconfirmed, name = "n") %>% arrange(n)

small_genera <- genus_n %>% filter(n < 5) %>% pull(genusconfirmed)

# extract_key_params() downstream pulls "genusconfirmedGymnoscopelus" and its
# interaction term specifically -- if Gymnoscopelus itself gets excluded here,
# those two rows will just be silently absent from the n5 model's half of
# smalln_comparison rather than erroring, which is easy to miss in a wide table.
if ("Gymnoscopelus" %in% small_genera) {
  warning("Gymnoscopelus is among the excluded small-n genera -- its rows ",
          "will be missing for the n<5-excluded model in smalln_comparison below.")
}

diet_n5 <- diet %>%
  filter(!genusconfirmed %in% small_genera) %>%
  mutate(genusconfirmed = droplevels(genusconfirmed))

message(sprintf("Small-n sensitivity: excluded %d genera (%s), retained %d of %d observations",
                length(small_genera), paste(small_genera, collapse = ", "),
                nrow(diet_n5), nrow(diet)))

best_model_n5 <- fit_or_load(paste0("models/", best_name, "_n5.rds"),
                             best_formula, data = diet_n5)

loo_n5 <- loo_or_load(paste0("models/loo_", best_name, "_n5.rds"), best_model_n5)

smalln_comparison <- bind_rows(
  extract_key_params(best_model,    "All genera"),
  extract_key_params(best_model_n5, "Genera with n≥5 only")
)

smalln_diagnostics <- tibble(
  Model = c("All genera", "Genera with n≥5 only"),
  n_obs = c(nobs(best_model), nobs(best_model_n5)),
  elpd_loo_per_n = c(
    loo_moderate$estimates["elpd_loo", "Estimate"] / nobs(best_model),
    loo_n5$estimates["elpd_loo", "Estimate"] / nobs(best_model_n5)
  ),
  Pareto_k_high = c(sum(loo_moderate$diagnostics$pareto_k > 0.7),
                    sum(loo_n5$diagnostics$pareto_k > 0.7))
)

# 
# SECTION 4: RESULTS — TABLES###################################################
# 

## 4-i: Designate best-supported trait model ------------------------------------
# The genus x size model (best_model) is designated in 3A-vii. The trait model
# is selected in 3B-ii; both are reported from fits on their own complete data.

best_trait_model <- ToothTypeModel_cc                # n = 333

## 4-ii: Model comparison -------------------------------------------------------
# Three comparison sets, each fit on identical rows so elpd differences are valid:
# (1) all 19 candidate models on the shared complete-case subset (3A-vi)
# (2) random-effects structures among the top genus x size models (3A-vii)
# (3) trait models vs. the selected genus x size model on the trait subset (3B-ii)
# elpd_diff and se_diff come from loo_compare() within each set.

var_labels <- c("genusconfirmed"       = "Genus",
                "sizez"                = "Size",
                "latmid"               = "Latitude",
                "season"               = "Season",
                "depthz"               = "Depth",
                "dist_to_shore_scaled" = "Distance to shore",
                "method_clean"         = "Method",
                "gillrakers_z"         = "Gill rakers",
                "jaw_class"            = "Jaw class",
                "raptorial_teeth"      = "Tooth type",
                " \\* "                = " \u00d7 ")

build_comparison_set <- function(models, loos, set_label) {
  cmp <- loo_compare(loos) %>%
    as.data.frame() %>%
    rownames_to_column("Model") %>%
    select(Model, elpd_diff, se_diff)
  
  imap_dfr(models, function(m, name) {
    l  <- loos[[name]]
    stopifnot(nrow(l$pointwise) == nobs(m))   # cached LOO matches the model
    f  <- deparse1(formula(m)$formula)
    r2 <- bayes_R2(m)
    tibble(
      Model          = name,
      Fixed_effects  = f %>%
        str_remove("^fish ~ ") %>%
        str_remove(" \\+ \\(.*\\)$") %>%
        str_replace_all(var_labels),
      Random_effects = str_extract(f, "\\(.*\\)"),
      LOOIC          = round(l$estimates["looic",    "Estimate"], 1),
      elpd_loo       = round(l$estimates["elpd_loo", "Estimate"], 1),
      p_loo          = round(l$estimates["p_loo",    "Estimate"], 1),
      Pareto_k_high  = sum(l$diagnostics$pareto_k > 0.7),
      R2             = sprintf("%.2f [%.2f, %.2f]",
                               r2[1, "Estimate"], r2[1, "Q2.5"], r2[1, "Q97.5"])
    )
  }) %>%
    left_join(cmp, by = "Model") %>%
    mutate(elpd_diff = round(elpd_diff, 1),
           se_diff   = round(se_diff, 1),
           Set       = set_label) %>%
    arrange(desc(elpd_diff))
}

supp_table2_final <- bind_rows(
  build_comparison_set(shared_models, loo_list_shared,
                       sprintf("Candidate models (n = %d)", nrow(diet_shared))),
  build_comparison_set(model_list_all[re_candidates], loo_list[re_candidates],
                       sprintf("Random-effects structure (n = %d)", nobs(best_model))),
  build_comparison_set(trait_models_shared, loo_list_traits,
                       sprintf("Traits vs. genus \u00d7 size (n = %d)", nrow(diet_traits)))
) %>%
  select(Set, Model, Fixed_effects, Random_effects, elpd_diff, se_diff,
         LOOIC, elpd_loo, p_loo, Pareto_k_high, R2)

k_note <- if (all(supp_table2_final$Pareto_k_high == 0)) {
  "All Pareto k values were < 0.7."
} else {
  "Observations with Pareto k > 0.7 are counted per model."
}

supp_table_loo <- flextable(supp_table2_final) %>%
  set_header_labels(
    Set = "Comparison set", Model = "Model", Fixed_effects = "Fixed effects",
    Random_effects = "Random effects", elpd_diff = "\u0394elpd", se_diff = "SE(\u0394elpd)",
    LOOIC = "LOOIC", elpd_loo = "elpd_loo", p_loo = "p_loo",
    Pareto_k_high = "Pareto k > 0.7 (n)", R2 = "Bayesian R2 [95% CrI]"
  ) %>%
  merge_v(j = "Set") %>%
  valign(j = "Set", valign = "top") %>%
  hline(i = which(!duplicated(supp_table2_final$Set))[-1] - 1,
        border = fp_border(color = "black", width = 1)) %>%
  align(align = "center", part = "all") %>%
  align(j = 1:4, align = "left", part = "all") %>%
  autofit() %>%
  set_caption(
    paste("Table S[X]. Model comparison by leave-one-out cross-validation.",
          "Models were selected in two stages: fixed effects were compared",
          "across all candidate models on observations complete for every",
          "candidate predictor (Candidate models), and random-effects",
          "structures were then compared among the top-ranked genus \u00d7",
          "size models on their common observations (Random-effects",
          "structure). Trait-based models were compared with the selected",
          "model on observations complete for all trait predictors.",
          "Each comparison set was fit on identical observations, so",
          "predictive performance is directly comparable within, but not",
          "between, sets. \u0394elpd is",
          "the difference in expected log pointwise predictive density from",
          "the top-ranked model in each set, with its standard error",
          "(SE(\u0394elpd)); differences smaller than ~2 SE indicate models",
          "are not distinguishable in predictive performance. p_loo is the",
          "effective number of parameters.", k_note,
          "Random effects: (1 | ref) = random intercept for study;",
          "(1 + sizez | ref) = correlated random intercept and size slope;",
          "(1 + sizez || ref) = uncorrelated random intercept and size slope.")
  ) %>%
  theme_booktabs()

## 4-iii: Diagnostics for the best-supported models ------------------------------
# Both models reported from their fits on their own complete data
# (genus x size: best_model, n = 271; tooth type: n = 333), not from the
# comparison subsets.

diag_models <- list(
  "Genus \u00d7 size" = list(
    model = best_model,
    loo   = loo_list[[best_name]]),
  "Tooth type" = list(
    model = best_trait_model,
    loo   = loo_or_load("models/loo_ToothTypeModel_cc.rds", best_trait_model))
)

model_diagnostics <- function(m, l) {
  stopifnot(nrow(l$pointwise) == nobs(m))   # cached LOO matches the model
  draws_sum <- posterior::as_draws_df(m) %>%
    posterior::subset_draws(variable = c("^b_", "^sd_", "^r_"), regex = TRUE) %>%
    posterior::summarise_draws("rhat", "ess_bulk", "ess_tail")
  np <- nuts_params(m)
  r2 <- bayes_R2(m)
  
  tibble(
    `Observations (n)`      = as.character(nobs(m)),
    `Studies (n)`           = as.character(n_distinct(m$data$ref)),
    `Max R-hat`             = sprintf("%.3f", max(draws_sum$rhat, na.rm = TRUE)),
    `Min bulk ESS`          = sprintf("%.0f", min(draws_sum$ess_bulk, na.rm = TRUE)),
    `Min tail ESS`          = sprintf("%.0f", min(draws_sum$ess_tail, na.rm = TRUE)),
    `Divergent transitions` = as.character(sum(np$Value[np$Parameter == "divergent__"])),
    `Max treedepth hits`    = as.character(sum(np$Value[np$Parameter == "treedepth__"] >=
                                                 brm_defaults$control$max_treedepth)),
    `elpd_loo (SE)`         = sprintf("%.1f (%.1f)",
                                      l$estimates["elpd_loo", "Estimate"],
                                      l$estimates["elpd_loo", "SE"]),
    `p_loo`                 = sprintf("%.1f", l$estimates["p_loo", "Estimate"]),
    `Pareto k > 0.7 (n)`    = as.character(sum(l$diagnostics$pareto_k > 0.7)),
    `Bayesian R2 [95% CrI]` = sprintf("%.2f [%.2f, %.2f]",
                                      r2[1, "Estimate"], r2[1, "Q2.5"], r2[1, "Q97.5"])
  )
}

diag_df <- imap_dfr(diag_models, function(x, name) {
  model_diagnostics(x$model, x$loo) %>% mutate(Model = name)
}) %>%
  pivot_longer(-Model, names_to = "Diagnostic", values_to = "Value") %>%
  pivot_wider(names_from = Model, values_from = Value)

supp_table_diag <- flextable(diag_df) %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  autofit() %>%
  set_caption(
    paste("Table S[X]. Convergence and predictive diagnostics for the",
          sprintf("best-supported genus \u00d7 size model (%s)", best_name),
          "and the best-supported trait model (ToothTypeModel_cc), each fit",
          "on its own complete data. R-hat and effective sample sizes (ESS)",
          "are reported across all population-level, group-level SD, and",
          "study-level parameters; values of R-hat < 1.01 and ESS > 400",
          "indicate adequate convergence. Divergent transitions and treedepth",
          "hits are summed across all 4 chains \u00d7 2,000 post-warmup draws.",
          "elpd_loo is the expected log pointwise predictive density from",
          "leave-one-out cross-validation; p_loo is the effective number of",
          "parameters. Values are not comparable between models, which are",
          "fit on different observations.")
  ) %>%
  theme_booktabs()

## 4-iv: Genus x size model parameters -------------------------------------------

# Genera with a single size value in the model data: their size interaction
# is not informed by the data and reflects the prior
slope_uninformed <- best_model$data %>%
  group_by(genusconfirmed) %>%
  summarise(n_sizes = n_distinct(sizez), .groups = "drop") %>%
  filter(n_sizes < 2) %>%
  pull(genusconfirmed) %>% as.character()

fe_tidy <- fixef(best_model, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>%
  rownames_to_column("Parameter") %>%
  rename(Est_Error = Est.Error) %>%
  mutate(
    OR        = round(exp(Estimate), 2),
    OR_lower  = round(exp(Q2.5),    2),
    OR_upper  = round(exp(Q97.5),   2),
    OR_CrI    = sprintf("%.2f [%.2f, %.2f]", OR, OR_lower, OR_upper),
    Estimate_fmt = sprintf("%.2f [%.2f, %.2f]",
                           round(Estimate, 2),
                           round(Q2.5,     2),
                           round(Q97.5,    2)),
    # Credible effect: 95% CrI excludes 1.0 on OR scale
    Credible  = ifelse(OR_lower > 1 | OR_upper < 1, "Yes", ""),
    # Clean parameter names for display in tables and figures.
    # Order matters: the interaction pattern must run before the plain
    # "genusconfirmed" pattern, or the latter strips the text the former
    # needs to match.
    Parameter_clean = Parameter %>%
      str_replace("^b_", "") %>%
      str_replace("^sizez:genusconfirmed", "Size \u00d7 ") %>%
      str_replace("^genusconfirmed", "") %>%
      str_replace("^method_cleanmolecular_isotope_inclusive$",
                  "Method: molecular/isotope (ref: visual)") %>%
      str_replace("^sizez$", "Size (ref: Diaphus)") %>%
      str_replace("^Intercept$", "Intercept (Diaphus, mean size)"),
    Parameter_clean = if_else(Parameter %in% paste0("sizez:genusconfirmed", slope_uninformed),
                              paste0(Parameter_clean, "\u2020"), Parameter_clean)
  )

# Study-level SDs (random intercept and random size slope), with 95% CrI
re_summary <- VarCorr(best_model)$ref$sd
re_df <- tibble(
  Parameter_clean = c("SD (intercept | study)", "SD (size slope | study)"),
  Estimate_fmt    = sprintf("%.2f [%.2f, %.2f]",
                            re_summary[c("Intercept", "sizez"), "Estimate"],
                            re_summary[c("Intercept", "sizez"), "Q2.5"],
                            re_summary[c("Intercept", "sizez"), "Q97.5"]),
  OR_CrI          = "\u2014",
  Credible        = ""
)

params_display <- bind_rows(
  fe_tidy %>% select(Parameter_clean, Estimate_fmt, OR_CrI, Credible),
  re_df
)

table1 <- flextable(params_display) %>%
  set_header_labels(
    Parameter_clean = "Parameter",
    Estimate_fmt    = "Posterior estimate [95% CrI]",
    OR_CrI          = "OR [95% CrI]",
    Credible        = "Credible effect"
  ) %>%
  bold(i = ~ Credible == "Yes") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  hline(i = nrow(fe_tidy), border = fp_border(color = "black", width = 1)) %>%
  autofit() %>%
  set_caption(
    paste("Table [X]. Posterior estimates from the best-supported model",
          "(Genus \u00d7 Size with uncorrelated study-level random intercepts",
          sprintf("and size slopes; %s). Fixed effects are reported as", best_name),
          "posterior mean log-odds with 95% credible intervals (CrI); odds",
          "ratios (ORs) are exponentiated values; for the intercept this is",
          "the baseline odds of piscivory. Genus ORs represent the",
          "multiplicative change in baseline piscivory odds relative to",
          "Diaphus at mean body size. Size \u00d7 genus interaction ORs",
          "indicate how the size effect for each genus differs",
          "multiplicatively from Diaphus. Study-level standard deviations are",
          "on the log-odds scale. Bold rows indicate effects whose 95% CrI",
          "excludes 1.0. \u2020 Genus represented by a single body size; its",
          "size interaction is not informed by the data and reflects the prior.")
  ) %>%
  theme_booktabs()

## 4-v: Tooth-type model parameters ----------------------------------------------
# Same structure as 4-iv, for the best-supported trait model.

fe_tidy_trait <- fixef(best_trait_model, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>%
  rownames_to_column("Parameter") %>%
  mutate(
    OR           = round(exp(Estimate), 2),
    OR_lower     = round(exp(Q2.5),     2),
    OR_upper     = round(exp(Q97.5),    2),
    OR_CrI       = sprintf("%.2f [%.2f, %.2f]", OR, OR_lower, OR_upper),
    Estimate_fmt = sprintf("%.2f [%.2f, %.2f]",
                           round(Estimate, 2), round(Q2.5, 2), round(Q97.5, 2)),
    Credible     = ifelse(OR_lower > 1 | OR_upper < 1, "Yes", ""),
    Parameter_clean = case_when(
      Parameter == "Intercept"                ~ "Intercept (villiform-only teeth)",
      Parameter == "raptorial_teethraptorial" ~ "Raptorial teeth (ref: villiform-only)",
      TRUE ~ Parameter
    )
  )

re_summary_trait <- VarCorr(best_trait_model)$ref$sd

params_display_trait <- bind_rows(
  fe_tidy_trait %>% select(Parameter_clean, Estimate_fmt, OR_CrI, Credible),
  tibble(
    Parameter_clean = "SD (intercept | study)",
    Estimate_fmt    = sprintf("%.2f [%.2f, %.2f]",
                              re_summary_trait["Intercept", "Estimate"],
                              re_summary_trait["Intercept", "Q2.5"],
                              re_summary_trait["Intercept", "Q97.5"]),
    OR_CrI          = "\u2014",
    Credible        = ""
  )
)

table_trait <- flextable(params_display_trait) %>%
  set_header_labels(
    Parameter_clean = "Parameter",
    Estimate_fmt    = "Posterior estimate [95% CrI]",
    OR_CrI          = "OR [95% CrI]",
    Credible        = "Credible effect"
  ) %>%
  bold(i = ~ Credible == "Yes") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  hline(i = nrow(fe_tidy_trait), border = fp_border(color = "black", width = 1)) %>%
  autofit() %>%
  set_caption(
    paste("Table [X]. Posterior estimates from the best-supported trait model",
          "(tooth type; ToothTypeModel_cc). Fixed effects are reported as",
          "posterior mean log-odds with 95% credible intervals (CrI); odds",
          "ratios (ORs) are exponentiated values; for the intercept this is",
          "the baseline odds of piscivory for species with villiform-only",
          "teeth. The raptorial-teeth OR is the multiplicative change in",
          "piscivory odds for species with raptorial teeth relative to",
          "villiform-only species. Bold rows indicate effects whose 95% CrI",
          "excludes 1.0.")
  ) %>%
  theme_booktabs()

## 4-vi: Sensitivity analyses ----------------------------------------------------
# All robustness checks of the best-supported model (Section 3C), in one table.

prior_or_rows <- sensitivity_df %>%
  rename(Variant = Prior, Parameter = Parameter_clean, Value = OR_fmt) %>%
  mutate(Check = "Prior specification", Metric = "OR [95% CrI]") %>%
  select(Check, Variant, Metric, Parameter, Value)

prior_loo_rows <- loo_sensitivity %>%
  pivot_longer(c(LOOIC, Pareto_k_high), names_to = "Metric", values_to = "Value") %>%
  mutate(Check = "Prior specification", Parameter = NA_character_, Value = as.character(Value)) %>%
  rename(Variant = Prior) %>%
  select(Check, Variant, Metric, Parameter, Value)

re_structure_rows <- re_structure_df %>%
  rename(Variant = Prior, Parameter = Parameter_clean, Value = OR_fmt) %>%
  mutate(Check = "Random-effects structure", Metric = "OR [95% CrI]") %>%
  select(Check, Variant, Metric, Parameter, Value)

weighting_rows <- weighting_diagnostics %>%
  pivot_longer(c(SD_intercept_study, SD_slope_study, Pareto_k_high),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Check = "Stomach-count weighting (attempted; not adopted)",
         Parameter = NA_character_,
         Value = as.character(round(as.numeric(Value), 2))) %>%
  rename(Variant = Model) %>%
  select(Check, Variant, Metric, Parameter, Value)

smalln_or_rows <- smalln_comparison %>%
  rename(Variant = Prior, Parameter = Parameter_clean, Value = OR_fmt) %>%
  mutate(Check = "Small-n genera (n<5 excluded)", Metric = "OR [95% CrI]") %>%
  select(Check, Variant, Metric, Parameter, Value)

smalln_diag_rows <- smalln_diagnostics %>%
  mutate(elpd_loo_per_n = round(elpd_loo_per_n, 4)) %>%
  pivot_longer(c(elpd_loo_per_n, Pareto_k_high), names_to = "Metric", values_to = "Value") %>%
  mutate(Check = "Small-n genera (n<5 excluded)", Parameter = NA_character_, Value = as.character(Value)) %>%
  rename(Variant = Model) %>%
  select(Check, Variant, Metric, Parameter, Value)

sensitivity_combined <- bind_rows(
  prior_or_rows, prior_loo_rows,
  re_structure_rows,
  weighting_rows,
  smalln_or_rows, smalln_diag_rows
) %>%
  mutate(Parameter = replace_na(Parameter, "\u2014"))

supp_table_sensitivity <- flextable(sensitivity_combined) %>%
  set_header_labels(Check = "Sensitivity check", Variant = "Model / specification",
                    Metric = "Metric", Parameter = "Parameter", Value = "Value") %>%
  merge_v(j = c("Check", "Variant")) %>%
  valign(j = c("Check", "Variant"), valign = "top") %>%
  hline(i = which(!duplicated(sensitivity_combined$Check))[-1] - 1,
        border = fp_border(color = "black", width = 1)) %>%
  align(align = "center", part = "all") %>%
  align(j = 1:2, align = "left", part = "all") %>%
  autofit() %>%
  set_caption(
    paste("Supplementary TableX. Sensitivity analyses for the best-supported",
          "model (Genus \u00d7 Size with study-level random intercepts and",
          "size slopes, uncorrelated). Random-effects structure: key estimates under each",
          "candidate random-effects structure (predictive comparison in the",
          "model comparison table).",
          "Prior specification: alternative weakly-informative priors,",
          "tested for consistency of key parameter estimates and predictive",
          "performance. Stomach-count weighting: two weighting schemes (raw",
          "and log-transformed stomach sample size) tested in response to",
          "reviewer comments on sample-size information loss; both produced",
          "numerically unstable fits (elevated between-study variance, high",
          "Pareto k) and are reported but not adopted. Small-n genera:",
          "primary model refit excluding genera with fewer than 5",
          "observations, confirming results are not driven by",
          "sparsely-sampled taxa.")
  ) %>%
  theme_booktabs()


## 4-vii: Identified fish prey ----------------------------------------------

# Fix spelling errors in the extracted prey names
prey_fixes <- c("argyopelecus" = "argyropelecus",
                "hemygimnus"   = "hemigymnus",
                "l\\.corodilus" = "l.crocodilus")

prey_groups <- list(
  "Lanternfishes"       = "myctoph|diaphus|lampanyct|diogenichthys|hygophum|ceratoscopelus|b\\.glaciale|c\\.maderensis|h\\.benoiti|h\\.hygomii|l\\.crocodilus|l\\.pusillus|m\\.punctatum|n\\.elongatus",
  "Other mesopelagic"   = "cyclothone|vinciguerria|argyropelecus|maurolicus|phosichth|chauliodus|stomias|paralepid|polyipnus|gonostomat|bregmaceros",
  "Coastal or neritic"  = "sardine|clupe|carangid|ariosoma|gnathophis|trichiur",
  "Unspecified"         = "larval fish|aulopiform|actinopterygii|osteichthyes|others"
)

prey_obs <- fish_rows %>%
  mutate(prey = str_squish(str_to_lower(coalesce(fishid, ""))),
         prey = str_replace_all(prey, prey_fixes)) %>%
  filter(!prey %in% c("", "unknown", "unidentified"))

# Prey group for each observation (an observation can fall in several groups)
prey_obs <- prey_obs %>%
  mutate(groups = map_chr(prey, function(p) {
    hits <- names(prey_groups)[map_lgl(prey_groups, ~ str_detect(p, .x))]
    if (length(hits) == 0) "Unclassified" else paste(hits, collapse = "; ")
  }))

stopifnot(!any(prey_obs$groups == "Unclassified"))   # every prey string assigned

# Counts for the Results text (studies per prey group)
cat("Observations with identified prey:", nrow(prey_obs),
    "| studies:", n_distinct(prey_obs$ref), "\n")
imap_dfr(prey_groups, ~ tibble(group = .y,
                               observations = sum(str_detect(prey_obs$prey, .x)),
                               studies = n_distinct(prey_obs$ref[str_detect(prey_obs$prey, .x)]))) %>%
  print()

# One row per study
format_ref <- function(r) {
  r <- str_squish(r)
  paste(str_to_title(str_remove(r, "\\s*\\d{4}[a-z]?$")), str_extract(r, "\\d{4}[a-z]?$"))
}

method_names <- c("gut" = "Gut contents", "sia" = "Stable isotopes",
                  "dna" = "DNA", "gut and sia" = "Gut contents + stable isotopes",
                  "gut and dna" = "Gut contents + DNA")

prey_table <- prey_obs %>%
  mutate(prey_tokens = str_split(prey, ",\\s*")) %>%
  group_by(ref) %>%
  summarise(
    Predator = paste(sort(unique(str_to_sentence(spconfirmed))), collapse = ", "),
    Method   = paste(unique(method_names[method]), collapse = "; "),
    Prey     = paste(setdiff(unique(unlist(prey_tokens)), c("unknown", "unidentified")),
                     collapse = ", "),
    Group    = paste(sort(unique(unlist(str_split(groups, "; ")))), collapse = "; "),
    .groups  = "drop"
  ) %>%
  mutate(Study = format_ref(ref),
         Prey  = str_to_sentence(Prey),
         Method = if_else(str_detect(Method, "^Stable isotopes$"),
                          paste0(Method, "†"), Method)) %>%
  arrange(Study) %>%
  select(Study, Predator, Method, Prey, Group)

supp_table_prey <- flextable(prey_table) %>%
  set_header_labels(Study = "Study", Predator = "Lanternfish species",
                    Method = "Method", Prey = "Fish prey identified",
                    Group = "Prey group") %>%
  italic(j = "Predator", part = "body") %>%
  valign(valign = "top", part = "body") %>%
  align(align = "left", part = "all") %>%
  width(j = c("Study", "Predator", "Method", "Prey", "Group"),
        width = c(1.1, 2.0, 1.3, 2.6, 1.3)) %>%
  fontsize(size = 9, part = "all") %>%
  set_caption(
    paste("Table S[X]. Fish prey identified in the compiled studies.",
          "Prey are listed as reported; observations in which fish prey were",
          "unidentified are omitted. Prey groups: lanternfishes (Myctophidae);",
          "other mesopelagic fishes (e.g., Gonostomatidae, Sternoptychidae,",
          "Phosichthyidae, Stomiidae, Paralepididae); coastal or neritic taxa",
          "(Clupeidae, Carangidae, Congridae, Trichiuridae); and unspecified.",
          "† Prey inferred from stable isotopes rather than observed in",
          "stomach contents.")
  ) %>%
  theme_booktabs()

## 4-viii: Write all tables to Word ----------------------------------------------
doc <- read_docx() %>%
  # Main-text tables (portrait)
  body_add_par("Table [X]", style = "heading 2") %>%
  body_add_flextable(table1) %>%
  body_add_break() %>%
  body_add_par("Table [X]", style = "heading 2") %>%
  body_add_flextable(table_trait) %>%
  body_end_section_portrait() %>%
  # Model comparison: 11 columns, so given its own landscape section
  body_add_par("Table S[X]", style = "heading 2") %>%
  body_add_flextable(supp_table_loo) %>%
  body_end_section_landscape() %>%
  # Remaining supplementary tables (portrait)
  body_add_par("Table S[X]", style = "heading 2") %>%
  body_add_flextable(supp_table_diag) %>%
  body_add_break() %>%
  body_add_par("Table S[X]", style = "heading 2") %>%
  body_add_flextable(supp_table_sensitivity)%>%
  # identified prey
  body_add_par("Table S[X]", style = "heading 2") %>%
  body_add_flextable(supp_table_prey) %>%
  body_end_section_landscape() %>%
  print(target = "tables/TableS_fish_prey.docx")

print(doc, target = "tables/piscivory_tables.docx")
message("Tables written to tables/piscivory_tables.docx")

# 
# SECTION 4B: RESULTS — FIGURES##################################################
# 

## Figure 1 — Global distribution map -----------------------------------------

diet_map <- diet %>%
  mutate(
    fish_bin = ifelse(fish == "yes", 1, 0),
    longmid  = round(longmid, 1),
    latmid   = round(latmid,  1)
  ) %>%
  group_by(longmid, latmid) %>%
  summarize(
    n_total        = n(),
    prop_piscivory = mean(fish_bin, na.rm = TRUE),
    .groups        = "drop"
  )

Figure1 <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    fill      = "grey80",
    color     = "white",
    linewidth = 0.3
  ) +
  geom_point(
    data = diet_map %>% arrange (desc(n_total)),
    aes(x = longmid, y = latmid,
        size = n_total,
        fill = prop_piscivory),
    shape = 21,
    color = "grey30"
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "Proportion\npiscivorous",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.50", "0.75", "1.0")
  ) +
  scale_size_continuous(
    name   = "Observations",
    range  = c(2, 12),
    breaks = c(1, 5, 10, 20, 36)
  ) +
  coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-90, 90)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.line       = element_blank(),
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 13),
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    panel.border    = element_rect(color = "grey60", fill = NA, linewidth = 0.5)
  ) +
  guides(
    fill = guide_colorbar(order = 1, barwidth = 1, barheight = 6),
    size = guide_legend(order = 2, override.aes = list(fill = "grey50"))
  )

ggsave("Figure1.png", width = 9, height = 7, dpi = 600)

## Figure 2 — Observations by genus (A) and season (B) ------------------------

# Panel A: by genus
n_genus      <- diet %>%
  count(genusconfirmed) %>%
  mutate(label = paste0(genusconfirmed, " (n=", n, ")"))
genus_labels <- setNames(n_genus$label, n_genus$genusconfirmed)

Figure2a <- ggplot(diet, aes(x = genusconfirmed, fill = fish)) +
  geom_bar(position = "fill") +
  scale_fill_manual(
    name   = "Piscivory",
    values = c("no" = "#721F81", "yes" = "#F0F921"),
    labels = c("no" = "No", "yes" = "Yes")
  ) +
  scale_x_discrete(labels = genus_labels) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Proportion of observations", tag = "A") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(face = "italic", size = 16, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    plot.tag     = element_text(size = 18, face = "bold")
  )


# Panel B: by season
season_levels <- c("spring", "summer", "fall", "winter", "multiple")

diet_season   <- diet %>%
  filter(season %in% season_levels) %>%
  mutate(season = factor(season, levels = season_levels))

n_season      <- diet_season %>% count(season)
season_labels <- n_season %>%
  mutate(
    label  = paste0(str_to_title(season), "\n(n=", n, ")"),
    season = as.character(season)
  ) %>%
  select(season, label) %>%
  deframe()

Figure2b <- ggplot(diet_season, aes(x = season, fill = fish)) +
  geom_bar(position = "fill") +
  scale_fill_manual(
    name   = "Piscivory",
    values = c("no" = "#721F81", "yes" = "#F0F921"),
    labels = c("no" = "No", "yes" = "Yes")
  ) +
  scale_x_discrete(labels = season_labels) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Proportion of observations", tag = "B") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    plot.tag     = element_text(size = 18, face = "bold")
  )

# Panel C: by identification method
method_levels <- c("gut", "sia", "dna")

diet_method <- diet %>%
  mutate(
    method_group = case_when(
      method == "gut"                     ~ "gut",
      method %in% c("sia", "gut and sia") ~ "sia",
      method %in% c("dna", "gut and dna") ~ "dna",
      TRUE                                ~ NA_character_
    ),
    method_group = factor(method_group, levels = method_levels)
  )

stopifnot(!anyNA(diet_method$method_group))   # every record assigned to a group

n_method      <- diet_method %>% count(method_group)
method_names  <- c(gut = "Visual", sia = "Stable\nisotopes", dna = "DNA")
method_labels <- n_method %>%
  mutate(label        = paste0(method_names[as.character(method_group)], "\n(n=", n, ")"),
         method_group = as.character(method_group)) %>%
  select(method_group, label) %>%
  deframe()

# Counts for the results text
diet_method %>%
  group_by(method_group) %>%
  summarise(n = n(), fish = sum(fish == "yes"), pct = round(100 * fish / n, 1)) %>%
  print()

Figure2c <- ggplot(diet_method, aes(x = method_group, fill = fish)) +
  geom_bar(position = "fill") +
  scale_fill_manual(
    name   = "Piscivory",
    values = c("no" = "#721F81", "yes" = "#F0F921"),
    labels = c("no" = "No", "yes" = "Yes")
  ) +
  scale_x_discrete(labels = method_labels) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Proportion of observations", tag = "C") +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    plot.tag     = element_text(size = 18, face = "bold")
  )

# Combine: A across the top, B and C side by side below, one shared legend
Figure2 <- Figure2a / (Figure2b | Figure2c + labs(y = NULL)) +
  plot_layout(heights = c(2, 1), guides = "collect")

Figure2
ggsave("Figure2.png", Figure2, width = 13, height = 10, dpi = 600)

## Figure 3 — Size ranges by genus --------------------------------------------

genus_summary <- diet %>%
  filter(!is.na(maxsize)) %>%
  group_by(genusconfirmed) %>%
  summarize(
    min_maxsize    = min(maxsize,    na.rm = TRUE),
    max_maxsize    = max(maxsize,    na.rm = TRUE),
    median_maxsize = median(maxsize, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  mutate(genusconfirmed = str_to_title(genusconfirmed))

Figure3 <- ggplot(genus_summary,
                  aes(x = reorder(genusconfirmed, max_maxsize))) +
  geom_linerange(
    aes(ymin = min_maxsize, ymax = max_maxsize),
    color     = "black",
    linewidth = 1.2,
    alpha     = 0.7
  ) +
  geom_point(aes(y = min_maxsize),    color = "black", size = 2) +
  geom_point(aes(y = max_maxsize),    color = "black", size = 2) +
  geom_point(aes(y = median_maxsize), color = "black", size = 3, shape = 18) +
  coord_flip() +
  labs(x = NULL, y = "Maximum standard length (mm)") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 13),
    axis.text.y  = element_text(face = "italic", size = 13),
    axis.title.x = element_text(size = 14)
  )
Figure3

ggsave("Figure3.png", width = 6, height = 6, dpi = 600)

## Figure 4 — Partial dependence heatmap --------------------------------------
# Compute conditional effects (shared by Figures 4 & 5) 
# Computed once here to avoid the expensive call being repeated for each figure.

ce_size_genus <- conditional_effects(
  best_model,
  effects = "sizez:genusconfirmed",
  prob    = 0.90
)

n_by_genus <- diet %>%
  mutate(genusconfirmed = str_to_title(genusconfirmed)) %>%
  count(genusconfirmed, name = "n")

df_heatmap <- ce_size_genus[[1]] %>%
  mutate(genusconfirmed = str_to_title(genusconfirmed)) %>%
  left_join(n_by_genus, by = "genusconfirmed")

genus_order <- df_heatmap %>%
  group_by(genusconfirmed) %>%
  summarize(mean_prob = mean(estimate__), .groups = "drop") %>%
  arrange(mean_prob) %>%
  left_join(n_by_genus, by = "genusconfirmed") %>%
  mutate(genus_label = paste0(genusconfirmed, " (n=", n, ")")) %>%
  pull(genus_label)

df_heatmap <- df_heatmap %>%
  mutate(
    genus_label = paste0(genusconfirmed, " (n=", n, ")"),
    genus_label = factor(genus_label, levels = genus_order)
  )

x_min       <- -1.5
x_max       <-  2.0
size_breaks <- c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
size_labels <- paste0(round(size_breaks * size_sd + size_mean, 0),
                      "\n(", size_breaks, " SD)")

df_heatmap_trimmed <- df_heatmap %>%
  filter(sizez >= x_min & sizez <= x_max)

Figure4 <- ggplot(df_heatmap_trimmed, aes(x = sizez, y = genus_label)) +
  geom_tile(aes(fill = estimate__)) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "Predicted\nprobability",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = c("0", "0.25", "0.50", "0.75", "1.0")
  ) +
  scale_x_continuous(
    breaks = size_breaks,
    labels = size_labels,
    limits = c(x_min, x_max),
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Body size (mm)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 12, hjust = 0.5),
    axis.text.y  = element_text(face = "italic", size = 13),
    axis.title.x = element_text(size = 14, margin = margin(t = 8)),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    panel.grid   = element_blank()
  )
Figure4

ggsave("Figure4.png", width = 8, height = 7, dpi = 600)

## Figure 5 — Predicted probability curves by genus ---------------------------

genus_range <- diet %>%
  group_by(genusconfirmed) %>%
  summarize(
    min_size  = min(maxsize, na.rm = TRUE),
    max_size  = max(maxsize, na.rm = TRUE),
    size_range = max_size - min_size,
    n         = n(),
    .groups   = "drop"
  ) %>%
  filter(n >= 3, size_range > 0)   # exclude genera with < 3 obs or no size variation

# Back-transform sizez to mm using global constants and filter to observed size ranges
ce_size_genus_main <- ce_size_genus[[1]] %>%
  mutate(size_actual = sizez * size_sd + size_mean) %>%
  inner_join(genus_range, by = "genusconfirmed") %>%
  filter(size_actual >= min_size & size_actual <= max_size)

Figure5 <- ggplot(ce_size_genus_main, aes(x = size_actual, y = estimate__)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__),
              fill = "black", alpha = 0.2) +
  facet_wrap(~genusconfirmed, scales = "free_x", ncol = 4) +
  labs(
    x = "Body size (mm)",
    y = "Predicted probability of piscivory"
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(n=3))+
  theme_classic(base_size = 14) +
  theme(
    strip.text       = element_text(face = "italic", size = 14),
    strip.background = element_blank()
  )
Figure5

ggsave("Figure5.png", width = 9, height = 9, dpi = 600)

## Figure6: genus-level piscivory grouped by tooth type --------------------------

# Species-level observed piscivory (same 333 records as the tooth-type model)
min_n <- 1   # omit species with fewer records; set to 1 to show all

species_obs <- diet %>%
  filter(!is.na(raptorial_teeth)) %>%
  group_by(sp_fixed) %>%
  summarise(n               = n(),
            prop_fish       = mean(fish == "yes"),
            n_tooth         = n_distinct(raptorial_teeth),
            raptorial_teeth = first(raptorial_teeth),
            .groups = "drop")

stopifnot(sum(species_obs$n) == nobs(best_trait_model))  # same records as the model
stopifnot(all(species_obs$n_tooth == 1))                 # one tooth type per species

species_plot <- species_obs %>%
  filter(n >= min_n) %>%
  mutate(tooth_group = factor(if_else(raptorial_teeth == "raptorial",
                                      "Raptorial", "Villiform-only"),
                              levels = c("Villiform-only", "Raptorial")))

# For the caption: species and records shown per tooth type
species_plot %>% count(tooth_group, wt = NULL, name = "species") %>%
  left_join(species_plot %>% group_by(tooth_group) %>%
              summarise(records = sum(n)), by = "tooth_group") %>% print()

# Tooth-type model predictions (population level)
tooth_nd <- tibble(raptorial_teeth = factor(c("villiform-only", "raptorial"),
                                            levels = levels(diet$raptorial_teeth)))
tooth_draws <- posterior_epred(best_trait_model, newdata = tooth_nd, re_formula = NA)

tooth_pred <- tooth_nd %>%
  mutate(est   = apply(tooth_draws, 2, median),
         lower = apply(tooth_draws, 2, quantile, 0.025),
         upper = apply(tooth_draws, 2, quantile, 0.975),
         tooth_group = factor(c("Villiform-only", "Raptorial"),
                              levels = c("Villiform-only", "Raptorial")))

Figure6 <- ggplot(species_plot, aes(x = tooth_group, y = prop_fish)) +
  geom_point(aes(size = n), colour = "grey40", alpha = 0.5,
             position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  geom_pointrange(data = tooth_pred, aes(y = est, ymin = lower, ymax = upper),
                  colour = "black", size = 1, linewidth = 1.2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_size_area(max_size = 5, name = "Records (n)") +
  labs(x = "Tooth type", y = "Probability of piscivory") +
  theme_classic(base_size = 13)


ggsave("Figure6.png", fig_traits, width = 6, height = 6, dpi = 600)

## Supplementary Figure 2 — Posterior predictive check ------------------------
# type = "bars" is appropriate for binary (Bernoulli) outcomes, showing the
# observed 0/1 counts against replicated datasets from the posterior.

SuppFig2 <- pp_check(best_model, type = "bars", ndraws = 100) +
  labs(
    x = "Fish presence (0 = absent, 1 = present)",
    y = "Count"
  ) +
  theme_classic()

ggsave("SuppFig2.png", width = 12, height = 10, dpi = 600)

## Supplementary Figure 3 — Trace plots ----------------------------------------
# Trace plots for key fixed-effect parameters of the best-fit model.
# Full convergence assessed via R-hat <= 1.01 and ESS > 400.

SuppFig3 <- mcmc_trace(
  as.array(best_model),
  pars = c("b_Intercept", "b_sizez"),
  facet_args = list(ncol = 1)
) +
  theme_classic()
SuppFig3

ggsave("SuppFig3.png", width = 18, height = 12, dpi = 600)

SuppFig4 <- mcmc_rank_hist(
  as.array(best_model),
  pars = c("b_Intercept", "b_sizez"),
  facet_args = list()
) +
  theme_classic()
SuppFig4

ggsave("SuppFig4.png", width = 18, height = 12, dpi = 600)
## Supplementary Figure 4 — Caterpillar plot of study-level random effects ----

SuppFig5 <- best_model %>%
  spread_draws(r_ref[ref, term]) %>%
  filter(term == "Intercept") %>%
  median_qi(.width = 0.95) %>%
  arrange(r_ref) %>%
  mutate(ref = factor(ref, levels = ref)) %>%
  ggplot(aes(y = ref, x = r_ref, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = "Random intercept estimate (log-odds)",
    y = "Study"
  ) +
  theme_classic()
SuppFig5
ggsave("SuppFig5.png", width = 14, height = 18, dpi = 600)

message("All figures written to figures/")

# 
# SECTION 5: VERIFICATION — HAND CALCULATION CHECK##############################
# 

# Manually verifies the predicted probability for Gymnoscopelus at mean size
# and mean + 1 SD against the conditional_effects() output, to confirm the
# OR table is consistent with the predicted probability figures.

fe_raw <- fixef(best_model)

intercept    <- fe_raw["Intercept",                          "Estimate"]
size_slope   <- fe_raw["sizez",                              "Estimate"]
genus_effect <- fe_raw["genusconfirmedGymnoscopelus",        "Estimate"]
interaction  <- fe_raw["sizez:genusconfirmedGymnoscopelus",  "Estimate"]

log_odds_mean    <- intercept + genus_effect
log_odds_plus1sd <- intercept + size_slope + genus_effect + interaction

prob_mean    <- plogis(log_odds_mean)
prob_plus1sd <- plogis(log_odds_plus1sd)

# Compare against conditional_effects output (reuse ce_size_genus from above)
ce_gymno <- ce_size_genus[[1]] %>%
  filter(genusconfirmed == "Gymnoscopelus")

ce_at_mean <- ce_gymno %>%
  mutate(diff = abs(sizez - 0)) %>%
  arrange(diff) %>%
  slice(1)

ce_at_plus1sd <- ce_gymno %>%
  mutate(diff = abs(sizez - 1)) %>%
  arrange(diff) %>%
  slice(1)

cat("\n--- Hand-calculation verification: Gymnoscopelus ---\n")
cat(sprintf("At mean size   | hand: %.3f | model: %.3f\n",
            prob_mean,    ce_at_mean$estimate__))
cat(sprintf("At mean + 1 SD | hand: %.3f | model: %.3f\n",
            prob_plus1sd, ce_at_plus1sd$estimate__))
cat("Values should agree to within rounding error.\n\n")

# 
# SECTION 6: CITATIONS##########################################################
# 

citation()
citation("brms")
citation("tidyverse")
citation("loo")
citation("tidybayes")
citation("bayesplot")
citation("sf")
citation("marmap")
citation("rnaturalearth")
citation("flextable")
citation("officer")
citation("car")

# 
# SECTION 7: SESSION INFO#######################################################
# 
# Records exact R and package versions for reproducibility.
sessionInfo()
summary(MultiModelA_rs)
