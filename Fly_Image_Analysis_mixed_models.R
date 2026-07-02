################################################################################
## Thermal plasticity erodes rapidly but independently of environmental
## predictability in experimental Drosophila populations
## ---------------------------------------------------------------------------
## Morphological mixed-model analyses.
## Author: A. C. Vinton   (contact: anna.vinton@maine.edu)
##
## MODEL (fit separately per trait and per assay):
##   Trait ~ TestTemp * Generation * EvoEnv * Sex + (1 | Population)
##     - TestTemp (23,25,27,29,31) and Generation (4,7,10): continuous
##     - EvoEnv (Constant/Predictable/Random) and Sex (F/M): categorical fixed
##     - Population (9 evolved populations): random intercept
##     - Type III SS, sum-to-zero contrasts (lmerTest; Satterthwaite df)
##
## TRAITS: body area, wing length, femur length, wing/body ratio, wing/femur ratio
## ASSAYS: PA (plasticity assay), CG (common garden)
##
## REQUIRED PACKAGES: dplyr, tidyr, lme4, lmerTest, emmeans
##   First-time install:
##     install.packages(c("dplyr","tidyr","lme4","lmerTest","emmeans"))
##   Exact versions used are recorded by sessionInfo() at the end of this script.
##
## HOW TO RUN: place the accompanying data file (see data_path below) in the same
##   folder as this script, set that folder as the working directory, then source
##   the file. All results are written to ./output/.
##
################################################################################

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)

set.seed(20260101)
options(contrasts = c("contr.sum", "contr.poly"))

proj_dir  <- "/Users/anna.vinton/Library/CloudStorage/Dropbox/Fly_eco_evo_project/Fly_morphology"
data_path <- file.path(proj_dir, "Fly_Image_Data.csv")
out_dir   <- file.path(proj_dir, "output")
if (!dir.exists(out_dir)) dir.create(out_dir)

raw <- read.csv(data_path, stringsAsFactors = FALSE, colClasses = "character")
cat("Raw rows:", nrow(raw), "\n")

# ---- Clean ----
fly_data <- raw %>%
  filter(Assay %in% c("PA", "CG"),
         !(Parent_Gen %in% c("0", ""))) %>%
  mutate(
    Incubator_Temp = case_when(
      Incubator_Temp %in% c("Flux", "FLux", "FlUX") ~ "Predictable",
      Incubator_Temp == "27"                        ~ "Constant",
      Incubator_Temp == "Random"                    ~ "Random",
      TRUE                                          ~ NA_character_
    ),
    WB_Temp    = ifelse(WB_Temp == "31`", "31", WB_Temp),
    Parent_Gen = ifelse(Parent_Gen == "77", "7", Parent_Gen),
    Femur_Length = suppressWarnings(as.numeric(Femur_Length)),
    Wing_Length  = suppressWarnings(as.numeric(Wing_Length)),
    Body_Area    = suppressWarnings(as.numeric(Body_Area)),
    # drop one implausible body-area value (data-entry error)
    Body_Area    = ifelse(Pop == "7" & Parent_Gen == "7" &
                          Incubator_Temp == "Predictable" & WB_Temp == "27" &
                          Sex == "M" & Replicate == "2", NA, Body_Area),
    Wing_Length  = ifelse(tolower(No_wing)  == "x", NA, Wing_Length),
    Femur_Length = ifelse(tolower(No_femur) == "x", NA, Femur_Length)
  ) %>%
  filter(Pop %in% as.character(1:9),
         !is.na(Incubator_Temp),
         WB_Temp %in% c("23", "25", "27", "29", "31"))

cat("Rows after cleaning:", nrow(fly_data), "\n")

# ---- One row per fly (body and appendages are recorded on separate rows) ----
fid <- c("Replicate", "Pop", "Assay", "Parent_Gen",
         "Incubator_Temp", "WB_Temp", "Sex")

body_data <- fly_data %>%
  filter(!is.na(Body_Area)) %>%
  group_by(across(all_of(fid))) %>%
  summarize(Body_Area = mean(Body_Area), .groups = "drop")

appendage_data <- fly_data %>%
  filter(!is.na(Wing_Length) | !is.na(Femur_Length)) %>%
  group_by(across(all_of(fid))) %>%
  summarize(Wing_Length  = mean(Wing_Length,  na.rm = TRUE),
            Femur_Length = mean(Femur_Length, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(c(Wing_Length, Femur_Length), ~ ifelse(is.nan(.), NA, .)))

fly <- full_join(body_data, appendage_data, by = fid)
cat("Unique flies:", nrow(fly),
    "| body:", sum(!is.na(fly$Body_Area)),
    "| wing:", sum(!is.na(fly$Wing_Length)),
    "| femur:", sum(!is.na(fly$Femur_Length)), "\n")

# ---- Derived variables ----
fly <- fly %>%
  mutate(
    log_body   = log(Body_Area),
    log_wing   = log(Wing_Length),
    log_femur  = log(Femur_Length),
    wing_body  = log_wing - 0.5 * log_body,   # log(wing / sqrt(body area))
    wing_femur = log_wing - log_femur,
    TestTemp   = as.numeric(WB_Temp),
    Generation = as.numeric(Parent_Gen),
    EvoEnv     = factor(Incubator_Temp, levels = c("Constant", "Predictable", "Random")),
    Sex        = factor(Sex, levels = c("F", "M")),
    Population = factor(Pop, levels = as.character(1:9))
  )

PA <- fly %>% filter(Assay == "PA")
CG <- fly %>% filter(Assay == "CG")
cat("PA flies:", nrow(PA), "| CG flies:", nrow(CG), "\n\n")

traits <- c("Body_Area"        = "log_body",
            "Wing_Length"      = "log_wing",
            "Femur_Length"     = "log_femur",
            "Wing_Body_Ratio"  = "wing_body",
            "Wing_Femur_Ratio" = "wing_femur")

# ---- Helpers ----
fit_model <- function(resp, dat, cross_assay = FALSE) {
  rhs <- if (cross_assay)
    "TestTemp * Generation * EvoEnv * Sex * Assay + (1 | Population)"
  else
    "TestTemp * Generation * EvoEnv * Sex + (1 | Population)"
  # singular fits expected: population variance ~ 0
  lmer(as.formula(paste(resp, "~", rhs)), data = dat, REML = TRUE)
}

get_icc <- function(m) {
  vc <- as.data.frame(VarCorr(m))
  pv <- vc$vcov[vc$grp == "Population"]
  rv <- vc$vcov[vc$grp == "Residual"]
  pv / (pv + rv)
}

skewness <- function(x) { x <- x[!is.na(x)]; mean((x - mean(x))^3) / sd(x)^3 }
exkurt   <- function(x) { x <- x[!is.na(x)]; mean((x - mean(x))^4) / sd(x)^4 - 3 }

diagnostics <- function(m) {
  r <- residuals(m)
  sw <- if (length(r) <= 5000) shapiro.test(r) else shapiro.test(sample(r, 5000))
  data.frame(Shapiro_W = unname(sw$statistic), Shapiro_p = sw$p.value,
             Skewness  = skewness(r), Excess_Kurtosis = exkurt(r))
}

term_row <- function(a, term) {
  if (term %in% rownames(a))
    sprintf("F(%g,%.0f)=%.2f, p=%.4g",
            a[term, "NumDF"], a[term, "DenDF"], a[term, "F value"], a[term, "Pr(>F)"])
  else "term absent"
}

# ---- Plasticity assay ----
cat("################ PLASTICITY ASSAY (PA) ################\n\n")
icc_summary <- list()

for (lab in names(traits)) {
  resp <- traits[[lab]]
  m <- fit_model(resp, PA)
  a <- anova(m, type = 3)
  icc <- get_icc(m)
  icc_summary[[paste0("PA_", lab)]] <- icc

  cat("---- PA:", lab, "(n =", nobs(m), ", ICC =", sprintf("%.1f%%", 100*icc), ") ----\n")
  cat("  Temp x Gen              :", term_row(a, "TestTemp:Generation"), "\n")
  cat("  Temp x Gen x Env        :", term_row(a, "TestTemp:Generation:EvoEnv"), "\n")
  cat("  Temp x Gen x Sex        :", term_row(a, "TestTemp:Generation:Sex"), "\n")
  cat("  EvoEnv (main)           :", term_row(a, "EvoEnv"), "\n\n")

  write.csv(cbind(Term = rownames(a), as.data.frame(a)),
            file.path(out_dir, paste0("ANOVA_PA_", lab, ".csv")), row.names = FALSE)
}

# ---- Common garden ----
cat("################ COMMON GARDEN (CG) ################\n\n")

for (lab in names(traits)) {
  resp <- traits[[lab]]
  m <- fit_model(resp, CG)
  a <- anova(m, type = 3)
  icc <- get_icc(m)
  icc_summary[[paste0("CG_", lab)]] <- icc

  cat("---- CG:", lab, "(n =", nobs(m), ", ICC =", sprintf("%.1f%%", 100*icc), ") ----\n")
  cat("  EvoEnv (main)           :", term_row(a, "EvoEnv"), "\n")
  cat("  Temp x Env              :", term_row(a, "TestTemp:EvoEnv"), "\n")
  cat("  Env x Sex               :", term_row(a, "EvoEnv:Sex"), "\n")
  cat("  Temp x Env x Sex        :", term_row(a, "TestTemp:EvoEnv:Sex"), "\n")
  cat("  Gen x Env x Sex         :", term_row(a, "Generation:EvoEnv:Sex"), "\n")
  cat("  Temp x Gen x Env x Sex  :", term_row(a, "TestTemp:Generation:EvoEnv:Sex"), "\n\n")

  write.csv(cbind(Term = rownames(a), as.data.frame(a)),
            file.path(out_dir, paste0("ANOVA_CG_", lab, ".csv")), row.names = FALSE)
}

# ---- Genetic-basis test (cross-assay; n.s. Temp:Gen:Assay = genetic) ----
cat("################ GENETIC-BASIS TEST (Temp x Gen x Assay) ################\n\n")

gb_rows <- list()
for (lab in names(traits)) {
  resp <- traits[[lab]]
  m <- fit_model(resp, fly, cross_assay = TRUE)
  a <- anova(m, type = 3)
  p_assay <- a["TestTemp:Generation:Assay", "Pr(>F)"]
  gb_rows[[lab]] <- data.frame(Trait = lab, p_TempGenAssay = p_assay)
  verdict <- if (p_assay >= 0.05) "n.s. -> equivalent in PA & CG (genetic)"
             else                  "significant -> differs between PA & CG"
  cat("  ", lab, ":", term_row(a, "TestTemp:Generation:Assay"), "|", verdict, "\n")
  write.csv(cbind(Term = rownames(a), as.data.frame(a)),
            file.path(out_dir, paste0("ANOVA_crossassay_", lab, ".csv")), row.names = FALSE)
}
write.csv(do.call(rbind, gb_rows), file.path(out_dir, "fig_genetic_pvals.csv"), row.names = FALSE)
cat("\n")

# ---- Figure inputs: slopes and means (written to ./output for the figures script) ----
cat("################ FIGURE INPUTS (slopes & means) ################\n\n")
tidy_tr <- function(x) { x <- as.data.frame(x); class(x) <- "data.frame"; names(x)[names(x) == "TestTemp.trend"] <- "Slope"; x }

# (a) Slopes by generation x assay  [slopes-by-generation bars; genetic-basis lines]
rows <- list()
for (lab in names(traits)) for (asy in c("PA", "CG")) {
  m <- fit_model(traits[[lab]], if (asy == "PA") PA else CG)
  s <- tidy_tr(summary(emtrends(m, ~ Generation, var = "TestTemp",
                                at = list(Generation = c(4, 7, 10)))))
  s$Trait <- lab; s$Assay <- asy
  rows[[paste(lab, asy)]] <- s
}
slopes_gen_assay <- do.call(rbind, rows)
write.csv(slopes_gen_assay, file.path(out_dir, "fig_slopes_gen_assay.csv"), row.names = FALSE)
cat("---- slopes by generation x assay ----\n"); print(slopes_gen_assay, row.names = FALSE)

# (b) Femur slopes by sex x generation x assay  [sex-specific femur figure]
rows <- list()
for (asy in c("PA", "CG")) {
  m <- fit_model("log_femur", if (asy == "PA") PA else CG)
  s <- tidy_tr(summary(emtrends(m, ~ Sex | Generation, var = "TestTemp",
                                at = list(Generation = c(4, 7, 10)))))
  s$Assay <- asy
  rows[[asy]] <- s
}
write.csv(do.call(rbind, rows), file.path(out_dir, "fig_femur_sex_assay.csv"), row.names = FALSE)

# (c) Slopes by environment x generation (PA, gen 4 vs 10)  [slope-by-environment figure]
rows <- list()
for (lab in c("Body_Area", "Wing_Length", "Femur_Length")) {
  m <- fit_model(traits[[lab]], PA)
  s <- tidy_tr(summary(emtrends(m, ~ Generation | EvoEnv, var = "TestTemp",
                                at = list(Generation = c(4, 10)))))
  s$Trait <- lab
  rows[[lab]] <- s
}
write.csv(do.call(rbind, rows), file.path(out_dir, "fig_slopes_env.csv"), row.names = FALSE)

# (d) CG mean morphology by environment (log scale)  [genetic divergence in means]
rows <- list()
for (lab in c("Body_Area", "Wing_Length", "Femur_Length")) {
  m <- fit_model(traits[[lab]], CG)
  e <- as.data.frame(summary(emmeans(m, ~ EvoEnv)))
  class(e) <- "data.frame"
  e$Trait <- lab
  rows[[lab]] <- e
}
write.csv(do.call(rbind, rows), file.path(out_dir, "fig_cg_means_env.csv"), row.names = FALSE)
cat("\n")

# ---- Diagnostics and summary exports ----
cat("################ MODEL DIAGNOSTICS (residuals) ################\n\n")
diag_rows <- list()
for (lab in names(traits)) for (asy in c("PA","CG")) {
  dat <- if (asy == "PA") PA else CG
  m <- fit_model(traits[[lab]], dat)
  d <- diagnostics(m); d$Trait <- lab; d$Assay <- asy
  diag_rows[[paste(asy, lab)]] <- d
}
diag_df <- do.call(rbind, diag_rows)[, c("Trait","Assay","Shapiro_W","Shapiro_p","Skewness","Excess_Kurtosis")]
print(diag_df, row.names = FALSE)
write.csv(diag_df, file.path(out_dir, "Diagnostics_summary.csv"), row.names = FALSE)

icc_df <- data.frame(Model = names(icc_summary),
                     ICC_percent = round(100 * unlist(icc_summary), 2))
write.csv(icc_df, file.path(out_dir, "ICC_summary.csv"), row.names = FALSE)
cat("\nICC summary:\n"); print(icc_df, row.names = FALSE)

cat("\nResults written to ./", out_dir, "/\n", sep = "")

cat("\n################ sessionInfo() ################\n")
print(sessionInfo())
