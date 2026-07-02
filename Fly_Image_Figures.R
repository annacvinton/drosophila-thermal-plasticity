################################################################################
## Thermal plasticity erodes rapidly but independently of environmental
## predictability in experimental Drosophila populations
## ---------------------------------------------------------------------------
## Publication figures. Author: A. C. Vinton   (contact: anna.vinton@maine.edu)
##
## Reaction-norm panels are drawn from the data (natural units). Slope and mean
## figures read the tables written by Fly_Image_Analysis_mixed_models.R, so RUN
## THAT SCRIPT FIRST (it fills ./output/ with the fig_*.csv inputs). Figures are
## written to ./figures/.
##
## REQUIRED PACKAGES: dplyr, tidyr, ggplot2, gridExtra
##   install.packages(c("dplyr","tidyr","ggplot2","gridExtra"))
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

proj_dir  <- "."
data_path <- file.path(proj_dir, "Fly_Image_Data.csv")
out_dir   <- file.path(proj_dir, "output")   # analysis tables are read from here
fig_dir   <- proj_dir                        # figures saved directly into Fly_morphology
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# ---- Load & clean (matches analysis script) ----
raw <- read.csv(data_path, stringsAsFactors = FALSE, colClasses = "character")
fly <- raw %>%
  filter(Assay %in% c("PA", "CG"), !(Parent_Gen %in% c("0", ""))) %>%
  mutate(
    Incubator_Temp = case_when(
      Incubator_Temp %in% c("Flux", "FLux", "FlUX") ~ "Predictable",
      Incubator_Temp == "27"                        ~ "Constant",
      Incubator_Temp == "Random"                    ~ "Random",
      TRUE                                          ~ NA_character_),
    WB_Temp    = ifelse(WB_Temp == "31`", "31", WB_Temp),
    Parent_Gen = ifelse(Parent_Gen == "77", "7", Parent_Gen),
    Femur_Length = suppressWarnings(as.numeric(Femur_Length)),
    Wing_Length  = suppressWarnings(as.numeric(Wing_Length)),
    Body_Area    = suppressWarnings(as.numeric(Body_Area)),
    Body_Area    = ifelse(Pop == "7" & Parent_Gen == "7" &
                          Incubator_Temp == "Predictable" & WB_Temp == "27" &
                          Sex == "M" & Replicate == "2", NA, Body_Area),
    Wing_Length  = ifelse(tolower(No_wing)  == "x", NA, Wing_Length),
    Femur_Length = ifelse(tolower(No_femur) == "x", NA, Femur_Length)
  ) %>%
  filter(Pop %in% as.character(1:9), !is.na(Incubator_Temp),
         WB_Temp %in% c("23", "25", "27", "29", "31"))

fid <- c("Replicate","Pop","Assay","Parent_Gen","Incubator_Temp","WB_Temp","Sex")
body_data <- fly %>% filter(!is.na(Body_Area)) %>%
  group_by(across(all_of(fid))) %>% summarize(Body_Area = mean(Body_Area), .groups = "drop")
app_data <- fly %>% filter(!is.na(Wing_Length) | !is.na(Femur_Length)) %>%
  group_by(across(all_of(fid))) %>%
  summarize(Wing_Length = mean(Wing_Length, na.rm = TRUE),
            Femur_Length = mean(Femur_Length, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(c(Wing_Length, Femur_Length), ~ ifelse(is.nan(.), NA, .)))

dat <- full_join(body_data, app_data, by = fid) %>%
  mutate(
    Wing_Length  = Wing_Length  * 0.05,     # ImageJ units -> mm  (1 unit = 0.05 mm)
    Femur_Length = Femur_Length * 0.05,
    Body_Area    = Body_Area    * 0.05^2,   # -> mm^2
    Temp = as.numeric(WB_Temp),
    Gen  = factor(Parent_Gen, levels = c("4", "7", "10")),
    Env  = factor(Incubator_Temp, levels = c("Constant", "Predictable", "Random")),
    Sex  = factor(Sex, levels = c("F", "M")),
    wb_nat = Wing_Length / sqrt(Body_Area),   # scale-free (unchanged by conversion)
    wf_nat = Wing_Length / Femur_Length)

# ---- Shared style ----
gen_colors   <- c("4" = "#377eb8", "7" = "#ff7f00", "10" = "#e41a1c")
gen2_colors  <- c("4" = "#377eb8", "10" = "#e41a1c")
sex_colors   <- c("F" = "#e41a1c", "M" = "#377eb8")
assay_colors <- c("PA" = "#ff7f00", "CG" = "#377eb8")
env_colors   <- c("Constant" = "#1f78b4", "Predictable" = "#33a02c", "Random" = "#b15928")
trait_lab    <- c(Wing_Length = "Wing length", Femur_Length = "Femur length", Body_Area = "Body area")
mm_labs      <- c(Wing_Length = "Wing length (mm)", Femur_Length = "Femur length (mm)",
                  Body_Area = "Body area (mm\u00b2)")
base_theme   <- theme_classic(base_size = 12) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey90"))

long_traits <- function(d) d %>%
  select(Temp, Gen, Env, Sex, Assay, Wing_Length, Femur_Length, Body_Area) %>%
  pivot_longer(c(Wing_Length, Femur_Length, Body_Area), names_to = "Trait", values_to = "val") %>%
  filter(!is.na(val)) %>%
  mutate(Trait = factor(mm_labs[Trait], levels = mm_labs))

# ---- Reaction-norm figures (data, natural units) ----
tsr <- long_traits(dat) %>% group_by(Trait, Temp) %>%
  summarize(m = mean(val), se = sd(val) / sqrt(n()), .groups = "drop")
fig_tsr <- ggplot(tsr, aes(Temp, m)) +
  geom_line() + geom_point() + geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.3) +
  facet_wrap(~ Trait, scales = "free_y") +
  labs(x = "Test temperature (\u00b0C)", y = NULL) + base_theme
ggsave(file.path(fig_dir, "Fig2_TSR.pdf"), fig_tsr, width = 11, height = 4)

rn <- long_traits(dat) %>% group_by(Trait, Gen, Temp) %>%
  summarize(m = mean(val), se = sd(val) / sqrt(n()), .groups = "drop")
fig_rn <- ggplot(rn, aes(Temp, m, color = Gen, group = Gen)) +
  geom_line() + geom_point() + geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.3, alpha = 0.6) +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_color_manual(values = gen_colors, name = "Generation") +
  labs(x = "Test temperature (\u00b0C)", y = NULL) + base_theme
ggsave(file.path(fig_dir, "Fig3_ReactionNorms.pdf"), fig_rn, width = 11, height = 4.5)

rr <- dat %>% select(Temp, Assay, wb_nat, wf_nat) %>%
  pivot_longer(c(wb_nat, wf_nat), names_to = "Ratio", values_to = "val") %>%
  filter(is.finite(val)) %>%
  mutate(Ratio = ifelse(Ratio == "wb_nat", "Wing / sqrt(body)", "Wing / femur")) %>%
  group_by(Ratio, Assay, Temp) %>% summarize(m = mean(val), se = sd(val) / sqrt(n()), .groups = "drop")
fig_ratios <- ggplot(rr, aes(Temp, m, color = Assay, group = Assay)) +
  geom_line() + geom_point() + geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.3, alpha = 0.6) +
  facet_wrap(~ Ratio, scales = "free_y") +
  scale_color_manual(values = assay_colors, labels = c("Common garden", "Plasticity assay"), name = NULL) +
  labs(x = "Test temperature (\u00b0C)", y = NULL) + base_theme
ggsave(file.path(fig_dir, "FigS1_Ratios.pdf"), fig_ratios, width = 9, height = 4.5)

# ---- Slope / mean figures (read analysis output) ----
need <- file.path(out_dir, c("fig_slopes_gen_assay.csv", "fig_femur_sex_assay.csv",
                             "fig_slopes_env.csv", "fig_cg_means_env.csv", "fig_genetic_pvals.csv"))
if (all(file.exists(need))) {

  sg <- read.csv(need[1]) %>% filter(Trait %in% names(trait_lab)) %>%
    mutate(Generation = factor(Generation, levels = c(4, 7, 10)),
           Assay = factor(Assay, levels = c("CG", "PA")),
           Trait = factor(trait_lab[Trait], levels = trait_lab))
  fig_slopes_gen <- ggplot(sg, aes(Generation, Slope)) +
    geom_col(fill = "grey65", width = 0.65) +
    geom_errorbar(aes(ymin = Slope - SE, ymax = Slope + SE), width = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    facet_grid(Trait ~ Assay, scales = "free_y") +
    labs(x = "Generation", y = "Temperature-size slope (log-units / \u00b0C)") +
    theme_classic(base_size = 12) + theme(strip.background = element_rect(fill = "grey90"))
  ggsave(file.path(fig_dir, "Fig4_Slopes_by_generation.pdf"), fig_slopes_gen, width = 8, height = 8)

  fs <- read.csv(need[2]) %>%
    mutate(Generation = factor(Generation, levels = c(4, 7, 10)),
           Assay = factor(Assay, levels = c("CG", "PA")))
  fig_femur_sex <- ggplot(fs, aes(Generation, Slope, color = Sex, group = Sex)) +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = Slope - SE, ymax = Slope + SE), width = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
    facet_wrap(~ Assay) +
    scale_color_manual(values = sex_colors, labels = c("Females", "Males"), name = NULL) +
    labs(x = "Generation", y = "Femur slope (log-units / \u00b0C)") + base_theme
  ggsave(file.path(fig_dir, "Fig5_Femur_by_sex.pdf"), fig_femur_sex, width = 8, height = 4.5)

  gp <- read.csv(need[5]) %>% filter(Trait %in% names(trait_lab)) %>%
    mutate(label = paste0("Temp\u00d7Gen\u00d7Assay\np = ", formatC(p_TempGenAssay, format = "f", digits = 3)),
           Trait = factor(trait_lab[Trait], levels = trait_lab))
  sg2 <- read.csv(need[1]) %>% filter(Trait %in% names(trait_lab)) %>%
    mutate(Trait = factor(trait_lab[Trait], levels = trait_lab),
           Assay = factor(Assay, levels = c("PA", "CG")), Gen = as.numeric(Generation))
  fig_gen_basis <- ggplot(sg2, aes(Gen, Slope, color = Assay, group = Assay)) +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = Slope - SE, ymax = Slope + SE), width = 0.3) +
    geom_text(data = gp, aes(x = 7, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.3, size = 3) +
    facet_wrap(~ Trait, scales = "free_y") +
    scale_color_manual(values = assay_colors, labels = c("Plasticity assay", "Common garden"), name = NULL) +
    scale_x_continuous(breaks = c(4, 7, 10)) +
    labs(x = "Generation", y = "Temperature-size slope (log-units / \u00b0C)") + base_theme
  ggsave(file.path(fig_dir, "Fig8_Genetic_basis.pdf"), fig_gen_basis, width = 11, height = 4.5)

  se_env <- read.csv(need[3]) %>%
    mutate(Generation = factor(Generation, levels = c(4, 10)),
           EvoEnv = factor(EvoEnv, levels = c("Constant", "Predictable", "Random")),
           Trait = factor(trait_lab[Trait], levels = trait_lab))
  fig_slopes_env <- ggplot(se_env, aes(EvoEnv, Slope, fill = Generation)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    geom_errorbar(aes(ymin = Slope - SE, ymax = Slope + SE), position = position_dodge(0.75), width = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    facet_wrap(~ Trait, scales = "free_y") +
    scale_fill_manual(values = gen2_colors, name = "Generation") +
    labs(x = "Evolutionary environment", y = "Temperature-size slope (log-units / \u00b0C)") +
    base_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(fig_dir, "Fig6_Slopes_by_environment.pdf"), fig_slopes_env, width = 11, height = 4.5)

  cm <- read.csv(need[4]) %>%
    mutate(scale    = ifelse(Trait == "Body_Area", 0.05^2, 0.05),   # -> mm / mm^2
           mean_nat = exp(emmean) * scale,
           lo       = exp(emmean - SE) * scale,
           hi       = exp(emmean + SE) * scale,
           EvoEnv   = factor(EvoEnv, levels = c("Constant", "Predictable", "Random")),
           Trait    = factor(mm_labs[Trait], levels = mm_labs))
  fig_cg_means <- ggplot(cm, aes(EvoEnv, mean_nat, color = EvoEnv)) +
    geom_point(size = 3, show.legend = FALSE) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, show.legend = FALSE) +
    facet_wrap(~ Trait, scales = "free_y") +
    scale_color_manual(values = env_colors) +
    labs(x = "Evolutionary environment", y = "Common-garden mean") +
    theme_classic(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey90"), axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(fig_dir, "Fig7_CG_means.pdf"), fig_cg_means, width = 11, height = 4.5)

  cat("Slope/mean figures written.\n")
} else {
  cat("NOTE: run Fly_Image_Analysis_mixed_models.R first (it writes the fig_*.csv inputs),\n",
      "then re-run this script for the slope and mean figures.\n", sep = "")
}

cat("Figures written to ./", fig_dir, "/\n", sep = "")
