# drosophila-thermal-plasticity

Data and R code for Vinton et al., Thermal plasticity erodes rapidly but independently of environmental predictability in experimental Drosophila populations.

Files


Fly_Image_Data.csv — per-fly morphology measurements (wing length, femur length, body area) by test temperature, generation, evolutionary environment, sex, and assay (PA = plasticity assay, CG = common garden). Raw images are not included; measurements were extracted from them.
Fly_Image_Analysis_mixed_models.R — mixed-model analyses; writes result tables to output/.
Fly_Image_Figures.R — generates all figures; run after the analysis script.


Requirements

R (developed under 4.5.2). Packages:

rinstall.packages(c("dplyr","tidyr","lme4","lmerTest","emmeans","ggplot2","gridExtra"))

Run

Set the working directory to this folder, then:

rsource("Fly_Image_Analysis_mixed_models.R")  # writes tables to output/
source("Fly_Image_Figures.R")                # writes figures

Released under the MIT License.
