# drosophila-thermal-plasticity

Data and code for Vinton AC, Vasquez T, Chen J, Ibbott B, Pillas JT, Koch D, Gascoigne SJL. *Environmental predictability shapes transgenerational, but not evolved, thermal plasticity in Drosophila.*

## Files

- `Fly_Image_Data.csv` — per-fly morphological measurements (wing length, femur length, body area) for nine experimentally evolved *Drosophila melanogaster* populations, by evolutionary regime, assay (PA = plasticity assay; CG = common garden), generation, test temperature, sex, and replicate. Measurements are in ImageJ units (1 unit = 0.05 mm; 1 unit² = 0.0025 mm²). Raw images are not included.
- `analysis_population_level.py` — reproduces every statistic, table (Tables 1–5, S1–S7) and figure (Figures 2–6, S1) in the manuscript from the raw data file.


## Key columns in `Fly_Image_Data.csv`

| Column | Meaning |
|---|---|
| `Pop` | Population 1–9 (1–3 random regime; 4–6 constant; 7–9 predictable) |
| `Incubator_Temp` | Evolutionary regime of the population: `27` = constant 27 °C, `Flux` = predictably fluctuating, `Random` = randomly fluctuating |
| `Assay` | `PA` = plasticity assay (reared directly from the evolutionary environment); `CG` = common garden (two generations at 27 °C) |
| `Parent_Gen` | Assay generation (4, 7, 10) |
| `WB_Temp` | Test (developmental) temperature, °C: 23, 25, 27, 29, 31 |
| `Sex` | F / M |
| `Replicate` | Replicate vial within a cell |
| `Wing_Length`, `Femur_Length`, `Body_Area` | Measurements (ImageJ units); blank if not measured |
| `No_wing`, `No_femur` | `x` marks an unusable wing or femur measurement |

## Requirements

Python 3 with pandas, numpy, scipy and matplotlib (tested with Python 3.12, pandas 3.0, NumPy 2.4, SciPy 1.17, Matplotlib 3.10).

## Run

```
python analysis_population_level.py Fly_Image_Data.csv
```

The script cleans the raw measurements to one row per fly, estimates each population's thermal reaction norm at generations 4, 7 and 10, tests for change in plasticity and for regime effects with the population as the unit of replication, compares mean morphology between the plasticity assay and the common garden, and writes all outputs to `results/` (CSV result tables, `tables.md`, `tables_supplement.md`, and `figures/` as PDF and PNG). It prints the headline statistics as a check.

## Data correction (26 October 2026)

Twenty-four cells in `Fly_Image_Data.csv` were corrected relative to the file originally deposited: 20 rows (population 7, plasticity assay, generation 10, females, replicates 4–5) had the regime recorded as `27` instead of `Flux`, and 4 rows (common garden, generation 4, 25 °C, females, replicates 4–5) had the population recorded as 9 instead of 1. No measurement values were changed. The script asserts that every population belongs to exactly one regime.

## License

MIT (see `LICENSE`).
