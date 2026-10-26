#!/usr/bin/env python3
"""
Population-level analysis of thermal plasticity in experimentally evolved
Drosophila melanogaster populations.

Reproduces every statistic, table, and figure (Figures 2-6, S1; Tables 1-5, S1)
reported in the manuscript from the raw measurement file Fly_Image_Data.csv.

Usage:  python analysis_population_level.py [path/to/Fly_Image_Data.csv]
Output: ./results/  (CSV tables, markdown tables, figures as PDF and PNG)

Requires: pandas, numpy, scipy, matplotlib.
"""
import sys, os, warnings
from itertools import combinations
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
warnings.filterwarnings("ignore")

RAW = sys.argv[1] if len(sys.argv) > 1 else "Fly_Image_Data.csv"
OUT = "results"; FIG = os.path.join(OUT, "figures"); os.makedirs(FIG, exist_ok=True)

REGIME = {"1": "Random", "2": "Random", "3": "Random", "4": "Constant", "5": "Constant",
          "6": "Constant", "7": "Predictable", "8": "Predictable", "9": "Predictable"}
REG = ["Constant", "Predictable", "Random"]
COL = {"Constant": "#0072B2", "Predictable": "#009E73", "Random": "#E69F00"}
TRAITS = {"Wing": "log_wing", "Femur": "log_femur", "Body": "log_body"}
SCALE = {"Wing": 0.05, "Femur": 0.05, "Body": 0.0025}      # pixel -> mm, pixel^2 -> mm^2
UNIT = {"Wing": "mm", "Femur": "mm", "Body": "mm²"}
TRL = {"Wing": "Wing length", "Femur": "Femur length", "Body": "Body area"}
ASL = {"PA": "Plasticity assay", "CG": "Common garden"}
GENS = [4.0, 7.0, 10.0]

# ----------------------------------------------------------------------------
# 1. Data cleaning: raw image measurements -> one row per fly
# ----------------------------------------------------------------------------
def load_and_clean(path):
    raw = pd.read_csv(path, dtype=str, keep_default_na=False)
    raw.columns = [c.strip() for c in raw.columns]
    d = raw[raw["Assay"].isin(["PA", "CG"])].copy()
    d["Parent_Gen"] = d["Parent_Gen"].replace({"77": "7"})           # typo
    d = d[d["Parent_Gen"].isin(["4", "7", "10"])]                    # assay generations
    d["Incubator_Temp"] = d["Incubator_Temp"].map(
        lambda x: "Predictable" if x.lower() == "flux" else "Constant" if x == "27"
        else "Random" if x == "Random" else None)
    d["WB_Temp"] = d["WB_Temp"].replace({"31`": "31"})               # typo
    for c in ["Femur_Length", "Wing_Length", "Body_Area"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    # one body-area measurement error (Pop 7, gen 7, 27 °C, male, replicate 2), as in the original analysis
    m = ((d.Pop == "7") & (d.Parent_Gen == "7") & (d.Incubator_Temp == "Predictable")
         & (d.WB_Temp == "27") & (d.Sex == "M") & (d.Replicate == "2"))
    d.loc[m, "Body_Area"] = np.nan
    d.loc[d["No_wing"].str.lower() == "x", "Wing_Length"] = np.nan
    d.loc[d["No_femur"].str.lower() == "x", "Femur_Length"] = np.nan
    d = d[d.Pop.isin(list(REGIME)) & d.Incubator_Temp.notna()
          & d.WB_Temp.isin(["23", "25", "27", "29", "31"])]
    # collapse duplicate images of the same fly to one row
    fid = ["Replicate", "Pop", "Assay", "Parent_Gen", "Incubator_Temp", "WB_Temp", "Sex"]
    body = d[d.Body_Area.notna()].groupby(fid, as_index=False)["Body_Area"].mean()
    app = d[d.Wing_Length.notna() | d.Femur_Length.notna()].groupby(fid, as_index=False).agg(
        Wing_Length=("Wing_Length", "mean"), Femur_Length=("Femur_Length", "mean"))
    fly = body.merge(app, on=fid, how="outer")
    fly["TestTemp"] = fly.WB_Temp.astype(float); fly["Generation"] = fly.Parent_Gen.astype(float)
    for tr, col in TRAITS.items():
        fly[col] = np.log(fly[{"Wing": "Wing_Length", "Femur": "Femur_Length", "Body": "Body_Area"}[tr]])
    fly["log_wb"] = fly.log_wing - 0.5 * fly.log_body      # log(wing length / sqrt(body area)): thermal allometry
    fly["log_wf"] = fly.log_wing - fly.log_femur           # log(wing length / femur length)
    fly["Regime"] = fly.Pop.map(REGIME)
    assert fly.groupby("Pop").Regime.nunique().max() == 1, "a population appears in more than one regime"
    return fly

fly = load_and_clean(RAW)
fly.to_csv(os.path.join(OUT, "flies_clean.csv"), index=False)
print(f"Flies: {len(fly)}  (PA {int((fly.Assay=='PA').sum())}, CG {int((fly.Assay=='CG').sum())}); "
      f"populations per regime: {fly.groupby('Regime').Pop.nunique().to_dict()}")

# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def ci95(x):
    x = np.asarray(x, float); n = len(x)
    return stats.t.ppf(0.975, n - 1) * x.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan

def fp(p): return "< 0.001" if p < 0.001 else f"{p:.3f}" if p < 0.01 else f"{p:.2f}"

def slope(sub, col, asy, pop, gen, min_n=5, min_temps=3):
    s = sub[(sub.Assay == asy) & (sub.Pop == pop) & (sub.Generation == gen) & sub[col].notna()]
    if len(s) < min_n or s.TestTemp.nunique() < min_temps: return np.nan, np.nan, len(s)
    r = stats.linregress(s.TestTemp, s[col]); return r.slope, r.stderr, len(s)

def rate_table(sub):
    """Two-stage estimate: per-population reaction-norm slope at each generation,
    then the rate of change of that slope across generations."""
    slopes, rates = [], []
    for tr, col in TRAITS.items():
        for asy in ["PA", "CG"]:
            for pop in REGIME:
                sl = []
                for g in GENS:
                    b, se, n = slope(sub, col, asy, pop, g); sl.append(b)
                    if not np.isnan(b): slopes.append(dict(Trait=tr, Assay=asy, Pop=pop, Regime=REGIME[pop], Gen=int(g), Slope=b, SE=se, n=n))
                if not any(np.isnan(sl)):
                    rates.append(dict(Trait=tr, Assay=asy, Pop=pop, Regime=REGIME[pop], S4=sl[0], S7=sl[1], S10=sl[2],
                                      Rate=stats.linregress([4, 7, 10], sl).slope))
    return pd.DataFrame(slopes), pd.DataFrame(rates)

def anova_by_regime(vals):
    v = vals.dropna(); g = [v[[p for p in v.index if REGIME[p] == r]].values for r in REG]
    if min(len(x) for x in g) < 2: return np.nan, np.nan
    return stats.f_oneway(*g)

def pairwise(vals):
    g = {r: vals[[p for p in vals.index if REGIME[p] == r]].values for r in REG}
    mse = np.mean([np.var(v, ddof=1) for v in g.values()]); hw = stats.t.ppf(0.975, 6) * np.sqrt(2 * mse / 3)
    return {f"{b[0]}−{a[0]}": (g[b].mean() - g[a].mean(), hw) for a, b in combinations(REG, 2)}

# ----------------------------------------------------------------------------
# 2. Plasticity change (Tables 1-3)
# ----------------------------------------------------------------------------
S, RP = rate_table(fly); RF = rate_table(fly[fly.Sex == "F"])[1]; RM = rate_table(fly[fly.Sex == "M"])[1]
S.to_csv(os.path.join(OUT, "pop_slopes.csv"), index=False)
for n, R in [("pooled", RP), ("female", RF), ("male", RM)]: R.to_csv(os.path.join(OUT, f"rates_{n}.csv"), index=False)

t1 = ["| Trait | Assay | Sex | n pops | Rate (log-units °C$^{-1}$ gen$^{-1}$) | 95% CI | t | p | Pops flattening | Mean change gen 4→10 |", "|---|---|---|---|---|---|---|---|---|---|"]
for tr in TRAITS:
    for asy in ["PA", "CG"]:
        for lab, R in [("Both", RP), ("F", RF), ("M", RM)]:
            w = R[(R.Trait == tr) & (R.Assay == asy)]; n = len(w)
            if n < 4: t1.append(f"| {tr} | {asy} | {lab} | {n} | — | — | — | — | — | — |"); continue
            t, p = stats.ttest_1samp(w.Rate, 0); m = w.Rate.mean(); c = ci95(w.Rate)
            pct = 100 * (w.S10.mean() - w.S4.mean()) / abs(w.S4.mean())
            t1.append(f"| {tr} | {asy} | {lab} | {n} | {m:+.5f} | {m-c:+.5f}, {m+c:+.5f} | {t:.2f} | {fp(p)} | {(w.Rate>0).sum()}/{n} | {pct:+.0f}% |")
            if lab == "Both": print(f"  {tr:5s} {asy}: rate {m:+.5f}  p = {p:.3f}  ({(w.Rate>0).sum()}/{n} flattening)")

t2 = ["| Trait | Assay | Constant | Predictable | Random | F₂,₆ | p | P − C (95% CI) | R − C (95% CI) | R − P (95% CI) |", "|---|---|---|---|---|---|---|---|---|---|"]
for tr in TRAITS:
    for asy in ["PA", "CG"]:
        w = RP[(RP.Trait == tr) & (RP.Assay == asy)].set_index("Pop").Rate
        F, p = anova_by_regime(w); pw = pairwise(w)
        mean = {r: w[[q for q in w.index if REGIME[q] == r]].mean() for r in REG}
        t2.append(f"| {tr} | {asy} | {mean['Constant']:+.5f} | {mean['Predictable']:+.5f} | {mean['Random']:+.5f} | {F:.2f} | {fp(p)} | "
                  + " | ".join(f"{d:+.5f} [{d-h:+.5f}, {d+h:+.5f}]" for d, h in pw.values()) + " |")

t3 = ["| Trait | CG rate | 95% CI | t₈ | p | Pops flattening in CG | PA − CG rate | 95% CI | t₈ | p | Pops PA steeper |", "|---|---|---|---|---|---|---|---|---|---|---|"]
for tr in TRAITS:
    pa = RP[(RP.Trait == tr) & (RP.Assay == "PA")].set_index("Pop").Rate; cg = RP[(RP.Trait == tr) & (RP.Assay == "CG")].set_index("Pop").Rate
    tc, pc = stats.ttest_1samp(cg, 0); d = (pa - cg.loc[pa.index]); td, pd_ = stats.ttest_1samp(d, 0)
    t3.append(f"| {tr} | {cg.mean():+.5f} | {cg.mean()-ci95(cg):+.5f}, {cg.mean()+ci95(cg):+.5f} | {tc:.2f} | {fp(pc)} | {(cg>0).sum()}/9 | "
              f"{d.mean():+.5f} | {d.mean()-ci95(d):+.5f}, {d.mean()+ci95(d):+.5f} | {td:.2f} | {fp(pd_)} | {(d>0).sum()}/9 |")

# ----------------------------------------------------------------------------
# 3. Mean morphology by regime, plasticity assay vs common garden (Tables 4-5)
# ----------------------------------------------------------------------------
def popmeans(sub, col, asy):
    m = sub[(sub.Assay == asy) & sub[col].notna()].groupby("Pop")[col].mean(); m.index = m.index.astype(str); return m

rows = []
for tr, col in TRAITS.items():
    for sex, lab in [(None, "Both"), ("F", "F"), ("M", "M")]:
        sub = fly if sex is None else fly[fly.Sex == sex]
        mPA, mCG = popmeans(sub, col, "PA"), popmeans(sub, col, "CG")
        FPA, pPA = anova_by_regime(mPA); FCG, pCG = anova_by_regime(mCG); Fi, pi = anova_by_regime(mPA - mCG)
        nat = lambda m, r: np.exp(m[[p for p in m.index if REGIME[p] == r]].mean()) * SCALE[tr]
        pw = pairwise(mPA)
        rows.append(dict(Trait=tr, Sex=lab, **{f"PA_{r[0]}": nat(mPA, r) for r in REG}, F_PA=FPA, p_PA=pPA,
                         **{f"CG_{r[0]}": nat(mCG, r) for r in REG}, F_CG=FCG, p_CG=pCG, F_int=Fi, p_int=pi,
                         **{f"PA_{k}_pct": 100 * (np.exp(d) - 1) for k, (d, h) in pw.items()},
                         **{f"PA_{k}_lo": 100 * (np.exp(d - h) - 1) for k, (d, h) in pw.items()},
                         **{f"PA_{k}_hi": 100 * (np.exp(d + h) - 1) for k, (d, h) in pw.items()}))
        if tr == "Body" and lab == "F":
            print(f"  Female body area PA: F = {FPA:.1f}, p = {pPA:.5f}; regime × assay F = {Fi:.1f}, p = {pi:.4f}")
RMs = pd.DataFrame(rows); RMs.to_csv(os.path.join(OUT, "regime_means_results.csv"), index=False)
t4 = ["| Trait | Sex | PA Constant | PA Predictable | PA Random | PA F₂,₆ | PA p | CG Constant | CG Predictable | CG Random | CG F₂,₆ | CG p | Regime × assay F₂,₆ | Regime × assay p |", "|" + "---|" * 14]
t5 = ["| Trait | Sex | Predictable − Constant | Random − Constant | Random − Predictable |", "|---|---|---|---|---|"]
for _, r in RMs.iterrows():
    t4.append(f"| {r.Trait} | {r.Sex} | {r.PA_C:.3f} | {r.PA_P:.3f} | {r.PA_R:.3f} | {r.F_PA:.1f} | {fp(r.p_PA)} | {r.CG_C:.3f} | {r.CG_P:.3f} | {r.CG_R:.3f} | {r.F_CG:.1f} | {fp(r.p_CG)} | {r.F_int:.1f} | {fp(r.p_int)} |")
    c = lambda k: f"{r[f'PA_{k}_pct']:+.1f}% [{r[f'PA_{k}_lo']:+.1f}, {r[f'PA_{k}_hi']:+.1f}]"
    t5.append(f"| {r.Trait} | {r.Sex} | {c('P−C')} | {c('R−C')} | {c('R−P')} |")

# per-generation check of the plasticity-assay regime effect, and sexual size dimorphism by regime
with open(os.path.join(OUT, "regime_means_per_generation.txt"), "w") as f:
    for tr, col in TRAITS.items():
        for g in GENS:
            m = fly[(fly.Assay == "PA") & (fly.Generation == g) & fly[col].notna()].groupby("Pop")[col].mean(); m.index = m.index.astype(str)
            F, p = anova_by_regime(m); f.write(f"{tr} PA gen {int(g)}: F(2,6) = {F:.2f}, p = {p:.3f}; " +
                                              ", ".join(f"{r} = {np.exp(m[[q for q in m.index if REGIME[q]==r]].mean())*SCALE[tr]:.3f}" for r in REG) + "\n")
    ssd = popmeans(fly[fly.Sex == "F"], "log_body", "PA") - popmeans(fly[fly.Sex == "M"], "log_body", "PA")
    F, p = anova_by_regime(ssd); f.write(f"Sexual size dimorphism (F − M log body area, PA) by regime: F(2,6) = {F:.2f}, p = {p:.3f}\n")

# baseline means (Table S1)
tS1 = ["| Test temperature (°C) | Wing length (mm) | Femur length (mm) | Body area (mm²) |", "|---|---|---|---|"]
for temp, g in fly.groupby("TestTemp"):
    tS1.append(f"| {int(temp)} | {g.Wing_Length.mean()*SCALE['Wing']:.2f} | {g.Femur_Length.mean()*SCALE['Femur']:.2f} | {g.Body_Area.mean()*SCALE['Body']:.2f} |")

with open(os.path.join(OUT, "tables.md"), "w") as f:
    for name, t in [("Table 1", t1), ("Table 2", t2), ("Table 3", t3), ("Table 4", t4), ("Table 5", t5), ("Table S1", tS1)]:
        f.write(f"## {name}\n\n" + "\n".join(t) + "\n\n")


# ----------------------------------------------------------------------------
# 3b. Supplementary tables (S2-S7)
# ----------------------------------------------------------------------------
RATIOS = {"Wing:body": "log_wb", "Wing:femur": "log_wf"}
tS2 = ["| Ratio | Assay | Sex | n pops | Slope gen 4 | Slope gen 10 | Negative slopes (gen 4, gen 10) | Rate of change | 95% CI | p | Pops flattening | Regime F₂,₆ | p |", "|" + "---|" * 13]
for name, col in RATIOS.items():
    for asy in ["PA", "CG"]:
        for sex, lab in [(None, "Both"), ("F", "F"), ("M", "M")]:
            sub = fly if sex is None else fly[fly.Sex == sex]; sl4, sl10, rt = {}, {}, {}
            for p in REGIME:
                v = [slope(sub, col, asy, p, g)[0] for g in GENS]
                if not any(np.isnan(v)): sl4[p], sl10[p], rt[p] = v[0], v[2], stats.linregress([4, 7, 10], v).slope
            n = len(rt)
            if n < 4: tS2.append(f"| {name} | {asy} | {lab} | {n} | — | — | — | — | — | — | — | — | — |"); continue
            r = pd.Series(rt); t, pv = stats.ttest_1samp(r, 0); F, pF = anova_by_regime(r)
            Fs, ps = ("—", "—") if np.isnan(F) else (f"{F:.2f}", fp(pF))
            tS2.append(f"| {name} | {asy} | {lab} | {n} | {np.mean(list(sl4.values())):+.4f} | {np.mean(list(sl10.values())):+.4f} | {sum(v<0 for v in sl4.values())}/{n}, {sum(v<0 for v in sl10.values())}/{n} | {r.mean():+.5f} | {r.mean()-ci95(r):+.5f}, {r.mean()+ci95(r):+.5f} | {fp(pv)} | {(r>0).sum()}/{n} | {Fs} | {ps} |")
            if lab == "Both": print(f"  {name:11s} {asy}: gen-4 slope {np.mean(list(sl4.values())):+.4f} ({sum(v<0 for v in sl4.values())}/{n} negative); rate p = {pv:.3f}")

tS3 = ["| Trait | Assay | Pop | Regime | Slope gen 4 (SE; n) | Slope gen 7 (SE; n) | Slope gen 10 (SE; n) | Rate of change |", "|---|---|---|---|---|---|---|---|"]
for (tr, asy, pop), g in S.groupby(["Trait", "Assay", "Pop"], sort=False):
    cells = {int(r.Gen): f"{r.Slope:+.4f} ({r.SE:.4f}; {int(r.n)})" for r in g.itertuples()}
    rate = RP[(RP.Trait == tr) & (RP.Assay == asy) & (RP.Pop == pop)].Rate
    tS3.append(f"| {tr} | {asy} | {pop} | {REGIME[pop]} | {cells.get(4,'—')} | {cells.get(7,'—')} | {cells.get(10,'—')} | {rate.iloc[0]:+.5f} |" if len(rate) else f"| {tr} | {asy} | {pop} | {REGIME[pop]} | {cells.get(4,'—')} | {cells.get(7,'—')} | {cells.get(10,'—')} | — |")

tS4 = ["| Trait | Assay | Pop | Regime | Female slopes gen 4 / 7 / 10 | Female rate | Male slopes gen 4 / 7 / 10 | Male rate |", "|---|---|---|---|---|---|---|---|"]
for tr in TRAITS:
    for asy in ["PA", "CG"]:
        for pop in REGIME:
            f = RF[(RF.Trait == tr) & (RF.Assay == asy) & (RF.Pop == pop)]; m = RM[(RM.Trait == tr) & (RM.Assay == asy) & (RM.Pop == pop)]
            fs = f"{f.S4.iloc[0]:+.4f} / {f.S7.iloc[0]:+.4f} / {f.S10.iloc[0]:+.4f}" if len(f) else "—"; fr = f"{f.Rate.iloc[0]:+.5f}" if len(f) else "—"
            ms = f"{m.S4.iloc[0]:+.4f} / {m.S7.iloc[0]:+.4f} / {m.S10.iloc[0]:+.4f}" if len(m) else "not estimable"; mr = f"{m.Rate.iloc[0]:+.5f}" if len(m) else "—"
            tS4.append(f"| {tr} | {asy} | {pop} | {REGIME[pop]} | {fs} | {fr} | {ms} | {mr} |")

tS5 = ["| Trait | Pop | Regime | PA both | PA female | PA male | CG both | CG female | CG male |", "|---|---|---|---|---|---|---|---|---|"]
for tr, col in TRAITS.items():
    for pop in REGIME:
        vals = []
        for asy in ["PA", "CG"]:
            for sex in [None, "F", "M"]:
                sub = fly if sex is None else fly[fly.Sex == sex]; m = popmeans(sub, col, asy)
                vals.append(f"{np.exp(m[pop])*SCALE[tr]:.3f}" if pop in m.index else "—")
        tS5.append(f"| {tr} ({UNIT[tr]}) | {pop} | {REGIME[pop]} | " + " | ".join(vals) + " |")

tS6 = ["| Trait | Generation | Constant | Predictable | Random | F₂,₆ | p |", "|---|---|---|---|---|---|---|"]
for tr, col in TRAITS.items():
    for g in GENS:
        m = fly[(fly.Assay == "PA") & (fly.Generation == g) & fly[col].notna()].groupby("Pop")[col].mean(); m.index = m.index.astype(str)
        F, p = anova_by_regime(m); nat = {r: np.exp(m[[q for q in m.index if REGIME[q] == r]].mean()) * SCALE[tr] for r in REG}
        tS6.append(f"| {tr} ({UNIT[tr]}) | {int(g)} | {nat['Constant']:.3f} | {nat['Predictable']:.3f} | {nat['Random']:.3f} | {F:.1f} | {fp(p)} |")

tS7 = ["| Pop | Regime | Assay | Generation | Flies (F / M) | With wing | With femur | With body area | Flies per test temperature (min–max) |", "|---|---|---|---|---|---|---|---|---|"]
for (pop, asy, g), d in fly.groupby(["Pop", "Assay", "Generation"]):
    per = d.groupby("TestTemp").size()
    tS7.append(f"| {pop} | {REGIME[pop]} | {asy} | {int(g)} | {len(d)} ({(d.Sex=='F').sum()} / {(d.Sex=='M').sum()}) | {d.Wing_Length.notna().sum()} | {d.Femur_Length.notna().sum()} | {d.Body_Area.notna().sum()} | {per.min()}–{per.max()} |")

with open(os.path.join(OUT, "tables_supplement.md"), "w") as f:
    for name, t in [("Table S2", tS2), ("Table S3", tS3), ("Table S4", tS4), ("Table S5", tS5), ("Table S6", tS6), ("Table S7", tS7)]:
        f.write(f"## {name}\n\n" + "\n".join(t) + "\n\n")
# ----------------------------------------------------------------------------
# 4. Figures
# ----------------------------------------------------------------------------
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
def save(fig, name): fig.savefig(f"{FIG}/{name}.pdf"); fig.savefig(f"{FIG}/{name}.png", dpi=300); plt.close(fig)
def ptxt(p): return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"
def column_headers(axes, top=1.22):
    for j, asy in enumerate(["PA", "CG"]):
        axes[0, j].text(0.5, top, ASL[asy], transform=axes[0, j].transAxes, ha="center", va="bottom", fontsize=10, fontweight="bold")
def letter(ax, L): ax.text(-0.16, 1.04, L, transform=ax.transAxes, fontsize=11, fontweight="bold")
regime_lines = [Line2D([0], [0], color=COL[r], marker="o", ms=5, lw=1.4, label=r) for r in COL]
regime_marks = [Line2D([0], [0], color=COL[r], marker="o", ms=6, lw=0, label=r) for r in COL]
rng = np.random.default_rng(1)

# Figure 2: temperature-size rule
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))
for ax, tr, raw in zip(axes, TRAITS, ["Wing_Length", "Femur_Length", "Body_Area"]):
    g = fly.groupby("TestTemp")[raw].agg(["mean", "sem"]) * SCALE[tr]
    ax.errorbar(g.index, g["mean"], yerr=g["sem"], fmt="o-", color="k", ms=4, capsize=3, lw=1.2)
    ax.set_xlabel("Test temperature (°C)"); ax.set_ylabel(f"{TRL[tr]} ({UNIT[tr]})"); ax.set_xticks([23, 25, 27, 29, 31])
    ax.text(-0.28, 1.05, "ABC"[list(TRAITS).index(tr)], transform=ax.transAxes, fontsize=11, fontweight="bold")
fig.tight_layout(); save(fig, "Fig2_TSR")

# Figure 3: per-population reaction-norm slopes across generations
fig, axes = plt.subplots(3, 2, figsize=(6.8, 7.8), sharex=True)
for i, tr in enumerate(TRAITS):
    ymin, ymax = np.inf, -np.inf
    for j, asy in enumerate(["PA", "CG"]):
        ax = axes[i, j]; d = S[(S.Trait == tr) & (S.Assay == asy)]
        for pop, g in d.groupby("Pop"):
            g = g.sort_values("Gen"); ax.errorbar(g.Gen, g.Slope, yerr=g.SE, color=COL[g.Regime.iloc[0]], lw=1, alpha=0.65, marker="o", ms=3.2, elinewidth=0.7, zorder=2)
        m = d.groupby("Gen").Slope.agg(["mean", ci95])
        ax.fill_between(m.index, m["mean"] - m["ci95"], m["mean"] + m["ci95"], color="k", alpha=0.10, zorder=3, lw=0)
        ax.plot(m.index, m["mean"], color="k", lw=2.2, marker="o", ms=5, zorder=4); ax.axhline(0, color="k", lw=0.7, ls="--", zorder=1)
        ax.set_xticks([4, 7, 10]); ax.set_title(TRL[tr], fontsize=9.5, loc="left"); letter(ax, "ABCDEF"[i * 2 + j])
        lo, hi = ax.get_ylim(); ymin, ymax = min(ymin, lo), max(ymax, hi)
    for j in range(2): axes[i, j].set_ylim(ymin, ymax)
    axes[i, 0].set_ylabel("Reaction-norm slope\n(log-units °C$^{-1}$)")
for j in range(2): axes[2, j].set_xlabel("Generation")
column_headers(axes)
fig.legend(handles=regime_lines + [Line2D([0], [0], color="k", lw=2.2, marker="o", ms=5, label="Mean of 9 populations ± 95% CI")], loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=[0, 0.04, 1, 0.98]); save(fig, "Fig3_population_trajectories")

# Figure 4: population-level rates by trait, assay and sex
groups = [("PA", RP), ("PA", RF), ("PA", RM), ("CG", RP), ("CG", RF), ("CG", RM)]
xpos = [0, 1, 2, 3.6, 4.6, 5.6]; xlab = ["Both", "♀", "♂", "Both", "♀", "♂"]
fig, axes = plt.subplots(3, 1, figsize=(6.2, 7.4), sharex=True)
for i, tr in enumerate(TRAITS):
    ax = axes[i]
    for x, (asy, R) in zip(xpos, groups):
        w = R[(R.Trait == tr) & (R.Assay == asy)]
        ax.scatter(x + rng.uniform(-0.18, 0.18, len(w)), w.Rate, c=[COL[r] for r in w.Regime], s=22, alpha=0.85, zorder=3, edgecolor="none")
        if len(w) >= 4: ax.errorbar(x, w.Rate.mean(), yerr=ci95(w.Rate), fmt="D", color="k", ms=5, capsize=4, elinewidth=1.4, zorder=4)
        else: ax.text(x, 0.96, f"n = {len(w)}", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7.5, color="0.35")
    ax.axhline(0, color="k", lw=0.7, ls="--", zorder=1); ax.axvline(2.8, color="0.8", lw=0.8)
    ax.set_ylabel("Rate of slope change\n(log-units °C$^{-1}$ gen$^{-1}$)"); ax.set_title(TRL[tr], fontsize=9.5, loc="left")
    ax.text(-0.13, 1.03, "ABC"[i], transform=ax.transAxes, fontsize=11, fontweight="bold")
axes[2].set_xticks(xpos); axes[2].set_xticklabels(xlab)
for lab, xc in [("Plasticity assay", 1.0), ("Common garden", 4.6)]: axes[2].text(xc, -0.22, lab, transform=axes[2].get_xaxis_transform(), ha="center", fontsize=9)
fig.legend(handles=regime_marks + [Line2D([0], [0], color="k", marker="D", ms=5, lw=0, label="Mean ± 95% CI")], loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.005))
fig.tight_layout(rect=[0, 0.05, 1, 1]); save(fig, "Fig4_population_rates")

# Figure 5: mean morphology by regime, plasticity assay vs common garden
panels = [("Wing length (mm)", fly, "log_wing", "Wing"), ("Femur length (mm)", fly, "log_femur", "Femur"),
          ("Body area (mm²)", fly, "log_body", "Body"), ("Female body area (mm²)", fly[fly.Sex == "F"], "log_body", "Body")]
rng5 = np.random.default_rng(3)
fig, axes = plt.subplots(4, 2, figsize=(6.9, 9.4), sharex=True)
for i, (lab, sub, col, tr) in enumerate(panels):
    ymin, ymax = np.inf, -np.inf
    for j, asy in enumerate(["PA", "CG"]):
        ax = axes[i, j]; m = np.exp(popmeans(sub, col, asy)) * SCALE[tr]
        for k, r in enumerate(REG):
            v = m[[p for p in m.index if REGIME[p] == r]].values
            ax.scatter([k] * 3 + rng5.uniform(-0.12, 0.12, 3), v, color=COL[r], s=28, zorder=3, edgecolor="k", linewidth=0.4)
            ax.errorbar(k + 0.32, v.mean(), yerr=ci95(v), fmt="D", color=COL[r], mec="k", ms=6, capsize=4, elinewidth=1.3, zorder=4)
        F, p = anova_by_regime(popmeans(sub, col, asy))
        ax.text(1.0, 1.03, f"F₂,₆ = {F:.1f}, {ptxt(p)}", transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
        ax.set_xticks(range(3)); ax.set_xticklabels(REG); ax.set_xlim(-0.5, 2.8)
        ax.set_title(lab, fontsize=9.5, loc="left"); letter(ax, "ABCDEFGH"[i * 2 + j])
        lo, hi = ax.get_ylim(); ymin, ymax = min(ymin, lo), max(ymax, hi)
    for j in range(2): axes[i, j].set_ylim(ymin, ymax)
    axes[i, 0].set_ylabel("Population mean")
for j in range(2): axes[3, j].set_xlabel("Evolutionary regime")
column_headers(axes)
fig.legend(handles=[Line2D([0], [0], color="k", marker="o", lw=0, ms=6, mfc="w", label="Population mean (n = 9)"), Line2D([0], [0], color="k", marker="D", lw=0, ms=6, label="Regime mean ± 95% CI (n = 3)")], loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.002))
fig.tight_layout(rect=[0, 0.03, 1, 0.985]); save(fig, "Fig5_regime_means_PA_vs_CG")

# Figure 6: persistence of the decline in the common garden
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7))
for i, tr in enumerate(TRAITS):
    ax = axes[i]; pa = RP[(RP.Trait == tr) & (RP.Assay == "PA")].set_index("Pop"); cg = RP[(RP.Trait == tr) & (RP.Assay == "CG")].set_index("Pop")
    x = pa.Rate; y = cg.loc[x.index].Rate; lim = [min(x.min(), y.min(), 0) - 0.0005, max(x.max(), y.max(), 0) + 0.0005]
    ax.plot(lim, lim, color="k", ls="--", lw=0.8, zorder=1); ax.axhline(0, color="0.7", lw=0.6); ax.axvline(0, color="0.7", lw=0.6)
    ax.scatter(x, y, c=[COL[r] for r in pa.Regime], s=30, zorder=3, edgecolor="k", linewidth=0.4)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal"); ax.set_title(TRL[tr], fontsize=9.5, loc="left"); ax.set_xlabel("Rate in plasticity assay")
    ax.text(-0.22, 1.05, "ABC"[i], transform=ax.transAxes, fontsize=11, fontweight="bold")
axes[0].set_ylabel("Rate in common garden\n(log-units °C$^{-1}$ gen$^{-1}$)")
fig.legend(handles=regime_marks + [Line2D([0], [0], color="k", ls="--", lw=0.8, label="1:1 (identical rates)")], loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.04))
fig.tight_layout(rect=[0, 0.08, 1, 1]); save(fig, "Fig6_PA_vs_CG")

# Figure S1: rate of plasticity change by regime
fig, axes = plt.subplots(3, 2, figsize=(6.4, 7.4), sharex=True)
for i, tr in enumerate(TRAITS):
    ymin, ymax = np.inf, -np.inf
    for j, asy in enumerate(["PA", "CG"]):
        ax = axes[i, j]; w = RP[(RP.Trait == tr) & (RP.Assay == asy)]; o = w.Rate.mean(); oc = ci95(w.Rate)
        ax.axhspan(o - oc, o + oc, color="k", alpha=0.07, zorder=0); ax.axhline(o, color="k", lw=1.2, zorder=1); ax.axhline(0, color="k", lw=0.7, ls="--", zorder=1)
        for k, r in enumerate(REG):
            v = w[w.Regime == r].Rate.values
            ax.scatter([k] * len(v) + rng.uniform(-0.12, 0.12, len(v)), v, color=COL[r], s=26, zorder=3, edgecolor="k", linewidth=0.4)
            ax.errorbar(k + 0.32, v.mean(), yerr=ci95(v), fmt="D", color=COL[r], mec="k", ms=6, capsize=4, elinewidth=1.4, zorder=4)
        ax.set_xticks(range(3)); ax.set_xticklabels(REG); ax.set_xlim(-0.5, 2.8); ax.set_title(TRL[tr], fontsize=9.5, loc="left"); letter(ax, "ABCDEF"[i * 2 + j])
        lo, hi = ax.get_ylim(); ymin, ymax = min(ymin, lo), max(ymax, hi)
    for j in range(2): axes[i, j].set_ylim(ymin, ymax)
    axes[i, 0].set_ylabel("Rate of slope change\n(log-units °C$^{-1}$ gen$^{-1}$)")
for j in range(2): axes[2, j].set_xlabel("Evolutionary regime")
column_headers(axes)
fig.legend(handles=[Line2D([0], [0], color="k", marker="o", lw=0, ms=6, mfc="w", label="Population (n = 9)"), Line2D([0], [0], color="k", marker="D", mec="k", ms=6, lw=0, label="Regime mean ± 95% CI (n = 3)"), Line2D([0], [0], color="k", lw=1.2, label="Mean of all 9 populations ± 95% CI")], loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.005))
fig.tight_layout(rect=[0, 0.04, 1, 0.98]); save(fig, "FigS1_regime_rate_contrasts")

print(f"Done. Outputs in ./{OUT}/ and ./{FIG}/")
