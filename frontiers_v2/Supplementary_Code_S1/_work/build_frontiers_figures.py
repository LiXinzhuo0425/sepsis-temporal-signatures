#!/usr/bin/env python3
"""Build six Frontiers-facing manuscript figures with Python only."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "_work" / "analysis"
REPO = ROOT / "_sources" / "repo_subset_v1.3.1"
INPUT = ROOT / "_work" / "input_snapshot" / "04_Data_S1.xlsx"
FIGDIR = ROOT / "final" / "Figures"
SRCDIR = ROOT / "final" / "Figure_Source_Data"
PREVIEW = ROOT / "_qa" / "figures"
for directory in [FIGDIR, SRCDIR, PREVIEW]:
    directory.mkdir(parents=True, exist_ok=True)

MM = 1 / 25.4
BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#7A5195"
SKY = "#56B4E9"
YELLOW = "#E69F00"
GREY = "#6E7781"
LIGHT = "#E9EEF3"
DARK = "#1F2933"
RED = "#B2182B"
COLORS = {"T24": BLUE, "T48": ORANGE}

SIG_ORDER = ["SIG001", "SIG002", "SIG003", "SIG004", "SIG022", "SIG023", "SIG033", "SIG034"]
SIG_NAMES = {
    "SIG001": "Sepsis MetaScore", "SIG002": "SeptiCyte LAB", "SIG003": "FAIM3:PLAC8", "SIG004": "sNIP",
    "SIG022": "Bacterial/Viral MetaScore", "SIG023": "Herberg two-gene DRS",
    "SIG033": "Lin mortality score", "SIG034": "Severe-or-Mild score",
}
SIG_SHORT = {
    "SIG001": "MetaScore", "SIG002": "SeptiCyte", "SIG003": "FAIM3:PLAC8", "SIG004": "sNIP",
    "SIG022": "B/V MetaScore", "SIG023": "Herberg DRS", "SIG033": "Lin score", "SIG034": "Severe/Mild",
}
SIG_COLORS = {
    "SIG001": "#0072B2", "SIG002": "#009E73", "SIG003": "#CC79A7", "SIG004": "#D55E00",
    "SIG022": "#56B4E9", "SIG023": "#E69F00", "SIG033": "#7A5195", "SIG034": "#6E7781",
}


def configure() -> None:
    sns.set_theme(style="white", context="paper")
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.8, "lines.linewidth": 2.0, "lines.markersize": 5,
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    })


def panel_label(ax, label: str) -> None:
    ax.text(-0.08, 1.04, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")


def clean(ax, grid: str | None = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(True, axis=grid, color="#D8DEE4", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)


def save(fig: plt.Figure, number: int) -> None:
    stem = f"Figure_{number}"
    fig.savefig(FIGDIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / f"{stem}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(PREVIEW / f"{stem}_preview.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def meta_data() -> pd.DataFrame:
    meta = pd.read_csv(REPO / "data/source_tables/stage3/03_05_signature_level_meta_analysis.csv")
    meta = meta[(meta["analysis_set"] == "PRIMARY_INDEPENDENT") & meta["time_window"].isin(["T1", "T2"])].copy()
    meta["landmark"] = meta["time_window"].map({"T1": "T24", "T2": "T48"})
    meta["signature_name"] = meta["signature_id"].map(SIG_NAMES)
    return meta


def figure1() -> None:
    cohorts = pd.read_excel(INPUT, sheet_name="S1_Cohorts")
    signatures = pd.read_excel(INPUT, sheet_name="S2_Signatures")
    cohorts.to_csv(SRCDIR / "Figure_1A_cohorts.csv", index=False)
    signatures.to_csv(SRCDIR / "Figure_1C_signatures.csv", index=False)

    fig = plt.figure(figsize=(180 * MM, 145 * MM), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.35], hspace=0.18, wspace=0.12)
    ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[0, 1]); ax_c = fig.add_subplot(gs[1, :])
    for ax in [ax_a, ax_b, ax_c]: ax.set_axis_off()

    panel_label(ax_a, "A")
    ax_a.set_title("Longitudinal cohort architecture", loc="left", fontweight="bold", pad=8)
    boxes = [
        (0.01, 0.67, 0.29, 0.18, "6 primary\ncohorts", BLUE),
        (0.355, 0.67, 0.29, 0.18, "302 baseline\npatients", BLUE),
        (0.70, 0.67, 0.29, 0.18, "264 with an early\npaired sample", BLUE),
        (0.355, 0.27, 0.29, 0.18, "T24\n172 patients", SKY),
        (0.70, 0.27, 0.29, 0.18, "T48\n146 patients", ORANGE),
        (0.01, 0.27, 0.29, 0.18, "Independent T24\nreplication · 47", GREEN),
    ]
    for x, y, w, h, text, color in boxes:
        ax_a.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02", facecolor=color, edgecolor="none", alpha=0.13))
        ax_a.text(x + w/2, y + h/2, text, ha="center", va="center", color=DARK, fontweight="bold", fontsize=6.8)
    for start, end in [((0.30, .76), (.355, .76)), ((.645, .76), (.70, .76)), ((.845, .67), (.845, .45)), ((.50, .67), (.50, .45))]:
        ax_a.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9, linewidth=1.4, color=GREY))
    ax_a.text(.03, .08, "Clinical context: GSE95233, GSE110487, GSE57065 and GSE54514", fontsize=7.2, color=GREY)

    panel_label(ax_b, "B")
    ax_b.set_title("Evidence layers", loc="left", fontweight="bold", pad=8)
    layers = [
        ("Temporal drift", "mean within-patient ΔZ"),
        ("Rank persistence", "baseline-to-follow-up ordering"),
        ("Clinical context", "source-defined course or severity"),
        ("Independent replication", "direction and transcriptome-matched null"),
        ("Formula architecture", "exact gene-level decomposition"),
    ]
    for i, (head, sub) in enumerate(layers):
        y = 0.82 - i * 0.17
        color = [BLUE, SKY, ORANGE, GREEN, PURPLE][i]
        ax_b.add_patch(FancyBboxPatch((.05, y-.09), .90, .14, boxstyle="round,pad=.008,rounding_size=.016", facecolor="white", edgecolor=color, linewidth=1.6))
        ax_b.text(.09, y+.006, head, fontweight="bold", color=color, va="center", fontsize=7.6)
        ax_b.text(.09, y-.045, sub, color=GREY, va="center", fontsize=6.7)

    panel_label(ax_c, "C")
    ax_c.set_title("Eight fixed published formulas", loc="left", fontweight="bold", pad=8)
    categories = [
        ("Sepsis / infection", ["SIG001", "SIG002", "SIG003", "SIG004"]),
        ("Pathogen class", ["SIG022", "SIG023"]),
        ("Outcome / severity", ["SIG033", "SIG034"]),
    ]
    formula = {
        "SIG001": "difference of geometric means", "SIG002": "four-gene linear contrast",
        "SIG003": "two-gene ratio", "SIG004": "three-gene normalized contrast",
        "SIG022": "difference of geometric means", "SIG023": "two-gene linear contrast",
        "SIG033": "seven-gene weighted sum", "SIG034": "four-module ratio",
    }
    x0 = [0.02, 0.52, 0.77]; widths = [0.47, 0.23, 0.21]
    for (cat, ids), x, width in zip(categories, x0, widths):
        ax_c.text(x, .90, cat, fontweight="bold", color=DARK)
        for j, sid in enumerate(ids):
            y = .76 - j * .19
            ax_c.add_patch(FancyBboxPatch((x, y-.11), width, .14, boxstyle="round,pad=.008,rounding_size=.012", facecolor=SIG_COLORS[sid], edgecolor="none", alpha=.12))
            ax_c.text(x+.015, y-.015, SIG_SHORT[sid], fontweight="bold", color=SIG_COLORS[sid], va="center")
            ax_c.text(x+.015, y-.073, formula[sid], color=GREY, fontsize=6.2, va="center")
    ax_c.text(.02, .02, "Published gene membership and algebra were retained; no formula was refitted in the longitudinal cohorts.", color=GREY, fontsize=7.2)
    save(fig, 1)


def figure2() -> None:
    meta = meta_data()
    cohort = pd.read_csv(REPO / "data/source_tables/stage3/03_04_cohort_signature_primary_effects.csv")
    cohort = cohort[cohort["time_window"].isin(["T1", "T2"])].copy()
    cohort["landmark"] = cohort["time_window"].map({"T1": "T24", "T2": "T48"})
    meta.to_csv(SRCDIR / "Figure_2AB_pooled_effects.csv", index=False)
    cohort.to_csv(SRCDIR / "Figure_2C_cohort_effects.csv", index=False)
    fig = plt.figure(figsize=(180 * MM, 165 * MM), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, .95])
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    y = np.arange(len(SIG_ORDER))[::-1]
    for ax, landmark, lab in zip(axes, ["T24", "T48"], ["A", "B"]):
        data = meta.set_index("signature_id").loc[SIG_ORDER]
        data = data[data["landmark"] == landmark] if "landmark" in data else data
        # Reindex after landmark filtering.
        data = meta[meta["landmark"] == landmark].set_index("signature_id").loc[SIG_ORDER]
        for pos, sid in zip(y, SIG_ORDER):
            row = data.loc[sid]
            ax.plot([row.prediction_lower, row.prediction_upper], [pos, pos], color="#C7CDD4", linewidth=2, zorder=1)
            ax.plot([row.ci95_lower, row.ci95_upper], [pos, pos], color=SIG_COLORS[sid], linewidth=3, zorder=2)
            ax.scatter(row.pooled_delta_z, pos, s=30, color=SIG_COLORS[sid], edgecolor="white", linewidth=.6, zorder=3)
        ax.axvline(0, color=DARK, linewidth=1)
        ax.set_yticks(y, [SIG_SHORT[s] for s in SIG_ORDER] if landmark == "T24" else [])
        ax.set_xlabel("Mean within-patient change, ΔZ")
        ax.set_title(landmark, fontweight="bold")
        ax.set_xlim(-2.35, 1.75)
        clean(ax)
        panel_label(ax, lab)
    axes[0].legend(handles=[Line2D([0],[0], color=GREY, lw=3, label="95% CI"), Line2D([0],[0], color="#C7CDD4", lw=2, label="95% prediction interval")], loc="lower left", frameon=False)

    ax = fig.add_subplot(gs[1, :]); panel_label(ax, "C")
    independent = {(r.signature_id, r.time_window): set(str(r.cohorts).split(";")) for r in meta.itertuples()}
    columns = []
    for landmark, window, datasets in [("T24", "T1", ["GSE236713", "GSE54514", "GSE57065", "GSE95233"]), ("T48", "T2", ["GSE54514", "GSE57065", "GSE95233", "GSE110487", "GSE8121"])]:
        for dataset in datasets:
            columns.append((landmark, window, dataset))
    matrix = np.full((len(SIG_ORDER), len(columns)), np.nan)
    for i, sid in enumerate(SIG_ORDER):
        for j, (landmark, window, dataset) in enumerate(columns):
            if dataset not in independent[(sid, window)]:
                continue
            row = cohort[(cohort.signature_id == sid) & (cohort.time_window == window) & (cohort.dataset == dataset)]
            if len(row): matrix[i, j] = row.iloc[0].mean_delta_z
    sns.heatmap(matrix, ax=ax, cmap=sns.diverging_palette(240, 10, as_cmap=True), center=0, vmin=-1.6, vmax=1.6,
                xticklabels=[f"{d}\n{l}" for l, _, d in columns], yticklabels=[SIG_SHORT[s] for s in SIG_ORDER],
                linewidths=.5, linecolor="white", cbar_kws={"label": "Cohort mean ΔZ", "shrink": .75})
    ax.set_title("Independent cohort effects", fontweight="bold", pad=5)
    ax.tick_params(axis="x", rotation=45)
    save(fig, 2)


def figure3() -> None:
    rank = pd.read_csv(AN / "02_rank_persistence_meta.csv")
    cohort = pd.read_csv(AN / "02_rank_persistence_cohort.csv")
    meta = meta_data()
    rank.to_csv(SRCDIR / "Figure_3AC_rank_meta.csv", index=False)
    cohort.to_csv(SRCDIR / "Figure_3B_rank_cohort.csv", index=False)
    fig = plt.figure(figsize=(180 * MM, 170 * MM), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05])
    ax = fig.add_subplot(gs[0, 0]); panel_label(ax, "A")
    base_y = np.arange(len(SIG_ORDER))[::-1]
    for offset, landmark in [(.13, "T24"), (-.13, "T48")]:
        data = rank[rank.time_window == landmark].set_index("signature_id").loc[SIG_ORDER]
        ypos = base_y + offset
        ax.errorbar(data.pooled_spearman_rho, ypos,
                    xerr=[data.pooled_spearman_rho-data.modified_hk_ci95_lower, data.modified_hk_ci95_upper-data.pooled_spearman_rho],
                    fmt="o", color=COLORS[landmark], ecolor=COLORS[landmark], elinewidth=2, capsize=2, label=landmark)
    ax.axvline(0, color=DARK, lw=1)
    ax.set_yticks(base_y, [SIG_SHORT[x] for x in SIG_ORDER])
    ax.set_xlim(-.3, 1.02); ax.set_xlabel("Pooled baseline-to-follow-up Spearman ρ")
    ax.legend(frameon=False, loc="lower left"); clean(ax)
    ax.set_title("Patient-rank persistence", fontweight="bold")

    ax = fig.add_subplot(gs[0, 1]); panel_label(ax, "B")
    columns = []
    for landmark, datasets in [("T24", ["GSE236713", "GSE54514", "GSE57065", "GSE95233"]), ("T48", ["GSE54514", "GSE57065", "GSE95233", "GSE110487", "GSE8121"])]:
        for dataset in datasets:
            columns.append((landmark, dataset))
    matrix = np.full((len(SIG_ORDER), len(columns)), np.nan)
    for i, sid in enumerate(SIG_ORDER):
        for j, (landmark, dataset) in enumerate(columns):
            row = cohort[(cohort.signature_id == sid) & (cohort.time_window == landmark) & (cohort.dataset == dataset) & (cohort.primary_independent_eligible == True)]
            if len(row): matrix[i, j] = row.iloc[0].spearman_rho
    sns.heatmap(matrix, ax=ax, cmap="YlGnBu", vmin=0, vmax=1, xticklabels=[f"{d}\n{l}" for l,d in columns],
                yticklabels=[SIG_SHORT[s] for s in SIG_ORDER], linewidths=.5, linecolor="white", cbar_kws={"label":"Spearman ρ", "shrink":.75})
    ax.tick_params(axis="x", rotation=55, labelsize=6); ax.set_title("Cohort-level rank correlations", fontweight="bold")

    ax = fig.add_subplot(gs[1, :]); panel_label(ax, "C")
    combined = rank.merge(meta[["signature_id","landmark","pooled_delta_z"]], left_on=["signature_id","time_window"], right_on=["signature_id","landmark"], validate="one_to_one")
    drift_rank_offsets = {
        "SIG001": (5, 3), "SIG002": (5, 3), "SIG003": (5, 3), "SIG004": (5, 3),
        "SIG022": (5, 5), "SIG023": (5, -11), "SIG033": (5, 3), "SIG034": (5, 11),
    }
    for sid in SIG_ORDER:
        group = combined[combined.signature_id == sid].set_index("time_window").loc[["T24","T48"]]
        ax.plot(group.pooled_delta_z.abs(), group.pooled_spearman_rho, color="#B8C0C8", lw=1.1, zorder=1)
        for landmark in ["T24","T48"]:
            row = group.loc[landmark]
            size = 25 + row.median_cohort_rank_displacement_pp * 3.2
            ax.scatter(abs(row.pooled_delta_z), row.pooled_spearman_rho, s=size, c=COLORS[landmark], alpha=.85, edgecolor="white", linewidth=.7, zorder=2)
        row = group.loc["T48"]
        ax.annotate(
            SIG_SHORT[sid],
            (abs(row.pooled_delta_z), row.pooled_spearman_rho),
            xytext=drift_rank_offsets[sid],
            textcoords="offset points",
            fontsize=6.4,
        )
    ax.set_xlabel("Absolute pooled temporal drift, |ΔZ|")
    ax.set_ylabel("Pooled rank persistence, Spearman ρ")
    ax.set_xlim(-.03, 1.02); ax.set_ylim(.32, .82); clean(ax, grid="both")
    ax.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=BLUE,label='T24'),Line2D([0],[0],marker='o',color='none',markerfacecolor=ORANGE,label='T48')],frameon=False, loc="lower right", title="Landmark")
    ax.text(.02, .34, "Point size: median cohort rank displacement", color=GREY, fontsize=7)
    combined.to_csv(SRCDIR / "Figure_3C_drift_rank_map.csv", index=False)
    save(fig, 3)


def clinical_panel(ax, data: pd.DataFrame, title: str, label: str, xlim=(-2.3, 2.6)) -> None:
    ax.text(-.08, 1.10, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")
    y = np.arange(4)[::-1]
    data = data.set_index("signature_id").loc[["SIG001","SIG002","SIG003","SIG004"]]
    for pos, sid in zip(y, ["SIG001","SIG002","SIG003","SIG004"]):
        row = data.loc[sid]
        ax.plot([row.primary_ci95_lower, row.primary_ci95_upper], [pos,pos], color=SIG_COLORS[sid], linewidth=2.3)
        ax.scatter(row.primary_beta_adverse_minus_favorable, pos, color=SIG_COLORS[sid], s=28, edgecolor="white", linewidth=.5, zorder=3)
        if np.isfinite(row.holm_adjusted_p_primary_four) and row.holm_adjusted_p_primary_four < .05:
            ax.text(xlim[1]-.08, pos, "*", ha="right", va="center", fontweight="bold", fontsize=9)
    ax.axvline(0, color=DARK, lw=1)
    ax.set_yticks(y, [SIG_SHORT[x] for x in ["SIG001","SIG002","SIG003","SIG004"]])
    ax.set_xlim(*xlim); ax.set_title(title, fontweight="bold", fontsize=7.5, pad=9)
    ax.set_xlabel("Adverse − favorable change, ΔZ")
    clean(ax)


def figure4() -> None:
    clinical = pd.read_csv(AN / "05_clinical_course_interactions.csv")
    days = pd.read_csv(AN / "05_gse95233_day_stratified_sensitivity.csv")
    clinical.to_csv(SRCDIR / "Figure_4_primary_clinical_contrasts.csv", index=False)
    days.to_csv(SRCDIR / "Figure_4_gse95233_day_sensitivity.csv", index=False)
    fig, axes = plt.subplots(2, 3, figsize=(180 * MM, 170 * MM), constrained_layout=True)
    clinical_panel(axes[0,0], clinical[(clinical.dataset=="GSE95233") & (clinical.contrast=="Admission to D2/D3")], "GSE95233, 28-day survival\n(day-adjusted)", "A")
    clinical_panel(axes[0,1], days[days.contrast=="Admission to D2"], "GSE95233, admission to D2", "B")
    clinical_panel(axes[0,2], days[days.contrast=="Admission to D3"], "GSE95233, admission to D3", "C")
    clinical_panel(axes[1,0], clinical[(clinical.dataset=="GSE110487") & (clinical.contrast=="T48")], "GSE110487, organ-function course", "D")
    # Two landmarks overlaid for severity and source-defined survival contexts.
    for ax, dataset, title, label in [(axes[1,1], "GSE57065", "GSE57065, baseline SAPS II category", "E"), (axes[1,2], "GSE54514", "GSE54514, source-defined survival", "F")]:
        ax.text(-.08, 1.10, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left"); y = np.arange(4)[::-1]
        for off, contrast in [(.13,"T24"),(-.13,"T48")]:
            d = clinical[(clinical.dataset==dataset)&(clinical.contrast==contrast)].set_index("signature_id").loc[["SIG001","SIG002","SIG003","SIG004"]]
            for pos, sid in zip(y+off, ["SIG001","SIG002","SIG003","SIG004"]):
                row=d.loc[sid]
                ax.plot([row.primary_ci95_lower,row.primary_ci95_upper],[pos,pos],color=COLORS[contrast],lw=2)
                ax.scatter(row.primary_beta_adverse_minus_favorable,pos,color=COLORS[contrast],s=22,edgecolor="white",linewidth=.4)
        ax.axvline(0,color=DARK,lw=1); ax.set_yticks(y,[SIG_SHORT[x] for x in ["SIG001","SIG002","SIG003","SIG004"]]); ax.set_xlim(-2.3,2.6)
        ax.set_xlabel("Adverse − favorable change, ΔZ"); ax.set_title(title,fontweight="bold",fontsize=7.3,pad=9); clean(ax)
        ax.legend(handles=[Line2D([0],[0],marker='o',color=BLUE,label='T24',lw=0),Line2D([0],[0],marker='o',color=ORANGE,label='T48',lw=0)],frameon=False,loc='lower right')
    save(fig, 4)


def figure5() -> None:
    effects = pd.read_csv(AN / "gse106878" / "GSE106878_directional_replication_effects.csv")
    effects = effects[effects.stratum == "all"].copy()
    pairs = pd.read_csv(AN / "gse106878" / "GSE106878_patient_paired_changes.csv")
    conc = pd.read_csv(AN / "03_cross_signature_concordance_meta.csv")
    ext_conc = pd.read_csv(AN / "07_gse106878_diagnostic_concordance.csv")
    null = pd.read_csv(AN / "08_formula_preserving_null_distribution.csv")
    bg = pd.read_csv(AN / "08_formula_preserving_background_summary.csv")
    for name, data in [("Figure_5A_replication.csv",effects),("Figure_5B_patient_changes.csv",pairs),("Figure_5C_primary_concordance.csv",conc),("Figure_5C_replication_concordance.csv",ext_conc),("Figure_5D_background_summary.csv",bg)]:
        data.to_csv(SRCDIR / name, index=False)
    fig = plt.figure(figsize=(180 * MM, 180 * MM), constrained_layout=True)
    gs = fig.add_gridspec(2,2)
    ax=fig.add_subplot(gs[0,0]); panel_label(ax,"A")
    data=effects.set_index("signature_id").loc[SIG_ORDER]; y=np.arange(8)[::-1]
    for pos,sid in zip(y,SIG_ORDER):
        r=data.loc[sid]; ax.plot([r.bootstrap_ci_low,r.bootstrap_ci_high],[pos,pos],color=SIG_COLORS[sid],lw=2.3); ax.scatter(r.mean_delta_z,pos,c=SIG_COLORS[sid],s=28,edgecolor="white",linewidth=.5)
    ax.axvline(0,color=DARK,lw=1); ax.set_yticks(y,[SIG_SHORT[s] for s in SIG_ORDER]); ax.set_xlabel("GSE106878 mean ΔZ"); ax.set_title("Independent T24 replication (n=47)",fontweight="bold"); clean(ax)

    ax=fig.add_subplot(gs[0,1]); panel_label(ax,"B")
    diag=pairs[pairs.signature_id.isin(["SIG001","SIG002","SIG003","SIG004"])].pivot(index="patient_id",columns="signature_id",values="delta_z")
    diag=diag.loc[diag.mean(axis=1).sort_values().index]
    sns.heatmap(diag,ax=ax,cmap=sns.diverging_palette(240,10,as_cmap=True),center=0,vmin=-3,vmax=3,yticklabels=False,xticklabels=[SIG_SHORT[s] for s in diag.columns],cbar_kws={"label":"Patient ΔZ","shrink":.7})
    ax.set_xlabel(""); ax.set_ylabel("47 paired patients"); ax.set_title("Individual diagnostic-score changes",fontweight="bold")
    ax.tick_params(axis="x", rotation=22, labelsize=6.5)

    ax=fig.add_subplot(gs[1,0]); panel_label(ax,"C")
    ids=["SIG001","SIG002","SIG003","SIG004"]; mat=np.eye(4)
    primary=conc[conc.time_window=="T48"]
    for r in primary.itertuples():
        i,j=ids.index(r.signature_a),ids.index(r.signature_b); mat[max(i,j),min(i,j)]=r.pooled_spearman_rho_delta_z
    for r in ext_conc.itertuples():
        i,j=ids.index(r.signature_a),ids.index(r.signature_b); mat[min(i,j),max(i,j)]=r.spearman_rho
    mask=np.eye(4,dtype=bool)
    sns.heatmap(mat,mask=mask,ax=ax,cmap=sns.diverging_palette(240,10,as_cmap=True),center=0,vmin=-.8,vmax=.8,annot=True,fmt=".2f",square=True,
                xticklabels=[SIG_SHORT[s] for s in ids],yticklabels=[SIG_SHORT[s] for s in ids],cbar_kws={"label":"Spearman ρ","shrink":.7})
    ax.tick_params(axis='x',rotation=35); ax.set_title("Change concordance\nupper: GSE106878; lower: pooled T48",fontweight="bold")

    ax=fig.add_subplot(gs[1,1]); panel_label(ax,"D")
    plot=null[null.status=="SUCCESS"].copy(); plot["name"]=plot.signature_id.map(SIG_SHORT)
    sns.violinplot(data=plot,x="pseudo_mean_delta_z",y="name",order=[SIG_SHORT[s] for s in SIG_ORDER],ax=ax,color=LIGHT,inner=None,cut=0,linewidth=.8)
    summ=bg.set_index("signature_id").loc[SIG_ORDER]
    ax.scatter(summ.observed_mean_delta_z,np.arange(8),c=[SIG_COLORS[s] for s in SIG_ORDER],s=28,edgecolor="white",linewidth=.5,zorder=4)
    ax.axvline(0,color=DARK,lw=1); ax.set_xlabel("Mean ΔZ under formula-preserving gene replacement"); ax.set_ylabel(""); ax.set_title("Baseline-matched temporal background",fontweight="bold"); clean(ax)
    ax.text(.98,.97,"Colored point: observed",transform=ax.transAxes,ha="right",va="top",fontsize=6.7,color=GREY)
    save(fig,5)


def figure6() -> None:
    contrib=pd.read_csv(REPO/"data/source_tables/stage4/04_07_gene_contribution_meta_analysis.csv")
    contrib=contrib[(contrib.analysis_set=="PRIMARY_INDEPENDENT")&(contrib.time_window=="T2")].copy()
    top=contrib.assign(abs_value=lambda d:d.pooled_contribution.abs()).sort_values(["signature_id","abs_value"],ascending=[True,False]).groupby("signature_id").head(2)
    arch=pd.read_csv(REPO/"data/source_tables/stage4/04_08_signature_drift_architecture.csv")
    recurrence=pd.read_csv(AN/"04_contribution_recurrence.csv"); recurrence=recurrence[recurrence.time_window=="T48"]
    rank=pd.read_csv(AN/"02_rank_persistence_meta.csv"); rank=rank[rank.time_window=="T48"]
    top.to_csv(SRCDIR/"Figure_6A_top_contributions.csv",index=False); arch.to_csv(SRCDIR/"Figure_6BD_architecture.csv",index=False); recurrence.to_csv(SRCDIR/"Figure_6C_recurrence.csv",index=False)
    fig=plt.figure(figsize=(180*MM,180*MM),constrained_layout=True); gs=fig.add_gridspec(2,2)
    ax=fig.add_subplot(gs[0,0]); panel_label(ax,"A")
    rows=[]
    for sid in SIG_ORDER:
        d=top[top.signature_id==sid].sort_values("abs_value",ascending=False)
        for j,r in enumerate(d.itertuples()): rows.append((sid,r.gene,r.pooled_contribution,r.ci95_lower,r.ci95_upper,j))
    ypos=np.arange(len(rows))[::-1]
    for pos,(sid,gene,val,lo,hi,j) in zip(ypos,rows):
        ax.plot([lo,hi],[pos,pos],color=SIG_COLORS[sid],lw=2); ax.scatter(val,pos,c=SIG_COLORS[sid],s=23,edgecolor="white",linewidth=.4)
    ax.axvline(0,color=DARK,lw=1); ax.set_yticks(ypos,[f"{SIG_SHORT[s]} · {g}" for s,g,*_ in rows]); ax.set_xlabel("Pooled gene contribution to ΔZ"); ax.set_title("Two largest T48 contributions",fontweight="bold"); clean(ax)

    ax=fig.add_subplot(gs[0,1]); panel_label(ax,"B")
    architecture_offsets={"SIG001":(3,3),"SIG002":(3,3),"SIG003":(-54,10),"SIG004":(3,3),"SIG022":(3,3),"SIG023":(-52,-12),"SIG033":(3,3),"SIG034":(3,3)}
    for r in arch.itertuples():
        ax.scatter(r.median_patient_dominance_ratio,r.median_patient_cancellation_index,s=35+70*r.median_patient_absolute_contribution_sum/arch.median_patient_absolute_contribution_sum.max(),c=SIG_COLORS[r.signature_id],edgecolor="white",linewidth=.6)
        ax.annotate(SIG_SHORT[r.signature_id],(r.median_patient_dominance_ratio,r.median_patient_cancellation_index),xytext=architecture_offsets[r.signature_id],textcoords="offset points",fontsize=6.5)
    ax.set_xlabel("Median dominance ratio"); ax.set_ylabel("Median cancellation index"); ax.set_xlim(0,.82); ax.set_ylim(-.08,.85); clean(ax,"both"); ax.set_title("Formula-level change architecture",fontweight="bold")

    ax=fig.add_subplot(gs[1,0]); panel_label(ax,"C")
    rec=recurrence.set_index("signature_id").loc[SIG_ORDER]; y=np.arange(8)[::-1]
    ax.barh(y,rec.cohort_mean_leading_gene_recurrence,color=[SIG_COLORS[s] for s in SIG_ORDER],alpha=.82)
    ax.set_yticks(y,[SIG_SHORT[s] for s in SIG_ORDER]); ax.set_xlim(0,1.08); ax.set_xlabel("Cohorts sharing the recurring leader")
    for pos,sid in zip(y,SIG_ORDER):
        r=rec.loc[sid]; ax.text(min(r.cohort_mean_leading_gene_recurrence+.025,1.0),pos,str(r.recurring_cohort_mean_leading_gene_set),va="center",fontsize=6.5)
    clean(ax); ax.set_title("T48 leading-contributor recurrence",fontweight="bold")

    ax=fig.add_subplot(gs[1,1]); panel_label(ax,"D")
    merged=arch.merge(rank[["signature_id","pooled_spearman_rho","median_cohort_rank_displacement_pp"]],on="signature_id")
    for r in merged.itertuples():
        size=30+100*abs(r.stage3_pooled_delta_z)
        ax.scatter(r.median_patient_cancellation_index,r.pooled_spearman_rho,s=size,c=SIG_COLORS[r.signature_id],edgecolor="white",linewidth=.6)
        ax.annotate(SIG_SHORT[r.signature_id],(r.median_patient_cancellation_index,r.pooled_spearman_rho),xytext=(3,3),textcoords="offset points",fontsize=6.7)
    ax.set_xlabel("Median cancellation index"); ax.set_ylabel("T48 pooled rank persistence, ρ"); ax.set_xlim(-.03,.85); ax.set_ylim(.35,.75); clean(ax,"both"); ax.set_title("Architecture and repeated-measurement properties",fontweight="bold")
    merged.to_csv(SRCDIR/"Figure_6D_architecture_rank.csv",index=False)
    save(fig,6)


def main() -> None:
    configure()
    for fn in [figure1,figure2,figure3,figure4,figure5,figure6]: fn()
    report={"figures":6,"backend":"Python matplotlib/seaborn","formats":["SVG","PDF","TIFF 600 dpi","PNG preview"],"status":"PASS"}
    (PREVIEW/"figure_build_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
