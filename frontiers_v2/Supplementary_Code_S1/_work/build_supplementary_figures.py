#!/usr/bin/env python3
"""Build Frontiers V2 Supplementary Figures S1-S6 with Python only."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "_work" / "analysis"
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
PINK = "#CC79A7"
GREY = "#6E7781"
LIGHT = "#E9EEF3"
DARK = "#1F2933"
RED = "#B2182B"

SIG_ORDER = ["SIG001", "SIG002", "SIG003", "SIG004", "SIG022", "SIG023", "SIG033", "SIG034"]
SIG_SHORT = {
    "SIG001": "MetaScore", "SIG002": "SeptiCyte", "SIG003": "FAIM3:PLAC8", "SIG004": "sNIP",
    "SIG022": "B/V MetaScore", "SIG023": "Herberg DRS", "SIG033": "Lin score", "SIG034": "Severe/Mild",
}
SIG_COLORS = {
    "SIG001": BLUE, "SIG002": GREEN, "SIG003": PINK, "SIG004": ORANGE,
    "SIG022": SKY, "SIG023": YELLOW, "SIG033": PURPLE, "SIG034": GREY,
}
TIME_COLORS = {"T24": BLUE, "T48": ORANGE}


def configure() -> None:
    sns.set_theme(style="white", context="paper")
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.8, "lines.linewidth": 1.7, "lines.markersize": 5,
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    })


def panel_label(ax: plt.Axes, label: str, x: float = -0.10, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")


def clean(ax: plt.Axes, grid: str | None = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(True, axis=grid, color="#D8DEE4", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)


def save(fig: plt.Figure, number: int) -> None:
    stem = f"Figure_S{number}"
    svg = FIGDIR / f"{stem}.svg"
    pdf = FIGDIR / f"{stem}.pdf"
    tiff = FIGDIR / f"{stem}.tiff"
    preview = PREVIEW / f"{stem}_preview.png"
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(tiff, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(tiff) as image:
        image.convert("RGB").save(tiff, compression="tiff_lzw", dpi=(600, 600))
    fig.savefig(preview, dpi=180, bbox_inches="tight")
    plt.close(fig)


def figure_s1() -> None:
    missing = pd.read_excel(INPUT, sheet_name="S22_LandmarkMissingness")
    scheduled = missing[missing["landmark_scheduled"].eq("YES")].copy()
    scheduled = scheduled.sort_values(["time_window", "dataset"])
    tipping = pd.read_excel(INPUT, sheet_name="S24_TippingPoint")
    tipping = tipping[tipping["signature_id"].astype(str).str.fullmatch(r"SIG\d+")].copy()
    tipping_sig = tipping[tipping["ci_excludes_zero"].eq("YES")].copy()
    scheduled.to_csv(SRCDIR / "Figure_S1A_landmark_missingness.csv", index=False)
    tipping.to_csv(SRCDIR / "Figure_S1B_tipping_point.csv", index=False)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(180 * MM, 105 * MM), gridspec_kw={"width_ratios": [1.12, 1]}, constrained_layout=True)
    y = np.arange(len(scheduled))
    observed = scheduled["observed_pairs_n"].to_numpy(float)
    absent = scheduled["missing_followup_n"].to_numpy(float)
    ax_a.barh(y, observed, color=BLUE, height=.68, label="Observed pair")
    ax_a.barh(y, absent, left=observed, color="#E6A57E", height=.68, label="Missing follow-up")
    labels = [f"{r.dataset}  {r.time_window}" for r in scheduled.itertuples()]
    ax_a.set_yticks(y, labels)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Baseline-score-eligible patients")
    ax_a.set_title("Scheduled landmark availability", fontweight="bold")
    for i, r in enumerate(scheduled.itertuples()):
        ax_a.text(float(r.baseline_score_eligible_n) + 1.5, i, f"{r.missing_percent:.0f}% missing", va="center", fontsize=6.2, color=GREY)
    ax_a.set_xlim(0, scheduled["baseline_score_eligible_n"].max() * 1.26)
    ax_a.legend(frameon=False, loc="lower right")
    clean(ax_a)
    panel_label(ax_a, "A")

    shown_sigs = [sid for sid in SIG_ORDER if sid in set(tipping_sig["signature_id"])]
    ymap = {sid: i for i, sid in enumerate(shown_sigs[::-1])}
    for landmark, offset in [("T24", .13), ("T48", -.13)]:
        sub = tipping_sig[tipping_sig["time_window"].eq(landmark)]
        ypos = [ymap[s] + offset for s in sub["signature_id"]]
        ax_b.scatter(sub["absolute_shift_for_ci_to_touch_zero"], ypos, s=32, color=TIME_COLORS[landmark], edgecolor="white", linewidth=.6, zorder=3, label=landmark)
        for x, yy in zip(sub["absolute_shift_for_ci_to_touch_zero"], ypos):
            ax_b.plot([0, x], [yy, yy], color=TIME_COLORS[landmark], alpha=.45, lw=1.6)
            ax_b.text(x + .11, yy, f"{x:.2f}", va="center", fontsize=6.2, color=TIME_COLORS[landmark])
    ax_b.set_yticks(np.arange(len(shown_sigs)), [SIG_SHORT[s] for s in shown_sigs[::-1]])
    ax_b.set_xlabel("Absolute common shift needed for the 95% CI\nto reach zero, ΔZ")
    ax_b.set_xlim(0, max(7.05, tipping_sig["absolute_shift_for_ci_to_touch_zero"].max() + .65))
    ax_b.set_title("Location-shift sensitivity", fontweight="bold")
    ax_b.legend(frameon=False, loc="lower right")
    clean(ax_b)
    panel_label(ax_b, "B")
    save(fig, 1)


def figure_s2() -> None:
    cohort = pd.read_csv(AN / "02_rank_persistence_cohort.csv")
    cohort = cohort[cohort["paired_n"].ge(10)].copy()
    cohort.to_csv(SRCDIR / "Figure_S2_rank_persistence_all_cohorts.csv", index=False)
    shown = cohort[cohort["primary_independent_eligible"].eq(True)].copy()

    fig, axes = plt.subplots(4, 2, figsize=(180 * MM, 225 * MM), constrained_layout=True, sharex=True)
    letters = list("ABCDEFGH")
    for ax, sid, letter in zip(axes.flat, SIG_ORDER, letters):
        sub = shown[shown["signature_id"].eq(sid)].copy()
        sub["label"] = sub["dataset"] + "  " + sub["time_window"]
        sub = sub.sort_values(["time_window", "dataset"], ascending=[False, False]).reset_index(drop=True)
        y = np.arange(len(sub))
        for i, r in sub.iterrows():
            color = TIME_COLORS[r["time_window"]]
            ax.plot([r["bootstrap_ci95_lower"], r["bootstrap_ci95_upper"]], [i, i], color=color, lw=2.1)
            ax.scatter(r["spearman_rho"], i, s=25, color=color, edgecolor="white", lw=.5, zorder=3)
            ax.text(1.05, i, f"n={int(r['paired_n'])}", ha="left", va="center", fontsize=5.7, color=GREY)
        ax.axvline(0, color=DARK, lw=.9)
        ax.set_yticks(y, sub["label"])
        ax.set_xlim(-.65, 1.20)
        ax.set_title(SIG_SHORT[sid], loc="left", fontweight="bold", color=SIG_COLORS[sid])
        ax.set_xlabel("Spearman ρ" if ax in axes[-1, :] else "")
        clean(ax)
        panel_label(ax, letter, x=-.13)
    fig.legend(handles=[Line2D([0], [0], marker="o", color=BLUE, label="T24", lw=2), Line2D([0], [0], marker="o", color=ORANGE, label="T48", lw=2)],
               loc="lower center", bbox_to_anchor=(.5, -.01), ncol=2, frameon=False)
    save(fig, 2)


def _clinical_order(dataset: str) -> list[str]:
    return {
        "GSE95233": ["Admission", "D2/D3"],
        "GSE110487": ["Baseline", "T48"],
        "GSE57065": ["Baseline", "T24", "T48"],
        "GSE54514": ["Baseline", "T24", "T48", "T72", "Day 5"],
    }[dataset]


def _clinical_colors(dataset: str) -> dict[str, str]:
    favorable = {"Survivor", "R", "SAPSII-Low"}
    statuses = {
        "GSE95233": ["Survivor", "Non-survivor"],
        "GSE110487": ["R", "NR"],
        "GSE57065": ["SAPSII-Low", "SAPSII-High"],
        "GSE54514": ["Survivor", "Non-survivor"],
    }[dataset]
    return {s: (BLUE if s in favorable else RED) for s in statuses}


def figure_s3() -> None:
    traj = pd.read_csv(AN / "05_clinical_course_trajectories.csv")
    traj.to_csv(SRCDIR / "Figure_S3_clinical_group_trajectories.csv", index=False)
    datasets = ["GSE95233", "GSE110487", "GSE57065", "GSE54514"]
    primary = SIG_ORDER[:4]
    fig, axes = plt.subplots(4, 4, figsize=(180 * MM, 225 * MM), constrained_layout=True)
    for col, dataset in enumerate(datasets):
        order = _clinical_order(dataset)
        colors = _clinical_colors(dataset)
        for row, sid in enumerate(primary):
            ax = axes[row, col]
            sub = traj[(traj["dataset"].eq(dataset)) & (traj["signature_id"].eq(sid))].copy()
            xmap = {t: i for i, t in enumerate(order)}
            for status, g in sub.groupby("clinical_status", sort=False):
                g = g.assign(x=g["display_time"].map(xmap)).sort_values("x")
                c = colors[status]
                ax.plot(g["x"], g["mean_score_z"], marker="o", color=c, label=status, lw=1.5, ms=3.8)
                ax.fill_between(g["x"].to_numpy(float), g["bootstrap_ci95_lower"].to_numpy(float), g["bootstrap_ci95_upper"].to_numpy(float), color=c, alpha=.12, linewidth=0)
                for _, r in g.iterrows():
                    direction = 1 if c == RED else -1
                    xoff = -6 if c == RED else 6
                    ax.annotate(f"n={int(r['n'])}", (r["x"], r["mean_score_z"]), xytext=(xoff, 7 * direction), textcoords="offset points",
                                ha="right" if c == RED else "left", va="center", fontsize=4.8, color=c)
            ax.axhline(0, color="#B8C0C8", lw=.7)
            ax.set_xticks(range(len(order)), order, rotation=35 if len(order) > 3 else 0, ha="right" if len(order) > 3 else "center")
            if row == 0:
                ax.set_title(dataset, fontweight="bold")
                panel_label(ax, list("ABCD")[col], x=-.10)
            if col == 0:
                ax.set_ylabel(f"{SIG_SHORT[sid]}\nScore z")
            else:
                ax.set_ylabel("")
            clean(ax, grid="y")
            if row == 0:
                ax.legend(frameon=False, fontsize=5.6, loc="best", handlelength=1.2, handletextpad=.3)
    save(fig, 3)


def figure_s4() -> None:
    mortality = pd.read_csv(AN / "07_gse106878_treatment_adjusted_mortality.csv")
    rank = pd.read_csv(AN / "07_gse106878_rank_persistence.csv")
    mortality.to_csv(SRCDIR / "Figure_S4A_gse106878_mortality.csv", index=False)
    rank.to_csv(SRCDIR / "Figure_S4B_gse106878_rank.csv", index=False)
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(180 * MM, 105 * MM), constrained_layout=True)
    y = np.arange(len(SIG_ORDER))[::-1]
    mort = mortality.set_index("signature_id").loc[SIG_ORDER]
    for pos, sid in zip(y, SIG_ORDER):
        r = mort.loc[sid]
        ax_a.plot([r["ci95_lower"], r["ci95_upper"]], [pos, pos], color=SIG_COLORS[sid], lw=2.3)
        ax_a.scatter(r["beta_adverse_minus_favorable"], pos, color=SIG_COLORS[sid], s=28, edgecolor="white", lw=.5, zorder=3)
    ax_a.axvline(0, color=DARK, lw=1)
    ax_a.set_yticks(y, [SIG_SHORT[s] for s in SIG_ORDER])
    ax_a.set_xlabel("Adjusted non-survivor minus survivor ΔZ")
    ax_a.set_title("Treatment-adjusted 28-day mortality context", fontweight="bold")
    clean(ax_a); panel_label(ax_a, "A")

    rr = rank.set_index("signature_id").loc[SIG_ORDER]
    for pos, sid in zip(y, SIG_ORDER):
        r = rr.loc[sid]
        ax_b.plot([r["rho_ci95_lower"], r["rho_ci95_upper"]], [pos, pos], color=SIG_COLORS[sid], lw=2.3)
        ax_b.scatter(r["spearman_rho"], pos, color=SIG_COLORS[sid], s=28, edgecolor="white", lw=.5, zorder=3)
    ax_b.axvline(0, color=DARK, lw=1)
    ax_b.set_yticks(y, [])
    ax_b.set_xlim(-.1, 1.02)
    ax_b.set_xlabel("Baseline-to-T24 Spearman ρ")
    ax_b.set_title("Patient-rank persistence", fontweight="bold")
    clean(ax_b); panel_label(ax_b, "B")
    save(fig, 4)


def figure_s5() -> None:
    null = pd.read_csv(AN / "08_formula_preserving_null_distribution.csv")
    summary = pd.read_csv(AN / "08_formula_preserving_background_summary.csv")
    null.to_csv(SRCDIR / "Figure_S5_formula_preserving_null_distributions.csv", index=False)
    summary.to_csv(SRCDIR / "Figure_S5_observed_effects_and_statistics.csv", index=False)
    fig, axes = plt.subplots(4, 2, figsize=(180 * MM, 190 * MM), constrained_layout=True)
    letters = list("ABCDEFGH")
    for ax, sid, letter in zip(axes.flat, SIG_ORDER, letters):
        g = null[(null["signature_id"].eq(sid)) & (null["status"].eq("SUCCESS"))]
        r = summary[summary["signature_id"].eq(sid)].iloc[0]
        sns.histplot(g["pseudo_mean_delta_z"], bins=30, stat="density", color="#BFC7CF", edgecolor="white", linewidth=.4, ax=ax)
        sns.kdeplot(g["pseudo_mean_delta_z"], color=GREY, lw=1.2, ax=ax, clip_on=False)
        ax.axvline(r["observed_mean_delta_z"], color=SIG_COLORS[sid], lw=2.5)
        ax.axvline(0, color=DARK, lw=.8, ls=":")
        ax.text(.97, .92, f"Pemp = {r['empirical_two_sided_p']:.3f}\nexact-bin match = {r['median_exact_bin_match_fraction']:.0%}",
                transform=ax.transAxes, ha="right", va="top", fontsize=6.3, color=DARK)
        ax.set_title(SIG_SHORT[sid], loc="left", fontweight="bold", color=SIG_COLORS[sid])
        ax.set_xlabel("Pseudo-formula mean ΔZ")
        ax.set_ylabel("Density")
        clean(ax, grid="y")
        panel_label(ax, letter, x=-.11)
        if sid == "SIG001":
            ax.legend(handles=[Line2D([0], [0], color=GREY, lw=7, alpha=.45, label="Matched pseudo-formulas"),
                               Line2D([0], [0], color=SIG_COLORS[sid], lw=2.5, label="Observed formula")],
                      loc="upper left", frameon=False, fontsize=5.8, handlelength=1.7)
    save(fig, 5)


def _short_pathway(name: str) -> str:
    return name.replace("HALLMARK_", "").replace("_", " ").title().replace("Tnfa", "TNFα").replace("Nfkb", "NF-κB")


def figure_s6() -> None:
    pathway = pd.read_excel(INPUT, sheet_name="S13_Pathway")
    pathway = pathway[(pathway["analysis_set"].eq("INDEPENDENT_ONLY")) & (pathway["tier"].eq("PRESET_PRIMARY")) & (pathway["time_label"].isin(["T24", "T48"]))].copy()
    cell = pd.read_excel(INPUT, sheet_name="S14_CellSource")
    cell = cell[cell["time_window"].isin(["T24", "T48"])].copy()
    pathway.to_csv(SRCDIR / "Figure_S6A_prespecified_pathway_correlations.csv", index=False)
    cell.to_csv(SRCDIR / "Figure_S6B_hpa_broad_lineage_contributions.csv", index=False)

    pathways = pathway["pathway"].drop_duplicates().tolist()
    fig = plt.figure(figsize=(180 * MM, 180 * MM), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, .9])
    heat_axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    for ax, landmark, letter in zip(heat_axes, ["T24", "T48"], ["A", "B"]):
        sub = pathway[pathway["time_label"].eq(landmark)]
        mat = sub.pivot(index="signature_id", columns="pathway", values="pooled_spearman_rho").reindex(index=SIG_ORDER, columns=pathways)
        cbar = landmark == "T48"
        sns.heatmap(mat, ax=ax, cmap=sns.diverging_palette(240, 10, as_cmap=True), center=0, vmin=-.65, vmax=.65,
                    xticklabels=[_short_pathway(p) for p in pathways], yticklabels=[SIG_SHORT[s] for s in SIG_ORDER] if landmark == "T24" else [],
                    linewidths=.45, linecolor="white", cbar=cbar, cbar_kws={"label": "Pooled Spearman ρ", "shrink": .75})
        for i, sid in enumerate(SIG_ORDER):
            for j, p in enumerate(pathways):
                row = sub[(sub["signature_id"].eq(sid)) & (sub["pathway"].eq(p))]
                if len(row) and row.iloc[0]["fdr_within_analysis_window_tier"] < .05:
                    ax.text(j + .5, i + .5, "*", ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        ax.set_title(f"Prespecified Hallmark correlations · {landmark}", fontweight="bold")
        ax.set_xlabel(""); ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=58, labelsize=5.7)
        if landmark == "T48":
            ax.text(.985, .02, "* FDR < 0.05", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8,
                    color=DARK, bbox={"facecolor": "white", "edgecolor": "none", "alpha": .82, "pad": 1.5})
        panel_label(ax, letter, x=-.12)

    lineage_order = ["GRANULOCYTE", "MONOCYTE", "DENDRITIC", "B_LYMPHOCYTE", "T_LYMPHOCYTE", "BROAD_OR_UNRESOLVED"]
    lineage_labels = {
        "GRANULOCYTE": "Granulocyte", "MONOCYTE": "Monocyte", "DENDRITIC": "Dendritic",
        "B_LYMPHOCYTE": "B lymphocyte", "T_LYMPHOCYTE": "T lymphocyte", "BROAD_OR_UNRESOLVED": "Broad/unresolved",
    }
    lineage_colors = {
        "GRANULOCYTE": BLUE, "MONOCYTE": ORANGE, "DENDRITIC": GREEN,
        "B_LYMPHOCYTE": PURPLE, "T_LYMPHOCYTE": SKY, "BROAD_OR_UNRESOLVED": "#A8ADB4",
    }
    bar_axes = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    for ax, landmark, letter in zip(bar_axes, ["T24", "T48"], ["C", "D"]):
        sub = cell[cell["time_window"].eq(landmark)]
        pivot = sub.pivot_table(index="signature_id", columns="broad_lineage", values="absolute_share", aggfunc="sum", fill_value=0).reindex(SIG_ORDER).fillna(0)
        y = np.arange(len(SIG_ORDER))
        left = np.zeros(len(SIG_ORDER))
        for lineage in lineage_order:
            vals = pivot[lineage].to_numpy() if lineage in pivot else np.zeros(len(SIG_ORDER))
            ax.barh(y, vals, left=left, color=lineage_colors[lineage], height=.65, label=lineage_labels[lineage])
            left += vals
        ax.set_yticks(y, [SIG_SHORT[s] for s in SIG_ORDER] if landmark == "T24" else [])
        ax.invert_yaxis(); ax.set_xlim(0, 1)
        ax.set_xlabel("Absolute contribution share")
        ax.set_title(f"HPA broad-lineage annotation · {landmark}", fontweight="bold")
        clean(ax)
        panel_label(ax, letter, x=-.12)
    handles = [Line2D([0], [0], color=lineage_colors[x], lw=7, label=lineage_labels[x]) for x in lineage_order]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, -.015), ncol=3, frameon=False)
    save(fig, 6)


def main() -> None:
    configure()
    figure_s1()
    figure_s2()
    figure_s3()
    figure_s4()
    figure_s5()
    figure_s6()
    print(f"Wrote Supplementary Figures S1-S6 to {FIGDIR}")


if __name__ == "__main__":
    main()
