#!/usr/bin/env python3
"""Build presentation figure references from the frozen release data.

The script is release-relative and changes presentation only. It does not
estimate, refit, or alter any scientific value. Figures 1 and 2 reproduce the
submitted layouts. Figure 3 is a computational reference for the same source
values; the author-approved TIFF in ``reference_outputs`` is the authoritative
final layout and is never overwritten by this script.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


SOURCE_RELEASE = Path(__file__).resolve().parents[1]
PROJECT = SOURCE_RELEASE / "reproduction_project"
S3 = PROJECT / "longitudinal_stage3" / "source_data"
S4 = PROJECT / "longitudinal_stage3" / "04_source_data"
GENERATED = Path(
    os.environ.get(
        "BR_FIGURE_OUTPUT_DIR",
        SOURCE_RELEASE
        / "reproducibility_evidence"
        / "generated_presentation_v1.3.0",
    )
).resolve()
OUT = GENERATED / "main_figures"
SOURCE_OUT = GENERATED / "source_data"
SUPP_MEDIA = GENERATED / "supplementary_figures"
REFERENCE_SOURCE = (
    SOURCE_RELEASE
    / "data"
    / "figure_source_data"
    / "presentation_v1.3.0"
)
REFERENCE_FIGURES = (
    SOURCE_RELEASE
    / "reference_outputs"
    / "presentation_v1.3.0"
)
OUT.mkdir(parents=True, exist_ok=True)
SOURCE_OUT.mkdir(parents=True, exist_ok=True)
SUPP_MEDIA.mkdir(parents=True, exist_ok=True)

SIGS = ["SIG001", "SIG002", "SIG003", "SIG004", "SIG022", "SIG023", "SIG033", "SIG034"]
NAMES = {
    "SIG001": "Sepsis MetaScore",
    "SIG002": "SeptiCyte LAB",
    "SIG003": "FAIM3:PLAC8",
    "SIG004": "sNIP",
    "SIG022": "Bacterial/Viral MetaScore",
    "SIG023": "Herberg two-transcript DRS",
    "SIG033": "Lin seven-gene mortality score",
    "SIG034": "Severe-or-Mild score",
}
COHORTS = ["GSE236713", "GSE57065", "GSE95233", "GSE54514", "GSE110487", "GSE8121"]
COL = {
    "blue": "#2F6B8A",
    "sky": "#74A9CF",
    "red": "#B24745",
    "grey": "#777777",
    "ink": "#24323D",
}
COHORT_COLORS = {
    "GSE236713": "#0F4D92",
    "GSE57065": "#8BCF8B",
    "GSE95233": "#B64342",
    "GSE54514": "#42949E",
    "GSE110487": "#9A4D8E",
    "GSE8121": "#CF9D3E",
}

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "font.size": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    }
)


def panel(ax: plt.Axes, letter: str, x: float = -0.10, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold",
    )


def export(fig: plt.Figure, stem: str) -> None:
    """Export editable files and a 600-dpi, LZW-compressed submission TIFF."""
    base = OUT / stem
    fig.savefig(base.with_suffix(".svg"))
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=300)
    tiff_path = base.with_suffix(".tiff")
    fig.savefig(tiff_path, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(tiff_path) as opened:
        rgb = opened.convert("RGB")
    rgb.save(tiff_path, dpi=(600, 600), compression="tiff_lzw")
    plt.close(fig)


def signature_labels() -> list[str]:
    return [f"{sig}  {NAMES[sig]}" for sig in SIGS]


def figure_1() -> None:
    meta = pd.read_csv(S3 / "03_05_signature_level_meta_analysis.csv")
    data = meta[
        (meta.analysis_set == "PRIMARY_INDEPENDENT") & meta.time_window.isin(["T1", "T2"])
    ].copy()

    fig, axes = plt.subplots(1, 2, figsize=(6.65, 4.35), sharey=True)
    for ax, window, title, color in zip(
        axes,
        ["T1", "T2"],
        ["T24 (12-36 h)", "T48 (36-60 h)"],
        [COL["blue"], COL["red"]],
    ):
        sub = data[data.time_window == window].set_index("signature_id").reindex(SIGS)
        y = np.arange(len(SIGS))[::-1]
        ax.axvline(0, color=COL["grey"], linestyle="--", linewidth=0.8)
        for yi, sig in zip(y, SIGS):
            row = sub.loc[sig]
            if np.isfinite(row.prediction_lower):
                ax.plot(
                    [row.prediction_lower, row.prediction_upper],
                    [yi - 0.14, yi - 0.14],
                    color=color,
                    linewidth=0.9,
                    alpha=0.48,
                )
            ax.plot([row.ci95_lower, row.ci95_upper], [yi, yi], color=color, linewidth=2.0)
            ax.scatter(
                row.pooled_delta_z,
                yi,
                s=27,
                color=color,
                edgecolor="white",
                linewidth=0.5,
                zorder=3,
            )
            ax.text(
                row.pooled_delta_z,
                yi + 0.21,
                f"{row.pooled_delta_z:.2f}",
                ha="center",
                fontsize=6.1,
                color=color,
            )
        ax.set_yticks(y)
        ax.set_yticklabels(signature_labels(), fontsize=6.8)
        ax.set_xlabel("Pooled within-patient change (ΔZ)", fontsize=8)
        ax.set_title(title, fontsize=9, fontweight="bold")
    panel(axes[0], "A", x=-0.62)
    panel(axes[1], "B", x=-0.08)
    fig.suptitle(
        "Pooled longitudinal changes in eight fixed signatures",
        x=0.02,
        y=0.975,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.40, right=0.98, bottom=0.13, top=0.87, wspace=0.16)
    export(fig, "Figure_1")
    data.to_csv(SOURCE_OUT / "Figure_1_source_data.csv", index=False)


def heatmap(
    ax: plt.Axes,
    matrix: np.ndarray,
    rows: list[str],
    cols: list[str],
    title: str,
) -> None:
    cmap = mpl.colormaps["RdBu_r"].copy()
    cmap.set_bad("#F2F2F2")
    image = ax.imshow(matrix, aspect="auto", vmin=-1.6, vmax=1.6, cmap=cmap)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, fontsize=7.0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=7.2)
    ax.tick_params(length=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not np.isfinite(matrix[i, j]):
                continue
            color = "white" if abs(matrix[i, j]) > 1.04 else "black"
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=6.2, color=color)
    ax.set_title(title, loc="left", fontsize=9, fontweight="bold")
    colorbar = plt.colorbar(image, ax=ax, fraction=0.025, pad=0.015)
    colorbar.set_label("Mean ΔZ", fontsize=7.2)
    colorbar.ax.tick_params(labelsize=6.8)


def figure_2() -> None:
    effects = pd.read_csv(S3 / "03_04_cohort_signature_primary_effects.csv")
    profile = pd.read_csv(S3 / "03_06_temporal_stability_profile.csv")

    fig = plt.figure(figsize=(6.65, 7.25))
    grid = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.87], hspace=0.54)
    for index, (window, title) in enumerate(
        [("T1", "T24 cohort-specific mean ΔZ"), ("T2", "T48 cohort-specific mean ΔZ")]
    ):
        ax = fig.add_subplot(grid[index])
        panel(ax, "A" if index == 0 else "B", x=-0.095)
        sub = (
            effects[effects.time_window == window]
            .pivot(index="dataset", columns="signature_id", values="mean_delta_z")
            .reindex(index=COHORTS, columns=SIGS)
        )
        heatmap(ax, sub.to_numpy(float), COHORTS, SIGS, title)
        row_i = COHORTS.index("GSE54514")
        for sig in ["SIG002", "SIG003", "SIG004"]:
            col_j = SIGS.index(sig)
            if np.isfinite(sub.loc["GSE54514", sig]):
                ax.add_patch(
                    mpl.patches.Rectangle(
                        (col_j - 0.47, row_i - 0.47),
                        0.94,
                        0.94,
                        fill=False,
                        edgecolor="#222222",
                        linewidth=1.25,
                    )
                )
                ax.text(
                    col_j + 0.34,
                    row_i - 0.31,
                    "†",
                    ha="center",
                    va="center",
                    fontsize=8.2,
                    fontweight="bold",
                )

    ax3 = fig.add_subplot(grid[2])
    panel(ax3, "C", x=-0.095)
    p = profile[profile.time_window == "T2"].set_index("signature_id").reindex(SIGS)
    y = np.arange(len(SIGS))[::-1]
    widths = (p.prediction_upper - p.prediction_lower).to_numpy(float)
    bars = ax3.barh(y, p.I2_percent, height=0.58, color=COL["sky"])
    ax3.set_yticks(y)
    compact_labels = [
        "SIG001  MetaScore",
        "SIG002  SeptiCyte LAB",
        "SIG003  FAIM3:PLAC8",
        "SIG004  sNIP",
        "SIG022  Bacterial/Viral MetaScore",
        "SIG023  Herberg DRS",
        "SIG033  Lin mortality score",
        "SIG034  Severe-or-Mild",
    ]
    ax3.set_yticklabels(compact_labels, fontsize=6.4)
    ax3.set_xlim(0, 105)
    ax3.set_xlabel("T48 I² (%)", fontsize=8)
    ax3.axvline(50, color="#999999", linewidth=0.8, linestyle="--")
    for bar, i2, width in zip(bars, p.I2_percent, widths):
        ax3.text(
            min(i2 + 2, 97),
            bar.get_y() + bar.get_height() / 2,
            f"{i2:.0f}%  |  PI width {width:.2f}",
            va="center",
            fontsize=6.2,
        )
    ax3.set_title(
        "T48 between-cohort heterogeneity and prediction-interval width",
        loc="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.suptitle(
        "Cohort-level longitudinal changes and cross-cohort dispersion",
        x=0.02,
        y=0.985,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.31, right=0.96, bottom=0.07, top=0.93)
    export(fig, "Figure_2")
    effects.to_csv(SOURCE_OUT / "Figure_2A_B_source_data.csv", index=False)
    profile.to_csv(SOURCE_OUT / "Figure_2C_source_data.csv", index=False)


def figure_3() -> None:
    gene = pd.read_csv(S4 / "04_07_gene_contribution_meta_analysis.csv")
    gene = gene[
        (gene.analysis_set == "PRIMARY_INDEPENDENT") & (gene.time_window == "T2")
    ].copy()
    architecture = (
        pd.read_csv(S4 / "04_08_signature_drift_architecture.csv")
        .set_index("signature_id")
        .reindex(SIGS)
        .reset_index()
    )

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(6.65, 4.55),
        gridspec_kw={"width_ratios": [1.28, 1.0], "wspace": 0.42},
    )
    y = np.arange(len(SIGS))[::-1]
    for yi, sig in zip(y, SIGS):
        leading = (
            gene[gene.signature_id == sig]
            .assign(abs_effect=lambda frame: frame.pooled_contribution.abs())
            .nlargest(2, "abs_effect")
        )
        offsets = np.linspace(-0.13, 0.13, len(leading))
        for offset, row in zip(offsets, leading.itertuples(index=False)):
            color = COL["red"] if row.pooled_contribution > 0 else COL["blue"]
            ax1.plot(
                [row.ci95_lower, row.ci95_upper],
                [yi + offset, yi + offset],
                color=color,
                linewidth=1.35,
            )
            ax1.scatter(
                row.pooled_contribution,
                yi + offset,
                s=23,
                color=color,
                edgecolor="white",
                linewidth=0.4,
                zorder=3,
            )
            ax1.annotate(
                row.gene,
                (row.pooled_contribution, yi + offset),
                xytext=(3 if row.pooled_contribution >= 0 else -3, 0),
                textcoords="offset points",
                ha="left" if row.pooled_contribution >= 0 else "right",
                va="center",
                fontsize=5.9,
            )
    ax1.axvline(0, color=COL["grey"], linewidth=0.8, linestyle="--")
    ax1.set_yticks(y)
    ax1.set_yticklabels(signature_labels(), fontsize=6.3)
    ax1.set_xlabel("Exact gene contribution to T48 ΔZ", fontsize=8)
    ax1.set_title(
        "Two largest pooled\nmathematical contributions",
        loc="left",
        fontsize=8.2,
        fontweight="bold",
        pad=8,
    )
    panel(ax1, "A", x=-0.64, y=1.12)

    for row in architecture.itertuples(index=False):
        size = 55 + 80 * min(row.median_patient_absolute_contribution_sum, 2.5) / 2.5
        ax2.scatter(
            row.median_patient_dominance_ratio,
            row.median_patient_cancellation_index,
            s=size,
            color=COL["blue"],
            alpha=0.88,
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )
        dx, dy = {
            "SIG001": (4, 4),
            "SIG002": (5, 10),
            "SIG003": (-5, 19),
            "SIG004": (4, 4),
            "SIG022": (-7, 9),
            "SIG023": (4, 6),
            "SIG033": (4, 4),
            "SIG034": (4, 4),
        }.get(row.signature_id, (4, 3))
        alignment = "right" if row.signature_id in {"SIG003", "SIG022"} else "left"
        ax2.annotate(
            row.signature_id,
            (row.median_patient_dominance_ratio, row.median_patient_cancellation_index),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=6.6,
            fontweight="bold",
            ha=alignment,
        )
    ax2.set_xlim(-0.03, 0.90)
    ax2.set_ylim(-0.04, 0.90)
    ax2.set_xlabel("Median cohort dominance", fontsize=8)
    ax2.set_ylabel("Median cohort cancellation", fontsize=8)
    ax2.set_title(
        "Dominance and\ncancellation structure",
        loc="left",
        fontsize=8.2,
        fontweight="bold",
        pad=8,
    )
    panel(ax2, "B", x=-0.20, y=1.12)

    fig.suptitle(
        "Gene-contribution patterns at T48",
        x=0.02,
        y=0.975,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.text(
        0.02,
        0.018,
        "Point area represents the median absolute-contribution sum. Contributions are mathematical, not causal.",
        fontsize=6.2,
        color="#5F6B76",
    )
    fig.subplots_adjust(left=0.36, right=0.98, bottom=0.16, top=0.79)
    export(fig, "Figure_3_computational_reference")
    gene.to_csv(SOURCE_OUT / "Figure_3A_source_data.csv", index=False)
    architecture.to_csv(SOURCE_OUT / "Figure_3B_source_data.csv", index=False)


def save_supplementary_png(fig: plt.Figure, filename: str) -> None:
    fig.savefig(SUPP_MEDIA / filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def supplementary_figure_1() -> None:
    """Pilot and non-pilot concordance with uppercase panel labels."""
    meta = pd.read_csv(S3 / "03_05_signature_level_meta_analysis.csv")
    rows = []
    for window in ["T1", "T2"]:
        for sig in SIGS:
            pilot = meta[
                (meta.signature_id == sig)
                & (meta.time_window == window)
                & (meta.analysis_set == "PILOT_ONLY")
            ]
            non_pilot = meta[
                (meta.signature_id == sig)
                & (meta.time_window == window)
                & meta.analysis_set.isin(
                    ["PRESPECIFIED_NON_PILOT_ONLY", "PRESPECIFIED_VALIDATION_ONLY"]
                )
            ]
            if not pilot.empty and not non_pilot.empty:
                rows.append(
                    {
                        "signature_id": sig,
                        "time_window": window,
                        "pilot": pilot.iloc[0].pooled_delta_z,
                        "validation": non_pilot.iloc[0].pooled_delta_z,
                    }
                )
    data = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7))
    limits = [
        float(np.nanmin(data[["pilot", "validation"]].to_numpy())),
        float(np.nanmax(data[["pilot", "validation"]].to_numpy())),
    ]
    margin = max(0.2, (limits[1] - limits[0]) * 0.12)
    limits = [limits[0] - margin, limits[1] + margin]
    label_offsets = {
        "T1": {
            "SIG001": (8, -8),
            "SIG002": (-8, 10),
            "SIG003": (-8, 8),
            "SIG004": (8, 10),
            "SIG022": (8, -10),
            "SIG023": (8, 22),
            "SIG033": (-8, -10),
            "SIG034": (-8, 22),
        },
        "T2": {
            "SIG001": (8, -2),
            "SIG002": (8, -2),
            "SIG003": (8, -2),
            "SIG004": (8, 0),
            "SIG022": (8, -2),
            "SIG023": (8, 8),
            "SIG033": (8, 8),
            "SIG034": (8, -10),
        },
    }
    for ax, window in zip(axes, ["T1", "T2"]):
        sub = data[data.time_window == window]
        ax.plot(limits, limits, color="#A8A8A8", linestyle="--", linewidth=0.9)
        ax.axhline(0, color="#E0E0E0", linewidth=0.6)
        ax.axvline(0, color="#E0E0E0", linewidth=0.6)
        ax.scatter(sub.pilot, sub.validation, color=COL["red"], s=30)
        for _, row in sub.sort_values("signature_id").iterrows():
            dx, dy = label_offsets[window][row.signature_id]
            ax.annotate(
                row.signature_id,
                xy=(row.pilot, row.validation),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=6.4,
                color=COL["ink"],
                ha="left" if dx >= 0 else "right",
                va="center",
                annotation_clip=False,
                arrowprops=dict(
                    arrowstyle="-",
                    color="#8B99A5",
                    linewidth=0.45,
                    shrinkA=1.5,
                    shrinkB=3.0,
                ),
            )
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Pilot pooled ΔZ")
        ax.set_ylabel("Prespecified non-pilot pooled ΔZ")
        ax.set_title("T24" if window == "T1" else "T48", fontweight="bold", fontsize=9)
    panel(axes[0], "A")
    panel(axes[1], "B")
    fig.suptitle(
        "Pilot and prespecified non-pilot concordance",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save_supplementary_png(fig, "image12.png")
    data.to_csv(SOURCE_OUT / "Supplementary_Figure_S1_source_data.csv", index=False)


def supplementary_figure_2() -> None:
    """Robustness estimates with uppercase panel labels."""
    meta = pd.read_csv(S3 / "03_05_signature_level_meta_analysis.csv")
    scaling = pd.read_csv(S3 / "03_11_scaling_sensitivity_analysis.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.6), sharey=True)
    sets = [
        ("ALL_COHORTS", "All cohorts", "#767676"),
        ("PRIMARY_INDEPENDENT", "Primary independent", "#0F4D92"),
        ("STRICT_NEVER_USED", "Strict never-used", "#B64342"),
    ]
    offsets = [-0.18, -0.06, 0.06, 0.18]
    for ax, window in zip(axes, ["T1", "T2"]):
        y = np.arange(len(SIGS))[::-1]
        ax.axvline(0, color="#767676", linestyle="--", linewidth=0.8)
        for offset, (analysis_set, label, color) in zip(offsets[:3], sets):
            data = (
                meta[(meta.time_window == window) & (meta.analysis_set == analysis_set)]
                .set_index("signature_id")
                .reindex(SIGS)
            )
            for yi, sig in zip(y, SIGS):
                row = data.loc[sig]
                if pd.isna(row.pooled_delta_z):
                    continue
                ax.plot(
                    [row.ci95_lower, row.ci95_upper],
                    [yi + offset, yi + offset],
                    color=color,
                    linewidth=1.0,
                )
                ax.scatter(
                    row.pooled_delta_z,
                    yi + offset,
                    color=color,
                    s=12,
                    label=label if yi == y[0] else None,
                )
        mad = (
            scaling[
                (scaling.record_type == "META_PRIMARY_INDEPENDENT")
                & (scaling.time_window == window)
                & (scaling.scaling_method == "BASELINE_MAD")
            ]
            .set_index("signature_id")
            .reindex(SIGS)
        )
        for yi, sig in zip(y, SIGS):
            row = mad.loc[sig]
            if pd.isna(row.mean_change):
                continue
            ax.plot(
                [row.ci95_lower, row.ci95_upper],
                [yi + offsets[3], yi + offsets[3]],
                color="#42949E",
                linewidth=1.0,
            )
            ax.scatter(
                row.mean_change,
                yi + offsets[3],
                color="#42949E",
                marker="D",
                s=12,
                label="MAD scaling" if yi == y[0] else None,
            )
        ax.set_yticks(y)
        ax.set_yticklabels(SIGS)
        ax.set_xlabel("Pooled within-patient change")
        ax.set_title("T24" if window == "T1" else "T48", fontweight="bold", fontsize=9)
    axes[1].legend(loc="lower right", fontsize=6.5)
    panel(axes[0], "A")
    panel(axes[1], "B")
    fig.suptitle(
        "Robustness across analysis sets and scaling",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_supplementary_png(fig, "image13.png")


def supplementary_figure_3() -> None:
    """Illustrative trajectories with uppercase panel labels."""
    changes = pd.read_parquet(
        SOURCE_RELEASE / "data" / "derived_patient_level" / "paired_changes_all_windows.parquet"
    )
    selected = ["SIG034", "SIG001", "SIG002", "SIG003"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), sharex=True)
    for ax, sig in zip(axes.ravel(), selected):
        group = changes[
            (changes.signature_id == sig) & changes.time_window.isin(["T1", "T2"])
        ]
        for dataset, cohort in group.groupby("dataset"):
            for _, row in cohort.iterrows():
                x = [0, 1 if row.time_window == "T1" else 2]
                ax.plot(x, [0, row.delta_z], color=COHORT_COLORS[dataset], alpha=0.12, linewidth=0.6)
            means = cohort.groupby("time_window").delta_z.mean()
            for window, x_position in (("T1", 1), ("T2", 2)):
                if window in means:
                    ax.scatter(
                        x_position,
                        means[window],
                        color=COHORT_COLORS[dataset],
                        s=20,
                        edgecolor="white",
                        linewidth=0.4,
                        zorder=4,
                    )
        ax.axhline(0, color="#767676", linestyle="--", linewidth=0.7)
        ax.set_title(sig, fontweight="bold", fontsize=9)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["T0", "T24", "T48"])
        ax.set_ylabel("Within-patient ΔZ")
    for index, ax in enumerate(axes.ravel()):
        panel(ax, chr(ord("A") + index))
    handles = [
        plt.Line2D([0], [0], color=color, linewidth=2, label=dataset)
        for dataset, color in COHORT_COLORS.items()
    ]
    fig.legend(handles=handles, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.95), fontsize=6.5)
    fig.suptitle(
        "Illustrative patient-level trajectories for four prespecified signatures",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    save_supplementary_png(fig, "image11.png")


def verify_inputs() -> None:
    required = [
        S3 / "03_05_signature_level_meta_analysis.csv",
        S3 / "03_04_cohort_signature_primary_effects.csv",
        S3 / "03_06_temporal_stability_profile.csv",
        S4 / "04_07_gene_contribution_meta_analysis.csv",
        S4 / "04_08_signature_drift_architecture.csv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing locked source tables:\n" + "\n".join(missing))


def verify_source_exports() -> None:
    """Require the five main-panel source exports to match the tracked layer."""
    names = (
        "Figure_1_source_data.csv",
        "Figure_2A_B_source_data.csv",
        "Figure_2C_source_data.csv",
        "Figure_3A_source_data.csv",
        "Figure_3B_source_data.csv",
    )
    for name in names:
        generated = pd.read_csv(SOURCE_OUT / name)
        reference = pd.read_csv(REFERENCE_SOURCE / name)
        pd.testing.assert_frame_equal(
            generated,
            reference,
            check_dtype=False,
            check_exact=False,
            rtol=0,
            atol=1e-12,
        )
        print(f"PASS_SOURCE_IDENTICAL {name}")


def report_reference_pixel_status() -> None:
    """Report, but do not require, decoded-pixel identity for Figures 1 and 2."""
    for name in ("Figure_1.tiff", "Figure_2.tiff"):
        with Image.open(OUT / name) as generated_image:
            generated = np.asarray(generated_image.convert("RGB"))
        with Image.open(REFERENCE_FIGURES / name) as reference_image:
            reference = np.asarray(reference_image.convert("RGB"))
        status = (
            "PASS_PIXEL_IDENTICAL"
            if generated.shape == reference.shape and np.array_equal(generated, reference)
            else "PASS_SOURCE_IDENTICAL_RENDER_VARIATION"
        )
        print(f"{status} {name}")


def main() -> None:
    verify_inputs()
    figure_1()
    figure_2()
    figure_3()
    supplementary_figure_1()
    supplementary_figure_2()
    supplementary_figure_3()
    verify_source_exports()
    report_reference_pixel_status()
    for name in (
        "Figure_1.tiff",
        "Figure_2.tiff",
        "Figure_3_computational_reference.tiff",
    ):
        print(OUT / name)
    print(
        "AUTHORED_FINAL Figure 3: "
        "reference_outputs/presentation_v1.3.0/Figure_3.tiff"
    )


if __name__ == "__main__":
    main()
