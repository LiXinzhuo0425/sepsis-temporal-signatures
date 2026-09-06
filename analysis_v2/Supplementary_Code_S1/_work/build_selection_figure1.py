#!/usr/bin/env python3
"""Build current Figure 1 with Python/matplotlib only.

Figure contract
---------------
Conclusion: dataset and signature eligibility were governed by dated,
pre-outcome rules before fixed formulas entered the longitudinal analysis.
Evidence logic: panel A separates the six primary cohort families from the
20 July 2026 public-data update; panel B shows the 34-candidate freeze and
the complete protein-coding A2 panel; panel C distinguishes two primary
longitudinal questions from complementary interpretation layers.
Export: flattened RGB LZW TIFF at 600 dpi plus PDF/SVG and PNG QA previews.
Review risks addressed: no PRISMA framing, no implication that GSE106878 was
pooled with the primary cohorts, and no conflation of arm sensitivity with
hydrocortisone efficacy.
"""

from __future__ import annotations

from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "final" / "Figures"
QA = ROOT / "_qa" / "selection_figure1"

NAVY = "#17365D"
BLUE = "#2C7FB8"
LIGHT_BLUE = "#E3F0F7"
TEAL = "#1B9E77"
LIGHT_TEAL = "#E1F3ED"
ORANGE = "#D95F02"
LIGHT_ORANGE = "#FCEBDE"
PURPLE = "#756BB1"
LIGHT_PURPLE = "#EEEAF6"
GREY = "#66727E"
LIGHT_GREY = "#F1F3F5"
DARK = "#1F2933"


def box(ax, xy, width, height, title, subtitle="", face=LIGHT_BLUE, edge=BLUE,
        title_color=DARK, subtitle_color=GREY, title_size=9.0, subtitle_size=8.0,
        lw=1.2, radius=0.02):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y), width, height,
        boxstyle=f"round,pad=0.009,rounding_size={radius}",
        linewidth=lw, edgecolor=edge, facecolor=face,
        transform=ax.transAxes, clip_on=False,
    )
    ax.add_patch(patch)
    title_size, subtitle_size = max(8.0, title_size), max(8.0, subtitle_size)
    title_height = len(title.splitlines()) * title_size * 1.12
    subtitle_height = len(subtitle.splitlines()) * subtitle_size * 1.12 if subtitle else 0
    total_height = title_height + subtitle_height + (4 if subtitle else 0)
    center = (x + width * 0.5, y + height * 0.5)
    title_artist = ax.annotate(title, xy=center, xycoords="axes fraction",
            xytext=(0, (total_height-title_height)/2), textcoords="offset points",
            ha="center", va="center", color=title_color,
            fontsize=title_size, fontweight="bold", linespacing=1.12)
    artists = [title_artist]
    if subtitle:
        artists.append(ax.annotate(subtitle, xy=center, xycoords="axes fraction",
                xytext=(0, -(total_height-subtitle_height)/2), textcoords="offset points",
                ha="center", va="center", color=subtitle_color,
                fontsize=subtitle_size, linespacing=1.12))
    patch._content_artists = artists
    return patch


def arrow(ax, start, end, color=GREY, lw=1.2):
    arr = FancyArrowPatch(start, end, transform=ax.transAxes,
                          arrowstyle="-|>", mutation_scale=9, linewidth=lw,
                          color=color, shrinkA=3, shrinkB=3)
    ax.add_patch(arr)


def panel_heading(ax, label, title):
    ax.text(0.0, 1.02, label, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=12, fontweight="bold", color=DARK)
    ax.text(0.045, 1.02, title, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=10.5, fontweight="bold", color=DARK)


def build():
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(7.5, 8.5), facecolor="white")
    grid = fig.add_gridspec(3, 1, height_ratios=[1.1, 1.1, 0.95], hspace=0.30,
                            left=0.045, right=0.975, top=0.95, bottom=0.035)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[2, 0])
    for ax in (ax_a, ax_b, ax_c):
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    panel_heading(ax_a, "A", "Cohort identification and roles")
    box(ax_a, (0.04, 0.76), 0.25, 0.18, "6 primary\ncohort families", face=LIGHT_BLUE, edge=BLUE)
    box(ax_a, (0.375, 0.76), 0.25, 0.18, "299 baseline\npatients", face=LIGHT_BLUE, edge=BLUE)
    box(ax_a, (0.71, 0.76), 0.25, 0.18, "264 with an early\npaired sample", face=LIGHT_BLUE, edge=BLUE)
    arrow(ax_a, (0.29, 0.85), (0.375, 0.85))
    arrow(ax_a, (0.625, 0.85), (0.71, 0.85))

    ax_a.text(0.04, 0.66, "PUBLIC-DATA UPDATE", transform=ax_a.transAxes,
              color=NAVY, fontsize=8, fontweight="bold")
    box(ax_a, (0.04, 0.30), 0.28, 0.22, "Search refreshed", "through 20 Jul 2026",
        face=LIGHT_GREY, edge=GREY)
    box(ax_a, (0.37, 0.30), 0.24, 0.22, "4 candidates", "formally adjudicated",
        face=LIGHT_GREY, edge=GREY)
    arrow(ax_a, (0.32, 0.41), (0.37, 0.41))
    box(ax_a, (0.66, 0.40), 0.30, 0.24, "GSE106878", "47 baseline-to-24-h pairs\nseparate replication",
        face=LIGHT_TEAL, edge=TEAL, title_color=TEAL)
    box(ax_a, (0.66, 0.075), 0.30, 0.28, "3 not analyzed", "E-MEXP-3850: no matrix\nPRJEB111201: raw reads only\nphs003608: controlled matrix",
        face=LIGHT_ORANGE, edge=ORANGE, title_color=ORANGE)
    arrow(ax_a, (0.61, 0.44), (0.66, 0.52), color=TEAL)
    arrow(ax_a, (0.61, 0.38), (0.66, 0.215), color=ORANGE)
    ax_a.text(0.04, 0.13, "The updated candidate set did not alter\nthe primary synthesis.",
              transform=ax_a.transAxes, fontsize=8, color=GREY, linespacing=1.3)

    panel_heading(ax_b, "B", "Signature eligibility frozen before outcomes")
    box(ax_b, (0.05, 0.68), 0.28, 0.27, "34 candidates", "frozen 15 Jul 2026\nbefore outcome inspection",
        face=LIGHT_BLUE, edge=BLUE, title_color=BLUE)
    box(ax_b, (0.39, 0.68), 0.26, 0.27, "Reproducibility", "A1 0  |  A2 10\nB 19  |  X 5",
        face=LIGHT_GREY, edge=GREY)
    arrow(ax_b, (0.33, 0.815), (0.39, 0.815))
    box(ax_b, (0.71, 0.68), 0.25, 0.27, "8 protein-coding\nA2 signatures", "included without refitting",
        face=LIGHT_TEAL, edge=TEAL, title_color=TEAL)
    box(ax_b, (0.71, 0.32), 0.25, 0.30, "2 lncRNA\nA2 signatures", "cross-platform mapping\nunresolved",
        face=LIGHT_ORANGE, edge=ORANGE, title_color=ORANGE)
    arrow(ax_b, (0.65, 0.83), (0.71, 0.815), color=TEAL)
    arrow(ax_b, (0.65, 0.73), (0.71, 0.47), color=ORANGE)
    box(ax_b, (0.05, 0.32), 0.58, 0.30,
        "Fixed publication-time rule required",
        "gene set + algebra + direction\nno feature selection, refitting or recalibration",
        face=LIGHT_PURPLE, edge=PURPLE, title_color=PURPLE)
    ax_b.text(0.05, 0.20,
              "Included: MetaScore | SeptiCyte LAB | FAIM3:PLAC8 | sNIP |\n"
              "B/V MetaScore | Herberg DRS | Lin | Severe-or-Mild",
              transform=ax_b.transAxes, fontsize=8, color=DARK, linespacing=1.35)
    ax_b.text(0.05, 0.015,
              "Selection did not use trajectory, effect magnitude, clinical association or statistical significance.",
              transform=ax_b.transAxes, fontsize=8, color=GREY)

    panel_heading(ax_c, "C", "Fixed-formula longitudinal evidence framework")
    ax_c.text(0.02, 0.86, "PRIMARY QUESTIONS", transform=ax_c.transAxes,
              color=NAVY, fontsize=8, fontweight="bold")
    box(ax_c, (0.02, 0.60), 0.45, 0.20, "Mean within-patient change", "baseline-standardized ΔZ at ~24 h and ~48 h",
        face=LIGHT_BLUE, edge=BLUE, title_color=BLUE, title_size=9.3)
    box(ax_c, (0.53, 0.60), 0.45, 0.20, "Patient-rank persistence", "Spearman ρ + percentile-rank displacement",
        face=LIGHT_BLUE, edge=BLUE, title_color=BLUE, title_size=9.3)
    ax_c.text(0.02, 0.48, "COMPLEMENTARY INTERPRETATION", transform=ax_c.transAxes,
              color=NAVY, fontsize=8, fontweight="bold")
    complementary = [
        (0.02, "Cross-signature\nconcordance", LIGHT_TEAL, TEAL),
        (0.27, "Cohort-specific\nclinical associations", LIGHT_ORANGE, ORANGE),
        (0.52, "Directional replication\n+ arm sensitivity", LIGHT_TEAL, TEAL),
        (0.77, "Exact gene-level\ncontributions", LIGHT_PURPLE, PURPLE),
    ]
    for x, title, face, edge in complementary:
        box(ax_c, (x, 0.21), 0.21, 0.19, title, face=face, edge=edge,
            title_color=edge, title_size=8.1)
    ax_c.text(0.50, 0.06,
              "Published gene membership and score algebra held fixed throughout",
              transform=ax_c.transAxes, ha="center", va="center",
              fontsize=8, color=GREY, fontweight="bold")

    FINAL.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fonts = [t.get_fontsize() for ax in (ax_a, ax_b, ax_c) for t in ax.texts]
    assert min(fonts) >= 8 and max(fonts) <= 12
    for ax in (ax_a, ax_b, ax_c):
        for patch in ax.patches:
            for artist in getattr(patch, "_content_artists", []):
                outer = patch.get_window_extent(renderer)
                inner = artist.get_window_extent(renderer)
                assert outer.contains(inner.x0, inner.y0) and outer.contains(inner.x1, inner.y1), artist.get_text()
    (QA / "Fig1_font_qa.json").write_text(json.dumps({"font":"Arial","min_pt":min(fonts),"max_pt":max(fonts),"text_elements":len(fonts),"box_containment":"passed"},indent=2))
    pdf = QA / "Fig1.pdf"
    svg = QA / "Fig1.svg"
    png = QA / "Fig1_preview.png"
    raw_tif = QA / "Fig1_raw.tif"
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    fig.savefig(png, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(raw_tif, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    with Image.open(raw_tif) as image:
        rgb = image.convert("RGB")
        rgb.save(FINAL / "Fig1.tif", compression="tiff_lzw", dpi=(600, 600))
        grey = rgb.convert("L")
        grey.thumbnail((1800, 1800))
        grey.save(QA / "Fig1_grayscale_preview.png")


if __name__ == "__main__":
    build()
