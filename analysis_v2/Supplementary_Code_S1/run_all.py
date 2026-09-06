#!/usr/bin/env python3
"""Run the frozen longitudinal analyses and current sensitivity in order."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "_work" / "reproduce_legacy_baseline.py",
    ROOT / "_work" / "run_v2_measurement_modules.py",
    ROOT / "_work" / "run_v2_clinical_anchoring.py",
    ROOT / "fetch_public_inputs.py",
    ROOT / "_work" / "reconstruct_gse106878.py",
    ROOT / "_work" / "run_gse106878_treatment_sensitivity.py",
    ROOT / "_work" / "run_gse106878_secondary_modules.py",
    ROOT / "_work" / "run_formula_preserving_background.py",
    ROOT / "_work" / "build_longitudinal_figures.py",
    ROOT / "_work" / "build_supplementary_figures.py",
    ROOT / "_work" / "build_selection_figure1.py",
]


def main() -> None:
    for script in SCRIPTS:
        print(f"\n[run] {script.relative_to(ROOT)}", flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
