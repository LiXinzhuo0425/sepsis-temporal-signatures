#!/usr/bin/env python3
"""Download and verify the two public GSE106878 reconstruction inputs."""

from __future__ import annotations

import hashlib
import shutil
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DESTINATION = ROOT / "_sources" / "gse106878_reconstruction" / "downloads"
FILES = {
    "GSE106878_series_matrix.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106878/matrix/GSE106878_series_matrix.txt.gz",
        "sha256": "8a03fda06ee787eb0570b6dfbf95e7b2e33ffac7f381772fac012913262a69f0",
    },
    "GPL10295_self_full.txt": {
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10295&targ=self&form=text&view=full",
        "sha256": "a0d1e3fcb7e5b50bba875461f22197d7641863eee9c4fecc5e8ae6a6a5d15582",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(name: str, specification: dict[str, str]) -> None:
    destination = DESTINATION / name
    if destination.exists() and sha256(destination) == specification["sha256"]:
        print(f"verified existing file: {name}")
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(specification["url"], headers={"User-Agent": "sepsis-longitudinal-reproducibility-package/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as output:
        shutil.copyfileobj(response, output)
    observed = sha256(temporary)
    if observed != specification["sha256"]:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"SHA-256 mismatch for {name}: {observed}")
    temporary.replace(destination)
    print(f"downloaded and verified: {name}")


def main() -> None:
    DESTINATION.mkdir(parents=True, exist_ok=True)
    for filename, specification in FILES.items():
        fetch(filename, specification)


if __name__ == "__main__":
    main()
