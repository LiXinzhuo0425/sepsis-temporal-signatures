#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_DIRECTORIES = {".git", "__pycache__"}
EXCLUDED_FILES = {".DS_Store"}
SELF_REFERENTIAL_ROOT_FILES = {"MANIFEST_SHA256.txt", "release_manifest.json"}
RELEASE_ID = "sepsis-temporal-signatures-v2.1.0"
RELEASE_DATE = "2026-09-06"
DEFAULT_VERSION_DOI = "10.5281/zenodo.21415496"
CONCEPT_DOI = "10.5281/zenodo.21415496"
REPOSITORY = "https://github.com/LiXinzhuo0425/sepsis-temporal-signatures"
VERSION_DOI = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_VERSION_DOI
PUBLICATION_STATUS = "prepared public release"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


files = []
for path in sorted(ROOT.rglob("*")):
    rel = path.relative_to(ROOT)
    if (
        not path.is_file()
        or path.name.endswith(".inspect.ndjson")
        or any(part in EXCLUDED_DIRECTORIES for part in rel.parts)
        or path.name in EXCLUDED_FILES
        or (len(rel.parts) == 1 and path.name in SELF_REFERENTIAL_ROOT_FILES)
    ):
        continue
    files.append(
        {
            "path": rel.as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": digest(path),
        }
    )

manifest = {
    "release_id": RELEASE_ID,
    "release_date": RELEASE_DATE,
    "doi": VERSION_DOI,
    "concept_doi": CONCEPT_DOI,
    "publication_status": PUBLICATION_STATUS,
    "repository": REPOSITORY,
    "files": files,
}
(ROOT / "release_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

# The JSON manifest cannot hash itself, but the outer SHA-256 manifest can and
# should protect it.  Generate that hash only after the JSON file is final.
manifest_entries = [
    *files,
    {
        "path": "release_manifest.json",
        "size_bytes": (ROOT / "release_manifest.json").stat().st_size,
        "sha256": digest(ROOT / "release_manifest.json"),
    },
]
(ROOT / "MANIFEST_SHA256.txt").write_text(
    "".join(f"{item['sha256']}  {item['path']}\n" for item in manifest_entries), encoding="utf-8"
)
print(
    json.dumps(
        {
            "release_id": manifest["release_id"],
            "content_files": len(files),
            "sha256_manifest_entries": len(manifest_entries),
        },
        indent=2,
    )
)
