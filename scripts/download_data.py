#!/usr/bin/env python
"""Fetch every input the pipeline needs into data/raw/.

    python scripts/download_data.py            # everything missing
    python scripts/download_data.py --check    # report sizes without downloading

Sources are public and unauthenticated:
  * Human Protein Atlas bulk TSVs (proteinatlas.org/about/download)
  * CELLxGENE Discover tumour scRNA-seq H5ADs (cellxgene.cziscience.com)

The HPA cell-line table is ~206 MB zipped / 1.07 GB extracted, so the whole set
is roughly 1.4 GB zipped and 2.5 GB on disk once extracted.
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
SCRNA = RAW / "scrna"

HPA_BASE = "https://www.proteinatlas.org/download"
HPA_FILES = {
    "proteinatlas.tsv.zip": f"{HPA_BASE}/proteinatlas.tsv.zip",
    "rna_single_cell_type.tsv.zip": f"{HPA_BASE}/tsv/rna_single_cell_type.tsv.zip",
    "subcellular_location.tsv.zip": f"{HPA_BASE}/tsv/subcellular_location.tsv.zip",
    "rna_celline.tsv.zip": f"{HPA_BASE}/tsv/rna_celline.tsv.zip",
}

# Tumour scRNA-seq with malignant cells annotated, chosen to span distinct
# lineages (epithelial, glial, haematological) while staying a manageable size.
SCRNA_FILES = {
    "breast_HTAPP895.h5ad": "https://datasets.cellxgene.cziscience.com/fe938180-0708-411a-b069-d4a5f49e101b.h5ad",
    "colorectal_kg146.h5ad": "https://datasets.cellxgene.cziscience.com/06409cab-8afb-4aa7-9e79-1c3c5dfe91fd.h5ad",
    "glioblastoma.h5ad": "https://datasets.cellxgene.cziscience.com/2e357b12-781b-4c90-9321-30206529cc4a.h5ad",
    "bladder_SCBO8.h5ad": "https://datasets.cellxgene.cziscience.com/59caba36-eb74-4d96-94cf-80a695d5dd44.h5ad",
    "dlbcl.h5ad": "https://datasets.cellxgene.cziscience.com/e9b97c4b-5000-488d-b2e4-505bfb327671.h5ad",
}

# The cell-line -> cancer-type mapping is parsed out of this help page, which
# lists every cell line beside its `cell_RNA_<cancer_type>` grouping.
HPA_HELP_URL = "https://www.proteinatlas.org/about/help/dataaccess"


def _report(path: Path) -> None:
    mb = path.stat().st_size / 1e6
    print(f"    {path.name}: {mb:,.1f} MB")


def download(url: str, dest: Path, force: bool = False) -> None:
    if dest.exists() and not force:
        print(f"  have {dest.name}")
        return
    print(f"  downloading {dest.name} ...", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as r, tmp.open("wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    tmp.replace(dest)
    _report(dest)


def unzip(path: Path) -> None:
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            target = path.parent / name
            if target.exists():
                print(f"  have {name}")
                continue
            print(f"  extracting {name} ...", flush=True)
            z.extract(name, path.parent)
            _report(target)


def build_celline_map(force: bool = False) -> None:
    """Parse cell line -> cancer type out of the HPA data-access help page."""
    out = RAW / "celline_cancer_map.tsv"
    if out.exists() and not force:
        print(f"  have {out.name}")
        return
    import re

    print("  fetching cell-line / cancer-type mapping ...", flush=True)
    with urllib.request.urlopen(HPA_HELP_URL) as r:
        html = r.read().decode("utf-8", errors="replace")
    text = re.sub(r"<[^>]*>", "\n", html)
    lines = [l.strip() for l in text.splitlines()]
    lines = [l for l in lines if l]

    rows, pending = [], None
    for l in lines:
        m = re.fullmatch(r"Cell line RNA - (.+) \[nTPM\]", l)
        if m:
            pending = m.group(1)
            continue
        if pending and l.startswith("cell_RNA_"):
            rows.append((pending, l[len("cell_RNA_"):]))
        pending = None

    seen, uniq = set(), []
    for line, ctype in rows:
        if line not in seen:
            seen.add(line)
            uniq.append((line, ctype))
    out.write_text("\n".join(f"{a}\t{b}" for a, b in uniq), encoding="utf-8")
    print(f"    {out.name}: {len(uniq)} cell lines")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="list what is present, download nothing")
    ap.add_argument("--force", action="store_true", help="re-download even if present")
    ap.add_argument("--skip-scrna", action="store_true", help="HPA tables only (~1.4 GB less)")
    args = ap.parse_args()

    RAW.mkdir(parents=True, exist_ok=True)
    SCRNA.mkdir(parents=True, exist_ok=True)

    if args.check:
        total = 0
        for d in (RAW, SCRNA):
            for p in sorted(d.glob("*")):
                if p.is_file():
                    total += p.stat().st_size
                    _report(p)
        print(f"  total: {total / 1e9:.2f} GB")
        return 0

    print("Human Protein Atlas tables:")
    for name, url in HPA_FILES.items():
        download(url, RAW / name, args.force)
        unzip(RAW / name)

    print("Cell-line annotation:")
    build_celline_map(args.force)

    if not args.skip_scrna:
        print("CELLxGENE tumour scRNA-seq:")
        for name, url in SCRNA_FILES.items():
            download(url, SCRNA / name, args.force)

    print("\nAll inputs present. Next: python run.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
