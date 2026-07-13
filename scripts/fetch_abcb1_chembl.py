"""Annotate table_drugs.csv with ChEMBL ABCB1 (P-gp) interaction data."""

from __future__ import annotations

import re
import time
from pathlib import Path

import pandas as pd
import requests

ABCB1_TARGET = "CHEMBL4302"
API_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
CSV_PATH = Path(__file__).resolve().parents[1] / "table_drugs.csv"
PAGE_SIZE = 1000


def classify_assay(description: str) -> str:
    text = (description or "").lower()
    if "substrate" in text:
        return "substrate"
    if any(word in text for word in ("inhibitor", "inhibition", "inhibit")):
        return "inhibitor"
    if "modulat" in text:
        return "modulator"
    return "interaction"


PRIORITY = {"inhibitor": 0, "substrate": 1, "modulator": 2, "interaction": 3}


def fetch_abcb1_annotations() -> dict[str, set[str]]:
    annotations: dict[str, set[str]] = {}
    offset = 0
    total = None

    while total is None or offset < total:
        params = {
            "target_chembl_id": ABCB1_TARGET,
            "limit": PAGE_SIZE,
            "offset": offset,
        }
        response = requests.get(API_URL, params=params, timeout=60)
        response.raise_for_status()
        payload = response.json()
        total = payload["page_meta"]["total_count"]

        for activity in payload["activities"]:
            mol_id = activity.get("molecule_chembl_id")
            if not mol_id:
                continue
            label = classify_assay(activity.get("assay_description", ""))
            annotations.setdefault(mol_id, set()).add(label)

        offset += PAGE_SIZE
        print(f"Fetched {min(offset, total)}/{total} ABCB1 activities")
        time.sleep(0.2)

    return annotations


def format_annotation(labels: set[str]) -> str:
    if not labels:
        return ""
    ordered = sorted(labels, key=lambda x: PRIORITY.get(x, 99))
    return "; ".join(ordered)


def main() -> None:
    print("Downloading ChEMBL ABCB1 annotations...")
    annotations = fetch_abcb1_annotations()
    print(f"Unique compounds with ABCB1 data in ChEMBL: {len(annotations)}")

    df = pd.read_csv(CSV_PATH, encoding="utf-8-sig")
    chembl_col = "ChEMBL ID"

    def lookup(chembl_id: object) -> str:
        if pd.isna(chembl_id):
            return ""
        key = str(chembl_id).strip()
        return format_annotation(annotations.get(key, set()))

    df["ABCB1_Pgp_CHEMBL"] = df[chembl_col].map(lookup)
    matched = (df["ABCB1_Pgp_CHEMBL"] != "").sum()
    print(f"Matched {matched}/{len(df)} drugs in local dataset")

    df.to_csv(CSV_PATH, index=False, encoding="utf-8-sig")
    print(f"Updated {CSV_PATH}")


if __name__ == "__main__":
    main()