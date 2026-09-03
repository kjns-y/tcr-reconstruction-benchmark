#!/usr/bin/env python3
"""Shared, non-destructive normalization for benchmark-layer receptor tables."""

from __future__ import annotations

import re
from typing import Any

MISSING = {"", "nan", "none", "null", "na", "n/a", "-"}


def normalize_text(value: Any) -> str:
    """Return stripped text, mapping common missing values to an empty string."""
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in MISSING else text


def normalize_barcode(value: Any) -> str:
    """Remove a terminal 10x gem-group suffix, e.g. AAAC-1 -> AAAC."""
    return re.sub(r"-[0-9]+$", "", normalize_text(value))


def normalize_gene(value: Any) -> str:
    """Remove allele suffixes and normalize case without changing gene family."""
    text = normalize_text(value).upper()
    if not text:
        return ""
    return re.sub(r"\*[0-9A-Z._-]+$", "", text)


def normalize_cdr3aa(value: Any) -> str:
    """Normalize amino-acid CDR3 strings for exact comparison."""
    return re.sub(r"\s+", "", normalize_text(value)).upper()


def infer_chain(*values: Any) -> str:
    """Infer TRA/TRB from one or more gene/chain labels."""
    text = " ".join(normalize_text(value).upper() for value in values)
    if re.search(r"\bTRA|TCRA|ALPHA", text):
        return "TRA"
    if re.search(r"\bTRB|TCRB|BETA", text):
        return "TRB"
    return ""

