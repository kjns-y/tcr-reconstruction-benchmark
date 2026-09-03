#!/usr/bin/env python3
"""Compatibility launcher for two packaging defects in YASIM-scTCR 1.0.1.

The tagged source imports ``src.yasim_sctcr`` and imports an absent, unused
``get_sample_data_path`` symbol. The shell wrappers put the read-only checkout
on PYTHONPATH; this launcher supplies only the absent symbol before delegating
to the original module frontend. No upstream or site-package file is changed.
"""

from __future__ import annotations

import runpy
from pathlib import Path

import yasim_sctcr._main as yasim_main


def _unused_get_sample_data_path(name: str) -> str:
    raise RuntimeError(
        "YASIM-scTCR requested get_sample_data_path unexpectedly; the 1.0.1 tag "
        f"does not implement it (requested: {name})."
    )


if not hasattr(yasim_main, "get_sample_data_path"):
    yasim_main.get_sample_data_path = _unused_get_sample_data_path  # type: ignore[attr-defined]

runpy.run_module("yasim_sctcr", run_name="__main__")

