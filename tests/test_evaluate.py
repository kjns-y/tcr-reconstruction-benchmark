#!/usr/bin/env python3
"""Small executable smoke test for multi-chain, cell-level metric behavior."""

from __future__ import annotations

import sys
import unittest
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/evaluation"))
from evaluate import Receptor, evaluate_maps  # noqa: E402
from normalization import normalize_barcode, normalize_cdr3aa, normalize_gene  # noqa: E402


def receptor_map():
    return defaultdict(lambda: defaultdict(set))


class EvaluationTests(unittest.TestCase):
    """Protect the paper-specific denominator and receptor-pair semantics."""

    def test_cell_denominators_and_multiple_receptors(self) -> None:
        truth = receptor_map()
        pred = receptor_map()
        truth["CELL1"]["TRA"].update(
            {Receptor("TRAV1", "TRAJ1", "CAVA"), Receptor("TRAV2", "TRAJ2", "CAVB")}
        )
        truth["CELL1"]["TRB"].add(Receptor("TRBV1", "TRBJ1", "CASA"))
        truth["CELL2"]["TRA"].add(Receptor("TRAV3", "TRAJ3", "CAVC"))
        pred["CELL1"]["TRA"].add(Receptor("TRAV2", "TRAJ2", "CAVB"))
        pred["CELL1"]["TRB"].add(Receptor("TRBV9", "TRBJ9", "WRONG"))
        # Prediction-only CELL3 belongs in Result cells but not All/True cells.
        pred["CELL3"]["TRA"].add(Receptor("TRAV3", "TRAJ3", "CAVC"))
        metrics, details = evaluate_maps(truth, pred, "TEST")
        row = metrics.query("chain == 'TRA' and metric == 'AsTCR'").iloc[0]
        self.assertEqual(int(row.all_cells), 2)
        self.assertEqual(int(row.result_cells), 2)
        self.assertEqual(int(row.true_cells), 1)
        self.assertEqual(float(row.sensitivity), 0.5)
        self.assertEqual(float(row.accuracy), 0.5)
        self.assertEqual(len(details), 6)

    def test_astcr_requires_one_matching_receptor_pair(self) -> None:
        truth = receptor_map()
        pred = receptor_map()
        truth["CELL1"]["TRA"].update(
            {Receptor("TRAV1", "TRAJ1", "CAVA"), Receptor("TRAV2", "TRAJ2", "CAVB")}
        )
        # Each individual field occurs in truth, but this recombinant tuple does not.
        pred["CELL1"]["TRA"].add(Receptor("TRAV1", "TRAJ2", "CAVB"))
        metrics, _ = evaluate_maps(truth, pred, "TEST")
        by_metric = metrics.query("chain == 'TRA'").set_index("metric")
        self.assertEqual(int(by_metric.loc["CDR3", "true_cells"]), 1)
        self.assertEqual(int(by_metric.loc["V", "true_cells"]), 1)
        self.assertEqual(int(by_metric.loc["J", "true_cells"]), 1)
        self.assertEqual(int(by_metric.loc["AsTCR", "true_cells"]), 0)

    def test_zero_result_cells_have_defined_zero_accuracy(self) -> None:
        truth = receptor_map()
        truth["CELL1"]["TRB"].add(Receptor("TRBV1", "TRBJ1", "CASA"))
        metrics, _ = evaluate_maps(truth, receptor_map(), "TEST")
        row = metrics.query("chain == 'TRB' and metric == 'AsTCR'").iloc[0]
        self.assertEqual(int(row.all_cells), 1)
        self.assertEqual(int(row.result_cells), 0)
        self.assertEqual(float(row.sensitivity), 0.0)
        self.assertEqual(float(row.accuracy), 0.0)

    def test_benchmark_layer_normalization(self) -> None:
        self.assertEqual(normalize_barcode("AAACCTGAGCAGACTG-1"), "AAACCTGAGCAGACTG")
        self.assertEqual(normalize_barcode("AAACCTGAGCAGACTG-12"), "AAACCTGAGCAGACTG")
        self.assertEqual(normalize_gene("trbv7-9*01"), "TRBV7-9")
        self.assertEqual(normalize_gene("TRBJ2-3*02"), "TRBJ2-3")
        self.assertEqual(normalize_cdr3aa(" cassl rayst " ), "CASSLRAYST")


if __name__ == "__main__":
    unittest.main(verbosity=2)
