# TCR Reconstruction Benchmark — minimal BC09 reproduction

## Project goal

This project is a deliberately small, inspectable reproduction of *Evaluation of T Cell Receptor Construction Methods from scRNA-Seq Data*. It follows one human 10x 5′ sample (`BC09_TUMOR1`) through the biological evidence chain:

```text
ordinary 10x 5′ scRNA-seq
        ↓
candidate reads overlapping TCR loci / TCR reference segments
        ↓
cell-aware assembly and receptor annotation
        ↓
V gene + J gene + CDR3 amino-acid sequence
        ↓
comparison with paired scTCR-seq ground truth
        ↓
cell-level accuracy and sensitivity
```

The first phase covers only TRUST4, MiXCR, and YASIM-scTCR. It does **not** expand to the other five paper methods, other datasets, mouse, 3′ libraries, pseudo-bulk, or multiple read-length experiments.

The objective is to reproduce the paper's trends and understand the mechanism. Current TRUST4/MiXCR versions and native single-cell modes may not reproduce the paper's historical numbers exactly. Parameters must be recorded; they must not be silently tuned to force agreement.

## Data

| GEO sample | Role | Expected public relation |
|---|---|---|
| `GSM3148575` / `BC09_TUMOR1` | ordinary 10x 5′ GEX reads used as reconstruction input | NCBI SRA experiment resolved at run time |
| `GSM3148580` / `BC09_TUMOR1_TCR` | paired scTCR-seq reference | official `filtered_contig_annotations.csv.gz` |

`scripts/download/resolve_accessions.py` reads the current GEO SOFT record, extracts its SRX relation, queries NCBI RunInfo, and writes exact SRRs to `data/metadata/*_accessions.tsv`. SRRs are never guessed.

The first ground truth is the official processed file:

```text
data/ground_truth/raw/GSM3148580_BC09_TUMOR1_filtered_contig_annotations.csv.gz
```

It is preserved unchanged. `build_ground_truth.py` keeps rows where `productive == True`, `full_length == True`, and `chain` is `TRA` or `TRB`, then writes the benchmark-layer table:

```text
barcode  chain  v_gene  j_gene  cdr3aa
```

Barcodes lose only a terminal 10x gem-group suffix (`-1`); V/J allele suffixes such as `*01` are removed only in processed/evaluation tables.

The resolved GEX experiment contains four runs: `SRR7191896` through `SRR7191899`. The
validated full download contains all eight paired FASTQs. Each run has exactly matching mate
counts; all R1 reads are 27 bp and all R2 reads are 99 bp, confirming the preview-derived read
roles on the complete dataset.

## Real-data pipeline

```text
GEX FASTQ
   |
   +---- TRUST4 ----+
   |                |
   +---- MiXCR -----+
                    |
                    v
              predictions
                    |
                    v
scTCR truth ---- evaluation
```

The real-data runners use native barcode-aware single-cell modes. They do not split the GEX library into thousands of per-cell FASTQs.

### Cell-level metric definitions

For one `chain + metric`:

- **All cells**: cells with at least one valid ground-truth receptor for that chain.
- **Result cells**: cells with at least one predicted receptor for that chain.
- **True cells**: result cells with at least one prediction matching at least one ground-truth receptor for the current metric.
- **Sensitivity**: `True cells / All cells`.
- **Accuracy**: `True cells / Result cells`.

Metrics are exact after normalization:

- `CDR3`: CDR3 amino-acid sequence matches.
- `V`: allele-stripped V gene matches.
- `J`: allele-stripped J gene matches.
- `AsTCR`: CDR3 **and** V **and** J match in the same predicted/truth receptor pair.

Both truth and prediction are stored as sets of receptors per `barcode + chain`; dual-alpha and other multi-chain cells are never collapsed to one record.

## Environment

The base environment is not modified. The declared isolated environment is `tcrbench`:

```bash
./scripts/setup/check_environment.sh
./scripts/setup/create_environment.sh
source scripts/setup/activate_tcrbench.sh
./scripts/setup/check_environment.sh
```

`envs/environment.yml` requests Python 3.9, common scientific packages, SRA tools, seqkit, pigz, samtools, TRUST4, ART, and YASIM-scTCR 1.0.1. The authors' README advertises PyPI version 1.0.1, but PyPI currently exposes only 0.1.0 and 1.0.0; the yml therefore pins the authors' immutable Git tag `1.0.1` and records this documentation/package-index difference.

The SRA Toolkit is constrained to `>=3.2`. The unconstrained initial solve selected 2.9.6,
whose bundled TLS stack rejected NCBI's current certificate; that failure is preserved in
`logs/preview_bc09_gex.log`.

MiXCR is intentionally separate because current distributions may require an academic license. No license bypass is attempted. Until both `mixcr -v` and `mixcr listPresets` succeed, MiXCR stages are marked:

```text
WAITING_FOR_MIXCR_LICENSE
```

See `scripts/setup/install_mixcr.sh`. It installs pinned MiXCR 4.7.0 from the official
`milaboratories` Conda channel into the isolated `tcrbench` environment and does not handle the
license itself. License activation follows MiLaboratories' official `mixcr activate-license`
workflow. The real-data runner verifies that the installed version actually provides
`10x-sc-5gex` before analysis.

Current validated server state: MiXCR 4.7.0 is installed, its license was activated through the
official command, and `10x-sc-5gex` is present. The credential is stored by MiXCR outside this
project and is never written to project logs or scripts.

> Current software versions may not produce exactly the paper's historical values. This phase aims to reproduce the trend, not force identical numbers.

## Reproduce the real-data minimum

All commands resolve the project root from their own location. Analysis scripts contain no server-specific absolute paths.

```bash
# 1. Inspect software only; installs and downloads nothing.
./run_pipeline.sh setup

# 2. Resolve exact accessions, download the small processed VDJ CSV, and build truth.
./run_pipeline.sh ground_truth

# 3. LARGE DATA — run only after reviewing metadata and the disk estimate.
./run_pipeline.sh download 16

# 4. Reconstruct and normalize predictions.
./run_pipeline.sh trust4 16
./run_pipeline.sh mixcr

# 5. Compute cell-level metrics and representative cells.
./run_pipeline.sh evaluate

# 6. Plot the available Figure-2-style benchmark.
./run_pipeline.sh plot
```

Calling `./run_pipeline.sh` with no stage only prints help. It never starts the SRA download.

Before the full download, a bounded preview can verify the actual read structure from the first
resolved SRR. It downloads at most 100,000 spots and defaults to 10,000:

```bash
./run_pipeline.sh preview 10000
```

The preview writes `BC09_fastq_preview_stats.tsv` plus a content-level structure inference. It
requires exact R1/R2 identifier agreement and uses observed lengths before assigning read roles.
When the ground-truth table exists, it also tests every possible 16 bp R1 window against known
BC09 T-cell barcodes; this verifies barcode coordinates independently of the FASTQ filename.

For the first 10,000 spots of `SRR7191896`, the measured layout is R1=27 bp and R2=99 bp.
Bases 1–16 of R1 match 7,694 reads to a BC09 ground-truth cell barcode, while the second-best
16 bp window matches only 7. Bases 17–26 are therefore treated as the 10 bp UMI; base 27 is a
non-UMI tail (91.94% T in the preview) and is intentionally excluded by TRUST4's
`bc:0:15,um:16:25` extraction coordinates.

### Controlled FASTQ download

`download_bc09_gex.sh`:

1. reads SRRs only from `GSM3148575_accessions.tsv`;
2. estimates required free space from NCBI `size_MB`, an 8× expansion factor, and a 20 GiB safety margin;
3. stops before download if free space is below that estimate;
4. uses `prefetch` and `fasterq-dump --split-files`;
5. compresses with pigz when available;
6. skips complete read pairs and refuses to overwrite partial output;
7. logs stdout/stderr.

For a multi-hour server download that must survive an SSH/Codex disconnect, source the environment
and use the explicit background wrapper after reviewing the same space gate:

```bash
./scripts/download/start_bc09_gex_background.sh 16
./scripts/download/status_bc09_gex.sh
```

It refuses to create a duplicate downloader and refuses to remove SRA lock files automatically.

During the recorded BC09 run, the server-to-NCBI connection for the last two accessions fell to
roughly 37--52 kB/s. As a documented recovery, the exact official NCBI S3 SRA objects were
downloaded on the client, checked against the server-reported `Content-Length`, uploaded under
temporary names, and accepted only after SRA Toolkit 3.4.1 reported metadata/READ/QUALITY MD5
checks as valid and each database as consistent. The incomplete server download was moved to a
timestamped archive rather than deleted. Transfer SHA-256 values and the archive path are recorded
in `logs/local_sra_transfer_manifest.tsv`; validation output is in
`logs/vdb_validate_SRR7191898_active.log`, `logs/vdb_validate_SRR7191899_active.log`, and the
corresponding `logs/SRR719189*.vdb_validate.log` transfer checks. After activation,
the normal `download_bc09_gex.sh` workflow found both SRA objects locally and performed
`fasterq-dump`, pigz compression, and the same full FASTQ QC. This recovery changes transport only,
not the SRA content or downstream analysis.

`check_fastq.sh` measures every file with `seqkit stats`, compares mate counts, and infers barcode/UMI versus cDNA roles from observed read lengths. TRUST4 refuses to run until that QC supports the 16 bp barcode + 10 bp UMI structure.

### TRUST4

The runner searches the active environment and `TRUST4_REFERENCE_DIR` for `human_IMGT+C.fa`; it does not hard-code a reference directory. For raw FASTQ, TRUST4's official documentation permits the IMGT sequence reference for both `-f` and `--ref`. The exact shell-escaped command is saved in `logs/trust4_bc09_command.txt`.

TRUST4 barcode AIRR/report output is parsed to:

```text
data/processed/BC09_trust4_prediction.tsv
```

Primary and secondary TRA/TRB chain fields are retained.

The runner also accepts `BC09_FASTQ_DIR`, `BC09_FASTQ_QC`, and `BC09_TRUST4_*`
output overrides. These are used for a bounded preview smoke test without mixing preview artifacts
with the full BC09 result; default paths remain the documented real-data paths.

### MiXCR

The runner verifies `10x-sc-5gex`, then uses:

```text
mixcr analyze 10x-sc-5gex --species hsa ...
```

Input-lane symlinks avoid concatenating or duplicating large FASTQs. A cell-split export explicitly requests cell ID, V gene, J gene, and CDR3 amino-acid sequence before parsing to:

```text
data/processed/BC09_mixcr_prediction.tsv
```

## Simulation

```text
YASIM-scTCR
    ↓
one shared 500-cell TCR truth
    ↓
2×   10×   50×   100×  (150 bp, otherwise fixed)
    ↓
ART 5′-oriented single-end TCR evidence
    ↓
deterministic 10x-like R1 barcode/UMI + R2 cDNA wrapper
    ↓
TRUST4 / MiXCR native single-cell reconstruction
    ↓
the same evaluate.py
    ↓
depth-performance curves
```

Only small receptor tables (`tcr_cache`, usage bias, insertion/deletion tables) are downloaded from the authors' Zenodo record. The example scRNA count matrix and transcriptome are not downloaded in this minimum.

YASIM-scTCR 1.0.1 uses `random.SystemRandom` internally and exposes no random-seed argument. This is recorded explicitly rather than claiming a false seed. Clean depth comparison is achieved by generating the receptor repertoire once and reusing the exact same per-cell TCR FASTA and truth table at all four depths.

ART is separately reproducible: `03_simulate_reads.sh` records base seed `20260829` and derives a
stable independent seed for every receptor from `SHA256(base_seed:receptor_basename)`. The full map
is saved in `logs/yasim_art_receptor_seeds.tsv`. Reusing the same receptor-specific random stream at
all depths avoids correlated fragment starts between different receptors while allowing higher
depth to extend the lower-depth evidence. The base can be overridden only explicitly with
`TCRBENCH_ART_SEED`; YASIM's unseedable repertoire generation remains a one-time, shared upstream
truth rather than being regenerated for each depth.

A repeated 2× run produced the identical multiset of 2,000 FASTQ records (identifier, sequence,
quality, and multiplicity), but not a byte-identical file because 16-thread YASIM assembly changes
record order. This content-level check is recorded in `logs/yasim_reproducibility_check.tsv`; no
claim of deterministic file ordering is made.

The tagged 1.0.1 wheel also contains two imports of `src.yasim_sctcr` and imports an absent `get_sample_data_path` symbol from `rearrange_tcr`. A normal installed wheel therefore cannot even show that subcommand's help. Simulation wrappers add the read-only author checkout root to `PYTHONPATH`; `yasim_sctcr_compat.py` injects only the absent, unused symbol before delegating to the original frontend. No upstream or site-package file is edited. Raw failures and compatibility help are both preserved in `logs/yasim_cli_help.txt`; remove the shim only after an upstream release fixes these imports.

The ART/YASIM output is single-end cDNA with the cell barcode encoded in the sequence identifier. `make_10x_fastqs.py` converts it to native-tool-compatible pairs without altering R2: R1 is the known 16 bp cell barcode plus a deterministic 10 bp synthetic UMI. This adapter is a benchmark interface, not a claim that YASIM natively emits raw 10x FASTQs.

Because these 500 deterministic simulation barcodes are synthetic, they are absent from MiXCR's
default 10x `737K-august-2016` whitelist. The simulation runner therefore supplies the exact shared
`simulation/truth/barcodes.txt` through MiXCR 4.7.0's documented
`--set-whitelist CELL=file:...` option. The real BC09 runner does not override the native 10x
whitelist.

For the one combined TCR depth TSV, the installed YASIM 3.2.1 interface is `python -m yasim art`. Its similarly named `python -m yasim_sc art` interface expects a **directory of per-cell depth files** and otherwise exits successfully after simulating zero records. Both help texts are captured so this easily missed API distinction remains visible.

```bash
./run_pipeline.sh simulate 16
python scripts/plotting/plot_depth_benchmark.py
```

Before simulation, `inspect_yasim_cli.sh` records top-level and subcommand help in `logs/yasim_cli_help.txt`. If README examples and installed CLI differ, the installed help controls.

## Validation gates

- Ground truth: reports unique cells, TRA/TRB receptors, cells by chain, and per-cell-chain receptor multiplicity.
- FASTQ: reports read count and length and aborts on mate-count mismatch.
- Simulated FASTQ: requires exactly `1000 × depth` read pairs at each depth, with fixed
  R1=26 bp and R2=150 bp; validation is saved in `simulation_fastq_stats.tsv`.
- TRUST4/MiXCR: reports cells with reconstructed TRA and TRB.
- Evaluation: raises immediately unless rates are in `[0,1]` and `true_cells <= all_cells` and `true_cells <= result_cells`.
- `python -m unittest discover -s tests -v`: tests prediction-only denominator handling,
  dual-alpha matching, zero-result behavior, and the requirement that AsTCR fields come from one
  matching receptor pair.
- Plotting consumes metrics TSVs only; figures never recompute biological matching.

## Outputs

```text
data/ground_truth/BC09_ground_truth.tsv
data/processed/BC09_{trust4,mixcr}_prediction.tsv
results/benchmark/BC09_{TRUST4,MiXCR}_metrics.tsv
results/benchmark/BC09_{TRUST4,MiXCR}_cell_details.tsv
results/benchmark/example_cells.tsv
results/figures/BC09_realdata_benchmark.{png,pdf}
results/benchmark/simulation_metrics.tsv
results/figures/simulation_depth_{sensitivity,accuracy}.pdf
```

## Validated minimum results

The completed real-data AsTCR benchmark is:

| method | chain | all cells | result cells | true cells | sensitivity | accuracy |
|---|---:|---:|---:|---:|---:|---:|
| TRUST4 | TRA | 5,849 | 1,127 | 701 | 0.1199 | 0.6220 |
| TRUST4 | TRB | 6,402 | 2,313 | 1,717 | 0.2682 | 0.7423 |
| MiXCR | TRA | 5,849 | 847 | 548 | 0.0937 | 0.6470 |
| MiXCR | TRB | 6,402 | 1,496 | 1,222 | 0.1909 | 0.8168 |

For this sample and these current native presets, TRUST4 returns receptors for more cells and has
higher sensitivity, while MiXCR's smaller result set has higher accuracy. Both tools reconstruct
TRB more readily than TRA. These are observed version-specific results, not parameters tuned to
match the historical paper exactly; CDR3-, V-, and J-specific values remain in the metrics TSVs.

The shared-repertoire simulation reproduces the core depth mechanism. At 2x, missing evidence
reduces result-cell coverage, especially for MiXCR; from 10x onward the curves largely plateau.
For AsTCR, TRUST4 sensitivity changes from 0.806/0.908 (TRA/TRB) at 2x to 0.820/0.940 at 10x,
while MiXCR changes from 0.778/0.864 to 0.836/0.946. The plateau is expected in this deliberately
simple 500-cell TCR-only minimum: once every reconstructable receptor has adequate spanning
evidence, additional reads do not fix method/reference/assembly limitations.

## Upstream repositories

These directories are read-only inputs and must not be edited:

```text
repo/Benchmarking_TCR_Construction/
repo/yasim-sctcr/
repo/TRUST4/  # official tool/reference checkout
```

All project-specific code lives under `scripts/`, `tests/`, and `notebooks/`.
