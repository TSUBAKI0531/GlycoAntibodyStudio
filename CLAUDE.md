# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install dependencies
pip install -r requirements.txt
# Optional heavy deps (conda recommended):
# conda install -c conda-forge openmm pdbfixer
# conda install -c bioconda hmmer muscle
# pip install anarci @ git+https://github.com/oxpig/ANARCI.git

# Run the app
streamlit run app.py

# Run with custom port
streamlit run app.py --server.port 8502
```

There are no tests or linting configs in this repository.

## Architecture

**Entry point:** `app.py` — Streamlit UI only. No business logic lives here. It drives three tabs (Analytics, 3D Viewer, Report) and calls into `modules/`.

**Core pipeline (`modules/engine.py`):**  
`GlycoEngine` is the central orchestrator. Instantiate it with a target SMILES; it pre-generates the glycan conformer ensemble at init time. Call `run_pipeline(cdr3_seq)` per CDR3 candidate. The pipeline runs four steps in sequence:

1. **CDR Grafting** — replaces CDR-H3 in the Trastuzumab VH scaffold (PDB: 1N8Z) using ANARCI IMGT numbering, with a position-based fallback when ANARCI is absent.
2. **Structure building** — uses OpenMM + PDBFixer + AMBER14 force field. Falls back to a placeholder helical PDB when OpenMM is not installed.
3. **Hydrated ensemble docking** — runs AutoDock Vina against each glycan conformer, including conserved waters. Falls back to `_simulate_docking()` (deterministic RNG, physically plausible values) when Vina is not in PATH.
4. **Scoring calibration** — H-bond weighting + desolvation penalty correction applied to raw Vina scores.

Results are `PipelineResult` dataclasses converted to dicts for DataFrame assembly; Pareto front identification runs in `app.py` via `_identify_pareto_front()`.

**`your_modules.py`** — two independent classes:
- `ComplexBuilder`: generates α-helix-approximated PDB strings from sequences for visualization (not physically accurate).
- `AntibodyGraftingEngine`: assembles full Fv sequences by grafting CDR sets onto Trastuzumab FR regions. `predict_cdrs()` is currently template-based (default / tn_antigen / sialyl_tn).

**`modules/utils.py`** — stateless utilities: SMILES validation (RDKit → regex fallback), CSV parsing with fuzzy column matching, CDR annotation via ANARCI, glycan conformer generation via RDKit ETKDGv3, Kyte-Doolittle hydrophobicity (BioPython → manual fallback). Exports `TRASTUZUMAB_VH` and `TRASTUZUMAB_VL` reference sequences.

**`modules/reports.py`** — `generate_pareto_chart()` (Plotly, interactive), `generate_pareto_chart_matplotlib()` (PNG for PDF embedding), `RankingReport.generate_pdf()` (FPDF2 → text fallback).

**`modules/visualization.py`** — `render_3d_structure()` with three-tier fallback: stmol → py3Dmol HTML embed → PDB text. CDR regions are highlighted by IMGT position ranges.

## Graceful degradation pattern

Every heavy optional dependency (OpenMM, ANARCI, RDKit, AutoDock Vina, stmol, fpdf2) is wrapped in `try/except ImportError` with a documented fallback. This is intentional — wet-lab users should be able to run the app even with a minimal install. Do not remove these fallbacks when editing the pipeline steps.

## Key data flow

```
target_smiles + CSV upload
  → GlycoEngine.__init__()  # conformer ensemble generated once
  → GlycoEngine.run_pipeline(cdr3_seq)  # per row in CSV
      → DockingResult (per conformer)
      → PipelineResult.to_dict()
  → pd.DataFrame  (st.session_state.results_df)
  → _identify_pareto_front()  # added as is_pareto column
  → Tab rendering + RankingReport.generate_pdf()
```

## Notes

- `EngineConfig` controls all numeric pipeline parameters; pass it to `GlycoEngine` to override defaults without touching call sites.
- `DockingResult` and `PipelineResult` are `@dataclass`; `CDRSet` validates amino acid characters on construction.
- `sample_candidates.csv` contains 10 test CDR3 sequences for manual testing.
- The `venv/` directory is git-ignored; do not commit it.
