# The manuscript's figures, benchmark and corpus validation

Everything the paper shows is regenerated from this folder.

- `figures/figN/` — one folder per figure: `panel_*.py` draws one panel each from `data/*.csv`,
  `merge.py` assembles the figure and writes `figN.svg`, `figN.pdf` and `figN.png`; `style.py` and
  `schematic_lib.py` (one level up) hold the shared colour code and the schematic drawing primitives.
  Figures 1 and S1 are drawn from `schematic.py` / `schematic_full.py` and carry no data.
  `make_data.py` in each folder is how `data/*.csv` was produced from the pipeline's run outputs;
  its `RUNS` paths are the authors' run folders and must be pointed at your own runs to re-derive
  the CSVs. The CSVs themselves are shipped, so `python figures/figN/merge.py` reproduces every
  figure without any run.
- `benchmark/` — the NIST SRM 1950 benchmark: `build_reference.py` harmonises the published lists
  (Bowden 2017 consensus as distributed with LipidQC, Quehenberger 2010, Godzien 2024) into
  `reference/srm1950_reference_lists.csv`; `score_nist.py` scores a run's `Results/` folder against
  it; `results/` holds the scores reported in the paper.
- `validation/` — the fifteen-study public corpus: how each deposit was screened, staged and scored
  (`SCOPE.md`, `deposits.py`, `score_agreement.py`, `results_table.py`), the per-study outcomes
  (`RESULTS_TABLE.md`, `validation_table.csv`) and the metadata used (`metadata/`). The staging
  paths inside the scripts are the authors' download area; the deposit accessions are in the tables.

Requirements beyond the pipeline: `pip install -e ".[paper]"` (matplotlib, scipy, cairosvg).
