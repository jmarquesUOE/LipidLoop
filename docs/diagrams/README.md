# Pipeline diagram — `.raw` to the delivered CSVs

**`pipeline_raw_to_csv.html` needs nothing installed** — double-click it and any browser opens it.
Self-contained: inline SVG, no CDN, no network. This is the one to use.

`pipeline_raw_to_csv.drawio` is the same diagram as an editable source. Draw.io does not need
installing either — <https://app.diagrams.net> opens it in a browser via File ▸ Open — but use it
only if you want to rearrange the boxes. Uncompressed XML, so it diffs in git.

Both are generated: `make_pipeline_html.py` and `make_pipeline_diagram.py`. Regenerate rather than
hand-edit, or the two will drift apart.

## What it shows

Seven bands, in the order `pipeline.py::run` executes them:

| band | |
|---|---|
| 1 | conversion, and the freshness check that makes a re-run cost 40 min rather than 4 h |
| 2 | identification per injection — artefact screen, library search, mass offset, purity |
| 3 | feature detection, RT alignment, linking, noise floor |
| 4 | the peak finder: join, score gate, naming, retention model, the adduct/in-source sweeps |
| 5 | post-filters: split peaks, adduct pairs, blanks, fatty acids, retention extension, presence |
| 6 | the six output tables per polarity |
| 7 | cross-polarity assignment and the report — needs BOTH polarities |

**Yellow diamonds are decisions**, and they are the point of the diagram. Each one is a place where
a row's fate changes, and most of them cost something to get wrong:

- `purity ≥ 75 %` decides molecular versus sum composition — but it scores the *winner* and does
  not gate the name, which is why `Chain Evidence` exists.
- `|z| > 3σ from its class model` — and a class with too few members has **no** model, so its rows
  pass by never taking the test (5 % of positive rows, 8 % of negative).
- `detected in ≥ 70 % of any one group` — per group, never globally, so a species present in only
  one arm survives.
- `class seen in both polarities` — **case A stops the run** rather than guessing.

## What it deliberately does not show

Per-injection detail (each of 152 files runs bands 1–2 independently), the QC report's internal
sections, and the standards path — standards are a facility measurement and are kept out of every
delivered table.

Regenerate with the script that produced it if the pipeline changes; the band order is taken from
the numbered comments in `src/lipidloop/pipeline.py`.
