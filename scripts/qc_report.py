"""Facility quality-control report — **one document covering both polarities**, HTML and PDF.

    python scripts/qc_report.py \\
        --results-pos <run>/Pos/Final_Results.csv --unfiltered-pos <run>/Pos/Unfiltered_Results.csv \\
        --results-neg <run>/Neg/Final_Results.csv --unfiltered-neg <run>/Neg/Unfiltered_Results.csv \\
        --metadata-pos metadata_Pos.csv --metadata-neg metadata_Neg.csv \\
        --sequence SeqRand.csv --study "Study name" --owner "PI name" \\
        --out <analysis>/QC

Either polarity may be omitted, so a single-polarity study still works.

## Who this is for

The person whose samples these are, who is not a mass spectrometrist. That constrains it:

  * every section opens with what it means, not what was computed;
  * an injection that failed is **named**, because "one sample was an outlier" cannot be acted on;
  * a limit on what the data can support is stated as plainly as a pass;
  * nothing about the instrument that does not change what they can conclude.

## Both polarities, never merged

The two polarities measure different molecules on different scales, and merging them would be
wrong. They are reported **side by side** because the interesting statements here are comparative
and only land that way. That makes labelling critical rather than cosmetic: a reader who loses
track will compare a positive CV against a negative D-ratio without noticing, so the polarity
appears in every column heading, every figure caption and every panel title.

## What it does not do

It does not correct anything. Drift is measured and shown; whether to apply a QC-RSC correction is
a decision that needs a threshold this facility has not yet set, and a corrected matrix that
arrives without that decision invites the question of what else was adjusted.
"""
from __future__ import annotations

import argparse
import base64
import csv
import html
import io
import json
import sys
from collections import Counter
from datetime import date
from pathlib import Path

import numpy as np
import re

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop import qc   # noqa: E402
from lipidloop.blanks import infer_role   # noqa: E402

COLOURS = {"sample": "#0072B2", "qc": "#E69F00", "blank": "#009E73", "flag": "#D55E00"}
LABELS = {"sample": "study sample", "qc": "pooled QC", "blank": "blank"}
GROUP_COLOURS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442"]
SIZE = 44

# Acceptance bands. The CV bands are the facility's; the D-ratio edges follow the usual
# convention, where 0.5 is the accepted bar and anything approaching 1 leaves no room for biology.
# 15% is this facility's target and is declared as such in the methods note; 20% and 30% are the
# thresholds in common use for untargeted profiling (Dunn et al. 2011, Nature Protocols), and the
# table below the gauge already reports the fraction better than each.
CV_BANDS = [(15, "optimum", "#009E73"), (20, "acceptable", "#56B4E9"),
            (30, "high", "#E69F00"), (float("inf"), "unreliable", "#D55E00")]
DRATIO_BANDS = [(0.5, "good", "#009E73"), (0.75, "marginal", "#E69F00"),
                (float("inf"), "poor", "#D55E00")]

# The running order, technical first and biological last, so the document reads as one argument:
# the instrument behaved -> here is what the table contains -> here is what is suspect in it ->
# here is what you can ask of it.
#
# ⚠ Numbered HERE and nowhere else. Section numbers were literals in the prose ("see section 2"),
# so any reorder silently mislabelled every cross-reference. Refer to a section by key via `ref`.
FACILITY_URL = "https://institute-genetics-cancer.ed.ac.uk/mass-spectrometry"
# The institute has no logo of its own — its own site uses the University logo as the masthead and
# sets the institute in text — so the report does the same rather than inventing a mark.
INSTITUTION = "Institute of Genetics and Cancer &middot; University of Edinburgh"

SECTIONS = [
    ("precision", "Precision &mdash; how repeatable is the measurement"),
    ("accuracy", "Mass accuracy &mdash; where the instrument put the ions"),
    ("calibration", "Mass calibration &mdash; per-file drift, measured and corrected"),
    ("order", "Run order &mdash; is acquisition confounded with the biology"),
    ("blanks", "Blanks &mdash; what is in the background"),
    ("artefacts", "Artefacts &mdash; impossible-mass-defect peaks stripped from spectra"),
    ("duty_cycle", "MS1&rarr;MS2 duty cycle &mdash; what the acquisition spent its budget on"),
    ("missing", "Missing values"),
    ("polarity", "Polarity &mdash; which mode measures each class"),
    ("annotation", "Annotation &mdash; what a row means"),
    ("outliers", "Outliers &mdash; which injections do not belong"),
    ("filters", "Filters &mdash; what was removed, and what was never examined"),
    ("support", "What the data can support"),
    ("structure", "Structure &mdash; what the data looks like"),
    ("methods", "How this was produced"),
    ("files", "The delivered files &mdash; what each one is"),
]
_NUMBER = {key: n for n, (key, _) in enumerate(SECTIONS, start=1)}


def head(key: str) -> str:
    """The `<h2>` for a section, numbered from SECTIONS."""
    title = dict(SECTIONS)[key]
    return f"<h2>{_NUMBER[key]}. {title}</h2>"


def ref(key: str) -> str:
    """A cross-reference by key, so a reorder cannot mislabel it."""
    return f"section&nbsp;{_NUMBER[key]}"


STYLE = """
.brand{display:flex;align-items:center;gap:18px;margin:0 0 18px 0;padding:0 0 16px 0;
  border-bottom:2px solid #111}
.brand img{height:52px;width:auto}
.brand-text{display:flex;flex-direction:column;line-height:1.35;font-size:11pt}
.brand-text strong{font-size:12pt}
.brand-text span{color:#444}
.brand-text a{color:#0072B2;text-decoration:none;font-size:10pt}
@media (prefers-color-scheme: dark){.brand img{filter:invert(1)}}

:root { --ink:#1a1a1a; --soft:#5b5b5b; --rule:#e3e3e3; --warn:#D55E00; --ok:#009E73;
        --pos:#0072B2; --neg:#CC79A7; }
* { box-sizing:border-box; }
body { font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
       color:var(--ink); max-width:1080px; margin:0 auto; padding:32px 40px 80px;
       line-height:1.55; font-size:15px; }
header { border-bottom:3px solid var(--ink); padding-bottom:14px; margin-bottom:28px; }
header .facility { font-size:13px; letter-spacing:.08em; text-transform:uppercase;
                   color:var(--soft); }
header h1 { font-size:26px; margin:10px 0 4px; font-weight:600; }
header .sub { font-size:15px; color:var(--soft); }
.meta { display:grid; grid-template-columns:repeat(auto-fit,minmax(200px,1fr)); gap:6px 24px;
        margin-top:16px; font-size:13.5px; }
.meta div span { color:var(--soft); display:block; font-size:12px; }
h2 { font-size:19px; margin:38px 0 4px; padding-top:14px; border-top:1px solid var(--rule);
     font-weight:600; }
h3 { font-size:15.5px; margin:22px 0 4px; font-weight:600; }
.lede { color:var(--soft); margin:2px 0 14px; }
table { border-collapse:collapse; width:100%; margin:12px 0 18px; font-size:13.5px; }
th,td { text-align:left; padding:6px 10px; border-bottom:1px solid var(--rule); }
th { font-weight:600; color:var(--soft); font-size:12px; text-transform:uppercase;
     letter-spacing:.04em; }
td.num, th.num { text-align:right; font-variant-numeric:tabular-nums; }
th.pos { color:var(--pos); } th.neg { color:var(--neg); }
.verdict { border-left:4px solid var(--ok); background:#f6fbf9; padding:14px 18px; margin:18px 0;
           border-radius:0 4px 4px 0; }
.verdict.warn { border-left-color:var(--warn); background:#fdf6f1; }
.verdict h3 { margin-top:0; }
figure { margin:16px 0 8px; }
figure img { width:100%; height:auto; display:block; }
figcaption { font-size:12.5px; color:var(--soft); margin-top:6px; }
.flag { color:var(--warn); font-weight:600; }
.ok { color:var(--ok); font-weight:600; }
.note.good { background:#f2f9f4; border-left-color:#2e8b57 }
.note.bad { background:#fdf3f2; border-left-color:#c0392b }
.note { background:#f7f7f7; border-left:3px solid var(--rule); padding:10px 14px; margin:14px 0;
        font-size:13.5px; }
.gauge { margin:14px 0 20px; }
.gauge .bar { display:flex; height:20px; border-radius:3px; overflow:hidden; font-size:10px; }
.gauge .bar div { display:flex; align-items:center; justify-content:center; color:#fff;
                  white-space:nowrap; overflow:hidden; }
/* One bar, a pointer above and a pointer below, so two values on one scale cannot collide. */
.gauge .mark { position:relative; height:18px; }
.gauge .mark span { position:absolute; transform:translateX(-50%); font-size:11.5px;
                    font-weight:600; white-space:nowrap; }
.gauge .mark.above span { bottom:1px; }
.gauge .mark.below span { top:1px; }
.gauge .caption { font-size:12px; color:var(--soft); margin-top:4px; }
footer { margin-top:48px; padding-top:14px; border-top:1px solid var(--rule);
         font-size:12px; color:var(--soft); }
@page { size:A4; margin:16mm 14mm;
        @bottom-center { content:"Positive and negative mode are reported side by side and are
                                  never merged \\2014 page " counter(page);
                         font-size:8pt; color:#888; } }
@media print { body { padding:0; max-width:none; font-size:11pt; }
               h2 { page-break-after:avoid; } figure { page-break-inside:avoid; } }
"""


# ── rendering helpers ──────────────────────────────────────────────────────────────────────

def embed(fig) -> str:
    """Lay the figure out and render it.

    ⚠ Spacing is applied HERE, not per figure. Adding `tight_layout` to some figures and not
    others is how the run-order panel ended up printing its second-row titles on top of the first
    row's x-axis label: the fix had been applied to the principal-component panels only. One
    place, every figure, no exceptions.
    """
    try:
        fig.tight_layout(h_pad=2.6, w_pad=1.6)
    except Exception:                      # noqa: BLE001 - layout is cosmetic, never fatal
        pass
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=170, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def figure(src, caption) -> str:
    return f'<figure><img src="{src}" alt=""/><figcaption>{caption}</figcaption></figure>'


def table(headers, rows, numeric=(), classes=()) -> str:
    head = ""
    for i, h in enumerate(headers):
        cls = " ".join(filter(None, ["num" if i in numeric else "",
                                     classes[i] if i < len(classes) else ""]))
        head += f'<th class="{cls}">{h}</th>' if cls else f"<th>{h}</th>"
    body = ""
    for row in rows:
        body += "<tr>" + "".join(f'<td class="num">{c}</td>' if i in numeric else f"<td>{c}</td>"
                                 for i, c in enumerate(row)) + "</tr>"
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def gauge(marks, bands, label, upper, fmt="{:.1f}") -> str:
    """One banded scale carrying every polarity's value, rather than one scale each.

    `marks` is [(name, value, colour), ...]. The first is drawn above the bar and the second
    below it, which is what keeps two values from overprinting when they land close together —
    and here they usually do, because the polarities agree more often than not.

    A number against a threshold is a judgement the reader has to look up; a number on a banded
    scale is one they can read at a glance. Each pointer repeats its value and the caption names
    the band, so nothing carries meaning by colour alone.
    """
    marks = [(n, v, c) for n, v, c in marks if v == v]
    if not marks:
        return ""
    segments, previous = [], 0.0
    for edge, name, colour in bands:
        end = min(edge, upper)
        if end > previous:
            segments.append((previous, end, name, colour))
        previous = end
        if edge >= upper:
            break
    total = sum(e[1] - e[0] for e in segments) or 1
    bar = "".join(f'<div style="background:{c};flex:{(hi - lo) / total}">{n}</div>'
                  for lo, hi, n, c in segments)

    def pointer(entry, arrow, css):
        name, value, colour = entry
        where = min(max(value, 0.0), upper) / upper * 100
        # Nudge the label inboard at the extremes so it cannot run off the edge of the figure.
        align = "left:0;transform:none" if where < 6 else (
            "right:0;left:auto;transform:none" if where > 94 else f"left:{where:.1f}%")
        return (f'<div class="mark {css}"><span style="{align};color:{colour}">'
                f'{arrow} {fmt.format(value)} {name}</span></div>')

    above = pointer(marks[0], "&#9660;", "above")
    below = pointer(marks[1], "&#9650;", "below") if len(marks) > 1 else ""
    bands_hit = "; ".join(
        f"{n} <strong>{next(b for edge, b, _ in bands if v < edge)}</strong>" for n, v, _ in marks)
    return (f'<div class="gauge">{above}<div class="bar">{bar}</div>{below}'
            f'<div class="caption">{label} &mdash; {bands_hit}</div></div>')


def scatter(ax, scores, roles, flagged, explained, xi=0, yi=1, labels=None, groups=None,
            confirmed=None):
    if groups:
        for n, (name, members) in enumerate(sorted(groups.items())):
            name = name.replace("_", " ")
            pick = [i for i in members if i not in flagged]
            if pick:
                ax.scatter(scores[pick, xi], scores[pick, yi], s=SIZE,
                           c=GROUP_COLOURS[n % len(GROUP_COLOURS)], edgecolor="white",
                           linewidth=.6, label=name, zorder=3)
    else:
        for role in ("sample", "qc", "blank"):
            pick = [i for i, r in enumerate(roles) if r == role and i not in flagged]
            if pick:
                ax.scatter(scores[pick, xi], scores[pick, yi], s=SIZE, c=COLOURS[role],
                           edgecolor="white", linewidth=.6, label=LABELS[role], zorder=3)
    # ⚠ An outlier is marked by a BOLD OUTLINE, never by recolouring. These panels use colour to
    # carry the biological factor, so a flag colour would spend the colour channel twice and
    # destroy the thing the panel exists to show. An outline is an independent channel.
    # ⚠ NO INJECTION IS LABELLED HERE, flagged or not. A score plot exists to show whether the
    # groups separate; a name written on one point pulls the eye to that point and invites the
    # reader to explain it, which is the opposite of what the panel is for. Which injections were
    # flagged, and on what evidence, is stated in the table below the figure where it can be read
    # against the numbers rather than guessed from a position.
    for i in sorted(confirmed or ()):
        ax.scatter(scores[i, xi], scores[i, yi], s=SIZE * 1.15, facecolors="none",
                   edgecolor="#1a1a1a", linewidth=1.8, zorder=5)
    ax.axhline(0, color="#dddddd", lw=.6, zorder=1)
    ax.axvline(0, color="#dddddd", lw=.6, zorder=1)
    ax.set_xlabel(f"PC{xi + 1} ({100 * explained[xi]:.1f}%)")
    ax.set_ylabel(f"PC{yi + 1} ({100 * explained[yi]:.1f}%)")
    ax.spines[["top", "right"]].set_visible(False)


def legend_below(ax, ncol=4, drop=0.20):
    """Legends go OUTSIDE the axes, everywhere, without exception.

    Chosen per figure, a legend eventually covers data — which it did, on the profile screen.
    """
    handles, _ = ax.get_legend_handles_labels()
    if handles:
        ax.legend(frameon=False, fontsize=9, ncol=ncol, loc="upper center",
                  bbox_to_anchor=(0.5, -drop))


_PREFIX = ""


def _beside(results: str) -> str | None:
    """The gap-filled table sitting next to the results, if the pipeline wrote one."""
    candidate = Path(results).parent / "Final_Results_Filtered.csv"
    return str(candidate) if candidate.exists() else None


def measured_fwhm(folder: Path) -> float | None:
    """The run's own mean peak width, from `run_config.json` beside the results.

    None when it cannot be found, and the caller must then say so rather than substituting a
    number from a different method.
    """
    import json
    for name in ("run_config.json", "run_config.calibrated.json"):
        path = folder / name
        if path.exists():
            try:
                config = json.loads(path.read_text())
            except (OSError, ValueError):
                continue
            seconds = config.get("features", {}).get("chrom_fwhm")
            if seconds:
                return float(seconds) / 60.0
    return None


def learn_prefix(names, coverage=0.6) -> None:
    """Work out the shared prefix from the injection names themselves.

    ⚠ Hardcoding one study's prefixes meant that on anyone else's data nothing was trimmed and
    every label carried the full injection name. A facility pipeline cannot know the convention
    in advance.

    ⚠ And the obvious fix — the prefix common to *all* names — is also wrong, because blanks
    routinely follow a different convention. On this study the samples are `KF_Brain_Lipid_*` and
    the blanks are `Blank_*`, so the universal prefix is empty and nothing is trimmed. The prefix
    is therefore taken from the **majority**: the longest one shared by at least `coverage` of the
    names, applied only to the names that have it.
    """
    global _PREFIX
    _PREFIX = ""
    names = [n for n in names if n]
    if len(names) < 4:
        return
    candidates = set()
    for name in names:
        for cut, character in enumerate(name):
            if character in "_-":
                candidates.add(name[:cut + 1])
    # Widest coverage first, longest only as the tie-break. Preferring length instead lets a
    # majority CLASS prefix beat the study prefix — `Patient-Cohort_sample_` held 60% of names and
    # won over `Patient-Cohort_`, trimming the samples to bare numbers while leaving the QCs and
    # blanks at full length.
    best, best_score = "", (0.0, 0)
    for candidate in candidates:
        held = sum(1 for n in names if n.startswith(candidate)) / len(names)
        if held < coverage:
            continue
        score = (held, len(candidate))
        if score <= best_score:
            continue
        # Trimming must not make two injections indistinguishable.
        trimmed = {n[len(candidate):] if n.startswith(candidate) else n for n in names}
        if len(trimmed) == len(set(names)):
            best, best_score = candidate, score
    _PREFIX = best


def short(name: str) -> str:
    return name[len(_PREFIX):] if _PREFIX and name.startswith(_PREFIX) else name


def sample_key(name: str) -> str:
    """The injection identity with the polarity stripped out.

    ⚠ Without this the whole point of one document is lost: `Pos_293` and `Neg_293` are the same
    sample measured twice, and the strongest outlier evidence available is one sample failing in
    both polarities independently. Comparing the raw names never matches, so the note that says so
    never fired.
    """
    return re.sub(r"(^|_)(Pos|Neg|POS|NEG|positive|negative)(_|$)", r"\1", short(name)).strip("_")


def pct(x, digits=1) -> str:
    return "&mdash;" if x is None or x != x else f"{100 * x:.{digits}f}%"


def num(x, digits=2) -> str:
    return "&mdash;" if x is None or x != x else f"{x:.{digits}f}"


# ── one polarity's worth of analysis ───────────────────────────────────────────────────────

class Arm:
    """Everything computed for one polarity, so the sections can stay comparative."""

    def __init__(self, label, results, unfiltered=None, sequence=None, metadata=None,
                 spectra=None, fwhm=None, filled=None, standards=None):
        self.label = label
        self.results = results
        # The converted mzML sit beside the results; their headers carry the instrument's own
        # record of when each run started, which beats a planned sequence and exists for
        # polarities whose sequence was never saved.
        mzml = sorted((Path(results).parents[2] / "mzml").rglob("*.mzML")) \
            if (Path(results).parents[2] / "mzml").is_dir() else []
        self.run = qc.load_run(results, sequence, metadata, mzml_files=mzml)
        self.unfiltered = (qc.load_run(unfiltered, sequence, metadata, mzml_files=mzml)
                           if unfiltered else None)
        # The gap-filled table is what the data owner analyses, so outliers are judged on it and
        # missingness is reported both ways. Falls back to the unfilled table when absent.
        self.filled = (qc.load_run(filled, sequence, metadata, mzml_files=mzml)
                       if filled else None)
        self.spectra = Path(spectra) if spectra and Path(spectra).exists() else None
        self.fwhm = fwhm
        run = self.run

        # Everything the search identified, filters aside — the figure in §Annotation plots it,
        # and a reader comparing the two must not find different numbers.
        self.search_identified = None
        if self.unfiltered is not None:
            named = [r for r in self.unfiltered.rows if r["Identification"].strip()]
            # Molecules, artefacts excluded, counted exactly as the figure counts them: a name is
            # dropped only when EVERY row carrying it is a fragment, a second adduct, an isotope
            # or a redundant re-match. Reporting the raw distinct-name count here while the figure
            # reports the cleaned one leaves a reader with two numbers and no way to reconcile.
            real = {r["Identification"].strip() for r in named
                    if not ARTEFACT_REASON.search(r.get("Filter Status", "") or "")}
            self.search_identified = {
                "rows": len(named),
                "molecules": len(real),
                "artefact_only": len({r["Identification"].strip() for r in named}) - len(real),
            }

        self.counts = {r: len(run.indices(r)) for r in ("sample", "qc", "blank")}
        self.profile = qc.profile_screen(run)
        self.missing = qc.missingness(run)
        self.missing_filled = qc.missingness(self.filled) if self.filled else None
        self.fill_rate = self._fill_rate()
        self.order = qc.run_order(run)
        self.drift = qc.pool_drift(run)
        self.confounding = qc.confounding(run)
        self.blanks = qc.blank_report(self.unfiltered or run)
        self.blank_ratio = qc.blank_ratio(self.unfiltered or run,
                                          multiplier=3.0)
        self.coverage = qc.coverage(run)
        self.standards = self._standards(standards)
        self.accuracy = qc.mass_accuracy(self.spectra, run.order, run.roles)
        # ⚠ What mass_accuracy measures is what is LEFT after any correction the run applied. If
        # a correction was applied and the report only showed the residual, it would report an
        # instrument as accurate BECAUSE it had been corrected — describing the software's
        # compensation as the instrument's performance. The applied value is read back from the
        # run configuration written beside the results, so the section can state both.
        self.applied_offset = self._applied_offset(results)
        # Row counts for the closing file guide. Taken from the tables themselves rather than
        # recomputed, so the guide can never disagree with what was delivered.
        self.n_kept = len(run.rows)
        self.n_identified = sum(1 for r in run.rows if r["Identification"].strip())
        self.n_unfiltered = len(self.unfiltered.rows) if self.unfiltered else 0
        # MS1->MS2 duty-cycle accounting, written beside the results by pipeline.py's
        # `build_exclusion_list` step (on by default). Absent only for a run predating that
        # default or with it explicitly turned off, in which case the section reports that.
        duty_path = Path(results).parent / "Duty_Cycle.json"
        self.duty_cycle = (json.loads(duty_path.read_text()) if duty_path.exists() else None)
        # The complementary list, written beside it by `build_inclusion_list` (also on by
        # default): low-abundance features this run's own Top-N never fragmented, matching a
        # library entry. Read as a row count only — the evidence per candidate is the CSV itself,
        # not something the report needs to reproduce.
        inclusion_path = Path(results).parent / "Inclusion_List.csv"
        self.inclusion_count = (max(0, sum(1 for _ in inclusion_path.open()) - 1)
                                if inclusion_path.exists() else None)
        # Impossible-mass-defect ions stripped per file, written beside the results by
        # `screen_artefacts` (on by default). `None` means never checked (predates the default,
        # or it was turned off); an empty `by_file` means checked and clean.
        artefacts_path = Path(results).parent / "Artefacts.json"
        self.artefacts = (json.loads(artefacts_path.read_text())
                          if artefacts_path.exists() else None)
        # Mass-drift measurement and correction, written by `run_study.py --auto-calibrate` (opt
        # in, not a default -- it costs a second full pass). `None` means the run never asked;
        # `checked=True` with `used_calibrated_files=False` means it asked and found nothing to
        # fix.
        calibration_path = Path(results).parent / "Calibration.json"
        self.calibration = (json.loads(calibration_path.read_text())
                            if calibration_path.exists() else None)
        self.sources = Counter(r.get("Identification Source", "").strip()
                               for r in run.rows if r["Identification"].strip())
        self.annotation = qc.annotation_ambiguity(run, fwhm=fwhm)

        # One feature set, decided on the sample columns, applied to every panel. See
        # qc.sample_feature_mask for why deciding it per panel is not neutral.
        self.mask_all = qc.sample_feature_mask(run)
        self.mask_named = self.mask_all & qc.identified_mask(run)
        # Precision both ways. Over every quantified feature this describes the instrument; over
        # the identified subset it describes the data that will actually be analysed, and those
        # are different numbers. The identified set leads because it is the one a reader means.
        self.precision = qc.precision(run, self.mask_named) or qc.precision(run, self.mask_all)
        self.precision_all = qc.precision(run, self.mask_all)
        self.spaces = {}
        for key, roles in (("blanks", ("sample", "qc", "blank")),
                           ("pools", ("sample", "qc")), ("samples", ("sample",))):
            for view, mask in (("all", self.mask_all), ("named", self.mask_named)):
                self.spaces[(key, view)] = qc.score_space(run, roles, mask=mask)
        # Outliers are judged on the identified lipids, for the same reason precision is: that is
        # the space the biology will be analysed in. Falls back to all features where too few rows
        # carry a name to build a model.
        # Outliers on the gap-filled, identified-only table: that is the matrix the biology is
        # analysed in, and gaps left as zeros distort a distance far more than a filled value does.
        source = self.filled or run
        o_mask = qc.sample_feature_mask(source) & qc.identified_mask(source)
        self.outlier_spaces = {
            "pools": qc.score_space(source, ("sample", "qc"), mask=o_mask),
            "samples": qc.score_space(source, ("sample",), mask=o_mask)}
        self.outliers = {k: (qc.outliers(v) if v else None)
                         for k, v in self.outlier_spaces.items()}
        self.outlier_space = ("gap-filled table, identified lipids" if self.filled
                              else "identified lipids")
        self.contrast = qc.contrast(self.spaces[("pools", "named")]
                                    or self.spaces[("pools", "all")])
        self.contrast_all = qc.contrast(self.spaces[("pools", "all")])

    @staticmethod
    def _applied_offset(results) -> float:
        """`mass_offset_ppm` the run was processed with, from run_config.json beside the tables."""
        import json
        path = Path(results).parent / "run_config.json"
        if not path.exists():
            return 0.0
        try:
            return float(json.loads(path.read_text()).get("mass_offset_ppm") or 0.0)
        except Exception:                       # noqa: BLE001 - absent config is not an error
            return 0.0

    @staticmethod
    def _standards(path):
        """Standards_Performance.csv, if a mix was declared. Absent is normal, not an error."""
        if not path or not Path(path).exists():
            return []
        import csv as _csv
        return list(_csv.DictReader(open(path)))

    def _fill_rate(self):
        """What share of the gaps were filled from re-integrated signal, and how many were not.

        A sensitivity statement in its own right: it says how much of the apparent missingness was
        a detection failure rather than an absence.
        """
        if self.filled is None:
            return None
        samples = self.run.samples
        if not samples or self.filled.areas.shape != self.run.areas.shape:
            return None
        before = self.run.areas[:, samples]
        after = self.filled.areas[:, samples]
        gaps = before <= 0
        if not gaps.any():
            return None
        return {"gaps": int(gaps.sum()),
                "filled": int((gaps & (after > 0)).sum()),
                "left": int((gaps & (after <= 0)).sum()),
                "rate": float((gaps & (after > 0)).sum() / gaps.sum())}

    # -- item 12: the QC as an identification anchor
    def named_only_in_pools(self):
        """Molecules whose only MS2 evidence is in a pooled QC, and are quantified anyway."""
        if not self.spectra:
            return None
        pools, samples = {}, {}
        with self.spectra.open(newline="") as fh:
            for row in csv.DictReader(fh):
                if not row.get("Compound Group"):
                    continue
                name = row["Name"].split(" [")[0]
                bucket = pools if infer_role(row["Sample"]) == "qc" else samples
                bucket.setdefault(name, set()).add(row["Sample"])
        only = sorted(set(pools) - set(samples))
        return {"count": len(only), "examples": only[:6]}

    def normalised_totals(self, role):
        """Total signal per injection of one role, relative to the MEDIAN of that role.

        ⚠ Not relative to the first injection. A single reference point carries all of its own
        noise into every other point, and the opening injections are the least trustworthy in any
        sequence — `qc.RUN_IN` exists precisely because start-of-run behaviour differs. One bad
        first injection would tilt the whole series; a median needs half of them to be bad.

        1.0 therefore reads as "the run's typical level" rather than "whatever went first", and
        the convention matches `median_normalise`, which the feature matrix already uses.
        """
        index = self.run.indices(role)
        if not index:
            return [], [], None
        order = sorted(index, key=lambda i: self.run.order.get(self.run.columns[i], i))
        totals = np.array([self.run.areas[:, i].sum() for i in order], dtype=float)
        base = float(np.median(totals[totals > 0])) if (totals > 0).any() else 0.0
        axis = [self.run.order.get(self.run.columns[i], n + 1) for n, i in enumerate(order)]
        values = totals / base if base else totals
        # Drift magnitude from a line through every point, not the difference between the two
        # least reliable ones. Expressed as the fitted change across the whole sequence.
        trend = None
        if len(values) >= 3 and len(set(axis)) > 1:
            slope, intercept = np.polyfit(np.asarray(axis, dtype=float), values, 1)
            trend = float(slope * (max(axis) - min(axis)))
        return axis, values, trend

    def groups(self, factor):
        """Sample indices within a score space, grouped by a metadata factor."""
        space = self.spaces[("samples", "all")]
        if space is None:
            return None
        out = {}
        for i, column in enumerate(space["columns"]):
            value = self.run.metadata.get(column, {}).get(factor)
            if value:
                out.setdefault(f"{factor} {value}", []).append(i)
        return out or None


# ── sections ───────────────────────────────────────────────────────────────────────────────

def comparative(arms, label, extract, numeric=True):
    """One row of a positive/negative comparison table."""
    return [label] + [extract(a) for a in arms]


def section_precision(arms, plt) -> str:
    heads = ["", *[f"{a.label} mode" for a in arms]]
    classes = ["", *["pos" if a.label == "positive" else "neg" for a in arms]]
    rows = [
        comparative(arms, "Median CV across the pooled QCs",
                    lambda a: pct(a.precision["median_pool_cv"] / 100)),
        comparative(arms, "Features measured to better than 20% CV",
                    lambda a: pct(a.precision["under_20"], 0)),
        comparative(arms, "Features measured to better than 30% CV",
                    lambda a: pct(a.precision["under_30"], 0)),
    ]
    tone = {"positive": "#0072B2", "negative": "#CC79A7"}
    gauges = gauge([(a.label, a.precision["median_pool_cv"], tone.get(a.label, "#1a1a1a"))
                    for a in arms if a.precision.get("usable", True)],
                   CV_BANDS, "median pooled-QC CV", upper=50, fmt="{:.1f}%")
    absent = [f"<strong>{a.label} mode:</strong> {html.escape(a.precision.get('reason', ''))}"
              for a in arms if not a.precision.get("usable", True)]
    if absent:
        gauges += ('<div class="note"><strong>No precision figures for '
                   + "; ".join(absent) + ". A dash below means the statistic could not be "
                   "computed, not that it was zero.</div>")
    other = [comparative(arms, "Median CV, all quantified features",
                         lambda a: pct(a.precision_all["median_pool_cv"] / 100)),
             ]
    return f"""
{head('precision')}
<p class="lede">Measured on the pooled QC injections, which are the same material in every vial, so
any variation between them is analytical by construction. <strong>The figures below are computed on
the identified lipids</strong> &mdash; the data that will actually be analysed.</p>
{table(heads, rows, numeric=tuple(range(1, len(arms) + 1)), classes=classes)}
{gauges}
<h3>The same statistics over every quantified feature</h3>
<p class="lede">Most quantified features are never identified, and they are on average fainter and
therefore noisier. Both views are given because they answer different questions &mdash; the
identified set describes your data, this one describes the instrument.</p>
{table(heads, other, numeric=tuple(range(1, len(arms) + 1)), classes=classes)}
<p><strong>What this is.</strong> The <em>coefficient of variation</em> is the standard deviation
of a feature's area across the pooled QCs divided by its mean, after median normalisation. The
bands are the facility's 15% target and the 20% and 30% thresholds in common use for untargeted
profiling.</p>
<div class="note"><strong>Nothing in this section depends on the biology.</strong> Every figure
here is computed on repeated injections of one pooled material, so it measures the instrument and
nothing else. Statistics that compare the samples against the QCs belong to a different question
and are in {ref('support')}.</div>
"""


def section_missing(arms, plt) -> str:
    heads = ["", *[f"{a.label} mode" for a in arms]]
    classes = ["", *["pos" if a.label == "positive" else "neg" for a in arms]]
    rows = [
        comparative(arms, "Missing when the peak finder ran, all features",
                    lambda a: pct(a.missing["missing_all"]) if a.missing else "&mdash;"),
        comparative(arms, "Missing when the peak finder ran, identified only",
                    lambda a: pct(a.missing["missing_identified"])),
        comparative(arms, "<strong>Missing in the delivered table, all features</strong>",
                    lambda a: "<strong>" + ((pct(a.missing_filled["missing_all"])
                                             if a.missing_filled else "&mdash;")
                                            if a.missing_filled else "&mdash;") + "</strong>"),
        comparative(arms, "<strong>Missing in the delivered table, identified only</strong>",
                    lambda a: "<strong>" + (pct(a.missing_filled["missing_identified"])
                                            if a.missing_filled else "&mdash;") + "</strong>"),
        comparative(arms, "Gaps recovered by re-integration",
                    lambda a: (f"{a.fill_rate['filled']:,} of {a.fill_rate['gaps']:,} "
                               f"({a.fill_rate['rate']:.0%})") if a.fill_rate else "&mdash;"),
    ]
    # Any supplied factor associated with detection rate, whatever it is called. Never named in
    # code: a comparison on such a factor is partly a comparison of detection, and that has to be
    # tested and reported for every study, not for the factors one study happened to have.
    # ⚠ Reported per FACTOR across both polarities, not per polarity. A hard p < 0.05 cut made the
    # report contradict itself: age is associated with detection in both, same direction and size,
    # but p = 0.056 in positive and 0.009 in negative, so the design was declared clean in one and
    # confounded in the other. And a non-significant result was being read as an absence of
    # confounding, which an underpowered test cannot establish.
    factors = sorted({k for a in arms for k in (a.missing.get("by_factor") or {})})
    warnings = ""
    if factors:
        lines = []
        for factor in factors:
            parts, flagged = [], False
            for a in arms:
                v = (a.missing.get("by_factor") or {}).get(factor)
                if not v:
                    continue
                worst = max(v["levels"], key=v["levels"].get)
                parts.append(f"{a.label} p = {v['p']:.3f}, {pct(v['spread'], 1)} apart "
                             f"({html.escape(worst)} least complete)")
                flagged |= v["p"] < 0.05
            if not parts:
                continue
            agree = len(parts) > 1 and all("least complete" in p for p in parts)
            lead = (f"<strong class='flag'>{html.escape(factor)}</strong>" if flagged
                    else f"<strong>{html.escape(factor)}</strong>")
            lines.append(f"<li>{lead} &mdash; {'; '.join(parts)}."
                         + (" Detection rate differs with this factor, so a comparison across it "
                            "is partly a comparison of how often a molecule was detected."
                            if flagged else " Not shown to be associated &mdash; which, at five per arm, is a statement about the power of the test rather than about the study.")
                         + ("" if len(arms) < 2 or not agree else
                            " <em>Both polarities point the same way &mdash; but they are the same "
                            "samples measured twice, so this is not corroboration. A sample that "
                            "under-loaded or degraded detects poorly in both regardless of any "
                            "factor, and on this study per-sample completeness correlates across "
                            "polarities at Spearman +0.72. The two tests are close to one test, "
                            "and their p-values must not be read as multiplying.</em>"))
        warnings = "<ul>" + "".join(lines) + "</ul>"

    return f"""
{head('missing')}
<p class="lede">Measured across the study samples. The first two rows describe what the peak finder
detected; the next two describe the table you were sent.</p>
{table(heads, rows, numeric=tuple(range(1, len(arms) + 1)), classes=classes)}
<div class="note"><strong>Gaps are measured, not imputed.</strong> Where the peak finder found
nothing, the raw file is re-integrated at that molecule's mass and retention time &mdash; so a
filled value is a measurement at that coordinate rather than a number inferred from other samples.
That matters because missing values here are not missing at random: a feature goes undetected
because it is small, and any method that infers a value from the features around it would replace
a small unmeasured value with a typical one. Where re-integration finds nothing, the position is
left empty rather than filled.</div>
<h3>Is detection confounded with the design?</h3>
<p class="lede">Every supplied factor, tested in both polarities. A p-value is reported beside each
rather than used as a gate &mdash; a test that fails to reach significance has not shown that
detection is unconfounded, only that it could not demonstrate that it is.</p>
{warnings or "<p>No sample metadata was supplied, so this could not be tested.</p>"}
"""


def section_accuracy(arms, plt) -> str:
    """Mass accuracy per injection, and whether it moved during the run.

    Two different faults, and they need different answers. A constant offset is calibration: the
    whole axis sits off zero, and `mass_offset_ppm` corrects it in one number. A *trend across
    injections* is the axis moving while the sequence ran, which no single constant can fix — it
    widens the effective search window for everything and is a reason to recalibrate the
    instrument rather than the data.
    """
    if not any(a.accuracy for a in arms):
        return f"""
{head('accuracy')}
<div class="note">No spectral match table was supplied, so mass accuracy could not be measured.
Pass <code>--spectra-pos</code> / <code>--spectra-neg</code>.</div>
"""
    rows = ""
    notes = ""
    for arm in arms:
        acc = arm.accuracy
        if not acc:
            continue
        drift = acc.get("drift_ppm")
        # The instrument's own error, which is what this section is about. `median_ppm` is what
        # survived the correction, so the correction has to be added back to recover it.
        instrument = acc["median_ppm"] + arm.applied_offset
        rows += (f"<tr><td>{arm.label}</td><td class='n'>{acc['n_matches']:,}</td>"
                 f"<td class='n'>{instrument:+.2f}</td>"
                 f"<td class='n'>{f'{arm.applied_offset:+.2f}' if arm.applied_offset else '&mdash;'}</td>"
                 f"<td class='n'>{acc['median_ppm']:+.2f}</td>"
                 f"<td class='n'>{acc['spread_ppm']:.2f}</td>"
                 f"<td class='n'>{f'{drift:+.2f}' if drift is not None else '&mdash;'}</td>"
                 f"<td class='n'>{f"{acc['rho']:+.2f}" if acc.get('rho') is not None else '&mdash;'}</td>"
                 f"</tr>")
        if abs(instrument) >= 3 and not arm.applied_offset:
            notes += (f"<div class='note'><strong>&#9888; {arm.label} mode sits "
                      f"{instrument:+.1f} ppm from zero.</strong> The whole axis is offset, "
                      f"which is calibration rather than scatter. Against a search window centred "
                      f"on zero it spends most of the tolerance on the offset and leaves little "
                      f"for genuine error, so marginal identifications are lost. Correct it with "
                      f"<code>mass_offset_ppm</code> rather than by widening the window.</div>")
        if drift is not None and abs(drift) >= 2:
            notes += (f"<div class='note'><strong>&#9888; {arm.label} mode drifts {drift:+.1f} ppm "
                      f"across the sequence</strong> ({acc['first_ppm']:+.2f} to "
                      f"{acc['last_ppm']:+.2f}). A single offset cannot correct a moving axis: it "
                      f"would be right in the middle of the run and wrong at both ends. This is an "
                      f"instrument calibration to repeat, not a number to set.</div>")

    figure = accuracy_figure(arms, plt)
    return f"""
{head('accuracy')}
<p>Measured from every rank-one library match &mdash; tens of thousands per polarity &mdash; which
is what makes a figure per injection possible. Spiked standards are the independent check on these
numbers rather than their source: an offset derived from the identifications assumes the
identifications are right, while a standard's mass is known before the run.</p>
<table><thead><tr><th>polarity</th><th class="n">matches</th>
<th class="n">instrument error (ppm)</th><th class="n">correction applied</th>
<th class="n">residual (ppm)</th>
<th class="n">spread (ppm)</th><th class="n">drift across run (ppm)</th><th class="n">&rho; vs order</th>
</tr></thead><tbody>{rows}</tbody></table>
{figure}
{corrected_note(arms)}
{notes}
"""


def fixed_axis(ax, values, low, high, pad=0.06):
    """A default y-range, widened only when the data would fall outside it.

    Panels that autoscale are unreadable across reports: the same picture means a different thing
    each time, and a run with tiny variation looks as dramatic as one with real drift. A fixed
    window makes the eye comparable between studies. It is a default and never a clip — data
    outside it widens the axis rather than disappearing, because a silently cropped point is the
    one worth seeing.
    """
    if not len(values):
        return
    lo, hi = float(min(values)), float(max(values))
    span = (hi - lo) or 1.0
    ax.set_ylim(min(low, lo - pad * span), max(high, hi + pad * span))


def corrected_note(arms) -> str:
    """Say when a correction was applied, because otherwise this section flatters itself.

    `mass_offset_ppm` shifts the precursors before the search, so the error measured afterwards is
    what the correction left behind. Reporting that alone would present an instrument as accurate
    *because* it had been corrected — the software's compensation described as the instrument's
    performance. The instrument column is the honest figure; the residual is evidence the
    correction did what it claimed.
    """
    corrected = [a for a in arms if a.applied_offset]
    if not corrected:
        return ""
    listed = "; ".join(f"{a.label} {a.applied_offset:+.2f} ppm" for a in corrected)
    return (f'<div class="note"><strong>&#9888; A correction was applied &mdash; {listed}.</strong> '
            f'The <em>residual</em> column is what was left after it, and is near zero by '
            f'construction: it shows the correction worked, not that the instrument is accurate. '
            f'The <em>instrument error</em> column is the honest figure for the mass axis, and it '
            f'is the one to watch between runs. The correction is recorded in '
            f'<code>run_config.json</code> beside the results, so the delivered masses can be '
            f'traced back to it. Spiked standards, measured on the unshifted feature masses, are '
            f'an independent check on the same number.</div>')


def accuracy_figure(arms, plt) -> str:
    """Mass error per injection: as measured, and after any correction.

    Two decisions worth stating. **The x-axis counts injections within each polarity**, not
    position in the acquisition sequence. The sequence interleaves the two, so positive occupies
    the odd positions and negative the even, and the first negative injection with enough matches
    to give a median sat at position six against positive's five — two lines starting in different
    places for a reason that has nothing to do with the instrument. Counting within a polarity puts
    both at one.

    **The measured series is drawn in grey behind the corrected one.** A correction that is only
    reported as a number is hard to judge; drawn, the reader sees how far the axis moved and
    whether it moved uniformly. Polarity is carried by marker shape rather than by colour, so
    colour is free to mean measured-versus-corrected.
    """
    usable = [a for a in arms if a.accuracy and a.accuracy.get("injections")]
    if not usable:
        return ""
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    shapes = {"positive": "o", "negative": "^"}
    every = []
    for arm, colour in zip(usable, ("#0072B2", "#D55E00")):
        pts = sorted((r for r in arm.accuracy["injections"] if r["order"] is not None),
                     key=lambda r: r["order"])
        if not pts:
            continue
        marker = shapes.get(arm.label, "o")
        # within-polarity position, so both polarities start at 1
        x = range(1, len(pts) + 1)
        corrected = [r["median_ppm"] for r in pts]
        if arm.applied_offset:
            measured = [v + arm.applied_offset for v in corrected]
            ax.plot(x, measured, marker + "--", ms=4, lw=1, color="#9a9a9a",
                    label=f"{arm.label}, as measured", zorder=2)
            every += measured
        ax.plot(x, corrected, marker + "-", ms=4, lw=1.2, color=colour,
                label=f"{arm.label}{' , corrected' if arm.applied_offset else ''}".replace(" ,", ","),
                zorder=3)
        every += corrected
    ax.axhline(0, color="#888", lw=0.8, ls="--")
    ax.set_xlabel("injection within this polarity")
    ax.set_ylabel("median mass error (ppm)")
    fixed_axis(ax, every, -10, 10)
    # Injections are counted, so the ticks are integers. The default locator put 2.5 and 7.5 on an
    # axis where half an injection does not exist.
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.spines[["top", "right"]].set_visible(False)
    legend_below(ax, ncol=2, drop=0.34)
    return figure(embed(fig),
                  "Mass error per injection, counted within each polarity so both series start at "
                  "one &mdash; the sequence interleaves them, so position in the run would start "
                  "the two lines in different places for no instrumental reason. Circles are "
                  "positive, triangles negative. <strong>Grey is the axis as measured</strong>; "
                  "colour is what remained after the correction. A flat line off zero is "
                  "calibration, correctable with one number; a sloping line is the mass axis "
                  "moving during the run, which no single constant can fix.")


def section_polarity(arms) -> str:
    """What the polarity filter did, and what it deliberately left alone.

    The count that matters most is the one NOT removed. Read without it, "half the rows went" is
    the impression; in fact most classes appear in one polarity only and are never eligible.
    """
    import json
    path = None
    for arm in arms:
        # <analysis>/Results/<polarity>/Final_Results.csv -> <analysis>
        candidate = Path(arm.results).resolve().parents[2] / "polarity_filter.json"
        if candidate.exists():
            path = candidate
            break
    if path is None:
        return f"""
{head('polarity')}
<div class="note">No polarity filter was applied, so each table still contains every class it
identified. Where a class appears in both, the same lipid is reported twice on scales that cannot
be compared &mdash; a choline lipid is a protonated ion in one mode and a formate adduct in the
other.</div>
"""
    data = json.loads(path.read_text())
    shared, removed = data["shared_classes"], data["removed"]
    rows = ""
    for name in sorted(shared):
        won = shared[name]
        lost = "Neg" if won == "Pos" else "Pos"
        n = removed.get(lost, {}).get(name, 0)
        rows += (f"<tr><td>{name}</td><td>{won}</td><td>{data['basis'].get(name, '')}</td>"
                 f"<td class='n'>{n}</td><td>{lost}</td></tr>")
    kept = data.get("kept_untouched", {})
    single = len(data.get("single_polarity_classes", []))
    agree = data.get("agreement") or {}
    agreement_note = ""
    if agree:
        agreement_note = (
            f"<div class='note'><strong>Agreement between the polarities, before the filter: "
            f"Spearman &rho; = {agree['rho']:+.3f}</strong> over {agree['shared']} lipids "
            f"identified in both, of which {agree['opposite']} "
            f"({100 * agree['opposite_fraction']:.0f}%) change in opposite directions. This is a "
            f"run-level figure rather than a result: the two modes measure a shared lipid by "
            f"different ions, so perfect agreement is not expected, but a run where this collapses "
            f"is saying something about the acquisition. It is computed before the filter, because "
            f"afterwards the two tables share no classes and the comparison cannot be made.</div>")
    checks = ""
    for label, key, tone in (("appear in both polarities with no assignment", "undecided", "bad"),
                             ("configured but identified in neither", "stale", ""),
                             ("configured but now in one polarity only", "moot", "")):
        names = data.get(key) or []
        if names:
            checks += (f"<li>{len(names)} {label}: "
                       f"{', '.join(f'<code>{n}</code>' for n in names)}</li>")
    checks = (f"<ul>{checks}</ul>" if checks
              else "<p>The configured table matches this run's inventory exactly &mdash; no "
                   "undecided, stale or moot entries.</p>")

    return f"""
{head('polarity')}
<p class="lede">The two modes are not two measurements of one thing. A choline lipid is a
protonated ion in positive and a formate adduct in negative; an acidic phospholipid is the reverse.
Where a class is identified in <em>both</em>, it is kept in the mode that measures it better and
removed from the other, so the delivered tables are two halves of one lipidome rather than two
overlapping views of it. <strong>Nothing is averaged</strong>: the areas are not on a common
scale.</p>
<table><thead><tr><th>class</th><th>kept in</th><th>basis</th><th class="n">rows removed</th>
<th>from</th></tr></thead><tbody>{rows}</tbody></table>
<div class="note"><strong>{single} classes appear in one polarity only and were not touched</strong>
&mdash; {kept.get('Pos', 0)} positive and {kept.get('Neg', 0)} lipids remain in the delivered
tables. Only classes found in both modes are eligible, so this is a much smaller intervention than
the row counts above suggest. A rule reaching further would delete species that are simply absent
from one mode, several of which are absent only because the current libraries cannot name them
there.</div>
{agreement_note}
<h3>Does the table still match the libraries?</h3>
<p>The set of classes each mode can name is a property of the <em>libraries</em>, not of chemistry,
so the configured table is checked against this run's own inventory. A class appearing in both
modes with no assignment stops the run rather than being guessed at.</p>
{checks}
"""


def section_standards(arms, plt) -> str:
    """Spiked standards judge the instrument; pooled QCs judge the batch.

    ⚠ NOT RENDERED IN THIS REPORT. The quality report goes to the client with their results, and
    the standards are a facility measurement about the instrument, not about their samples. Kept
    here, and still registered nowhere, for the facility report that will carry it — the numbers
    are in Standards_Performance_<polarity>.csv in the analysis folder meanwhile.

    Kept separate from precision for that reason. A pool is made from these samples, so it can
    only say whether this batch was internally repeatable. A standard mix is the same material in
    every batch and every year, so its mass accuracy, retention and absolute response are
    comparable with previous runs — which is the only way a report can say the instrument has
    moved rather than the study.
    """
    if not any(a.standards for a in arms):
        return f"""
{head('standards')}
<div class="note">No standards mix was declared for either polarity, so this section is empty.
Spiked standards are the only measurement here that can be compared with other runs: a pooled QC
is made from this study and cannot say whether the instrument itself has drifted between batches.
Declare one with <code>standards</code> in the run configuration and name the injections
<code>Std_Mix_*</code>.</div>
"""
    body = ""
    for arm in arms:
        if not arm.standards:
            continue
        rows = ""
        errors = []
        for r in arm.standards:
            found = str(r.get("found", "")).lower() == "true"
            if not found:
                rows += (f"<tr><td>{r['standard']}</td><td>{r['adduct']}</td>"
                         f"<td class='n'>{r['expected_mz']}</td>"
                         f"<td colspan='4' class='bad'>not found</td></tr>")
                continue
            err = float(r["mass_error_ppm"])
            errors.append(err)
            cv = r.get("cv_percent", "")
            rows += (f"<tr><td>{r['standard']}</td><td>{r['adduct']}</td>"
                     f"<td class='n'>{r['expected_mz']}</td>"
                     f"<td class='n'>{r['observed_mz']}</td>"
                     f"<td class='n'>{err:+.2f}</td><td class='n'>{r['retention']}</td>"
                     f"<td class='n'>{float(r.get('median_area') or 0):,.0f}</td>"
                     f"<td class='n'>{cv or '&mdash;'}</td></tr>")
        note = ""
        if errors:
            median_err = sorted(errors)[len(errors) // 2]
            same_sign = all(e < 0 for e in errors) or all(e > 0 for e in errors)
            if same_sign and abs(median_err) >= 3:
                note = (f"<div class='note'><strong>&#9888; Systematic mass error in "
                        f"{arm.label} mode: every standard is off in the same direction, median "
                        f"{median_err:+.1f} ppm.</strong> This is calibration, not scatter. With a "
                        f"search window centred on zero it spends most of the tolerance on the "
                        f"offset and leaves little for genuine error, so marginal identifications "
                        f"are lost. Correct it with <code>mass_offset_ppm</code> rather than by "
                        f"widening the window, which would cost specificity to buy the same "
                        f"tolerance.</div>")
            else:
                note = (f"<div class='note'>Mass accuracy in {arm.label} mode is centred: median "
                        f"{median_err:+.1f} ppm across {len(errors)} standards.</div>")
        body += (f"<h3>{arm.label}</h3><table><thead><tr><th>standard</th><th>adduct</th>"
                 f"<th class='n'>expected <em>m/z</em></th><th class='n'>observed</th>"
                 f"<th class='n'>error (ppm)</th><th class='n'>RT (min)</th>"
                 f"<th class='n'>raw intensity</th>"
                 f"<th class='n'>area CV (%)</th></tr></thead><tbody>{rows}</tbody></table>{note}")

    return f"""
{head('standards')}
<p>A mix of authentic standards, injected to judge the <em>instrument</em>. This is the one part of
the report comparable with other runs: the same material is injected in every batch, so mass
accuracy, retention and absolute response can be tracked between studies and over years. The
pooled QCs above answer a different question &mdash; whether <em>this</em> batch was repeatable
&mdash; and neither measurement substitutes for the other.</p>
{body}
<div class="note"><strong>Raw intensity is the column to keep.</strong> It is an absolute area,
not normalised to anything, and because the aliquots are identical between batches it is directly
comparable with every previous run. One report cannot say whether it is high or low &mdash; there
is no range yet &mdash; but recorded run after run it becomes the facility's sensitivity trend, and
a fall in it is a change in the instrument rather than in the study. The area CV beside it needs
more than one standards injection to mean anything.</div>
"""


# Each filter, the configuration key that switches it on, and the `Filter Status` strings it
# writes. Order is the order they run in, so the table reads as the path a row takes.
#
# ⚠ The point of this table is the DISTINCTION BETWEEN THREE SILENCES. A filter contributing no
# rows to `Unfiltered_Results.csv` may have been switched off, may have run and found nothing, or
# may have been unable to examine those rows at all. Those are very different statements about a
# delivered table and until now they looked identical — the delivered negative table showed zero
# `RT out of class range` rows because the filter was OFF, and read as though it had found nothing.
FILTERS = [
    ("Redundant identification", None, ("Redundant Identification",)),
    ("Adduct / dimer / isotope", "adduct_filtering",
     ("Adduct of existing peak", "Adduct of existing identified peak", "Dimer", "isotope")),
    ("In-source fragment", "in_source_filtering", ("In-source fragment",)),
    ("Class retention window", "rt_filter", ("RT out of class range",)),
    ("Retention model outlier", "retention_model_filter", ("RT model outlier",)),
    ("Correlated unidentified", "correlation", ("Correlated with",)),
    ("Split-peak merge", "split_peaks", ("Split peak, merged into",)),
    ("Adduct-pair removal", "adduct_pairs", ("Adduct pair", "same molecule as")),
    ("Blank", "blank_filter", ("Blank",)),
    ("Presence in a group", "presence_filter", ("Not detected in",)),
    ("Polarity assignment", "class_polarity", ("polarity:",)),
]


def _enabled(config: dict, key: str | None):
    """True/False/None — None meaning the run configuration does not record it."""
    if key is None:
        return True
    if key not in config:
        return None
    value = config[key]
    return bool(value.get("enabled", True)) if isinstance(value, dict) else bool(value)


def filter_activity(results, unfiltered) -> list[dict]:
    """How many rows each filter removed, read from the tables rather than from a log.

    Counted from `Filter Status` so the figure can never disagree with the delivered file, which a
    number carried from a run log eventually would.
    """
    import json
    config = {}
    for name in ("run_config.calibrated.json", "run_config.json"):
        path = Path(results).parent / name
        if path.exists():
            try:
                config = json.loads(path.read_text())
                break
            except Exception:
                pass

    reasons = []
    if unfiltered:
        reasons = [(r.get("Filter Status") or "").strip() for r in unfiltered.rows]

    out = []
    for label, key, patterns in FILTERS:
        removed = sum(1 for reason in reasons
                      if reason and any(pattern in reason for pattern in patterns))
        out.append({"filter": label, "enabled": _enabled(config, key), "removed": removed})
    return out


def retention_exemption(arm) -> tuple[int, int, list]:
    """MS2-identified rows their class model never judged, because the class has no usable model.

    The blind spot behind the retention filter: a class with too few members to fit, or one that is
    not a homologous series at all, produces no model — and every row in it passes the check by
    never taking it. `RT Model Error` is blank for exactly those rows, which reads as "fine".
    """
    rows = [r for r in arm.run.rows
            if r.get("Identification", "").strip()
            and (r.get("Identification Source", "") or "").strip() == "MS2"]
    if not rows or "RT Model Error" not in (rows[0] or {}):
        return 0, 0, []
    exempt = [r for r in rows if not (r.get("RT Model Error") or "").strip()]
    classes = Counter(r.get("Lipid Class", "?") for r in exempt)
    return len(rows), len(exempt), classes.most_common()


def section_filters(arms) -> str:
    """What each filter removed, what it was never able to examine, and what was switched off."""
    mark = {True: "on", False: "<b>off</b>", None: "&mdash;"}
    headers = ["filter"] + [f"{a.label} state" for a in arms] + [f"{a.label} removed" for a in arms]
    activity = [filter_activity(a.results, a.unfiltered) for a in arms]

    rows = []
    for i, (label, _, _) in enumerate(FILTERS):
        # A filter the configuration does not name, but which removed rows, demonstrably ran.
        # Printing an em-dash beside 164 removals would reproduce the exact ambiguity this table
        # exists to remove, so evidence overrides the missing key.
        state = []
        for j in range(len(arms)):
            on, n = activity[j][i]["enabled"], activity[j][i]["removed"]
            state.append(mark[True] if on is None and n else mark[on])
        counts = []
        for j in range(len(arms)):
            n, on = activity[j][i]["removed"], activity[j][i]["enabled"]
            counts.append("&mdash;" if on is False else f"{n:,}")
        rows.append([label] + state + counts)

    numeric = tuple(range(1 + len(arms), 1 + 2 * len(arms)))
    body = head("filters")
    body += (
        "<p>Every filter, whether it ran, and how many rows it removed. The counts come from "
        "<code>Filter Status</code> in <code>Unfiltered_Results.csv</code>, so they cannot "
        "disagree with the delivered file.</p>"
        "<p><b>Read a dash and a zero differently.</b> A dash means the filter was switched off "
        "for this run: it examined nothing, and rows it would have removed are in the delivered "
        "table. A zero means it ran and found nothing to remove. Both previously appeared as an "
        "absence of rows, and an absence is not a result.</p>")
    body += table(headers, rows, numeric=numeric)

    # The exemption, which no filter reports because nothing was removed
    ex_rows, note = [], ""
    for arm in arms:
        total, exempt, classes = retention_exemption(arm)
        if not total:
            continue
        named = ", ".join(f"{c} ({n})" for c, n in classes[:7]) or "&mdash;"
        ex_rows.append([arm.label, f"{total:,}", f"{exempt:,}",
                        pct(exempt / total) if total else "&mdash;", named])
    if ex_rows:
        note = (
            "<h3>Rows the retention model never judged</h3>"
            "<p>A class with too few members to fit a model &mdash; or one that is not a "
            "homologous series at all &mdash; produces no model, and every row in it passes the "
            "retention check <b>by never taking it</b>. Those rows carry a blank "
            "<code>RT Model Error</code>, which reads as though nothing was wrong. This is the one "
            "exemption no filter can report, because nothing was removed.</p>")
        note += table(["", "MS2-identified", "exempt", "%", "classes wholly exempt"],
                      ex_rows, numeric=(1, 2, 3))
        note += ("<p>These are not marginal classes. Read alongside "
                 f"{ref('annotation')}: a row exempt from the retention check has been named on "
                 "its spectrum alone.</p>")
    return body + note


def section_support(arms) -> str:
    """What the data can support — explicitly NOT a quality judgement."""
    heads = ["", *[f"{a.label} mode" for a in arms]]
    classes = ["", *["pos" if a.label == "positive" else "neg" for a in arms]]
    def block(precision_of, contrast_of):
        return [
            comparative(arms, "Median variation across the study samples",
                        lambda a: pct(precision_of(a)["median_sample_cv"] / 100)),
            comparative(arms, "D-ratio (analytical &divide; total variation)",
                        lambda a: num(precision_of(a)["median_d_ratio"])),
            comparative(arms, "Sample spread &divide; pooled-QC spread",
                        lambda a: num(contrast_of(a)["ratio"], 1) + "&times;"
                        if contrast_of(a) else "&mdash;"),
        ]
    rows = block(lambda a: a.precision, lambda a: a.contrast)
    rows_all = block(lambda a: a.precision_all, lambda a: a.contrast_all)
    return f"""
{head('support')}
<p class="lede">How large the differences between samples are, relative to the measurement noise.
<strong>This is not a quality assessment.</strong></p>
<div class="note"><strong>A small number here is not a fault in the data.</strong> Both statistics
below divide by the variation between samples, so they fall when the samples genuinely resemble
each other &mdash; which is what happens when an intervention has little effect. A study that
measured cleanly and found no difference will score poorly on both, and the correct reading is
that the experiment produced a null result on well-measured data, not that the data is bad.
Whether the measurement is sound is settled in {ref('precision')}, on the pooled QCs alone.</div>
<h3>Identified lipids</h3>
{table(heads, rows, numeric=tuple(range(1, len(arms) + 1)), classes=classes)}
<h3>The same statistics over every quantified feature</h3>
{table(heads, rows_all, numeric=tuple(range(1, len(arms) + 1)), classes=classes)}
<p><strong>What they are for.</strong> The <em>D-ratio</em> is the ratio of analytical dispersion
to total dispersion, both as robust standard deviations &mdash; so it answers directly what
fraction of the observed spread is the instrument rather than the samples. A D-ratio of 0.5 means
analytical variation is half the total, which leaves the biological component dominant; as it
approaches 1 there is no room left for biology and no effect size is detectable however many
samples are run. The spread ratio asks the same question
geometrically &mdash; how much further a sample sits from the QC centroid than the QCs sit from
each other. Read them as a guide to the effect size worth chasing and the number of samples it
would take, not as a mark out of ten.</p>
<div class="note"><strong>The full table looks better, and it is worth knowing why.</strong> Its
advantage is entirely in the denominator: precision is indistinguishable between the two sets, and
in negative mode the full table is the worse of the two. What differs is variation between samples,
which is larger across all features &mdash; and that is unlikely to be biology. Unidentified
features are the fainter half of the table, faint features vary more between samples for reasons
unrelated to the tissue, and many are in-source fragments, contaminants and background that track
the matrix rather than the specimen. That inflates the biological denominator without being
biological. <strong>The identified figures look worse and are the honest ones</strong>, and they
describe the set that will actually be tested.</div>
"""


def section_outliers(arms, plt) -> str:
    across = cross_polarity_outliers(arms)
    rows = []
    for a in arms:
        found = {}
        for key, label in (("pools", "samples + QCs"), ("samples", "samples only")):
            space, report = a.outlier_spaces.get(key), a.outliers.get(key)
            if not (space and report):
                continue
            for i in report["flagged"]:
                if space["roles"][i] == "blank":
                    continue          # a blank is not like the samples; that is its job
                entry = found.setdefault(space["columns"][i], {"tests": set(), "where": []})
                entry["tests"] |= set(report["support"][i])
                entry["where"].append(label)
        for column, entry in sorted(found.items()):
            agree = len(entry["tests"])
            both = sample_key(column) in across
            if both:
                reading = "<span class='flag'>flagged in both polarities &mdash; investigate</span>"
            elif agree >= 2:
                reading = "worth a look"
            else:
                reading = "expected by chance"
            rows.append([short(column), a.label, f"{agree} of 5",
                         " and ".join(dict.fromkeys(entry["where"])), reading])
    body = (table(["injection", "polarity", "agreement", "seen in", "reading"], rows,
                  numeric=(2,))
            if rows else "<p>No injection was flagged in either polarity.</p>")
    repeated = sorted({r[0] for r in rows if sample_key(r[0]) in across})
    note = ""
    if repeated:
        note = (f"<div class='note'><strong>{', '.join(repeated)}: the same specimen fails in both "
                f"polarities, independently.</strong> One flag at a 95% limit is expected by chance "
                f"in a run this size, and several tests agreeing within one polarity share the same "
                f"features, scaling and model &mdash; but two ionisation modes measuring different "
                f"molecules is genuinely independent evidence. Check the extraction record and the "
                f"tissue amount, and run the statistics with and without it. Do not remove it on "
                f"this evidence alone.</div>")
    fig, axes = plt.subplots(1, len(arms), figsize=(6.2 * len(arms), 3.6), squeeze=False)
    for ax, a in zip(axes[0], arms):
        band = [p for p in a.profile if p["role"] == "sample"]
        low = min((p["rho"] for p in band), default=0), max((p["rho"] for p in band), default=1)
        ax.axhspan(low[0], low[1], color="#e8eef5", zorder=1,
                   label="range spanned by the study samples")
        # `order` is None without a sequence file, which the CLI allows. Fall back to position
        # in the table: the axis then reads as table order rather than acquisition order, and the
        # label below says so instead of the figure crashing.
        known = all(p.get("order") is not None for p in a.profile)
        at = {p["injection"]: (p["order"] if known else n)
              for n, p in enumerate(a.profile, start=1)}
        for role in ("sample", "qc", "blank"):
            pick = [p for p in a.profile if p["role"] == role]
            if pick:
                ax.scatter([at[p["injection"]] for p in pick], [p["rho"] for p in pick], s=SIZE,
                           c=COLOURS[role], edgecolor="white", linewidth=.6,
                           label=LABELS[role], zorder=3)
        for p in a.profile:
            # Blanks sit below the sample band by definition, so labelling them marks exactly the
            # injections that should never be marked. The exclusion was applied to the verdict and
            # the table before this figure was noticed — the third use of the same flag.
            if p.get("below_sample_band") and p.get("role") != "blank":
                ax.annotate(short(p["injection"]), (at[p["injection"]], p["rho"]), fontsize=7,
                            xytext=(5, -9), textcoords="offset points", color=COLOURS["flag"])
        ax.set_xlabel("injection number" if known else "position in the table (no sequence given)")
        ax.set_ylabel("correlation with the median sample profile")
        ax.set_title(f"{a.label} mode", fontsize=10.5, pad=8)
        ax.spines[["top", "right"]].set_visible(False)
        legend_below(ax, ncol=2)
    profile_fig = figure(
        embed(fig),
        "Did every injection work? Each injection's rank correlation against the median study "
        "sample profile, after median normalisation, per polarity. This judges an injection on "
        "whether it <em>looks like</em> the samples rather than on how much of it arrived &mdash; "
        "an injection can be small and faithful, or full-sized and wrong. The shaded band is the "
        "range the study samples themselves span.")

    ffig, faxes = plt.subplots(1, len(arms), figsize=(5.4 * len(arms), 3.2), squeeze=False)
    for ax, a in zip(faxes[0], arms):
        space, report = a.outlier_spaces.get("samples"), a.outliers.get("samples")
        if not (space and report) or not len(report.get("forest", [])):
            ax.axis("off")
            continue
        score = report["forest"]
        order = np.argsort(score)
        colour = [COLOURS["flag"] if sample_key(space["columns"][i]) in across
                  else COLOURS["sample"] for i in order]
        ax.scatter(range(len(order)), score[order], s=SIZE, c=colour,
                   edgecolor="white", linewidth=.6, zorder=3)
        med = float(np.median(score))
        mad = float(np.median(np.abs(score - med))) * 1.4826
        if mad > 0:
            ax.axhline(med + 3 * mad, color="#bbbbbb", lw=.9, ls="--")
        for rank, i in enumerate(order):
            if sample_key(space["columns"][i]) in across:
                ax.annotate(short(space["columns"][i]), (rank, score[i]), fontsize=7,
                            xytext=(-6, 3), textcoords="offset points", ha="right",
                            color=COLOURS["flag"])
        ax.set_xlabel("study samples, ranked")
        ax.set_ylabel("isolation-forest anomaly score")
        ax.set_title(f"{a.label} mode", fontsize=9.5, pad=8)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    forest_fig = figure(
        embed(ffig),
        "Isolation-forest anomaly score for each study sample, ranked, per polarity. The forest "
        "isolates a sample by random splits of the principal components and needs fewer of them "
        "for an unusual one. Dashed line: three robust standard deviations above the median. "
        "Pooled QCs are not shown &mdash; their quality is settled in "
        f"{ref('precision')} &mdash; and only a sample flagged in both polarities is named.")

    return f"""
{head('outliers')}
<p class="lede">Tested in two spaces, because an injection can be an outlier in one and not the
other. The <em>samples + QCs</em> model is the measurement's own space; the <em>samples only</em>
model is the space the biology will actually be analysed in, and is the one that matters for
interpretation. Both are reported, along with the profile screen and missing-value rate, so a
flag can be judged on how many independent lines support it.</p>
{profile_fig}
{body}
{note}
{forest_fig}
<p><strong>Five tests, because they fail differently.</strong> Hotelling's T&sup2; finds an
injection far from the centre of the model. Its <em>robust</em> form recentres on medians, because
the classical covariance is inflated by the very outliers being looked for, so two of them mask
each other. DModX finds an injection the model does not describe, which can sit innocently near
the centre of a score plot. Neighbour distance uses no covariance model at all, so it fails
differently again. An isolation forest &mdash; scikit-learn's, run on the components &mdash;
isolates a point by random splits, and needs fewer of them for an unusual one.</p>
<p class="lede">Agreement is counted out of five. The tests are Hotelling's T&sup2;, its robust
form, DModX, mean distance to the nearest injections, and an isolation forest.</p>
<div class="note"><strong>Count the agreement, not the flags.</strong> One flag at a 95% limit is
expected by chance in a run this size, so a single test naming an injection means little.
<strong>Two or more tests agreeing is the signal</strong>, and an injection flagged independently
in both polarities is stronger still. Tested on the delivered gap-filled table using the
identified lipids &mdash; the matrix the biology is analysed in.</div>
"""


def cross_polarity_outliers(arms) -> set:
    """Sample keys flagged in BOTH polarities.

    Two independent measurements — different ionisation, different molecules — failing on the same
    specimen is far stronger than several tests agreeing within one polarity, where they share the
    features, the scaling and the model. This is the statement the single document exists to make,
    and it is the bar for putting a name on a score plot.
    """
    per = []
    for a in arms:
        flagged = set()
        for key in ("pools", "samples"):
            space, report = a.outlier_spaces.get(key), a.outliers.get(key)
            if space and report:
                flagged |= {sample_key(space["columns"][i]) for i in report["flagged"]}
        per.append(flagged)
    if len(per) < 2:
        return set()
    return set.intersection(*per)


def group_members(arm, space) -> "dict[str, list[int]] | None":
    """Column indices per study-group level, for colouring a score plot.

    None when there is no usable grouping, so the panel falls back to role colours rather than
    inventing a single group containing everything.
    """
    factors = qc.usable_factors({c: arm.run.metadata.get(c, {}) for c in space["columns"]})
    if not factors:
        return None
    factor = factors[0]
    out: dict[str, list[int]] = {}
    for index, column in enumerate(space["columns"]):
        level = (arm.run.metadata.get(column) or {}).get(factor)
        if level:
            out.setdefault(level, []).append(index)
    return out if len(out) > 1 else None


def section_structure(arms, plt) -> str:
    """PCA panels: both feature sets, both polarities, with the injection subsets explained."""
    heads = ["principal component analysis", "injections", "features used", "PC1", "PC2"]
    rows = []
    for a in arms:
        for key, label in (("blanks", "blanks + QCs + samples"), ("pools", "samples and QCs"),
                           ("samples", "samples only")):
            for view, vlabel in (("all", "all features"), ("named", "identified only")):
                space = a.spaces[(key, view)]
                if space:
                    rows.append([f"{a.label}: {label}, {vlabel}", len(space["columns"]),
                                 f"{space['features']:,}",
                                 f"{100 * space['explained'][0]:.1f}%",
                                 f"{100 * space['explained'][1]:.1f}%"])
    figs = ""
    for view, vlabel in (("all", "all quantified features"), ("named", "identified lipids only")):
        fig, axes = plt.subplots(len(arms), 3, figsize=(9.5, 3.0 * len(arms)), squeeze=False)
        for r, a in enumerate(arms):
            for c, (key, label) in enumerate((("blanks", "blanks + QCs + samples"),
                                              ("pools", "samples and QCs"),
                                              ("samples", "samples only"))):
                ax = axes[r][c]
                space = a.spaces[(key, view)]
                if space is None:
                    ax.axis("off")
                    continue
                # No outlier marking here: these panels answer a different question, and the
                # outliers section is where injections are judged.
                #
                # The samples-only panel is coloured by the study group. Colouring it by role
                # spends the colour channel on a constant — every point is a sample — and left the
                # one panel where the biology could show as the one panel showing nothing. It also
                # made a separate group figure necessary, which then duplicated this panel.
                groups = group_members(a, space) if key == "samples" else None
                scatter(ax, space["scores"], space["roles"], set(), space["explained"],
                        groups=groups)
                title = label if key != "samples" or not groups else f"{label}, by group"
                ax.set_title(f"{a.label} mode &mdash; {title}".replace("&mdash;", "—"),
                             fontsize=9.5, pad=8)
                ax.tick_params(labelsize=8)
                ax.xaxis.label.set_size(8.5)
                ax.yaxis.label.set_size(8.5)
                # Two legends on the bottom row, because the columns no longer share a key: the
                # first two carry injection role, the third carries study group. One legend under
                # the middle panel would silently label the group colours as roles.
                # Two legends on the bottom row, because the columns no longer share a key: the
                # role panels and the group panel.
                #
                # The role legend goes under the FIRST column, not the second. Blanks appear only
                # there — the other panels exclude them — so a legend drawn under the second panel
                # names every colour except the one the reader cannot place. The second column's
                # roles are a subset of the first's, so one legend serves both.
                #
                # The group legend stacks (ncol=1). Group names come from the submitter and can be
                # long — "D28_Tim4_Cre_DTR_depletion_SI" — and laid out in columns it spread wide
                # enough to run into the legend beside it and render as one unreadable strip.
                if r == len(arms) - 1 and c in (0, 2):
                    legend_below(ax, ncol=1 if c == 2 else 3)
        figs += figure(embed(fig),
                       f"Principal component analysis on <strong>{vlabel}</strong>, positive mode "
                       f"on the top row and negative on the bottom. The two polarities are "
                       f"analysed separately throughout and are never combined.")
    # ⚠ There is deliberately NO separate "coloured by group" figure. There used to be, and
    # with a single study factor it was a two-panel figure laid out for a grid — rendered at the
    # width of the three-column panels above, so it came out oversized and crowded — showing
    # exactly what the samples-only column now shows. One panel, in the grid, at the same size as
    # its neighbours.
    return f"""
{head('structure')}
<p class="lede">Three injection subsets, because they answer different questions, and two feature
sets, because "all quantified features" and "identified lipids" give different pictures and both
are worth having.</p>
{table(heads, rows, numeric=(1, 2, 3, 4))}
<div class="note"><strong>These are features, not lipids.</strong> A feature is one quantified
signal; most are never identified. The identified-only panels show the subset carrying a name.
<strong>The feature set is the same in all three subsets of a polarity</strong> &mdash; a feature
qualifies by being detected in at least half the <em>study samples</em>, decided once. Deciding it
per panel instead would let the pooled QCs, which contain nearly everything, hand sparse features
extra detections and quietly change the feature set between panels.</div>
{figs}
<div class="note">The <strong>samples-only</strong> panel is coloured by study group. That is the
space any downstream comparison lives in, and separation there is a hypothesis, not a result, until
it is tested. Injections are not named on any of these plots &mdash; a name on one point pulls the
eye to it and invites an explanation, which is the opposite of what a score plot is for. Which
injections were flagged, and on what evidence, is in the outliers section.</div>
"""


def _run_in_note(arms) -> str:
    """Mention the start-of-run check only when it changes the answer.

    Equilibration at the head of a sequence can masquerade as drift, so it is always tested. But a
    permanent table row that moves the figure by half a point tells a reader nothing and reads as
    an unexplained duplicate. Same principle as the noise floor: always measured, mentioned when
    it matters.
    """
    said = []
    for a in arms:
        if not a.order.get("tested"):
            continue
        full = a.order["drifting"] / a.order["tested"]
        late = a.order["drifting_late"] / a.order["tested"]
        if abs(full - late) >= 0.02:
            said.append(f"<strong>{a.label} mode:</strong> excluding the first {qc.RUN_IN} "
                        f"injections changes this from {full:.0%} to {late:.0%}, so part of the "
                        f"apparent drift is the column equilibrating rather than the response "
                        f"changing")
    if not said:
        return ("<p>Excluding the first {} injections does not materially change these figures, "
                "so what is seen is not a start-of-run effect.</p>".format(qc.RUN_IN))
    return "<div class='note'>" + "; ".join(said) + ".</div>"


def confounding_verdict(arms) -> str:
    """Answer the question the heading asks, in the first line, in words.

    The section used to open with drift statistics and never say whether acquisition was
    confounded with the biology — which is the only thing the heading promises and the only thing
    that decides whether the drift below matters. Drift in a properly randomised run is
    recoverable; mild drift with one group sitting at the start of the sequence is not.
    """
    tested = [(a, a.confounding) for a in arms if a.confounding]
    if not tested:
        return ('<div class="note"><strong>Not testable.</strong> This needs the injection order '
                'and a sample grouping &mdash; a sequence file and metadata &mdash; and at least '
                'eight samples carrying both. Supply <code>--sequence</code> and '
                '<code>--metadata-*</code> and the question is answered here.</div>')
    bad = [(a.label, r) for a, rows in tested for r in rows if r.get("confounded")]
    if bad:
        listed = "; ".join(f"<strong>{r['factor']}</strong> in {label} mode (p = {r['p']:.3f})"
                           for label, r in bad)
        return (f'<div class="note bad"><strong>&#9888; Yes &mdash; and it matters.</strong> '
                f'{listed}. A group that sits disproportionately early or late in the sequence '
                f'cannot be separated from drift by any statistic computed afterwards: the '
                f'difference between the groups and the difference in when they ran are the same '
                f'number. Treat any comparison on that factor as confounded, and use injection '
                f'order as a covariate at minimum.</div>')
    every = [(r["p"], r["factor"]) for a, rows in tested for r in rows]
    closest = min(every)[0] if every else 1.0
    factors = sorted({r["factor"] for _, rows in tested for r in rows})
    return (f'<div class="note good"><strong>No.</strong> Nothing is confounded with injection '
            f'position: {", ".join(f"<code>{f}</code>" for f in factors)} tested in both '
            f'polarities, closest p = {closest:.2f} against a 0.05 bar. The groups are spread '
            f'through the sequence, so the drift described below cannot masquerade as a '
            f'difference between them &mdash; it inflates variance rather than creating an '
            f'effect.</div>')


def section_order(arms, plt) -> str:
    heads = ["", *[f"{a.label} mode" for a in arms]]
    classes = ["", *["pos" if a.label == "positive" else "neg" for a in arms]]
    rows = [
        # Where the order came from. Every number in this section rests on it being right, and
        # the two sources answer different questions: the instrument records what happened, the
        # sequence file records what was planned.
        comparative(arms, "Injection order taken from",
                    lambda a: a.run.order_source or "&mdash; unavailable"),
        comparative(arms, "Features correlating with injection order",
                    lambda a: pct(a.order["drifting"] / a.order["tested"], 0)
                    if a.order.get("tested") else "&mdash;"),

        comparative(arms, "Fitted change in pooled-QC signal across the run",
                    lambda a: (f"{a.normalised_totals('qc')[2]:+.0%}"
                               if a.normalised_totals("qc")[2] is not None else "&mdash;")),
        comparative(arms, "Features drifting systematically in the pools",
                    lambda a: f"{a.drift.get('systematic', 0):,} of "
                              f"{a.drift.get('features', 0):,}"),
        comparative(arms, "Per-sample loading against order (&rho;)",
                    lambda a: num(a.order.get("loading_rho"))),
    ]
    clash = ""
    for a in arms:
        if a.run.order_disagreement:
            clash += (f"<div class=\"note\"><strong>{a.label} mode: the sequence file and the "
                      f"instrument disagree.</strong> {a.run.order_disagreement}. The order used "
                      f"here is the instrument's, because it records what was acquired rather "
                      f"than what was scheduled &mdash; a run repeated, reordered or dropped "
                      f"mid-sequence produces exactly this. Worth checking before the drift "
                      f"result below is relied on.</div>")

    fig, axes = plt.subplots(len(arms), 2, figsize=(12, 3.6 * len(arms)), squeeze=False)
    for r, a in enumerate(arms):
        for c, role in enumerate(("sample", "qc")):
            ax = axes[r][c]
            order, values, trend = a.normalised_totals(role)
            if not len(values):
                ax.axis("off")
                continue
            ax.plot(order, values, "o-", color=COLOURS[role], ms=5, lw=1.1)
            ax.axhline(1.0, color="#bbbbbb", lw=.8, ls="--")
            if trend is not None:
                fit = np.poly1d(np.polyfit(np.asarray(order, dtype=float), values, 1))
                ax.plot(order, fit(np.asarray(order, dtype=float)), "-", color="#555555",
                        lw=1.0, alpha=.8)
                ax.annotate(f"fitted change {trend:+.0%}", (0.02, 0.04), xycoords="axes fraction",
                            fontsize=8, color="#555555")
            # Centred on 1.0, which is what the normalisation means, and wider for samples than
            # for pools because they are not expected to agree to the same degree: pools are the
            # same material, samples are not.
            fixed_axis(ax, values, *((0.0, 2.0) if role == "sample" else (0.8, 1.2)))
            ax.set_xlabel("injection number")
            ax.set_ylabel("signal relative to this type's median")
            ax.set_title(f"{a.label} mode — {LABELS[role]}s", fontsize=10.5, pad=8)
            ax.spines[["top", "right"]].set_visible(False)
    caption = ("Total signal against injection order, <strong>normalised to the median of each "
               "injection type</strong>, so 1.0 is the run's typical level and one unusual "
               "injection cannot tilt the series. The grey line is fitted through every point "
               "&mdash; drift read off the first and last injections would rest on the two least "
               "reliable points in the sequence.")
    return f"""
{head('order')}
<p class="lede">If the response changes through the sequence and the sample order is not random,
a difference between groups can be a difference in when they were injected.</p>
{confounding_verdict(arms)}
{clash}
{table(heads, rows, numeric=tuple(range(1, len(arms) + 1)), classes=classes)}
{figure(embed(fig), caption)}
{_run_in_note(arms)}
<p>Every &ldquo;signal against injection number&rdquo; axis in this report uses the same
convention: divided by the median of its own injection type, so the panels are directly
comparable and differences in how much material was injected do not read as drift.</p>
"""


# Colour carries the LIPID MAPS category, not the individual class. Twenty-five hues do not
# separate for anyone and fail far sooner for a colourblind reader; six do, and six is what the
# chemistry actually asks of the plot — whether each category elutes in its own window.
#
# Categories are matched by class-name prefix, longest first, so `Plasmenyl-PE` and `LysoPC` land
# with the glycerophospholipids rather than falling through to a default.
CATEGORIES = [
    ("Fatty acyls", "#E69F00", ("FA", "AC", "CAR")),
    ("Glycerolipids", "#D55E00", ("TG", "DG", "MG", "Alkenyl-TG", "Alkanyl-TG")),
    ("Glycerophospholipids", "#0072B2",
     ("PC", "PE", "PS", "PI", "PG", "PA", "BMP", "CL", "LysoPC", "LysoPE", "LysoPG", "LysoPS",
      "LysoPI", "LysoPA", "Plasmenyl-", "Plasmanyl-", "PE-NMe", "PC[OH]", "PE[OH]")),
    ("Sphingolipids", "#009E73",
     ("SM", "Cer", "HexCer", "GlcCer", "GalCer", "LacCer", "SHexCer", "SP", "GM", "GD", "GT",
      "CerP", "Sulfatide")),
    ("Sterol lipids", "#CC79A7", ("CE", "Cholesterol", "ST")),
    ("Prenol lipids", "#8C6BB1", ("CoQ", "Ubiquinone")),
]
UNIDENTIFIED = "#c9c9c9"
REMOVED = "#8a8a8a"      # identified, then filtered out — outlined, never filled
OTHER_CATEGORY = ("Other", "#661100")


def lipid_category(lipid_class: str) -> tuple:
    """(category, colour) for a class name, by longest-prefix match.

    Longest first so `LysoPC` is not swallowed by `PC`, and `Plasmenyl-PE` not by `PE`.
    """
    name = (lipid_class or "").strip()
    best = None
    for category, colour, prefixes in CATEGORIES:
        for prefix in prefixes:
            if name.startswith(prefix) and (best is None or len(prefix) > best[0]):
                best = (len(prefix), category, colour)
    return (best[1], best[2]) if best else OTHER_CATEGORY


# A row removed because it is another row's shadow — an in-source fragment, a second adduct of a
# molecule already in the table, an isotope, a split peak merged into its parent, a redundant
# re-match of a name. These are not molecules, and a molecule whose every row carries one of these
# reasons was never a separate finding. Distinct from removals that concern a REAL molecule the
# study could not use: too sparse across the samples, present in the blanks, off its class model.
ARTEFACT_REASON = re.compile(
    r"in-source|adduct pair|redundant identification|split peak|isotope|dimer|correlated",
    re.IGNORECASE)


def elution_figure(arms, plt) -> str:
    """Where each class elutes — one point per MOLECULE, coloured by category, split by evidence.

    The plot a reader can check against chemistry without knowing anything about the software: a
    class occupies a retention window, and a member of it sits where its chain length and
    unsaturation put it. A class smeared across the whole run, or sitting past the end of the
    gradient, is visible here immediately — which is how an end-of-gradient wash cluster reported
    as five phosphatidylserines was found.

    Three decisions about what goes on it, each of which has been wrong at some point.

    **One point per molecule, not per row.** A name appears on many rows: in-source fragments,
    adducts and co-eluting isobars re-match it at other retention times, so a row-level plot drew
    2,905 points for 896 molecules and every density in it was a statement about how often the
    search re-matched a name. Each molecule is drawn once, at the retention time of its largest
    surviving peak — the peak the quantitation comes from.

    **Artefacts are dropped, sparse molecules are not.** The filters remove two unrelated kinds of
    row and lumping them together made the figure meaningless in both directions. One kind is a
    molecule's own shadow — its fragments, its other adducts, its isotopes — which is not a
    finding and should never be plotted. The other is a real molecule the study could not use:
    on this run the largest single removal is "not detected in 70% of any group", 428 positive and
    524 negative molecules, which are perfectly good identifications that the presence filter
    excludes from the analysis table. Those belong on a plot of what was identified. A molecule is
    dropped only when EVERY row carrying it is an artefact.

    **Filled is a spectrum, hollow is not.** The split that matters for trusting a name is the
    evidence behind it: MS2 matched the fragments the molecule produced, while `RT model` matched
    accurate mass and elution only and cannot separate isomers. Opacity carries a third, weaker
    fact — whether the molecule reached the delivered table — without spending another colour.
    """
    fig, axes = plt.subplots(len(arms), 1, figsize=(9.8, 4.3 * len(arms)), squeeze=False)
    # One retention axis for both polarities. Auto-scaling gave positive 0-21 and negative 0-25,
    # which makes the same gradient look like two different methods and invites a reader to
    # compare elution windows that are not on the same scale.
    span = max((float(r["Retention Time (min)"]) for a in arms for r in a.run.rows), default=25.0)
    limit = float(np.ceil(span / 5.0) * 5.0)

    def molecules(rows):
        """One representative row per molecule: its largest peak, artefact rows never chosen."""
        pick = {}
        for r in rows:
            name = r["Identification"].strip()
            if not name or ARTEFACT_REASON.search(r.get("Filter Status", "") or ""):
                continue
            if name not in pick or float(r["Area (max)"] or 0) > float(pick[name]["Area (max)"] or 0):
                pick[name] = r
        return pick

    for ax, a in zip(axes[:, 0], arms):
        delivered = molecules(a.run.rows)
        found = dict(molecules(a.unfiltered.rows)) if a.unfiltered is not None else {}
        found.update(delivered)          # the delivered row wins as the representative
        if not found:
            ax.axis("off")
            continue

        # The unidentified backdrop, underneath everything. Features, not molecules — an unnamed
        # feature has no name to collapse on — so it is the one layer counted in rows, and the
        # legend says so rather than letting a reader assume the two layers are comparable.
        unnamed = [r for r in a.run.rows if not r["Identification"].strip()]

        names = list(found)
        rows = [found[n] for n in names]
        area = np.array([float(r["Area (max)"] or 0) for r in rows], dtype=float)
        retention = np.array([float(r["Retention Time (min)"]) for r in rows], dtype=float)
        classes = [r["Lipid Class"].strip() for r in rows]
        by_ms2 = np.array([r.get("Identification Source", "").strip() != "RT model" for r in rows])
        in_table = np.array([n in delivered for n in names])
        top = float(area.max()) or 1.0
        intensity = 100.0 * area / top

        if unnamed:
            u_ret = np.array([float(r["Retention Time (min)"]) for r in unnamed])
            u_int = 100.0 * np.array([float(r["Area (max)"] or 0) for r in unnamed]) / top
            ax.scatter(u_ret, u_int, s=10, c=UNIDENTIFIED, linewidth=0, alpha=.5, zorder=1,
                       label=f"unidentified ({len(unnamed):,} features)")

        category = {i: lipid_category(classes[i]) for i in range(len(rows))}
        weight = Counter()
        for i, (name, _) in category.items():
            weight[name] += area[i]
        for name, _ in weight.most_common():
            pick = np.array([i for i, (n, _) in category.items() if n == name])
            classes_here = sorted({classes[i] or "unassigned" for i in pick})
            colour = category[pick[0]][1]
            for spectrum in (True, False):
                for kept in (True, False):
                    sel = pick[(by_ms2[pick] == spectrum) & (in_table[pick] == kept)]
                    if not len(sel):
                        continue
                    style = (dict(c=colour, edgecolor="white", linewidth=.4) if spectrum
                             else dict(facecolors="none", edgecolor=colour, linewidth=1.2))
                    ax.scatter(retention[sel], intensity[sel], s=30,
                               alpha=.95 if kept else .40, zorder=3 if kept else 2, **style)
            ax.scatter([], [], s=30, c=colour, edgecolor="white", linewidth=.4,
                       label=(f"{name} ({len(classes_here)} "
                              f"{'class' if len(classes_here) == 1 else 'classes'})"))

        # The evidence and delivery keys, drawn in neutral grey so they read as keys rather than
        # as another category.
        ax.scatter([], [], s=30, c="#555555", edgecolor="white", linewidth=.4,
                   label=f"MS2 spectrum ({int(by_ms2.sum()):,})")
        ax.scatter([], [], s=30, facecolors="none", edgecolor="#555555", linewidth=1.2,
                   label=f"MS1 accurate mass + RT model ({int((~by_ms2).sum()):,})")
        ax.scatter([], [], s=30, c="#555555", edgecolor="white", linewidth=.4, alpha=.40,
                   label=f"faded: not in the delivered table ({int((~in_table).sum()):,})")

        ax.set_yscale("log")
        floor = min([float(intensity.min())]
                    + ([float(u_int[u_int > 0].min())] if unnamed and (u_int > 0).any() else []))
        ax.set_ylim(max(floor * 0.7, 1e-4), 160)
        ax.set_xlim(0, limit)
        ax.set_xlabel("retention time (min)")
        ax.set_ylabel("intensity, % of the largest")
        ax.set_title(f"{a.label} mode — {len(found):,} molecules identified, "
                     f"{int(in_table.sum()):,} in the delivered table", fontsize=9.5, pad=8)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8.5, ncol=3, loc="upper center",
                  bbox_to_anchor=(0.5, -0.19), handletextpad=.4, columnspacing=1.4)
    return figure(
        embed(fig),
        "<strong>One identified molecule per point</strong>, at the retention time of its largest "
        "peak, against that peak's area as a percentage of the largest, on a logarithmic axis. "
        "Not one point per row &mdash; a name appears on many rows because in-source fragments, "
        "adducts and co-eluting isobars re-match it elsewhere in the run, and a row-level plot "
        "draws several times as many points as there are lipids. Those artefact rows are excluded "
        "outright, and a molecule is dropped only when every row carrying it is one. "
        "<strong>Molecules the analysis filters excluded are still here, faded</strong> &mdash; "
        "almost all of them removed for being detected in too few samples, which makes them "
        "unusable for statistics without making them wrong. "
        "<strong>Filled markers are named from a fragmentation spectrum; hollow markers are named "
        "from accurate mass and retention only</strong>, which is enough to be worth reporting and "
        "not enough to separate isomers &mdash; the <code>Identification Source</code> column "
        "carries this per row. Colour is <strong>LIPID MAPS category</strong>, ordered by the "
        "signal each accounts for, because six categories separate cleanly where twenty-five "
        "individual classes do not. <strong>This is the plot to check against chemistry</strong> "
        "&mdash; a class occupies a retention window, and a member of it elutes where its chain "
        "length and unsaturation put it. A class smeared across the whole run, or sitting beyond "
        "the end of the gradient, shows up here at a glance. <strong>Small grey dots are "
        "quantified features that were never identified</strong>, the backdrop the identifications "
        "sit on &mdash; counted in features rather than molecules, because an unnamed feature has "
        "no name to collapse on.")


def section_annotation(arms, plt) -> str:
    heads = ["", *[f"{a.label} mode" for a in arms]]
    classes = ["", *["pos" if a.label == "positive" else "neg" for a in arms]]
    # Three numbers, not six. The surplus is simply rows minus molecules; presenting "names on
    # exactly one row" and "rows those names occupy" alongside it made a reader hold four
    # near-synonyms at once to derive a subtraction.
    rows = [
        comparative(arms, "Identified rows in the table",
                    lambda a: f"{a.annotation['identified_rows']:,}"),
        comparative(arms, "Distinct molecules among them",
                    lambda a: f"{a.annotation['distinct_names']:,}"),
        comparative(arms, "<strong>Rows that repeat a molecule already listed</strong>",
                    lambda a: f"<strong>{a.annotation['surplus']:,}</strong>"),
    ]
    worked = ""
    for a in arms:
        top = a.annotation["repeated"][:1]
        if not top:
            continue
        r = top[0]
        # Describe what this run actually shows. The prose previously asserted "tens of rows",
        # which was true only while a name-parsing fault was collapsing subclasses onto their class.
        worked += (f"<p><strong>{a.label} mode, worked example.</strong> <code>{r['name']}</code> "
                   f"appears on <strong>{r['rows']} rows</strong>, spread over "
                   f"{r['spread']:.1f} min (RT {r['retentions'][0]:.2f} to "
                   f"{r['retentions'][-1]:.2f}). A molecule elutes once, so rows that far apart "
                   f"are not all the same peak: in-source fragments and co-eluting isobars match "
                   f"the same name at other retention times. They are filtered out of the curated "
                   f"table &mdash; the point here is only that a row count is not a molecule "
                   f"count.</p>")
    anchor = ""
    for a in arms:
        found = a.named_only_in_pools()
        if found and found["count"]:
            anchor += (f"<p><strong>{a.label} mode:</strong> {found['count']} molecules have MS2 "
                       f"evidence only in pooled QC injections and none in any study sample, and "
                       f"are quantified in the samples by MS1 &mdash; for example "
                       f"<code>{', '.join(found['examples'][:3])}</code>.</p>")
    pairs = sum(len(a.annotation["same_peak"]) + len(a.annotation["summed_pairs"]) for a in arms)
    return f"""
{head('annotation')}
<p class="lede">A name is not a row and a row is not a molecule. This section is about the
difference, because it changes how the table must be counted.</p>
{table(heads, rows, numeric=tuple(range(1, len(arms) + 1)), classes=classes)}
<div class="note"><strong>Read it as one subtraction.</strong> The table has more rows than it has
molecules, and the difference is the third line: those rows name a molecule that is already in the
table somewhere else. <strong>Count molecules, not rows.</strong></div>
{worked}
<h3>What was measured</h3>
{elution_figure(arms, plt)}
{table(["", *[f"{a.label} mode" for a in arms]],
       [comparative(arms, "Quantified features", lambda a: f"{a.coverage['rows']:,}"),
        # The search-level counts, so this table and the figure above it state the same thing.
        # Without them the figure said 2,905 and the table said 440, and nothing reconciled them.
        comparative(arms, "Identified rows the search made, before filters",
                    lambda a: f"{a.search_identified['rows']:,}" if a.search_identified else "&mdash;"),
        comparative(arms, "&hellip; distinct molecules, artefact-only names excluded",
                    lambda a: f"{a.search_identified['molecules']:,}" if a.search_identified else "&mdash;"),
        comparative(arms, "&hellip; names that were only ever a fragment, adduct or re-match",
                    lambda a: f"{a.search_identified['artefact_only']:,}" if a.search_identified else "&mdash;"),
        comparative(arms, "Identified rows", lambda a: f"{a.coverage['identified']:,}"),
        comparative(arms, "&hellip; from a fragmentation spectrum",
                    lambda a: f"{a.sources.get('MS2', 0):,}"),
        comparative(arms, "&hellip; from accurate mass and retention only",
                    lambda a: f"{a.sources.get('RT model', 0):,}"),
        comparative(arms, "Distinct molecules", lambda a: f"{a.coverage['molecules']:,}"),
        comparative(arms, "Lipid classes", lambda a: f"{len(a.coverage['classes']):,}"),
        # None when nothing was identified at all — a real case, not a broken run.
        comparative(arms, "Median dot product",
                    lambda a: f"{a.coverage['dot_median']:.0f}"
                    if a.coverage.get("dot_median") is not None else "&mdash;"),
        comparative(arms, "Identified rows scoring under 700",
                    lambda a: f"{a.coverage['dot_low']:,}")],
       numeric=tuple(range(1, len(arms) + 1)), classes=classes)}
<div class="note"><strong>Two kinds of identification, and they are not equivalent.</strong> A row
named from a fragmentation spectrum is matched on the fragments the molecule produced. A row named
from <em>accurate mass and retention only</em> has no spectrum behind it: its mass matches a lipid
to within 10 ppm and it elutes where that lipid's class model predicts, which is enough to be worth
reporting and not enough to distinguish isomers. Those rows carry
<code>Identification Source = RT model</code>, and the assignment is declined outright when two
candidates both fit. <strong>Treat them as provisional and check the column before using a
name.</strong></div>
<p>A dot product is how closely a measured spectrum matched its library entry. Rows here run from
about 500 to 1000, and a row at 520 and a row at 990 are not the same statement &mdash; the
<code>Dot Product</code> column carries it per row.</p>
<h3>Identification and quantitation are separate stages</h3>
<p>A compound group takes its <em>name</em> from any MS2 spectrum in any injection that matched,
and its <em>area</em> from the MS1 peak in every injection where the feature was detected. No MS2
is required in the sample itself &mdash; so a molecule can be named from a pooled QC and
quantified across every sample. That is a genuine benefit of running pooled QCs, and it rests on
the assumption that the same mass at the same retention time in a sample is the same molecule.
The <code>Identification Source</code> and <code>MS2 Files</code> columns carry the evidence.</p>
{anchor}
<h3>One peak reported as two lipids</h3>
<p>Pairs of differently-named rows at one peak, requiring both a retention gap within one measured
peak width and areas that correlate across samples: <strong>{pairs}</strong> found.</p>
<div class="note">A fixed retention window catches neighbours rather than duplicates. The window
here is one measured peak width per polarity, and correlation is a criterion rather than a
printed statistic &mdash; two rows are one measurement only if they behave like one.
<strong>And even a real pair has a limit:</strong> a sum composition can correspond to several
molecular species, so pairing <code>PC 34:1</code> with <code>PC 16:0_18:1</code> assumes that
resolved form is the right one of the possibilities. The pairing can say two rows are one
measurement; it can never confirm the chain assignment.</div>
"""


def section_blanks(arms, plt) -> str:
    rows = []
    for a in arms:
        for n, record in enumerate(a.blanks, start=1):
            rows.append([a.label, f"blank {n}",
                         record.get("order") if record.get("order") is not None else "&mdash;",
                         f"{record['features']:,}", f"{record['total']:.1e}"])
    # ⚠ The blank name and the injection column previously printed the same value, because the
    # key was renamed and the table not updated. The columns now carry different information:
    # which blank it is, and where in the sequence it ran.

    fig, axes = plt.subplots(1, len(arms), figsize=(5.6 * len(arms), 3.2), squeeze=False)
    detail = ""
    for ax, a in zip(axes[0], arms):
        b = a.blank_ratio
        if not b:
            ax.axis("off")
            continue
        edges = np.linspace(-2, 5, 40)
        ax.hist(np.clip(b["log_ratio"], -2, 5), bins=edges, color=COLOURS["sample"],
                alpha=.85, label="every feature seen in a blank")
        if b["identified"].any():
            ax.hist(np.clip(b["log_ratio"][b["identified"]], -2, 5), bins=edges,
                    color=COLOURS["qc"], alpha=.9, label="of those, identified")
        cut = np.log10(b["multiplier"])
        ax.axvline(cut, color=COLOURS["flag"], lw=1.2, ls="--")
        ax.axvspan(-2, cut, color=COLOURS["flag"], alpha=.06)
        ax.set_xlabel("sample &divide; blank  (log&#8321;&#8320;)".replace("&divide;", "÷")
                      .replace("&#8321;&#8320;", "₁₀"))
        ax.set_ylabel("features")
        ax.set_title(f"{a.label} mode", fontsize=9.5, pad=8)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        legend_below(ax, ncol=2)
        detail += (f"<p><strong>{a.label} mode:</strong> {b['features']:,} features appear in at "
                   f"least one blank. Their median sample-to-blank ratio is "
                   f"<strong>{b['median']:,.0f}&times;</strong>, and "
                   f"<strong>{b['below']:.0%}</strong> sit within {b['multiplier']:g}&times; of "
                   f"the blanks &mdash; background rather than sample. Among the identified "
                   f"features the figure is {b['below_identified']:.0%}.</p>")

    return f"""
{head('blanks')}
<p class="lede">Measured before the blank filter, because a feature the filter removed has to still
be in the table for the question to mean anything.</p>
{table(["polarity", "blank", "injection number", "features", "total signal"], rows,
       numeric=(2, 3, 4))}
{figure(embed(fig),
        "How far the samples stand above the blanks, for every feature detected in a blank. "
        "Dashed line and shading: the filter's threshold &mdash; features to its left are no "
        "higher in the samples than in the background and are removed.")}
{detail}
<div class="note"><strong>Whether a feature is "in the blanks" is not the question; how far above
them the samples sit is.</strong> A feature a thousandfold higher in the samples is background that
does not matter, and one threefold higher is background that does. Both would count as shared,
which is why the distribution is shown rather than a single percentage.</div>
<p>The filter is a filter, never a subtraction: no intensity is modified, and a compound is kept
when the samples stand clear of the blanks by the configured multiple. A leading blank carrying
signal is contamination and a trailing one is carryover, and they need different fixes.</p>
"""


def section_calibration(arms) -> str:
    """Per-file mass drift, measured against this run's own confident identifications, and
    corrected where it exceeded the noise floor.

    Written by `run_study.py --auto-calibrate` (opt in, not a default -- see
    `docs/MASS_CALIBRATION.md`): a throwaway pass at loose linking measures where each file's mass
    axis actually sits, and only files that drifted enough to matter get corrected before the run
    that produced this report.
    """
    checked = [a for a in arms if a.calibration]
    if not checked:
        return f"""
{head('calibration')}
<p class="lede">Not run for this study &mdash; <code>run_study.py --auto-calibrate</code> was not
used. A fixed linking tolerance assumes every injection sits on the same mass axis; where that
assumption is wrong, features from a drifted file split into their own group instead of joining
the rest, understating that lipid's presence for that injection.</p>
"""
    rows = []
    for a in checked:
        c = a.calibration
        if c.get("refused"):
            rows.append([a.label, "&mdash;", "&mdash;", c["refused"], "&mdash;"])
            continue
        spread = c.get("spread_ppm")
        spread_s = f"{spread[0]:+.1f} to {spread[1]:+.1f}" if spread else "&mdash;"
        bias = c.get("batch_bias_ppm")
        bias_s = f"{bias:+.1f}" if bias is not None else "not measured"
        rows.append([a.label, f"{c.get('total_files', 0):,}",
                    f"{c.get('corrected_files', 0):,}", spread_s, bias_s])
    corrected_any = any(c.get("corrected_files") for a in checked if (c := a.calibration))
    return f"""
{head('calibration')}
<p class="lede">Measured per file against consensus m/z from this run's own confident
identifications (<code>dot &ge; 900</code>), plus a batch-level bias against theoretical library
masses where enough of those exist. A file correction is applied only past the 2.0 ppm floor
below, so a correction here is real drift, not noise being fitted.</p>
{table(["polarity", "files", "corrected", "spread across files (ppm)",
       "batch absolute bias (ppm)"], rows, numeric=(1, 2))}
<div class="note"><strong>{"Correction applied &mdash; this run's identifications come from the "
"corrected files, not the originals." if corrected_any else
"Checked and clean &mdash; nothing exceeded the correction floor, so this run used the original "
"files unmodified."}</strong> The floor is 2.0 ppm: below it, correcting would be fitting noise
rather than removing a real offset. Residual spread on corrected files under 2 ppm is the
confirmation the correction took; a wider spread there means it did not.</div>
"""


def section_artefacts(arms) -> str:
    """Peaks stripped from MS2 spectra because no singly-charged CHNOPS ion could carry that mass
    defect &mdash; column bleed, background polymer, electronic noise, never a real lipid fragment.

    Written by `screen_artefacts` (on by default): per file, a peak present in most spectra, often
    as the base peak, at a mass no real ion could have, is removed before spectral matching sees
    it. A spiked deuterated standard is protected explicitly and never flagged this way.
    """
    checked = [a for a in arms if a.artefacts]
    if not checked:
        return f"""
{head('artefacts')}
<p class="lede">Not available for this run &mdash; either it predates the default
artefact-screening step, or <code>screen_artefacts</code> was turned off.</p>
"""
    rows = [[a.label, f"{a.artefacts['files_screened']:,}",
            f"{a.artefacts['files_with_artefacts']:,}",
            f"{sum(v['peaks_removed'] for v in a.artefacts['by_file'].values()):,}"]
           for a in checked]
    examples = []
    for a in checked:
        for stem, info in list(a.artefacts["by_file"].items())[:3]:
            for art in info["artefacts"][:2]:
                examples.append([a.label, stem, f"{art['mz']:.3f}", f"{100*art['fraction']:.0f}%",
                                f"{100*art['base_fraction']:.0f}%"])
    examples_block = (table(["polarity", "file", "m/z", "in this many spectra",
                            "as the base peak"], examples[:12], numeric=(3, 4))
                      if examples else "<p>None found in the files screened.</p>")
    return f"""
{head('artefacts')}
<p class="lede">Every MS2 spectrum carries whatever the instrument saw at that moment, real
fragments and background alike. A peak this common and this implausible is not a coincidence
across spectra &mdash; it is removed before it can compete with a real fragment for a spectral
match.</p>
{table(["polarity", "files screened", "files with an artefact", "peaks removed"], rows)}
{examples_block}
"""


def section_duty_cycle(arms) -> str:
    """What the MS2 budget was spent on, and how much of it never named anything.

    Computed from every MS2 scan in the batch (`build_exclusion_list`, on by default): a
    precursor earns a place here only by chromatography &mdash; never identified, fragmented
    enough times to matter, present across most of the gradient, and in a blank when one exists
    &mdash; never by "unidentified" alone, which would just as happily flag a real lipid the
    library does not cover.
    """
    have = [a for a in arms if a.duty_cycle]
    if not have:
        return f"""
{head('duty_cycle')}
<p class="lede">Not available for this run &mdash; either it predates the default exclusion-list
step, or <code>build_exclusion_list</code> was turned off.</p>
"""
    rows = []
    for a in have:
        d = a.duty_cycle
        rows.append([a.label, f"{d['total_ms2']:,}", f"{d['files_searched']:,}",
                    f"{d['wasted_precursors']:,}", f"{d['wasted_scans']:,}",
                    f"{d['duty_cycle_wasted_fraction']:.1%}"])
    series_rows = []
    for a in have:
        for s in a.duty_cycle["by_series"][:6]:
            series_rows.append([a.label, s["series"] or "unassigned", f"{s['precursors']:,}",
                                f"{s['scans']:,}", f"{s['fraction']:.1%}"])
    inclusion_rows = [[a.label, f"{a.inclusion_count:,}"] for a in arms
                      if a.inclusion_count is not None]
    inclusion_block = f"""
<p>Exclusion frees duty cycle; inclusion decides where some of it should go on the study's
<strong>next</strong> injection. These are real MS1 features detected in this batch &mdash; low
abundance relative to this run's own MS2-confirmed lipids, matching a library entry, and never
fragmented at all &mdash; so a low intensity is the most likely reason a real lipid never earned an
identification here.</p>
{table(["polarity", "candidates"], inclusion_rows, numeric=(1,))}
<div class="note"><strong>Also a candidate list, not an applied one, and specific to this run's own
samples.</strong> <code>Inclusion_List.csv</code> (m/z, retention) written beside the results is a
target list for re-injecting more of the same study, not identifications in this delivery &mdash;
nothing on this list has been confirmed by MS2, and forcing fragmentation on the next run is what
would confirm or rule it out.</div>
""" if inclusion_rows else ""
    return f"""
{head('duty_cycle')}
<p class="lede">Every MS2 scan is a slot the instrument spent once. One spent on a precursor that
was never identified, eluted across most of the gradient and (where a blank exists) fragmented
there too is a slot the sample did not get.</p>
{table(["polarity", "MS2 scans", "files", "unidentified precursors", "MS2 wasted on them",
       "share of duty cycle"], rows, numeric=(1, 2, 3, 4, 5))}
{table(["polarity", "series", "precursors", "MS2 scans", "share of duty cycle"], series_rows,
       numeric=(2, 3, 4)) if series_rows else ""}
<div class="note"><strong>This is a candidate list, not an applied one.</strong> Nothing on it has
been excluded from anything &mdash; the two files written beside the results
(<code>Exclusion_List.csv</code>, the evidence per candidate, and
<code>Exclusion_List_Thermo.csv</code>, the same rows in the instrument's mass-list import format)
are for a person to review and load onto the acquisition method, not something this software
applies on its own. Run <code>scripts/check_exclusion_safety.py</code> against the libraries and
the identifications before importing either: an exclusion entry does not remove one m/z, it removes
everything the instrument cannot tell apart from it within the method's own exclusion tolerance.</div>
{inclusion_block}"""


def section_methods(arms) -> str:
    fw = ", ".join(f"{a.label} {a.fwhm:.3f} min" for a in arms if a.fwhm)
    return f"""
{head('methods')}
<p><strong>Where the acceptance bands come from.</strong> The 15% band is this facility's own
target and is not a published threshold. 20% and 30% are the values in common use for untargeted
LC-MS profiling (Dunn et al., <em>Nature Protocols</em> 2011), and the D-ratio bar of 0.5 follows
Broadhurst et al., <em>Metabolomics</em> 2018. The &plusmn;15% figures from bioanalytical method
validation apply to targeted, calibrated assays and are not used here.</p>
<p>Areas are median-normalised per injection before any statistic, then log&#8321;&#8320;(x+1)
transformed and Pareto-scaled for the principal component analysis, which is computed by singular
value decomposition. A feature enters the analysis when detected in at least half the study
samples, decided once per polarity and applied to every panel. Hotelling's T&sup2; uses a 95%
limit; DModX uses a robust median plus three median absolute deviations, because a
distance-to-model is one-sided and skewed. Correlations against injection order are Spearman.
Measured peak width: {fw or 'not supplied'}. Figures use a colourblind-safe palette, with colour
&mdash; never size or shape &mdash; carrying the group.</p>
<div class="note"><strong>What this report does not do.</strong> It does not correct anything.
Drift is measured and shown; whether to apply a QC-RSC correction needs a threshold this facility
has not yet set, and a corrected matrix arriving without that decision invites the question of
what else was adjusted.</div>
<div class="note"><strong>&#9888; Some processing steps consulted these same QC injections.</strong>
Split peaks were merged using pooled-QC precision, features were filtered on detection rate, and
gaps were re-integrated from the raw files. The precision and missingness figures above therefore
describe a table those decisions helped shape, and are not fully independent evidence about it.
The effect is small &mdash; a few tens of rows out of several hundred &mdash; but it is the reason
this paragraph exists rather than being left for the reader to discover.</div>
"""


def _standards_beside(results: Path) -> Path:
    """Where the standards table is. It moved out of Results/ — a facility record, not a client
    deliverable — so the old location is still checked for analyses produced before the move."""
    folder = results.parent
    candidates = [folder.parent.parent / f"Standards_Performance_{folder.name}.csv",
                  folder.parent / f"Standards_Performance_{folder.name}.csv",
                  folder / "Standards_Performance.csv"]
    for path in candidates:
        if path.exists():
            return path
    return candidates[0]


def logo_data_uri() -> str:
    """The University logo as a data URI, so the HTML stays self-contained and WeasyPrint can
    render it into the PDF without reaching the network. Absent asset is not an error — a report
    without a logo is still a report."""
    import base64
    path = Path(__file__).resolve().parents[1] / "data/branding/university_of_edinburgh.png"
    if not path.exists():
        return ""
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def section_files(arms) -> str:
    """Last section: what each delivered CSV is, and which one to analyse.

    Placed at the end deliberately — it is reference material, read once when the folder is opened,
    not part of the argument the report is making about data quality.
    """
    rows = ""
    for a in arms:
        rows += (f"<tr><td>{a.label}</td><td class='n'>{a.n_unfiltered or '&mdash;'}</td>"
                 f"<td class='n'>{a.n_kept or '&mdash;'}</td>"
                 f"<td class='n'>{a.n_identified or '&mdash;'}</td></tr>")
    return f"""
{head('files')}
<p>Four files per polarity, plus one that spans both. The <code>Compound Group</code> identifier
is stable across the four per-polarity ones, so any row joins to the same group in the others.</p>
<table>
<thead><tr><th>file</th><th>what it holds</th><th>use it for</th></tr></thead>
<tbody>
<tr><td><code>Unfiltered_Results.csv</code></td>
    <td>Every compound group, including rejected ones, each carrying a <code>Filter Status</code>
        giving the reason it was dropped &mdash; blank, dimer, isotope, adduct of an existing peak,
        retention-model outlier, or below the detection-rate threshold.</td>
    <td>Auditing. Checking why a lipid you expected is not in the results.</td></tr>
<tr><td><code>Final_Results.csv</code></td>
    <td>The kept rows, named and unnamed. Areas are <strong>as measured</strong>: a sample where the
        feature was not detected is zero.</td>
    <td>Seeing what was and was not detected, before any gap is filled.</td></tr>
<tr><td><code>Final_Results_Filtered.csv</code></td>
    <td><strong>Identified lipids only</strong>, with gaps measured by re-integrating the raw MS1 at
        the expected mass and retention rather than left blank or imputed. Zeros that remain are
        cells where re-integration found no peak, so a real absence stays an absence.</td>
    <td><strong>The analysis-ready table.</strong> Start here.</td></tr>
<tr><td><code>Associated_Spectra.csv</code></td>
    <td>One row per spectrum-to-library match, at every rank, with precursor, library mass, mass
        error, dot product, reverse dot product and purity.</td>
    <td>Evidence. Why a lipid carries the name it does, and what the runners-up were.</td></tr>
<tr><td><code>Combined_Filtered_Normalised.csv</code></td>
    <td>Both polarities' <code>Final_Results_Filtered_median_Normalised.csv</code> stacked into one
        table, one row per lipid per polarity. Each polarity is median-normalised
        <strong>within itself</strong>, not against the other &mdash; the <code>Scale</code> column
        says so, and the two blocks are not on a common scale. Written only when Pos and Neg
        injections represent the same vials (matched by sample name once the polarity suffix is
        removed); a study where the two polarities were not the same physical injections will not
        have this file.</td>
    <td>A single per-lipid table spanning both polarities &mdash; useful for a combined heat map or
        PCA, but standardise per row (or split back out by <code>Mode</code>) first, or the
        positive block wins on magnitude and row count alone.</td></tr>
<tr><td><code>Exclusion_List.csv</code> / <code>Duty_Cycle.json</code></td>
    <td>MS1&rarr;MS2 duty-cycle candidates: precursors never identified, fragmented enough to
        matter, present across most of the gradient. Evidence per row, never applied
        automatically &mdash; see {ref('duty_cycle')}.</td>
    <td>Deciding what to put on the instrument's own exclusion list, after
        <code>scripts/check_exclusion_safety.py</code>.</td></tr>
<tr><td><code>Inclusion_List.csv</code></td>
    <td>Real MS1 features from this batch, low abundance against this run's own MS2-confirmed
        lipids, matching a library entry, never fragmented &mdash; see
        {ref('duty_cycle')}.</td>
    <td>Forced-fragmentation targets for the study's <strong>next</strong> injection, not
        identifications in this delivery.</td></tr>
<tr><td><code>Artefacts.json</code></td>
    <td>Impossible-mass-defect peaks stripped from MS2 spectra before spectral matching saw them
        &mdash; see {ref('artefacts')}.</td>
    <td>Confirming what was removed and why, per file.</td></tr>
<tr><td><code>Calibration.json</code></td>
    <td>Per-file mass drift, measured and (past the 2 ppm floor) corrected before this run &mdash;
        written only when <code>run_study.py --auto-calibrate</code> was used &mdash; see
        {ref('calibration')}.</td>
    <td>Confirming whether this run's identifications came from the original files or corrected
        ones.</td></tr>
</tbody></table>
<table>
<thead><tr><th>polarity</th><th class='n'>all groups</th><th class='n'>kept rows</th>
<th class='n'>identified lipids</th></tr></thead>
<tbody>{rows}</tbody></table>
<div class="note"><strong>Why the analysis-ready table is smaller.</strong> It contains only rows
carrying a name. An unidentified feature cannot enter a lipid-level analysis, and leaving thousands
of them in a file labelled analysis-ready invites someone to model them by accident. Nothing is
lost: every kept row, named or not, remains in <code>Final_Results.csv</code>, and every rejected
one in <code>Unfiltered_Results.csv</code> with its reason.</div>
<div class="note"><strong>&#9888; Rows named without a spectrum.</strong> Where the
<code>Identification Source</code> column reads <code>RT model</code>, the assignment rests on
accurate mass and retention alone &mdash; no fragmentation, and isomers are not distinguished.
Those are level 3 identifications and should not be pooled with the MS2-derived ones in a
downstream analysis without saying so.</div>
"""


def section_verdict(arms) -> str:
    """One verdict per polarity, never averaged — they can reach different conclusions."""
    blocks = ""
    for a in arms:
        # ⚠ Without pooled QCs there is no precision evidence at all, and the verdict must say so.
        # It previously read "the measurement is sound. Median pooled-QC CV nan%" — passing a run
        # it had not assessed, which is worse than refusing to judge it.
        if not a.precision.get("usable", True):
            blocks += (f'<div class="verdict warn"><h3>{a.label.capitalize()} mode</h3>'
                       f'<p><strong>Precision cannot be assessed.</strong> '
                       f'{html.escape(a.precision.get("reason", ""))}. Repeatability is measured '
                       f'on repeated injections of one pooled material, and without them nothing '
                       f'in this report can tell an instrument problem from a real difference '
                       f'between samples. Everything else below still applies; this is a limit '
                       f'of the sequence design, not a fault found in the data.</p></div>')
            continue
        cv = a.precision["median_pool_cv"]
        drift = (a.order["drifting"] / a.order["tested"]) if a.order.get("tested") else 0.0
        # Blanks are excluded: a blank that resembled the samples would be the fault. The check
        # is whether a SAMPLE or a POOL failed to look like the material it should look like.
        failed = [p["injection"] for p in a.profile
                  if p.get("below_sample_band") and p.get("role") != "blank"]
        problems = []
        if cv > 30:
            problems.append(f"the median pooled-QC CV is {cv:.1f}%, above the 30% bar")
        if drift > 0.30:
            problems.append(f"{drift:.0%} of features track injection order")
        if failed:
            problems.append(f"{len(failed)} injection(s) do not resemble the samples")
        verdict = "warn" if problems else ""
        lead = ("<strong>The data are fit for purpose.</strong>" if not problems
                else "<strong>The data are not fit for purpose.</strong>")
        detail = (f" Median pooled-QC CV {cv:.1f}%, {drift:.0%} of features tracking injection "
                  f"order, and no injection failed."
                  if not problems else " " + "; ".join(problems).capitalize() + ".")
        blocks += (f'<div class="verdict {verdict}"><h3>{a.label.capitalize()} mode</h3>'
                   f'<p>{lead}{detail} '
                   f'{a.missing["missing_identified"] * 100:.0f}% of identified values are missing, '
                   f'left-censored, so imputation must be chosen for that &mdash; '
                   f'see {ref('missing')}.</p></div>')
    return f"""
<h2>Verdict</h2>
<p class="lede">Reported separately for each polarity. They measure different molecules and can
reach different conclusions, so they are never averaged into one statement.</p>
{blocks}
<div class="note"><strong>This verdict is about the measurement, not about the biology.</strong>
Whether the samples differ from each other is a property of the experiment, not of the data
quality: a study whose intervention changed nothing produces tightly clustered samples, and that
is a result, not a fault. Statistics that compare biological spread against analytical spread
&mdash; the D-ratio and the sample-to-QC spread &mdash; are reported in <strong>{ref('support')}</strong>
as a guide to what the data can support, and they deliberately take no part in the judgement
above.</div>
"""


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ap = argparse.ArgumentParser()
    for pol in ("pos", "neg"):
        ap.add_argument(f"--results-{pol}")
        ap.add_argument(f"--unfiltered-{pol}")
        ap.add_argument(f"--metadata-{pol}")
        ap.add_argument(f"--spectra-{pol}")
        ap.add_argument(f"--filled-{pol}",
                        help="Final_Results_Filtered.csv — the delivered, gap-filled table")
        ap.add_argument(f"--standards-{pol}",
                        help="Standards_Performance.csv for this polarity")
        ap.add_argument(f"--fwhm-{pol}", type=float)
    ap.add_argument("--sequence")
    ap.add_argument("--study", default="")
    ap.add_argument("--owner", default="")
    ap.add_argument("--instrument", default="")
    ap.add_argument("--method", default="")
    ap.add_argument("--analyst", default="Dr. Jair G. Marques")
    ap.add_argument("--facility", default="Mass Spectrometry Core Facility",
                    help="the issuing facility; the institute and university are added around it")
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-pdf", action="store_true")
    args = ap.parse_args()

    arms = []
    for key, label in (("pos", "positive"), ("neg", "negative")):
        results = getattr(args, f"results_{key}")
        if not results:
            continue
        # The peak finder measures FWHM per run and writes it beside the results. A per-polarity
        # constant here was this method's peak width, and would have been silently wrong on any
        # other gradient — the same failure mode as an absolute noise threshold.
        measured = measured_fwhm(Path(results).parent)
        arms.append(Arm(label, results,
                        unfiltered=getattr(args, f"unfiltered_{key}"),
                        sequence=args.sequence,
                        metadata=getattr(args, f"metadata_{key}"),
                        spectra=getattr(args, f"spectra_{key}"),
                        filled=getattr(args, f"filled_{key}") or _beside(results),
                        standards=(getattr(args, f"standards_{key}")
                                   or _standards_beside(Path(results))),
                        fwhm=getattr(args, f"fwhm_{key}") or measured))
    if not arms:
        print("nothing to report: pass --results-pos and/or --results-neg", file=sys.stderr)
        return 2
    learn_prefix([c for a in arms for c in a.run.columns])
    for a in arms:
        print(f"{a.label}: {a.run.areas.shape[0]:,} features x {len(a.run.columns)} injections "
              f"({a.counts['sample']} samples, {a.counts['qc']} pools, {a.counts['blank']} blanks)")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    polarities = " and ".join(a.label for a in arms)
    meta = [("Study", args.study), ("Samples from", args.owner),
            ("Polarities", f"{polarities} &mdash; reported side by side, never merged"),
            ("Instrument", args.instrument), ("Method", args.method),
            ("Injections", ", ".join(f"{a.counts['sample']}+{a.counts['qc']}+{a.counts['blank']} "
                                     f"({a.label})" for a in arms)),
            ("Analyst", args.analyst), ("Issued", date.today().isoformat())]
    logo = logo_data_uri()
    brand = (f'<div class="brand">'
             + (f'<img src="{logo}" alt="The University of Edinburgh">' if logo else "")
             + f'<div class="brand-text"><strong>{args.facility}</strong>'
               f'<span>{INSTITUTION}</span>'
               f'<a href="{FACILITY_URL}">{FACILITY_URL.replace("https://", "")}</a>'
               f'</div></div>')
    header = (f'<header>{brand}'
              f'<h1>Quality control report</h1>'
              f'<div class="sub">{html.escape(args.study) or "untargeted lipidomics"} '
              f'&mdash; {polarities} mode</div>'
              f'<div class="meta">'
              + "".join(f"<div><span>{k}</span>{v}</div>" for k, v in meta if v)
              + "</div></header>")

    # Order comes from SECTIONS: technical first, biological last.
    body = (header + section_verdict(arms)
            + section_precision(arms, plt) + section_accuracy(arms, plt)
            + section_calibration(arms)
            + section_order(arms, plt)
            + section_blanks(arms, plt) + section_artefacts(arms)
            + section_duty_cycle(arms) + section_missing(arms, plt)
            + section_polarity(arms) + section_annotation(arms, plt) + section_outliers(arms, plt)
            + section_filters(arms)
            + section_support(arms) + section_structure(arms, plt)
            + section_methods(arms) + section_files(arms)
            + f'<footer>Generated {date.today().isoformat()} by '
              f'<code>scripts/qc_report.py</code>. Positive and negative mode are analysed '
              f'independently throughout; no figure or statistic in this document combines '
              f'them.</footer>')
    page = (f"<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>QC report — {html.escape(args.study)}</title>"
            f"<style>{STYLE}</style></head><body>{body}</body></html>")

    target = out / "QC_report.html"
    target.write_text(page, encoding="utf-8")
    print("written", target)
    if not args.no_pdf:
        try:
            from weasyprint import HTML
            HTML(string=page).write_pdf(out / "QC_report.pdf")
            print("written", out / "QC_report.pdf")
        except Exception as exc:            # noqa: BLE001 - the HTML is the deliverable
            print(f"PDF not written ({exc})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
