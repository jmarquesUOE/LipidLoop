# Known limitations

- **The decoy estimate covers chain-assignment error only.** It measures errors within a correct
  class and sum composition on spectrum-supported rows. Rows that are another lipid's isotope,
  adduct or in-source ion arise before spectral matching and no spectral decoy represents them; on
  ground-truth data they were 4 to 8 % of delivered rows. The co-elution step recognises many of
  them, not all.
- **Retention-model names carry no error estimate.** A row named by mass and retention without a
  spectrum (MSI level 3) cannot be assigned to a decoy. These rows are labelled in every table.
- **Classes without a retention model pass unjudged.** A class with too few members, or a poor fit,
  gets no model, so its rows are never tested on retention.
- **Isotopes of identified rows.** The isotope rule removes an M+1 or M+2 trace only from an
  unidentified peak; an identified row that is the M+2 isotope of the species with one more double
  bond is left to the co-elution step. Coenzyme Q10H2 lies within the MS1 tolerance of the CoQ10 M+2
  isotope and is separable only by retention: run `scripts/check_ubiquinol_rt.py` on any table that
  reports it.
- **Feature finding is not fully deterministic.** Reruns of the same commit on the same files can
  differ by one or two rows.
- **Ion-trap MS2.** The default fragment window (0.01 Da) suits Orbitrap and TOF MS2. For ion-trap
  MS2 the run measures the axis and widens the window when the configured value is implausible;
  check the decision in the configuration.
- **Ether lipids.** A plasmanyl and a plasmenyl name supported by a single head-group fragment
  carry no discriminating information; such ties are broken by library order.
- **Scope.** Validated on data-dependent reversed-phase lipidomics from five vendors. HILIC,
  data-independent and MS1-only acquisitions are out of scope.
- **Installed copies.** The pipeline reads its tables from `data/` beside the source, so install
  from a checkout (`pip install -e .`); a copy that is not a git checkout records the code version
  as `unknown`.
