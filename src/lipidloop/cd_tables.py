"""Read Compound Discoverer's exported peak tables.

Only needed to validate the peak finder against a reference run: with CD's own features going
in, any difference in `Final_Results.csv` is the peak finder's rather than the feature
detector's. The pipeline proper builds the same structures from `features.py`.

The aligned export is one CSV holding three nested tables, distinguished by how many empty
leading columns a row has:

    Checked,Name,Molecular Weight,RT [min],Area (Max.),MS2,Area: QC_01.raw (F1),...   <- group
    ,Checked,Molecular Weight,RT [min],FWHM [min],Max. # MI,# Adducts,Area,Study File ID
    ,,Checked,Ion,Charge,Molecular Weight,m/z,RT [min],FWHM [min],# MI,Area,Parent Area [%],...
"""
from __future__ import annotations

import csv
from pathlib import Path

from .peaks import UNALIGNED_PPM, Compound, CompoundGroup, Feature, Sample, ppm_diff


def _f(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _i(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def read_aligned(path: str | Path, min_feature_count: int = 1
                 ) -> tuple[list[CompoundGroup], list[Sample]]:
    """Parse the aligned export into compound groups, compounds and features."""
    groups: list[CompoundGroup] = []
    samples: list[Sample] = []
    by_id: dict[str, Sample] = {}

    group_header: list[str] = []
    compound_header: list[str] = []
    feature_header: list[str] = []
    area_columns: list[int] = []

    current_group: CompoundGroup | None = None
    current_compound: Compound | None = None

    with Path(path).open(newline="") as fh:
        for row in csv.reader(fh):
            if not row or not any(row):
                continue
            depth = 0
            while depth < len(row) and row[depth] == "":
                depth += 1

            if "Checked" in row:
                if depth == 0:
                    group_header = row
                    area_columns = [i for i, c in enumerate(row) if c.startswith("Area:")]
                    for i in area_columns:
                        name = row[i].split("Area:", 1)[1].strip()
                        sample = Sample(file=name)
                        if "(" in name:
                            sample.cd_id = name[name.rfind("(") + 1:name.rfind(")")]
                        samples.append(sample)
                        by_id[sample.cd_id] = sample
                elif depth == 1:
                    compound_header = row
                else:
                    feature_header = row
                continue

            if depth == 0:
                if current_group is not None:
                    if current_compound is not None:
                        current_group.compounds.append(current_compound)
                        current_compound = None
                    if len(current_group.compounds) >= min_feature_count:
                        groups.append(current_group)
                current_group = _parse_group(row, group_header, area_columns)
            elif depth == 1:
                if current_compound is not None and current_group is not None:
                    current_group.compounds.append(current_compound)
                current_compound = _parse_compound(row, compound_header, by_id)
            else:
                if current_compound is not None:
                    feature = _parse_feature(row, feature_header, by_id)
                    if feature is not None:
                        current_compound.features.append(feature)

    if current_group is not None:
        if current_compound is not None:
            current_group.compounds.append(current_compound)
        if len(current_group.compounds) >= min_feature_count:
            groups.append(current_group)

    return groups, samples


def _column(header: list[str], *names: str) -> int:
    for name in names:
        for i, c in enumerate(header):
            if c.strip() == name:
                return i
    for name in names:
        for i, c in enumerate(header):
            if name in c:
                return i
    return -1


def _parse_group(row, header, area_columns) -> CompoundGroup:
    name_i = _column(header, "Name")
    mw_i = _column(header, "Molecular Weight")
    rt_i = _column(header, "RT [min]")
    area_i = _column(header, "Area (Max.)")
    ms2_i = _column(header, "MS2")
    areas = [_f(row[i]) if i < len(row) else 0.0 for i in area_columns]
    return CompoundGroup(
        name=row[name_i] if 0 <= name_i < len(row) else "",
        mw=_f(row[mw_i]) if 0 <= mw_i < len(row) else 0.0,
        retention=_f(row[rt_i]) if 0 <= rt_i < len(row) else 0.0,
        max_area=_f(row[area_i]) if 0 <= area_i < len(row) else 0.0,
        has_ms2=bool(row[ms2_i].strip()) if 0 <= ms2_i < len(row) else False,
        areas=areas)


def _parse_compound(row, header, by_id) -> Compound:
    file_id = row[_column(header, "Study File ID")] if _column(header, "Study File ID") >= 0 else ""
    sample = by_id.get(file_id.strip()) or Sample(file=file_id.strip(), cd_id=file_id.strip())
    return Compound(
        mw=_f(row[_column(header, "Molecular Weight")]),
        retention=_f(row[_column(header, "RT [min]")]),
        fwhm=_f(row[_column(header, "FWHM [min]")]),
        max_mi=_i(row[_column(header, "Max. # MI")]),
        n_adducts=_i(row[_column(header, "# Adducts")]),
        area=_f(row[_column(header, "Area")]),
        sample=sample)


def _parse_feature(row, header, by_id) -> Feature | None:
    file_i = _column(header, "Study File ID")
    if file_i < 0 or file_i >= len(row):
        return None
    sample = by_id.get(row[file_i].strip()) or Sample(file=row[file_i].strip(),
                                                      cd_id=row[file_i].strip())
    charge = _i(row[_column(header, "Charge")], 1)
    adduct = row[_column(header, "Ion")]
    if "]-" in adduct and charge > 0:
        charge = -charge
    return Feature(
        adduct=adduct,
        charge=charge,
        mw=_f(row[_column(header, "Molecular Weight")]),
        mass=_f(row[_column(header, "m/z")]),
        retention=_f(row[_column(header, "RT [min]")]),
        fwhm=_f(row[_column(header, "FWHM [min]")]),
        mi=_i(row[_column(header, "# MI")]),
        area=_f(row[_column(header, "Area")]),
        parent_area_percent=_f(row[_column(header, "Parent Area [%]")]),
        sample=sample)


def read_unaligned(path: str | Path, groups: list[CompoundGroup], avg_fwhm: float) -> int:
    """Attach each feature's *measured* retention time from the unaligned export.

    That measured-versus-aligned difference is the whole basis of the retention-time
    correction: without it an MS2's retention time cannot be compared to a feature's.
    """
    index: dict[int, list[CompoundGroup]] = {}
    bucket = 0.1
    for group in groups:
        index.setdefault(int(group.retention / bucket), []).append(group)

    matched = 0
    with Path(path).open(newline="") as fh:
        header: list[str] = []
        for row in csv.reader(fh):
            if not row or not any(row):
                continue
            if "Checked" in row:
                header = row
                continue
            if not header:
                continue
            mass = _f(row[_column(header, "m/z")])
            retention = _f(row[_column(header, "RT [min]")])
            area = _f(row[_column(header, "Area")])
            file_id = row[_column(header, "Study File ID")].strip()
            adduct = row[_column(header, "Ion")]
            polarity = "+" if "]+" in adduct else "-"

            window = max(avg_fwhm * 3.0, 1e-6)
            for b in range(int((retention - window) / bucket), int((retention + window) / bucket) + 1):
                done = False
                for group in index.get(b, ()):
                    if group.quant_ion is None:
                        continue
                    if abs(group.retention - retention) >= 2.0:
                        continue
                    for compound in group.compounds:
                        for feature in compound.features:
                            if (ppm_diff(mass, feature.mass) < UNALIGNED_PPM
                                    and abs(retention - feature.retention) < 1.0
                                    and polarity == feature.polarity
                                    and feature.area and area / feature.area > 0.95
                                    and feature.sample.cd_id == file_id):
                                feature.real_retention = retention
                                matched += 1
                                done = True
                                break
                        if done:
                            break
                    if done:
                        break
                if done:
                    break
    return matched
