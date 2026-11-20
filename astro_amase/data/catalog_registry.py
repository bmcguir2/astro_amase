"""Resolve the same source-qualified spectroscopy in candidate search, fitting and saved plots."""

import csv
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import shutil

import numpy as np


class MeasuredPartition:
    """Delegate to the repaired catalog within its measured temperature interval."""
    def __init__(self, partition, bounds, label):
        self.partition, self.bounds, self.label = partition, bounds, label

    def _check(self, temperature):
        lo, hi = self.bounds
        if not np.all((lo <= np.asarray(temperature)) & (np.asarray(temperature) <= hi)):
            raise ValueError(f"{self.label} has verified partition data from {lo:g} to {hi:g} K. "
                             "Choose a supported temperature or exclude this catalog.")

    def q(self, temperature):
        self._check(temperature)
        return self.partition.q(temperature)

    def qrot(self, temperature):
        self._check(temperature)
        return self.partition.qrot(temperature)

    def qvib(self, temperature):
        self._check(temperature)
        return self.partition.qvib(temperature)


@dataclass(frozen=True)
class CatalogEntry:
    source: str
    tag: int
    label: str
    path: Path
    smiles: str
    upload: object = None
    temperature_bounds: tuple = ()

    def load(self):
        if self.upload is not None:
            return self.upload.molecule()
        with self.path.open("rb") as handle:
            molecule = pickle.load(handle)
        if self.temperature_bounds:
            molecule.qpart = MeasuredPartition(molecule.qpart, self.temperature_bounds, self.label)
        return molecule


class CatalogRegistry:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.by_tag = {}
        self.by_label = {}
        self.uploads = []
        record_file = self.directory / "build_record.json"
        self.build = json.loads(record_file.read_text()) if record_file.is_file() else {}
        self.band = self.build.get("validated_observing_band_MHz")
        for source, filename, label_field, tag_field, folder in (
            ("CDMS", "all_cdms_final_official.csv", "mol", "tag", "cdms_pkl"),
            ("JPL", "all_jpl_final_official.csv", "name", "save tag", "jpl_pkl"),
        ):
            with (self.directory / filename).open(encoding="utf-8-sig", newline="") as handle:
                for row in csv.DictReader(handle):
                    tag = int(row[tag_field])
                    stem = f"{tag:06d}" if source == "CDMS" else str(tag)
                    member = f"{folder}/{stem}.pkl"
                    limits = self.build.get("files", {}).get(member, {}).get("validated_temperature_K", ())
                    self._register(CatalogEntry(source, tag, row[label_field],
                                                self.directory / member, row["smiles"], temperature_bounds=tuple(limits)))
        base = self.directory / "uploaded_catalogs"
        if base.exists():
            from .upload_catalog import read_uploads
            for item in read_uploads(base):
                if item.label in self.by_label:
                    raise ValueError(f"Uploaded catalog label is already active: {item.label}")
                limits = (float(item.partition.temps[0]), float(item.partition.temps[-1]))
                self._register(CatalogEntry("LOCAL", item.tag, item.label, item.path, item.smiles, item, limits))
                self.uploads.append(item)

    def _register(self, entry):
        key = entry.source, entry.tag
        if key in self.by_tag:
            raise ValueError(f"Repeated active catalog source and tag: {key}")
        self.by_tag[key] = entry
        # Preserve the upstream CDMS-first fallback for old analyses saved by label.
        self.by_label.setdefault(entry.label, {}).setdefault(entry.source, entry)

    def entry(self, label, source=None):
        matches = self.by_label.get(label, {})
        if source in matches:
            return matches[source]
        if source is None and matches:
            return next(iter(matches.values()))
        raise KeyError(f"No active catalog for {label!r} ({source or 'any source'}).")

    @staticmethod
    def _same_transitions(upload, frequency, energy, degeneracy):
        freq, _, _, upper, gu, _, _ = upload.columns
        order = np.argsort(frequency, kind="stable")
        predictions = frequency[order]
        left = np.searchsorted(predictions, freq)
        right = np.searchsorted(predictions, freq, side="right")
        if np.any(left == right):
            return False
        # Equal frequencies can occur for physically distinct transitions.
        for i, (lo, hi) in enumerate(zip(left, right)):
            rows = order[lo:hi]
            if not np.any((degeneracy[rows] == gu[i]) & np.isclose(energy[rows], upper[i], rtol=1e-8, atol=.001)):
                return False
        return True

    def check_upload_overlap(self, arrays, upload, preceding):
        """Catch the same predictions reintroduced in text form under a different label."""
        for previous in preceding:
            freq, _, _, upper, gu, _, _ = previous.columns
            if self._same_transitions(upload, freq, upper, gu) or self._same_transitions(
                    previous, upload.columns[0], upload.columns[3], upload.columns[4]):
                raise ValueError(f"Uploads {upload.label!r} and {previous.label!r} contain the same "
                                 "transitions. Select one representation.")
        freq = upload.columns[0]
        native_freq, native_tags, native_sources = arrays[0], arrays[2], arrays[3]
        candidates = None
        for i in np.unique(np.linspace(0, len(freq)-1, min(8,len(freq)), dtype=int)):
            first = np.searchsorted(native_freq, freq[i], side="left")
            last = np.searchsorted(native_freq, freq[i], side="right")
            at_frequency = set(zip(native_sources[first:last], native_tags[first:last]))
            candidates = at_frequency if candidates is None else candidates & at_frequency
            if not candidates:
                return
        for key in sorted(candidates):
            entry = self.by_tag[key]
            native = entry.load().catalog
            if self._same_transitions(upload, native.frequency, native.eup, native.gup):
                raise ValueError(f"All transitions in the upload {upload.label!r} already occur in "
                                 f"{entry.source}:{entry.tag:06d} ({entry.label}). Select one representation.")

    def index(self, velocity=False):
        filename = "vlsr_database_logint.pkl.gz" if velocity else "transitions_database.pkl.gz"
        with gzip.open(self.directory / filename, "rb") as handle:
            arrays = list(pickle.load(handle))
        # The established velocity-tracer list stays curated. An explicitly requested
        # uploaded tracer is available through the full transition index fallback.
        if velocity or not self.uploads:
            return tuple(arrays)
        for i, item in enumerate(self.uploads):
            self.check_upload_overlap(arrays, item, self.uploads[:i])
        extras = [upload.index_columns() for upload in self.uploads]
        order = None
        for i, base in enumerate(arrays):
            joined = np.concatenate([base, *(extra[i] for extra in extras)])
            if order is None:
                order = np.argsort(joined, kind="stable")
            arrays[i] = joined[order]
        return tuple(arrays)

    def validate_structure_mode(self, enabled, model):
        if enabled:
            for item in self.uploads:
                if item.smiles:
                    try:
                        model.embed_smiles(item.smiles)
                        continue
                    except Exception:
                        pass
                raise ValueError(f"{item.label} has no supported structural embedding. "
                                 "Set consider_structure=False to run spectroscopy-only analysis.")

    def validate_temperature(self, temperature):
        for item in self.uploads:
            item.validate_temperature(temperature)
        for entry in self.by_tag.values():
            if entry.temperature_bounds:
                low, high = entry.temperature_bounds
                if not low <= temperature <= high:
                    raise ValueError(f"{entry.label} has validated partition data from {low:g} to {high:g} K. "
                                     "Choose a temperature in this range or exclude this catalog.")

    def temperature_interval(self, low, high):
        for entry in self.by_tag.values():
            if entry.temperature_bounds:
                low = max(low, entry.temperature_bounds[0])
                high = min(high, entry.temperature_bounds[1])
        if low >= high:
            raise ValueError("The active catalogs have no shared temperature interval for this fit.")
        return low, high

    def validate_band(self, frequencies):
        if self.band and (np.min(frequencies) < self.band[0] or np.max(frequencies) > self.band[1]):
            lo, hi = self.band
            raise ValueError(f"This catalog selection was verified only between {lo:g} and {hi:g} MHz. "
                             "Use the general core for other observing bands.")

    def scope_record(self):
        files = {}
        for name in ("active_catalog_manifest.json", "build_record.json"):
            path = self.directory / name
            if path.is_file():
                files[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        for item in self.uploads:
            for path in item.source_files:
                name = path.relative_to(self.directory).as_posix()
                files[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        return {"format": "astroamase-catalog-scope-v1", "files": files,
                "uploads_without_structural_prior": [u.label for u in self.uploads if not u.smiles],
                "missing_structure_policy": "Use declared elements and spectroscopic scores. "
                                            "A SMILES-based isotope-parent prior requires an identified structure."}

    def save_scope(self, output):
        record = self.scope_record()
        if record["files"]:
            dest = Path(output)
            (dest / "catalog_scope.json").write_text(json.dumps(record, indent=2) + "\n")
            for name in record["files"]:
                path = dest / name
                path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(self.directory / name, path)

    def check_saved_scope(self, stored_parameters):
        path = Path(stored_parameters).parent / "catalog_scope.json"
        if path.is_file() and json.loads(path.read_text()) != self.scope_record():
            raise ValueError("The active catalogs differ from the ones used for this saved fit. "
                             "Restore the recorded database and uploaded catalog files.")
