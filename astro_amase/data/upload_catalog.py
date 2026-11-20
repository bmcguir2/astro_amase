"""Validated text-only local catalogs. Uploaded metadata never selects filesystem paths."""
import csv
from dataclasses import dataclass
import fcntl
import io
import json
import math
import os
from pathlib import Path
import re
from tempfile import TemporaryDirectory

import numpy as np
from rdkit import Chem
from scipy.constants import c, h, k

from ..utils.molsim_classes import Catalog, Molecule

CSV_SCHEMA = "astroamase-local-csv-v1"
CAT_SCHEMA = "astroamase-local-spcat-v1"
COLUMNS = ["frequency_MHz", "uncertainty_MHz", "A_s-1", "E_upper_K", "g_upper", "qn_upper", "qn_lower"]
MAX_META = 64 * 1024
MAX_LINES = 100_000
MAX_BYTES = 16 * 1024**2
MAX_CATALOGS = 32
MAX_SCOPE_BYTES = 64 * 1024**2
MHZ_TO_K = h / k * 1e6
CM_TO_K = h * c * 100 / k
# Same SPCAT coefficients as molsim_classes.Catalog._set_sijmu_aij.
A_COEFFICIENT = 2.40251e4 * 1.16395e-20
TEXT_FIELD = re.compile(r"[\w ()\[\]+,.'=:#^–−-]+", re.UNICODE)


def _keys(value, names):
    if not isinstance(value, dict) or set(value) != set(names.split()):
        raise ValueError(f"Expected exactly these metadata keys: {names}")


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate metadata key {key!r}")
        result[key] = value
    return result


QN_FIELD = re.compile(r"[\w ()\[\]+,./'=:#^–−-]+", re.UNICODE)
NUMBER_FIELD = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[Ee][+-]?[0-9]+)?")


def _string(value, name, limit=96, display=False, qn=False):
    if (not isinstance(value, str) or not 0 < len(value) <= limit or
            value.strip() != value or not value.isprintable() or
            (display and not (QN_FIELD if qn else TEXT_FIELD).fullmatch(value))):
        raise ValueError(f"Invalid {name}")
    return value


def _int(value, name, lo, hi):
    if type(value) is not int or not lo <= value <= hi:
        raise ValueError(f"{name} must be an integer between {lo} and {hi}")
    return value


def _number(value, name, lo=0, hi=1e12, *, include_zero=False):
    if (type(value) not in (float, int, str) or
            (isinstance(value, str) and not NUMBER_FIELD.fullmatch(value.strip()))):
        raise ValueError(f"Invalid {name}")
    try:
        number = float(value)
    except (ValueError, OverflowError) as err:
        raise ValueError(f"Invalid {name}") from err
    if not math.isfinite(number) or not lo <= number <= hi or (number == 0 and not include_zero):
        raise ValueError(f"{name} outside supported range")
    return number


def _read(path, limit):
    if path.is_symlink() or not path.is_file() or path.stat().st_size > limit:
        raise ValueError(f"Expected a regular text file of at most {limit} bytes: {path.name}")
    with path.open("rb") as stream:
        data = stream.read(limit + 1)
    if len(data) > limit:
        raise ValueError(f"Input too large: {path.name}")
    return data


@dataclass
class TabulatedQ:
    temps: np.ndarray
    vals: np.ndarray

    def q(self, T):
        t = np.asarray(T, dtype=float)
        if not (np.isfinite(t) & (self.temps[0] <= t) & (t <= self.temps[-1])).all():
            raise ValueError(f"Temperature must be within the uploaded Q table [{self.temps[0]}, {self.temps[-1]}] K")
        return np.interp(t, self.temps, self.vals)

    qrot = q

    def qvib(self, T):
        return np.ones_like(self.q(T))


def _metadata(data):
    m = json.loads(data.decode("utf-8-sig"), object_pairs_hook=_unique_object)
    if not isinstance(m, dict) or m.get("schema") not in (CSV_SCHEMA, CAT_SCHEMA):
        raise ValueError("Unsupported local catalog schema")
    base = ("schema tag label formula state isotope_count vibrational_state hyperfine_resolved "
            "elements origin version degeneracy_convention partition")
    _keys(m, base + (" smiles" if "smiles" in m else "") +
          (" spcat_tag intensity_convention" if m["schema"] == CAT_SCHEMA else ""))
    m.setdefault("smiles", None)
    _int(m["tag"], "tag", 1, 2147483647)
    _int(m["isotope_count"], "isotope_count", 0, 1024)
    _int(m["vibrational_state"], "vibrational_state", 0, 100000)
    if type(m["hyperfine_resolved"]) is not bool:
        raise ValueError("hyperfine_resolved must be a boolean")
    for name in ("label", "formula", "state"):
        _string(m[name], name, display=True)
    if not m["label"][0].isalnum() or m["label"].casefold() in {"nan", "none", "unidentified", "null", "na", "n/a"}:
        raise ValueError("Use a distinct readable catalog label starting with a letter or digit")
    for name in ("origin", "version", "degeneracy_convention"):
        _string(m[name], name, 2048)
    if m["schema"] == CAT_SCHEMA:
        _int(m["spcat_tag"], "spcat_tag", 1, 999999)
        if m["intensity_convention"] != "log10_nm2_MHz_at_300_K":
            raise ValueError("SPCAT requires log10 integrated intensity in nm^2 MHz at 300 K")
    periodic = Chem.GetPeriodicTable()
    symbols = {periodic.GetElementSymbol(i) for i in range(1, 119)}
    elements = m["elements"]
    if (not isinstance(elements, list) or not elements or
            any(not isinstance(x, str) or x not in symbols for x in elements) or
            len(elements) != len(set(elements))):
        raise ValueError("elements must contain distinct chemical element symbols")
    if m["smiles"] is not None:
        s = _string(m["smiles"], "smiles", limit=192)
        options = Chem.SmilesParserParams()
        options.parseName = False
        options.allowCXSMILES = False
        molecule = Chem.MolFromSmiles(s, options)
        if molecule is None or {a.GetSymbol() for a in Chem.AddHs(molecule).GetAtoms()} != set(elements):
            raise ValueError("Invalid SMILES or declared elements disagree with SMILES (including hydrogen)")
        m["smiles"] = Chem.MolToSmiles(molecule)
    p = m["partition"]
    _keys(p, "convention interpolation source degeneracy_convention temperature_K Q")
    if (p["convention"] != "effective_lte" or p["interpolation"] != "linear_T_Q" or
            p["degeneracy_convention"] != m["degeneracy_convention"]):
        raise ValueError("Require effective LTE Q, linear T/Q interpolation and matching upper-state degeneracy conventions")
    _string(p["source"], "partition source", 2048)
    for field in ("temperature_K", "Q"):
        if not isinstance(p[field], list) or not 2 <= len(p[field]) <= 1024:
            raise ValueError(f"Invalid partition {field} table length")
        if any(type(x) not in (int, float) for x in p[field]):
            raise ValueError(f"Invalid partition {field} number")
    q = TabulatedQ(np.array([_number(x, "temperature_K", hi=1e6) for x in p["temperature_K"]]),
                   np.array([_number(x, "Q", hi=1e30) for x in p["Q"]]))
    if len(q.temps) != len(q.vals) or np.any(np.diff(q.temps) <= 0) or np.any(np.diff(q.vals) < 0):
        raise ValueError("Q must be nondecreasing on a strictly increasing temperature grid")
    q.q(300)  # Plot line rankings and SPCAT intensities use this reference temperature.
    return m, q


def _csv_rows(data):
    reader = csv.reader(io.StringIO(data.decode("utf-8-sig"), newline=""), strict=True)
    if next(reader, None) != COLUMNS:
        raise ValueError("Expected transitions.csv columns " + ", ".join(COLUMNS))
    for values in reader:
        if len(values) != len(COLUMNS):
            raise ValueError(f"Invalid CSV field count at line {reader.line_num}")
        f, err, a, eu, gu, up, low = values
        if not re.fullmatch(r"[0-9]{1,9}", gu):
            raise ValueError(f"Invalid g_upper at line {reader.line_num}")
        yield f, err, a, eu, int(gu), up, low


def _spcat_rows(data, meta, q300):
    # Restricted decimal SPCAT layout, with 1–6 occupied QN slots per state.
    # Preserve the actual two-character QN tokens, including encoded values.
    first_format = None
    for n, row in enumerate(io.StringIO(data.decode("ascii")), 1):
        row = row.rstrip("\r\n")
        if not 69 <= len(row) <= 79 or not row.isprintable() or "*" in row:
            raise ValueError(f"Unsupported SPCAT record at line {n}")
        row = row.ljust(79)
        freq = _number(row[:13], "frequency", hi=1e9)
        err = _number(row[13:21], "uncertainty", include_zero=True)
        logi = _number(row[21:29], "SPCAT log intensity", -300, 100, include_zero=True)
        _int(int(row[29:31]), "SPCAT degrees of freedom", 0, 3)
        el = _number(row[31:41], "SPCAT E_lower", hi=1e7, include_zero=True) * CM_TO_K
        gu = _int(int(row[41:44]), "decimal SPCAT g_upper", 1, 999)
        if abs(int(row[44:51])) != meta["spcat_tag"]:
            raise ValueError(f"SPCAT source tag differs at line {n}")
        fmt = _int(int(row[51:55]), "SPCAT QNFMT", 1, 9999)
        if first_format is not None and first_format != fmt:
            raise ValueError(f"Mixed SPCAT quantum number formats at line {n}")
        first_format = fmt
        count = fmt % 10
        if not 1 <= count <= 6:
            raise ValueError(f"Unsupported SPCAT QN count at line {n}")
        labels = []
        for start in (55, 67):
            tokens = [row[j:j + 2].strip() for j in range(start, start + 12, 2)]
            if (any(not re.fullmatch(r"[A-Za-z0-9+-]{1,2}", t) for t in tokens[:count])
                    or any(tokens[count:])):
                raise ValueError(f"Unsupported SPCAT quantum numbers at line {n}")
            labels.append(",".join(tokens[:count]))
        eu = el + freq * MHZ_TO_K
        log_a = (math.log(A_COEFFICIENT) + math.log(10) * logi + 2 * math.log(freq)
                 + math.log(q300) - math.log(gu) + el / 300 - math.log(-math.expm1(-freq * MHZ_TO_K / 300)))
        if not -700 <= log_a <= math.log(1e12):
            raise ValueError(f"Unrepresentable SPCAT Einstein A at line {n}")
        yield freq, err, math.exp(log_a), eu, gu, *labels


@dataclass
class UploadedCatalog:
    path: Path
    source_files: list
    metadata: dict
    partition: TabulatedQ
    columns: tuple

    @property
    def tag(self):
        return self.metadata["tag"]

    @property
    def label(self):
        return self.metadata["label"]

    @property
    def smiles(self):
        return self.metadata["smiles"]

    @property
    def elements(self):
        return self.metadata["elements"]

    @property
    def structure_key(self):
        return self.smiles or f"LOCAL:{self.tag}"

    @property
    def hyperfine_resolved(self):
        return self.metadata["hyperfine_resolved"]

    def validate_temperature(self, T):
        self.partition.q(T)

    def index_columns(self):
        frequency, err = self.columns[:2]
        vals = (self.tag, "LOCAL", self.structure_key, self.label, self.metadata["isotope_count"],
                self.metadata["vibrational_state"], self.metadata["formula"] + " " + self.metadata["state"])
        return (frequency, err, *(np.repeat(v, len(frequency)) for v in vals))

    def molecule(self):
        freq, err, a, eu, gu, up, low = self.columns
        delta = freq * MHZ_TO_K
        el = np.maximum(0, eu - delta)
        # Invert the upstream SPCAT formula in log space to preserve weak-line rankings.
        logi = ((np.log(a) + np.log(gu) - np.log(A_COEFFICIENT) - 2 * np.log(freq)
                 - np.log(self.partition.q(300)) - el / 300 + np.log(-np.expm1(-delta / 300))) / np.log(10))
        empty = np.full(len(freq), None, dtype=object)
        qns = {f"qn{i}{side}": empty for side in ("up", "low") for i in range(2, 9)}
        cat = Catalog(catid=self.tag, molecule=self.label, frequency=freq, freq_err=err,
                      aij=a, gup=gu, eup=eu, elow=el, logint=logi, qnup_str=up,
                      qnlow_str=low, qn1up=up, qn1low=low, source="LOCAL",
                      version=self.metadata["version"], refs=self.metadata["origin"], **qns)
        return Molecule(id=self.tag, name=self.label, formula=self.metadata["formula"],
                        catalog=cat, qpart=self.partition)


def _load(path):
    if path.is_symlink() or not path.is_dir():
        raise ValueError(f"Expected a local catalog directory: {path}")
    meta_file = path / "metadata.json"
    meta, q = _metadata(_read(meta_file, MAX_META))
    is_csv = meta["schema"] == CSV_SCHEMA
    cat_file = path / ("transitions.csv" if is_csv else "transitions.cat")
    if set(path.iterdir()) != {meta_file, cat_file}:
        raise ValueError("A local catalog must contain exactly its metadata and transition files")
    content = _read(cat_file, MAX_BYTES)
    rows = _csv_rows(content) if is_csv else _spcat_rows(content, meta, q.q(300))
    seen = set()
    parsed = []
    for i, (f, err, a, eu, gu, up, low) in enumerate(rows, 1):
        if i > MAX_LINES:
            raise ValueError("Too many transitions")
        f = _number(f, "frequency_MHz", lo=1e-6, hi=1e9)
        err = _number(err, "uncertainty_MHz", hi=f, include_zero=True)
        eu = _number(eu, "E_upper_K", hi=1e8)
        if eu + 1e-10 * max(eu, 1) < f * MHZ_TO_K:
            raise ValueError(f"E_upper_K implies a negative lower energy at transition {i}")
        a = _number(a, "A_s-1", lo=1e-300)
        gu = _int(gu, "g_upper", 1, 999999999)
        key = (f, _string(up, "upper QNs", display=True, qn=True), _string(low, "lower QNs", display=True, qn=True))
        if key in seen:
            raise ValueError(f"Duplicate transition at row {i}")
        seen.add(key)
        parsed.append((f, err, a, eu, gu, up, low))
    if not parsed:
        raise ValueError("Transition file is empty")
    parsed.sort(key=lambda r: r[0])
    return UploadedCatalog(meta_file, [meta_file, cat_file], meta, q,
                           tuple(np.asarray(col) for col in zip(*parsed)))


def read_uploads(directory):
    directory = Path(directory)
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("Expected an uploaded_catalogs directory")
    paths = sorted(directory.iterdir())
    if len(paths) > MAX_CATALOGS:
        raise ValueError("Too many local catalogs")
    uploads = []
    tags, labels = set(), set()
    total = 0
    for path in paths:
        if not re.fullmatch(r"[1-9][0-9]{0,9}", path.name):
            raise ValueError("Uploaded catalog folder must be its numeric LOCAL tag")
        item = _load(path)
        if path.name != str(item.tag) or item.tag in tags or item.label.casefold() in labels:
            raise ValueError("Conflicting local tag or label")
        tags.add(item.tag)
        labels.add(item.label.casefold())
        total += sum(f.stat().st_size for f in item.source_files)
        if total > MAX_SCOPE_BYTES:
            raise ValueError("Combined uploaded catalogs exceed supported size")
        uploads.append(item)
    return uploads


def stage_upload(run_catalog_dir, metadata_path, transitions_path):
    """Validate two local text files and atomically add a LOCAL tag to a run catalog directory."""
    directory = Path(run_catalog_dir)
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("Choose an existing run catalog directory")
    meta_data = _read(Path(metadata_path), MAX_META)
    meta, _ = _metadata(meta_data)
    filename = "transitions.csv" if meta["schema"] == CSV_SCHEMA else "transitions.cat"
    if Path(transitions_path).suffix.lower() != Path(filename).suffix:
        raise ValueError("Transition file extension must match its declared format")
    with TemporaryDirectory(dir=directory) as tmp:
        package = Path(tmp) / str(meta["tag"])
        package.mkdir()
        (package / "metadata.json").write_bytes(meta_data)
        (package / filename).write_bytes(_read(Path(transitions_path), MAX_BYTES))
        item = _load(package)
        # Serialize competing additions using the trusted run directory as the lock.
        fd = os.open(directory, os.O_RDONLY)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX)
            base = directory / "uploaded_catalogs"
            if base.is_symlink():
                raise ValueError("uploaded_catalogs cannot be a symlink")
            base.mkdir(exist_ok=True)
            active = read_uploads(base)
            if any(u.tag == item.tag or u.label.casefold() == item.label.casefold() for u in active):
                raise ValueError("Local catalog tag or label is already active")
            for name, key in (("all_cdms_final_official.csv", "mol"), ("all_jpl_final_official.csv", "name")):
                table = directory / name
                if table.is_file():
                    with table.open(encoding="utf-8-sig", newline="") as stream:
                        if any(row[key].casefold() == item.label.casefold() for row in csv.DictReader(stream)):
                            raise ValueError("Catalog label is already active in the native catalog")
            if len(active) >= MAX_CATALOGS or sum(p.stat().st_size for u in [*active, item] for p in u.source_files) > MAX_SCOPE_BYTES:
                raise ValueError("Combined uploaded catalogs exceed supported limits")
            target = base / package.name
            package.rename(target)
        finally:
            os.close(fd)
    return target
