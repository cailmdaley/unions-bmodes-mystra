"""Compare rerun defaults with the pre-wiring snapshot (run on Slurm).

Floating arrays must agree to max(abs(new-old))/max(abs(old)) <= 1e-12.
Explicitly audited, unchanged computational paths may agree to 1e-4 across CPU
architectures (PASS-arch); see work/notes/decision-wiring-proof.md.
Archive timestamps are ignored; explicitly renamed members are paired by name.
No comparison is made against a recomputed reference.
"""
import argparse
import csv
import hashlib
import io
import json
import re
from pathlib import Path
import tarfile
import zlib

import numpy as np


# No science changes on these paths: only option parsing / config selection /
# neutral archive filenames. The audit records the source diffs and CPU evidence.
ARCH_UNCHANGED = {
    "submitted/pure_eb/pure_eb_semianalytic_data.npz",
    "submitted/pure_eb/pure_eb_pte_per_cut.npz",
    "submitted/pure_eb/pure_eb_sweep.tar.gz",
    "submitted/pure_eb/pure_eb_ptes.tar.gz",
    "submitted/glass_mocks/pure_cov_mockmatched.npz",
}
ARCH_REASON = ("Intel n23 versus AMD n31-33; producing numerical path unchanged "
               "(flag/config/name plumbing only); source audit: "
               "work/notes/decision-wiring-proof.md; float arrays <=1e-4")


# Renames required because the old names encoded a now-variable decision.
def old_name(name):
    name = re.sub(r"(pseudo_cl(?:_cov)?_sp_v[^/]+)_fine\.", r"\1_leak_corr_nbins96.", name)
    if "/cosebis_harmonic_scan/" in name and "_leak_corr." not in name:
        name = re.sub(r"(harmonic_modes_sp_v[^/.]+)(\.)", r"\1_leak_corr\2", name)
    return (name.replace("cov_integration.", "cov_integration_g.")
            .replace("_integration/", "_g_minsep=0.5_maxsep=300.0_nbins=1000_masked/")
            .replace("_integration_processed.txt", "_g_minsep=0.5_maxsep=300.0_nbins=1000_masked_processed.txt")
            .replace("fiducial_counterpart", "sp_v1p4p6p3"))


def digest(data):
    return hashlib.sha256(data).hexdigest()


def array_delta(old, new, tolerance=1e-12):
    if old.shape != new.shape or old.dtype != new.dtype:
        raise ValueError(f"array shape/dtype changed: {old.shape}/{old.dtype} -> {new.shape}/{new.dtype}")
    if old.dtype.kind in "fc":
        if not np.array_equal(np.isnan(old), np.isnan(new)) or not np.array_equal(np.isinf(old), np.isinf(new)):
            raise ValueError("nonfinite entries changed")
        finite = np.isfinite(old)
        if not np.array_equal(old[~finite], new[~finite], equal_nan=True):
            raise ValueError("nonfinite entries changed")
        scale = float(np.max(np.abs(old[finite]), initial=0))
        delta = float(np.max(np.abs(new[finite] - old[finite]), initial=0))
        ratio = delta / scale if scale else (0.0 if not delta else float("inf"))
        if ratio > tolerance:
            raise ValueError(f"float-array relative delta {ratio:.17g} exceeds {tolerance:g}")
        return ratio
    if not np.array_equal(old, new):
        raise ValueError("nonfloating array changed")
    return 0.0


def pdf_payload(data):
    """Exact Matplotlib PDF objects, without writing-time serialization noise.

    Drawing commands, coordinates, text, font programs and trailer stay exact.
    Only dates, generated font-subset names and compression/xref bookkeeping
    differ in the observed files; no rendered-pixel tolerance is used.
    """
    def metadata(value):
        value = re.sub(rb"/[A-Z]{6}\+", b"/SUBSET+", value)
        return re.sub(rb"/(?:CreationDate|ModDate)\s*\([^)]*\)", b"", value)

    objects = []
    for match in re.finditer(rb"(?m)^(\d+ \d+) obj\n(.*?)\nendobj", data, re.S):
        body = match[2]
        if b"stream\n" in body:
            header, tail = body.split(b"stream\n", 1)
            stream, end = tail.rsplit(b"\nendstream", 1)
            if b"/FlateDecode" not in header:
                raise ValueError("unsupported PDF stream encoding")
            stream = zlib.decompress(stream)
            header = re.sub(rb"/Length (\d+)(\s+\d+\s+R)?",
                            lambda m: m[0] if m[2] else b"/Length #", header)
            body = metadata(header) + b"stream\n" + metadata(stream) + b"\nendstream" + metadata(end)
        else:
            body = metadata(body)
        objects.append((match[1], body))
    trailer = re.search(rb"trailer\s*(<<.*?>>)", data, re.S)
    if not objects or trailer is None:
        raise ValueError("unsupported PDF object layout")
    return data.split(b"\n", 2)[:2], objects, metadata(trailer[1])


def compare(old, new, name, tolerance=1e-12):
    if old == new:
        return "byte_identical", 0.0
    if name.endswith(".json"):
        a, b = json.loads(old), json.loads(new)
        af, bf = [], []
        def walk(x, y, path=""):
            if isinstance(x, dict) and isinstance(y, dict):
                if set(x) != set(y):
                    raise ValueError(f"JSON keys changed at {path}")
                for k in x:
                    if k == "generated":  # wall-clock timestamp written by the evidence emitters
                        continue
                    walk(x[k], y[k], path + "/" + k)
            elif isinstance(x, list) and isinstance(y, list):
                if len(x) != len(y):
                    raise ValueError(f"JSON list length changed at {path}")
                for i, (xx, yy) in enumerate(zip(x, y)):
                    walk(xx, yy, path + f"/{i}")
            elif isinstance(x, float) and isinstance(y, float):
                af.append(x); bf.append(y)
            elif x != y:
                raise ValueError(f"nonfloating JSON value changed at {path}: {x!r} -> {y!r}")
        walk(a, b)
        ratio = array_delta(np.asarray(af), np.asarray(bf), tolerance)
        return "json_float_payload_within_1e-12" if ratio else "json_payload_identical", ratio
    if name.endswith(".npz"):
        with np.load(io.BytesIO(old), allow_pickle=False) as a, np.load(io.BytesIO(new), allow_pickle=False) as b:
            if set(a.files) != set(b.files):
                raise ValueError("NPZ keys changed")
            ratios = []
            for k in a.files:
                if k == "provenance" and a[k].shape == b[k].shape == ():
                    # This scalar is JSON metadata, including the same generated
                    # wall-clock field handled in standalone JSON evidence.
                    _, r = compare(str(a[k]).encode(), str(b[k]).encode(), "provenance.json", tolerance)
                else:
                    r = array_delta(a[k], b[k], tolerance)
                ratios.append(r)
            ratio = max(ratios, default=0.0)
        return "float_arrays_within_1e-12" if ratio else "array_payload_identical", ratio
    if name.endswith(".tar.gz"):
        with tarfile.open(fileobj=io.BytesIO(old), mode="r:gz") as a, tarfile.open(fileobj=io.BytesIO(new), mode="r:gz") as b:
            am = {m.name: m for m in a.getmembers() if m.isfile()}
            bm = {old_name(m.name): m for m in b.getmembers() if m.isfile()}
            if set(am) != set(bm):
                raise ValueError(f"archive members changed: {set(am) ^ set(bm)}")
            ratio = 0.0
            for k in am:
                _, r = compare(a.extractfile(am[k]).read(), b.extractfile(bm[k]).read(), k, tolerance)
                ratio = max(ratio, r)
        return "archive_float_arrays_within_1e-12" if ratio else "archive_payload_identical", ratio
    if name.endswith(".txt"):
        # TreeCorr text has a main table and a covariance block of a different
        # width. Compare each array separately, retaining every metadata line.
        def text_blocks(data):
            headers, blocks, lines = [], [], []
            for line in data.splitlines():
                if line.lstrip().startswith(b"#"):
                    if lines:
                        blocks.append(np.loadtxt(io.BytesIO(b"\n".join(lines)), ndmin=2))
                        lines = []
                    headers.append(line)
                elif line.strip():
                    lines.append(line)
            if lines:
                blocks.append(np.loadtxt(io.BytesIO(b"\n".join(lines)), ndmin=2))
            return headers, blocks
        try:
            ah, a = text_blocks(old)
            bh, b = text_blocks(new)
        except ValueError:
            pass
        else:
            if ah != bh or len(a) != len(b):
                raise ValueError("text metadata or block count changed")
            if a:
                ratios = []
                for x, y in zip(a, b):
                    if x.shape != y.shape:
                        raise ValueError("text array shape changed")
                    # Do not let pair counts (~1e15) hide shear changes (~1e-5).
                    ratios.extend(array_delta(xx, yy, tolerance) for xx, yy in zip(x.T, y.T))
                return "float_arrays_within_1e-12", max(ratios)
    if name.endswith(".pdf"):
        if pdf_payload(old) != pdf_payload(new):
            raise ValueError("PDF content changed beyond dates/font subset identifiers")
        return "pdf_payload_identical", 0.0
    if name.endswith(".fits"):
        from astropy.io import fits
        with fits.open(io.BytesIO(old)) as a, fits.open(io.BytesIO(new)) as b:
            if len(a) != len(b):
                raise ValueError("FITS extensions changed")
            ratio = 0.0
            for ah, bh in zip(a, b):
                # DATE is the writing timestamp, not an analysis result.
                for h in (ah.header, bh.header):
                    h.remove("DATE", ignore_missing=True, remove_all=True)
                if str(ah.header) != str(bh.header):
                    raise ValueError("FITS metadata changed")
                if ah.data is not None:
                    if ah.data.dtype.names:
                        for k in ah.data.dtype.names:
                            ratio = max(ratio, array_delta(ah.data[k], bh.data[k], tolerance))
                    else:
                        ratio = max(ratio, array_delta(ah.data, bh.data, tolerance))
        return "float_arrays_within_1e-12" if ratio else "array_payload_identical", ratio
    raise ValueError("non-array bytes changed")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--before", type=Path, default=Path("_scratch/decision-wiring-before"))
    ap.add_argument("--results", type=Path, default=Path("results"))
    ap.add_argument("--out", type=Path, default=Path("work/notes/decision-wiring-identity.csv"))
    ap.add_argument("--rerun-targets", type=Path,
                    help="Explicit universe/output IDs; require every rerun in the CSV even if its manifest is unchanged")
    args = ap.parse_args()
    required = set(args.rerun_targets.read_text().split()) if args.rerun_targets else set()
    seen = set()
    rows = []
    failed = False
    for universe in ("baseline", "submitted"):
        for manifest in sorted((args.results / universe).rglob(".*.manifest.json")):
            rel = manifest.relative_to(args.results)
            old_rel = Path(old_name(str(rel)))
            ref_manifest = args.before / old_rel
            stem = manifest.name[1:-len(".manifest.json")]
            target = universe + "/" + ".".join((*rel.parent.parts[1:], stem))
            unchanged = ref_manifest.exists() and ref_manifest.read_bytes() == manifest.read_bytes()
            protected = "/mock_xi_scan/" in str(rel) or manifest.name == ".mock_xi_reporting.manifest.json"
            if protected and (not unchanged or target in required):
                raise ValueError(f"protected mock measurement changed or requested: {manifest}")
            if unchanged and target not in required:
                continue
            seen.add(target)
            # Manifest and artifact share a stem; the manifest records the format.
            matches = [p for p in manifest.parent.glob(stem + ".*") if not p.name.startswith(".")]
            if len(matches) != 1:
                raise ValueError(f"cannot identify artifact for {manifest}: {matches}")
            artifact = matches[0]
            artifact_rel = artifact.relative_to(args.results)
            reference = args.before / old_name(str(artifact_rel))
            row = dict(universe=universe, output=str(artifact_rel.relative_to(universe)),
                       original_sha256="", current_sha256="", max_relative_delta="", verdict="", detail="")
            new = artifact.read_bytes()
            row["current_sha256"] = digest(new)
            if not reference.exists():
                row["verdict"] = "new_output"
            else:
                old = reference.read_bytes()
                row["original_sha256"] = digest(old)
                try:
                    verdict, ratio = compare(old, new, str(artifact_rel))
                    row.update(verdict=verdict, max_relative_delta=f"{ratio:.17g}")
                except ValueError as exc:
                    if str(artifact_rel) in ARCH_UNCHANGED:
                        try:
                            _, ratio = compare(old, new, str(artifact_rel), tolerance=1e-4)
                            row.update(verdict="PASS-arch", max_relative_delta=f"{ratio:.17g}",
                                       detail=ARCH_REASON)
                        except ValueError as arch_exc:
                            row.update(verdict="FAIL", detail=str(arch_exc))
                            failed = True
                    else:
                        row.update(verdict="FAIL", detail=str(exc))
                        failed = True
            rows.append(row)
    if required - seen:
        raise ValueError(f"rerun targets missing from identity proof: {sorted(required - seen)}")
    with args.out.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=["universe", "output", "original_sha256", "current_sha256", "max_relative_delta", "verdict", "detail"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    from collections import Counter
    print(dict(Counter(row["verdict"] for row in rows)))
    for row in rows:
        if row["verdict"] == "FAIL":
            print(row["universe"], row["output"], row["detail"])
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
