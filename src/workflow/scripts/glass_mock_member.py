"""One GLASS realization, generated and measured, packed as one member.

Generates the catalog in a scratch directory (``glass_mock_generate.py``),
measures it with this analysis's recipes, and writes a tar.gz of the three
measurements under the names the mock consumers read; the catalog is deleted.

- ``xi_fine/gg_glass_mock_<id>_nbins=1000.fits``: TreeCorr ξ± on the integration
  grid (``run_glass_mock_2pcf.compute_one``);
- ``xi_reporting/gg_glass_mock_<id>_nbins=20.fits``: the same call on the
  20-bin reporting grid;
- ``cl/cl_glass_mock_<id>_4096.npy``: catalog-based NaMaster bandpowers on the
  data's 32 √ℓ-spaced bands, 8 ≤ ℓ ≤ 2048, rows ℓ_eff, EE, EB, BE, BB.

    python glass_mock_member.py --number N --summary SUMMARY --suite SUITE --fiducial-catalog CAT \
        --nz NZ --shell-cls C --sampling S --out-file member.tar.gz
"""

import argparse
import os
import tarfile
import tempfile
from pathlib import Path

import fitsio
import numpy as np

from decision_values import INTEGRATION_GRID, REPORTING_BINNING
from glass_mock_generate import footprint_mask, realization
from run_glass_mock_2pcf import compute_one


def mock_cl(catalog):
    """Decoupled EE, EB, BE, BB bandpowers of an unweighted mock shear catalog."""
    import pymaster as nmt

    d = fitsio.read(catalog, columns=["RA", "Dec", "e1", "e2"])
    lmin, lmax, nbins = 8, 2048, 32
    ells = np.arange(lmin, lmax + 1)
    edges = np.linspace(np.sqrt(lmin), np.sqrt(lmax), nbins + 1) ** 2
    bpws = np.digitize(ells.astype(float), edges) - 1
    bpws[0], bpws[-1] = 0, nbins - 1
    b = nmt.NmtBin(ells=ells, bpws=bpws, lmax=lmax - 1)
    field = nmt.NmtFieldCatalog(
        positions=[d["RA"] % 360, d["Dec"]],
        weights=np.ones(len(d)),
        field=[d["e1"], -d["e2"]],
        lmax=lmax - 1,
        lmax_mask=lmax - 1,
        spin=2,
        lonlat=True,
    )
    wsp = nmt.NmtWorkspace.from_fields(field, field, b)
    cl = wsp.decouple_cell(nmt.compute_coupled_cell(field, field))
    return np.concatenate([b.get_effective_ells()[None], cl])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--number", type=int, required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--suite", required=True)
    p.add_argument("--fiducial-catalog", required=True)
    p.add_argument("--nz", required=True)
    p.add_argument("--shell-cls", required=True)
    p.add_argument("--sampling", required=True)
    p.add_argument("--out-file", required=True)
    a = p.parse_args(argv)
    mid = f"{a.number:05d}"
    threads = int(os.environ.get("OMP_NUM_THREADS", "1") or 1)

    with tempfile.TemporaryDirectory(dir=os.environ.get("SLURM_TMPDIR")) as tmp:
        tmp = Path(tmp)
        catalog = tmp / f"unions_glass_sim_{mid}_4096.fits"
        mask = footprint_mask(a.summary, a.suite, a.fiducial_catalog)
        realization(a.number, mask, a.nz, a.shell_cls, a.sampling, catalog)
        names = []
        for sub, grid in (("xi_reporting", REPORTING_BINNING["bins_20"]), ("xi_fine", INTEGRATION_GRID["bins_1000"])):
            path = tmp / sub / f"gg_glass_mock_{mid}_nbins={grid.nbins}.fits"
            path.parent.mkdir()
            compute_one(catalog, path, min_sep=grid.min_sep, max_sep=grid.max_sep, nbins=grid.nbins, num_threads=threads)
            names.append(path)
        path = tmp / "cl" / f"cl_glass_mock_{mid}_4096.npy"
        path.parent.mkdir()
        np.save(path, mock_cl(catalog))
        names.append(path)
        catalog.unlink()
        with tarfile.open(a.out_file, "w:gz") as tar:
            for path in names:
                tar.add(path, arcname=str(path.relative_to(tmp)))


if __name__ == "__main__":
    main()
