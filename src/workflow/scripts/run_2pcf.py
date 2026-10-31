"""TreeCorr ξ±(θ) two-point correlation function for one catalog version.

Dual-mode. Under Snakemake (``script:`` directive) the injected ``snakemake``
object supplies the parameters; as a standalone CLI (argparse) the same compute
runs from explicit flags. The CLI form is what the lightcone/ASTRA recipe calls,
so the measurement is driven directly (no nested Snakemake) with lc handling
orchestration:

    python run_2pcf.py \
        --ver SP_v1.4.6.3_leak_corr \
        --min-sep 1.0 --max-sep 250.0 --nbins 20 --npatch 1 \
        --cat-config /path/to/cosmo_val/cat_config.yaml \
        --out <output_dir>

The measurement itself is unchanged — ``CosmologyValidation.calculate_2pcf``
does the TreeCorr work and writes the ``.txt`` dump plus ξ+/ξ- FITS files into
``output_dir``. ``output_dir`` is passed explicitly (rather than via the
``COSMO_VAL`` env hook) so lc can point each run at its own ``{output}`` tree.
"""

import argparse
import functools
import os
import shutil
import tempfile
from pathlib import Path

import treecorr
from sp_validation.cosmo_val import CosmologyValidation


def _cap_treecorr_threads():
    """Cap TreeCorr's OpenMP pool at the cores this task actually owns.

    TreeCorr resolves its default ``num_threads=None`` through
    ``multiprocessing.cpu_count()`` (treecorr/util.py) and then calls
    ``omp_set_num_threads()`` itself, so it sizes from the NODE's core count
    and *overrides* ``OMP_NUM_THREADS``. Under lc's SLURM venue the task slot
    count already equals the allocated cores, so a task that grabs the whole
    node oversubscribes by construction: 16 slots x 48 threads produced a load
    average of 338 on a 16-core allocation, and everything on the node was
    time-sliced to a fraction of a core.

    ``set_max_omp_threads`` is a ceiling applied inside ``set_omp_threads``
    regardless of what any caller passes down, so it holds even though the
    TreeCorr call itself is inside pinned sp_validation code. Read from
    ``OMP_NUM_THREADS`` rather than a new knob: that is the variable the
    runner already sets, and this makes TreeCorr honour it like every other
    library. Default 1 = one core per task slot.
    """
    treecorr.set_max_omp_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))


def _skip_psf_catalogue(cv, ver):
    """Load only the shear catalogue for this version's ξ± measurement.

    ``calculate_2pcf`` opens its data through ``temporarily_read_data()``,
    whose ``psf=True`` default also reads the star/PSF catalogue for the
    PSF-leakage statistics. The ξ± measurement itself never touches it: it
    takes RA/Dec, the calibrated shears and the weights from ``dat_shear``
    alone.

    Under ``--out-file`` the leakage side products are discarded with the
    scratch directory and nothing downstream consumes them, so that read is
    pure cost — a second large FITS file, and a dependency on a catalogue
    the declared inputs have no reason to carry. Bound per instance, so the
    legacy directory mode keeps reading it.
    """
    runner = cv.results[ver]
    runner.temporarily_read_data = functools.partial(
        runner.temporarily_read_data, psf=False
    )


def run_2pcf(
    ver,
    min_sep,
    max_sep,
    nbins,
    npatch,
    cat_config,
    output_dir,
    save_fits=True,
    read_psf=True,
):
    """Measure ξ±(θ) for ``ver`` and write it under ``output_dir``.

    Parameters mirror the TreeCorr reporting/integration grids: ``min_sep`` /
    ``max_sep`` in arcmin, ``nbins`` logarithmic bins, ``npatch`` spatial
    patches (1 for the paper fiducial). ``cat_config`` is an absolute path to
    the catalog configuration; ``output_dir`` overrides
    ``cat_config['paths']['output']`` so products land where lc expects.
    ``read_psf=False`` skips the PSF-leakage catalogue read the ξ± itself
    does not need (see :func:`_skip_psf_catalogue`).
    """
    _cap_treecorr_threads()
    cv = CosmologyValidation(
        versions=[ver],
        catalog_config=cat_config,
        output_dir=output_dir,
    )
    if not read_psf:
        _skip_psf_catalogue(cv, ver)
    return cv.calculate_2pcf(
        ver=ver,
        npatch=npatch,
        save_fits=save_fits,
        min_sep=min_sep,
        max_sep=max_sep,
        nbins=nbins,
    )


def _from_snakemake(smk):
    p = smk.params
    run_2pcf(
        ver=p["ver"],
        min_sep=float(p["min_sep"]),
        max_sep=float(p["max_sep"]),
        nbins=int(p["nbins"]),
        npatch=int(p["npatch"]),
        # cat_config / output_dir were previously resolved via an os.chdir into
        # the cosmo_val dir + the COSMO_VAL env var; expose them as optional
        # params so the rule can pass them explicitly, falling back to the
        # class defaults (./cat_config.yaml, COSMO_VAL env) otherwise.
        cat_config=p.get("cat_config", "./cat_config.yaml"),
        output_dir=p.get("output_dir", None),
        save_fits=True,
    )


def _from_cli(argv=None):
    ap = argparse.ArgumentParser(
        description="TreeCorr ξ± 2PCF for one catalog version."
    )
    ap.add_argument(
        "--ver",
        required=True,
        help="Catalog version key in cat_config, e.g. SP_v1.4.6.3_leak_corr",
    )
    ap.add_argument(
        "--min-sep", type=float, required=True, help="Min separation [arcmin]"
    )
    ap.add_argument(
        "--max-sep", type=float, required=True, help="Max separation [arcmin]"
    )
    ap.add_argument("--nbins", type=int, required=True, help="Number of log bins")
    ap.add_argument(
        "--npatch", type=int, default=1, help="TreeCorr patch count (paper fiducial: 1)"
    )
    ap.add_argument(
        "--cat-config", required=True, help="Absolute path to cat_config.yaml"
    )
    out_group = ap.add_mutually_exclusive_group(required=True)
    out_group.add_argument("--out", help="Output directory (lc {output})")
    out_group.add_argument(
        "--out-file",
        help="Write the TreeCorr text dump to exactly this path "
        "(no directory, no siblings — FITS/patches/leakage side "
        "products are confined to a scratch dir and discarded).",
    )
    ap.add_argument("--no-fits", action="store_true", help="Skip ξ+/ξ- FITS export")
    a = ap.parse_args(argv)

    if a.out_file is not None:
        scratch = tempfile.mkdtemp(prefix="run_2pcf_")
        try:
            run_2pcf(
                ver=a.ver,
                min_sep=a.min_sep,
                max_sep=a.max_sep,
                nbins=a.nbins,
                npatch=a.npatch,
                cat_config=a.cat_config,
                output_dir=scratch,
                save_fits=not a.no_fits,
                read_psf=False,
            )
            produced = (
                Path(scratch)
                / f"{a.ver}_xi_minsep={a.min_sep}_maxsep={a.max_sep}"
                f"_nbins={a.nbins}_npatch={a.npatch}.txt"
            )
            Path(a.out_file).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(produced, a.out_file)
        finally:
            shutil.rmtree(scratch, ignore_errors=True)
        return

    run_2pcf(
        ver=a.ver,
        min_sep=a.min_sep,
        max_sep=a.max_sep,
        nbins=a.nbins,
        npatch=a.npatch,
        cat_config=a.cat_config,
        output_dir=a.out,
        save_fits=not a.no_fits,
    )


if __name__ == "__main__":
    # A membership test, not a caught NameError: the exception form leaves
    # every CLI-path traceback prefixed with "During handling of the above
    # exception" and a NameError that is not the failure, which reads as
    # the bug and is not one.
    if "snakemake" in globals():
        _from_snakemake(snakemake)  # noqa: F821
    else:
        _from_cli()
