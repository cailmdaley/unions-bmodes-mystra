"""Run treecorr on GLASS mock galaxy catalogs.

Computes fine-binned ξ±(θ) for mock validation. Each output is a treecorr
GGCorrelation FITS file readable via gg.read(path).

Mock catalogs have columns: RA, Dec, e1, e2, w.
No response correction or PSF leakage subtraction (these are Gaussian mocks).

Vendored from sp_validation `workflow/scripts/run_glass_mock_2pcf.py`. The
upstream file is a bare snakemake script running one mock per job; the body is
unchanged apart from being wrapped in `compute_one` and given a CLI. The lc
recipes call the single-mock form (`--mock-id` + `--out-file`), one sweep
member per realisation, which is the same fan-out snakemake did.
"""

import argparse
import os
from pathlib import Path

import treecorr
from astropy.io import fits

from decision_values import INTEGRATION_GRID

# v1.4.8 = v1.4.6.3 + the faint/bright star-halo cuts (masks.txt in
# /n17data/UNIONS/WL/v1.4.x/v1.4.8/).  mask_r_nside131072.hsp is a bit-packed
# BOOLEAN HealSparse map, True where any structural bit (1|2|4|8|64|1024)
# masks the sky; validated 2026-08-23 against the six per-bit _n{bit}.hsp
# maps (99.74% pixel agreement) and the real catalogs: keeps 85.4% of
# v1.4.6.3 rows vs the real v1.4.8's 87.8% (residual = pixelized mask vs
# per-object flags; bits 4/8/64/1024 confirmed already absent from
# v1.4.6.3).  Out-of-coverage returns the sentinel False = "no hole" = keep.


def apply_v148_mask(ra, dec, mask_path):
    """Return the coordinates' keep mask for the v1.4.8 footprint."""
    import healsparse as hsp
    import numpy as np

    survey_mask = hsp.HealSparseMap.read(mask_path)
    mask_values = survey_mask.get_values_pos(ra, dec, lonlat=True)
    if mask_values.dtype != np.bool_:
        raise TypeError(
            f"{mask_path} expected boolean union mask, got {mask_values.dtype}"
        )
    return ~mask_values


def _default_num_threads():
    """TreeCorr's own default is every core on the node, which oversubscribes
    a task that was given a slice of one.  Honour OMP_NUM_THREADS instead."""
    try:
        return max(1, int(os.environ.get("OMP_NUM_THREADS", "1")))
    except ValueError:
        return 1


def compute_one(
    catalog_path,
    output_path,
    *,
    min_sep,
    max_sep,
    nbins,
    mask_file=None,
    num_threads=None,
):
    """Run treecorr on one mock catalog and write the GGCorrelation FITS."""

    with fits.open(catalog_path) as hdul:
        data = hdul["SOURCE_CATALOGUE"].data
        ra = data["RA"]
        dec = data["Dec"]
        e1 = data["e1"]
        e2 = data["e2"]
        w = data["w"]

    if mask_file:
        keep = apply_v148_mask(ra, dec, mask_file)
        n_before = len(ra)
        ra, dec, e1, e2, w = (x[keep] for x in (ra, dec, e1, e2, w))
        print(
            f"Applied v1.4.8 mask {mask_file}: kept {len(ra):,}/{n_before:,} "
            f"galaxies ({len(ra) / n_before:.3%})"
        )
    else:
        print(f"Loaded {len(ra)} galaxies from {catalog_path}")

    cat = treecorr.Catalog(
        ra=ra,
        dec=dec,
        g1=e1,
        g2=e2,
        w=w,
        ra_units="degrees",
        dec_units="degrees",
    )

    gg = treecorr.GGCorrelation(
        min_sep=min_sep,
        max_sep=max_sep,
        nbins=nbins,
        sep_units="arcmin",
    )

    if num_threads is None:
        num_threads = _default_num_threads()
    print(
        f"Running treecorr: {nbins} bins, [{min_sep}, {max_sep}] arcmin, "
        f"{num_threads} threads..."
    )
    gg.process(cat, num_threads=num_threads)
    gg.write(output_path)
    print(f"Saved to {output_path}")
    return output_path


def _parser(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog-dir",
        required=True,
        help="Directory holding unions_glass_sim_<id>_4096.fits realizations.",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory for gg_glass_mock_<id>_nbins=<nbins>.fits (batch mode).",
    )
    parser.add_argument(
        "--mock-id",
        help="Single realization to process; requires --out-file.",
    )
    parser.add_argument(
        "--out-file",
        help="Exact FITS path for the --mock-id realization.",
    )
    parser.add_argument(
        "--n-mocks",
        type=int,
        default=350,
        help="Realizations 1..N to process (default: the full 350-mock suite).",
    )
    parser.add_argument(
        "--mock-ids", nargs="*", help="Explicit five-digit IDs, overriding --n-mocks."
    )
    parser.add_argument(
        "--integration-grid", type=INTEGRATION_GRID.__getitem__, required=True
    )
    parser.add_argument(
        "--num-threads",
        type=int,
        default=None,
        help="TreeCorr OpenMP threads (default: $OMP_NUM_THREADS, else 1).",
    )
    parser.add_argument(
        "--mask-file",
        default=None,
        help="HealSparse v1.4.8 union mask; omit for the unmasked footprint.",
    )
    args = parser.parse_args(argv)
    single = args.mock_id is not None
    if single != (args.out_file is not None):
        parser.error("--mock-id and --out-file must be given together")
    if not single and not args.output_dir:
        parser.error("batch mode needs --output-dir")
    return args


def _standalone_main(argv=None):
    args = _parser(argv)
    grid = args.integration_grid
    catalog_dir = Path(args.catalog_dir)

    if args.mock_id is not None:
        mock_id = str(args.mock_id).zfill(5)
        out_path = Path(args.out_file)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        compute_one(
            catalog_dir / f"unions_glass_sim_{mock_id}_4096.fits",
            out_path,
            min_sep=grid.min_sep,
            max_sep=grid.max_sep,
            nbins=grid.nbins,
            mask_file=args.mask_file,
            num_threads=args.num_threads,
        )
        return


    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ids = args.mock_ids or [f"{i:05d}" for i in range(1, args.n_mocks + 1)]
    for raw_id in ids:
        mock_id = str(raw_id).zfill(5)
        catalog_path = catalog_dir / f"unions_glass_sim_{mock_id}_4096.fits"
        output_path = output_dir / f"gg_glass_mock_{mock_id}_nbins={grid.nbins}.fits"
        if output_path.is_file():
            continue
        if not catalog_path.is_file():
            print(f"Skipping {mock_id}: missing catalog {catalog_path}")
            continue
        compute_one(
            catalog_path,
            output_path,
            min_sep=grid.min_sep,
            max_sep=grid.max_sep,
            nbins=grid.nbins,
            mask_file=args.mask_file,
            num_threads=args.num_threads,
        )


def _snakemake_main(smk):
    compute_one(
        smk.input.catalog,
        smk.output.gg,
        min_sep=smk.params.min_sep,
        max_sep=smk.params.max_sep,
        nbins=smk.params.nbins,
        mask_file=getattr(smk.params, "mask_file", None),
    )


if "snakemake" in globals():
    _snakemake_main(snakemake)  # noqa: F821
elif __name__ == "__main__":
    _standalone_main()
