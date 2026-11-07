"""Generate Cl covariances (no data vector).

The CLI runs the measurement from explicit flags. The primitive's native
``pseudo_cl_cov_{ver}.fits`` is produced in a private scratch directory and
moved to the single ``--out-file`` path the caller gives (the fresh scratch
also means the primitives' skip-if-exists never collides across nbins runs):

    python generate_cl_cov.py \
        --ver SP_v1.4.6.3_leak_corr \
        --cat-config /path/to/cosmo_val/cat_config.yaml \
        --out-file <path.fits> \
        --binning powspace --nbins 32 --power 0.5

Supports two binning modes:
- Linear binning with configurable nbins for COSEBIS
- Power-space binning with configurable nbins for standard C_ell analysis

See generate_cl.py for data vector generation.
"""

import argparse
import json
import os
import re
import shutil
import tempfile

import yaml

from decision_values import FIDUCIAL_BLIND
from astropy.io import fits

from sp_validation.cosmo_val import CosmologyValidation


def generate_cl_cov(
    version: str,
    output_dir: str,
    cat_config: str,
    nside: int = 1024,
    npatch: int = 1,
    blind: str = None,
    cosmo_params: dict = None,
    binning: str = "powspace",
    nbins: int = 32,
    power: float = 0.5,
):
    """Generate a Cl covariance into ``output_dir``.

    Parameters
    ----------
    version : str
        Catalog version (e.g., "SP_v1.4.6_leak_corr")
    output_dir : str
        Directory the covariance FITS file is written into. The primitive writes
        its native ``pseudo_cl_cov_{version}.fits`` here; callers that need a
        tagged filename rename it themselves.
    cat_config : str
        Path to catalog configuration YAML
    nside : int
        HEALPix nside for map-based estimation
    npatch : int
        Number of jackknife patches
    blind : str, optional
        Blind identifier (A, B, or C) to override n(z) path
    cosmo_params : dict, optional
        Cosmological parameters. Keys: Omega_m, sigma_8, n_s, h, Omega_b.
        If None, uses Planck 2018 defaults.
    binning : str
        Binning mode: "linear", "logspace", or "powspace"
    nbins : int
        Number of ell bins
    power : float
        Power for powspace binning (0.5 = sqrt spacing)

    Returns
    -------
    str
        Path to the primitive's native ``pseudo_cl_cov_{version}.fits`` product.
    """
    os.makedirs(output_dir, exist_ok=True)

    blind_str = f" blind={blind}" if blind else ""
    if binning == "linear":
        ell_step = max(1, (2048 - 2) // nbins)
        bin_str = f"nbins={nbins} (ell_step={ell_step})"
    elif binning == "logspace":
        bin_str = f"nbins={nbins} (geomspace)"
    else:
        bin_str = f"nbins={nbins}, power={power}"

    print(f"\n{'=' * 60}")
    print(f"Generating Cl covariance for {version}{blind_str}")
    print(f"Binning: {binning} ({bin_str})")
    if cosmo_params:
        print(
            f"Cosmology: Om={cosmo_params.get('Omega_m')}, s8={cosmo_params.get('sigma_8')}"
        )
    else:
        print("Cosmology: Planck 2018 defaults")
    print(f"{'=' * 60}\n")

    # Remap config cosmology param names to get_cosmo() expected names
    if cosmo_params:
        cosmo_params = {
            "Omega_m": cosmo_params.get("Omega_m"),
            "Omega_b": cosmo_params.get("Omega_b"),
            "h": cosmo_params.get("h"),
            "sig8": cosmo_params.get("sigma_8"),
            "ns": cosmo_params.get("n_s"),
        }

    # CosmologyValidation.calculate_pseudo_cl_eb_cov reads the theory n(z) from
    # cat_config's `redshift_path` and ignores `blind`, so the blind is applied
    # here, to a copy of the config.
    if blind:
        with open(cat_config) as f:
            cc = yaml.safe_load(f)
        entry = cc[version.removesuffix("_leak_corr")]["shear"]
        entry["redshift_path"] = re.sub(r"_[ABC]\.txt$", f"_{blind}.txt", entry["redshift_path"])
        cat_config = os.path.join(output_dir, "cat_config_blind.yaml")
        with open(cat_config, "w") as f:
            yaml.safe_dump(cc, f)
        print(f"Theory n(z): {entry['redshift_path']}")

    # Build CV kwargs based on binning mode
    cv_kwargs = dict(
        versions=[version],
        catalog_config=cat_config,
        output_dir=output_dir,
        binning=binning,
        nside=nside,
        cell_method="catalog",
        nrandom_cell=100,
        blind=blind,
        cosmo_params=cosmo_params,
        npatch=npatch,
        theta_min=1.0,
        theta_max=250.0,
        nbins=20,
    )
    if binning == "linear":
        cv_kwargs["ell_step"] = ell_step
    else:
        cv_kwargs["n_ell_bins"] = nbins
        cv_kwargs["power"] = power

    cv = CosmologyValidation(**cv_kwargs)

    # Calculate covariance only
    print("Calculating covariance...")
    cv.calculate_pseudo_cl_eb_cov()

    # Report on the native product
    src_cov = os.path.join(output_dir, f"pseudo_cl_cov_{version}.fits")
    if os.path.exists(src_cov):
        with fits.open(src_cov) as hdul:
            # CV outputs covariance blocks as COVAR_XX_YY extensions
            cov = hdul["COVAR_BB_BB"].data
            n_ell = int(len(cov) ** 0.5)
            print(f"Generated covariance matrix: {n_ell}x{n_ell}")
    return src_cov


def _from_snakemake(smk):
    p = smk.params
    output_cov = smk.output.cl_cov
    src_cov = generate_cl_cov(
        version=p["version"],
        output_dir=os.path.dirname(output_cov),
        cat_config=p["cat_config"],
        nside=int(p["nside"]),
        npatch=int(p["npatch"]),
        blind=p.get("blind", None),
        cosmo_params=p.get("cosmo_params", None),
        binning=p["binning"],
        nbins=int(p["nbins"]),
        power=float(p.get("power", 0.5)),
    )
    # Rename the native product to the tagged output filename.
    if os.path.exists(src_cov) and src_cov != output_cov:
        os.rename(src_cov, output_cov)
        print(f"Saved to: {output_cov}")


def _from_cli(argv=None):
    ap = argparse.ArgumentParser(
        description="Cl covariance for one catalog version."
    )
    ap.add_argument(
        "--ver",
        required=True,
        help="Catalog version key in cat_config, e.g. SP_v1.4.6.3_leak_corr",
    )
    ap.add_argument(
        "--cat-config", required=True, help="Absolute path to cat_config.yaml"
    )
    ap.add_argument(
        "--out-file", required=True, help="Output FITS file path (lc {output})"
    )
    ap.add_argument(
        "--nside", type=int, default=1024, help="HEALPix nside for map estimation"
    )
    ap.add_argument(
        "--npatch",
        type=int,
        default=1,
        help="Jackknife patch count (paper fiducial: 1)",
    )
    ap.add_argument(
        "--binning",
        choices=["linear", "logspace", "powspace"],
        default="powspace",
        help="Ell binning mode",
    )
    ap.add_argument("--nbins", type=int, required=True,
                    help="Number of ell bins")
    ap.add_argument(
        "--power",
        type=float,
        default=0.5,
        help="Power for powspace binning (0.5 = sqrt spacing)",
    )
    ap.add_argument("--blind", type=FIDUCIAL_BLIND.__getitem__)
    ap.add_argument(
        "--cosmo-json",
        default=None,
        help="Path to a Planck18-style cosmology JSON; omit for Planck18 defaults",
    )
    a = ap.parse_args(argv)

    cosmo_params = None
    if a.cosmo_json:
        with open(a.cosmo_json) as f:
            cosmo_params = json.load(f)

    scratch = tempfile.mkdtemp(prefix="cl_cov_")
    src_cov = generate_cl_cov(
        version=a.ver,
        output_dir=scratch,
        cat_config=a.cat_config,
        nside=a.nside,
        npatch=a.npatch,
        blind=a.blind,
        cosmo_params=cosmo_params,
        binning=a.binning,
        nbins=a.nbins,
        power=a.power,
    )
    shutil.move(src_cov, a.out_file)
    print(f"Saved to: {a.out_file}")
    shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    if "snakemake" in globals():
        _from_snakemake(snakemake)  # noqa: F821
    else:
        _from_cli()
