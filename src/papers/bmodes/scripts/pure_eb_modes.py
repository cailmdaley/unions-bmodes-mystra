"""Pure E/B data vector from the measured ξ±.

The Schneider (2022) decomposition of the reporting and integration grids, in
the geometry ``--integration-method`` implies: pointwise on the reporting grid
under ``adaptive_dqags``, on the fine grid then averaged into the reporting bins
under ``fixed_operator``, with the weights ``--rebinning-weights`` selects.

    python pure_eb_modes.py \\
        --xi-reporting <xi 20-bin .txt> --xi-integration <xi 1000-bin .txt> \\
        --integration-method fixed_operator --integration-limits fine_grid_extent \\
        --rebinning-weights npairs --out-file <pure_eb_data.npz>
"""

import argparse

import numpy as np
from pure_eb_transform import (
    add_integration_limits_argument,
    add_integration_method_argument,
    add_rebinning_weights_argument,
    data_vector,
    load_xi,
)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--xi-reporting", required=True)
    ap.add_argument("--xi-integration", required=True)
    ap.add_argument("--min-sep", type=float, default=1.0)
    ap.add_argument("--max-sep", type=float, default=250.0)
    ap.add_argument("--nbins", type=int, default=20)
    ap.add_argument("--nbins-int", type=int, default=1000)
    ap.add_argument("--out-file", required=True, help="Output npz (lc {output})")
    add_integration_method_argument(ap)
    add_integration_limits_argument(ap)
    add_rebinning_weights_argument(ap)
    ap.add_argument(
        "--geometry",
        choices=["pointwise", "averaged"],
        default=None,
        help="Override the geometry the integration method defaults to.",
    )
    ap.add_argument(
        "--local-from-eval-grid",
        action="store_true",
        help="Pointwise fixed_operator only: take the local 0.5*(xi_p +/- xi_m) "
        "term from the reporting-grid xi at the requested abscissa rather than "
        "from the snapped fine node. The integral term still comes from the "
        "fine grid at the snapped node.",
    )
    a = ap.parse_args(argv)
    gg = load_xi(a.xi_reporting, a.nbins)
    gg_int = load_xi(a.xi_integration, a.nbins_int)
    xip_E, xim_E, xip_B, xim_B, xip_amb, xim_amb = data_vector(
        gg,
        gg_int,
        min_sep=a.min_sep,
        max_sep=a.max_sep,
        nbins=a.nbins,
        integration_method=a.integration_method,
        rebinning_weights=a.rebinning_weights,
        integration_limits=a.integration_limits,
        geometry=a.geometry,
        local_from_eval_grid=a.local_from_eval_grid,
    )
    np.savez(
        a.out_file,
        theta=gg["meanr"],
        theta_int=gg_int["meanr"],
        xip_total=gg["xip"],
        xim_total=gg["xim"],
        xip_E=xip_E,
        xim_E=xim_E,
        xip_B=xip_B,
        xim_B=xim_B,
        xip_amb=xip_amb,
        xim_amb=xim_amb,
    )
    print(f"Saved to {a.out_file}")


if __name__ == "__main__":
    main()
