#!/usr/bin/env python

import sys

import matplotlib.pyplot as plt
import numpy as np


def get_cov(filename):

    data = np.loadtxt(filename)
    ndata = int(np.max(data[:, 0])) + 1

    print("Dimension of cov: %dx%d" % (ndata, ndata))

    cov_g = np.zeros((ndata, ndata))
    cov_ng = np.zeros((ndata, ndata))
    for i in range(0, data.shape[0]):
        cov_g[int(data[i, 0]), int(data[i, 1])] = data[i, 8]
        cov_g[int(data[i, 1]), int(data[i, 0])] = data[i, 8]
        cov_ng[int(data[i, 0]), int(data[i, 1])] = data[i, 9]
        cov_ng[int(data[i, 1]), int(data[i, 0])] = data[i, 9]

    return cov_g, cov_ng, ndata


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print(
            "Usage: python cosmocov_process.py <input_file> <output_stub>"
            " [<plot_path>]"
        )
        sys.exit(1)

    covfile = sys.argv[1]
    output_base = sys.argv[2]
    # Optional explicit plot path; a .png here is rendered at 300 dpi (the
    # project figure convention) instead of the legacy 2000-dpi PDF.
    plot_arg = sys.argv[3] if len(sys.argv) == 4 else None

    c_g, c_ng, ndata = get_cov(covfile)

    cov = c_ng + c_g
    cov_g = c_g

    b = np.sort(np.linalg.eigvals(cov))
    print("min+max eigenvalues cov: %e, %e" % (np.min(b), np.max(b)))
    if np.min(b) <= 0.0:
        print("non-positive eigenvalue encountered! Covariance Invalid!")
        exit()

    print("Covariance is positive definite!")

    np.savetxt(str(output_base) + ".txt", cov)
    print("covmat saved as %s" % (str(output_base) + ".txt"))

    np.savetxt(str(output_base) + "_g.txt", cov_g)
    print("Gaussian covmat saved as %s" % (str(output_base) + "_g.txt"))

    cmap = "seismic"

    pp_norm = np.zeros((ndata, ndata))
    for i in range(ndata):
        for j in range(ndata):
            pp_norm[i][j] = cov[i][j] / np.sqrt(cov[i][i] * cov[j][j])

    print("Plotting correlation matrix ...")

    plot_path = plot_arg if plot_arg else str(output_base) + "_plot.pdf"
    plot_dpi = 300 if str(plot_path).lower().endswith(".png") else 2000
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    extent = (0, ndata, ndata, 0)
    im3 = ax.imshow(pp_norm, cmap=cmap, vmin=-1, vmax=1, extent=extent)

    plt.axvline(x=int(ndata / 2), color="black", linewidth=1.0)
    plt.axhline(y=int(ndata / 2), color="black", linewidth=1.0)

    fig.colorbar(im3, orientation="vertical")

    ax.text(int(ndata / 4), ndata + 5, r"$\xi_+^{ij}(\theta)$", fontsize=12)
    ax.text(3 * int(ndata / 4), ndata + 5, r"$\xi_-^{ij}(\theta)$", fontsize=12)
    ax.text(-9, int(ndata / 4), r"$\xi_+^{ij}(\theta)$", fontsize=12)
    ax.text(-9, 3 * int(ndata / 4), r"$\xi_-^{ij}(\theta)$", fontsize=12)

    plt.savefig(plot_path, dpi=plot_dpi)
    plt.close()
    print("Plot saved as %s" % (plot_path))
