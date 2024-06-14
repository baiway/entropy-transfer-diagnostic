import argparse
from pathlib import Path
import json
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


def plot_Ts(slice, plot_path):
    """Plots the entropy transfer, zonal phi2 and turbulent phi2 as a function of
    the poloidal angle theta."""
    fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=600)

    t = slice["t"]
    theta = slice["theta"]
    Ts = slice["entropy_transfer_3D"]
    T = t[-1] - t[0]
    Ts_by_theta = np.reciprocal(T) * Ts.integrate(coord="t").sum(dim=["kys", "kxs", "kxt"])
    Ts_by_theta /= abs(Ts_by_theta).max() # normalise


    phi2 = slice["phi_t"].sel(ri=0) ** 2 + slice["phi_t"].sel(ri=1) ** 2
    phi2_by_mode_by_theta = np.reciprocal(T) * phi2.integrate(coord=["t"])

    phi2_by_theta = phi2_by_mode_by_theta.sum(dim=["kx", "ky"])
    phi2_zonal_by_theta = phi2_by_mode_by_theta.sel(ky=0).sum(dim="kx")
    phi2_turb_by_theta = phi2_by_theta - phi2_zonal_by_theta

    # normalise
    phi2_zonal_by_theta /= phi2_by_theta
    phi2_turb_by_theta /= phi2_by_theta
    
    ax.plot(theta, Ts_by_theta, label=r"$\langle T_s \rangle_{t,x,y}$")
    ax.plot(theta, phi2_turb_by_theta, label="turbulent " + 
            r"$\langle|\hat{\phi}_0|^2\rangle_{t,x,y}$")
    ax.plot(theta, phi2_zonal_by_theta, label="zonal " + 
            r"$\langle|\hat{\phi}_0|^2\rangle_{t,x,y}$")
    ax.set_xlabel(r"$\theta$")
    ax.set_xticks(np.linspace(-np.pi, np.pi, 5),
        [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.legend()

    # Save frame and return filepath
    fig.savefig(plot_path)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot entropy transfer function.")
    parser.add_argument("input_file", type=str, help="Full path to JSON input file.")
    args = parser.parse_args()
    input_file = args.input_file
    
    with open(input_file, "r") as file:
        data = json.load(file)

    GS2_output_file = data["GS2_output_file"]
    tsat = data["tsat"]
    T = data["T"]

    ds = xr.open_dataset(GS2_output_file, engine="netcdf4")
    slice = ds.sel(t=slice(tsat, tsat+T))

    plot_path = Path(GS2_output_file).parent / "Ts_poloidal_structure.png"

    plot_Ts(slice, plot_path)
