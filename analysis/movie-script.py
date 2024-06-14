import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import multiprocessing
import imageio.v2 as iio
from matplotlib.colors import LogNorm


def configure_matplotlib():
    plt.rcParams["text.usetex"] = True
    plt.rcParams["xtick.labelsize"] = 8
    plt.rcParams["ytick.labelsize"] = 8
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Computer Modern"


# TODO change this to two subplots
def plot_poloidal_structure(ax, slice, tstart, tend, filename=None):
    """Plots the poloidal structure of the entropy transfer function 
    and the zonal and turbulent components of phi2 (the square of the 
    electrostatic potential). Here quantities are integrated over x, y and t."""
    theta = slice["theta"]
    Ts = slice["entropy_transfer_3D"]
    T = tend - tstart
    Ts_by_theta = np.reciprocal(T) * Ts.integrate(coord="t").sum(dim=["kys","kxs","kxt"])

    # Use this scalar to normalise the turbulent and zonal phi2_by_theta
    phi2_total = np.reciprocal(T) * slice["phi2"].integrate(coord=["t"])

    # Integrate over t and sum over kx and ky to obtain turbulent and zonal
    # phi2 as a function of theta
    phi2_t = slice["phi_t"].isel(ri=0)**2 + slice["phi_t"].isel(ri=1)**2
    phi2_by_mode_by_theta = np.reciprocal(T) * phi2_t.integrate(coord=["t"])
    phi2_by_theta = phi2_by_mode_by_theta.sum(dim=["kx", "ky"])
    phi2_zonal_by_theta = phi2_by_mode_by_theta.sel(ky=0).sum(dim="kx")
    phi2_turb_by_theta = phi2_by_theta - phi2_zonal_by_theta
    
    # normalise
    phi2_zonal_by_theta /= phi2_total
    phi2_turb_by_theta /= phi2_total
    Ts_by_theta /= abs(Ts_by_theta).max()
    
    ax.plot(theta, Ts_by_theta, label=r"$\langle T_s \rangle_{t,x,y}$")
    ax.plot(theta, phi2_turb_by_theta, label="turbulent " + 
              r"$\langle|\hat{\phi}_0|^2\rangle_{t,x,y}$")
    ax.plot(theta, phi2_zonal_by_theta, label="zonal " + 
              r"$\langle|\hat{\phi}_0|^2\rangle_{t,x,y}$")
    ax.set_xlabel(r"$\theta$")
    ax.set_xticks(np.linspace(-np.pi, np.pi, 5), 
            [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    
    ax.grid()
    ax.set_axisbelow(True)
    ax.legend()
        

# TODO change this to three subplots
def plot_transfer_trace(ax, slice, tstart, tend, filename=None):
    """Plots the time trace of the entropy transfer function and the zonal 
    and turbulent components of phi2 (the square of the electrostatic 
    potential). Here quantities are integrated over x, y and theta."""
    t = slice["t"]
    Ts = slice["entropy_transfer_3D"]
    Ts_over_t = np.reciprocal(2*np.pi) * Ts.integrate(coord="theta").sum(dim=["kys", "kxs", "kxt"])
    phi2_total = slice["phi2"]
    phi2_zonal = slice["phi2_by_mode"].sel(ky=0).sum(dim="kx")
    phi2_turb = phi2_total - phi2_zonal

    # normalise
    Ts_over_t /= abs(Ts_over_t).max()
    phi2_turb /= phi2_total
    phi2_zonal /= phi2_total

    # compute time-average of entropy transfer (appears net negative)
    avg_Ts_over_t= Ts_over_t.mean(dim="t")
    
    ax.plot(t, Ts_over_t, label=r"$\langle T_s \rangle_V$")
    ax.plot(t, phi2_turb, label=r"$\langle|\hat{\phi}_0|^2\rangle_{V,\mathrm{turb}}$")
    ax.plot(t, phi2_zonal, label=r"$\langle|\hat{\phi}_0|^2\rangle_{V,\mathrm{zonal}}$")
    ax.axhline(avg_Ts_over_t, linestyle="dashed", alpha=0.5)
    ax.set_xlabel(r"$t$" + " / " + r"$(L_\mathrm{ref}/v_\mathrm{th})$")
    ax.legend()


def get_boundaries(centres):
    """Given a 1D array of N cell centres, returns an array of 
    (N+1) cell boundaries"""
    # First check if input array is uniformly spaced
    diffs = np.diff(centres)
    if not np.allclose(diffs, diffs[0]):
        raise ValueError(f"Input array is not uniformly spaced. Differences: {diffs}")
    edges = np.zeros(centres.size + 1)
    edges[1:-1] = (centres[:-1] + centres[1:]) / 2.0
    edges[0] = centres[0] - (centres[1] - centres[0]) / 2.0
    edges[-1] = centres[-1] + (centres[-1] - centres[-2]) / 2.0

    return edges


def plot_turbulence_spectrum(ax, ds, tstart, tend):
    """Plots the (kx, ky) spectrum of phi2 (the square of the electrostatic
    potential) integrated over t and theta."""
    T = tend - tstart
    # 'roll' and 'sort' to change structure of kx 
    # from [0.0, 0.08, ..., 1.89, -1.89, ..., -0.08] 
    # to [-1.89, ..., -0.08, 0.0, 0.08, ..., 1.89]
    ds = ds.roll(kx=len(ds["kx"] // 2), roll_coords=True)
    ds = ds.sortby("kx")

    integrated_phi2_by_mode = (1 / T) * ds["phi2_by_mode"].integrate(coord="t")
    kx_edges = get_boundaries(ds["kx"].values)
    ky_edges = get_boundaries(ds["ky"].values)

    c = ax.pcolormesh(kx_edges, ky_edges, integrated_phi2_by_mode, 
                   cmap="Reds", shading="auto", norm=LogNorm())

    fig = ax.figure
    fig.colorbar(c, ax=ax, label=r"$\langle|\hat{\phi}_0|^2\rangle_{t,\theta}$" + 
                " / " + r"$\phi_\mathrm{ref}^2$")

    ax.set_xlabel(r"$k_x\rho_\mathrm{ref}$")
    ax.set_ylabel(r"$k_y\rho_\mathrm{ref}$")


def plot_phi2_trace(ax, ds, tstart, tend):
    """Plots phi2 (the square of the electrostatic potential) over time
    and highlights the specified time range."""
    # Create a nice plot to display the selected time slice
    t = ds["t"].values
    phi2_total = ds["phi2"]
    phi2_zonal = ds["phi2_by_mode"].sel(ky=0).sum(dim="kx")
    phi2_turb = phi2_total - phi2_zonal

    ax.plot(t, phi2_total, label="total")
    ax.plot(t, phi2_turb, label="turbulent")
    ax.plot(t, phi2_zonal, label="zonal")
    ax.set_xlabel(r"$t$" + " / " + r"$(L_\mathrm{ref}/v_\mathrm{th})$")
    ax.set_ylabel(r"$\langle|\hat{\phi}_0|^2\rangle_V$" + " / " + 
                      r"$\phi_\mathrm{ref}^2$")
    ax.set_yscale("log")

    # Highlight saturated region on the plot
    ax.axvspan(tstart, tend, color="green", alpha=0.2) 
    ax.legend()


def draw_frame(filename, times, frames_dir):
    """Draw plots for a single time slice. Return the path to the frame
    so it can be stitched together with the others to form a movie."""
    ds = xr.open_dataset(filename, engine="netcdf4")
    start, stop = times
    t = ds["t"]
    drop_indices = t.where((t < start) | (t > stop), drop=True).values
    slice = ds.drop_sel(t=drop_indices)
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=600)

    plot_phi2_trace(axes[0, 0], ds, start, stop)
    plot_turbulence_spectrum(axes[0, 1], slice, start, stop)
    plot_transfer_trace(axes[1, 0], slice, start, stop)
    plot_poloidal_structure(axes[1, 1], slice, start, stop)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.10)
    fig.text(0.05, 0.02, "Quantities in the bottom row are normalised to " + r"$\pm 1$"
             + " by dividing by their maximum values.", ha="left", va="center")

    # Save frame and return filepath
    outfile = frames_dir / f"t-{int(np.round(start))}.png"
    plt.savefig(outfile)
    return outfile

def main(filename, frames_dir, from_file, movie_name, nproc, fps, keep_frames):
    """Produces a movie showing the (kx, ky) spectrum of phi2, the time trace
    of phi2 and T_s, and the poloidal structure of phi2 and T_s at different
    time slices."""
    # Validate inputs and set some sensible defaults
    if filename is not None and not filename.is_file():
        raise FileNotFoundError(f"File does not exist: {filename}")

    if frames_dir is None:
        frames_dir = Path.cwd() / "frames"
    
    if frames_dir is not None and not frames_dir.is_dir():
        if from_file:
            raise FileNotFoundError(f"Directory does not exist: {frames_dir}")
        else:
            print(f"'frames_dir' does not yet exist. Creating: {frames_dir}")
            frames_dir.mkdir()
            
    if movie_name is None:
        if filename is None:
            movie_name = Path.cwd() / "plots.mp4"
        else:
            basename = str(filename.name).split(".out.nc")[0]
            movie_name = Path.cwd() / f"{basename}-plots.mp4"
    
    if not from_file:
        # Load dataset
        ds = xr.open_dataset(filename, engine="netcdf4")

        # Plot phi2 over time to determine rough saturation region
        plt.plot(ds["t"], ds["phi2"])
        plt.yscale("log")
        plt.show(block=False)
        tsat = float(input("Saturation time: "))
        T = float(input("Window width: "))
        dt = float(input("Step size: "))
        plt.close()

        # Generate (tstart, tend) pairs
        t = ds["t"]
        tlast = t.max().item()
        start_times = np.arange(tsat, tlast - T, dt)
        times = [(t.sel(t=start, method="nearest").item(),
                  t.sel(t=start+T, method="nearest").item()) for start in start_times]
        times.append((times[-1][1], t[-1].item())) # add final time slice
    
        # Set font, fontsize, etc.
        configure_matplotlib()

        # Generate frames (parallelised)    
        with multiprocessing.Pool(processes=nproc) as pool:
            pool.starmap(draw_frame, [(filename, t, frames_dir) for t in times])

    with iio.get_writer(movie_name, fps=fps) as writer:
        for file in sorted(frames_dir.glob("t-*.png")):
            frame = iio.imread(file)
            writer.append_data(frame)
            if not keep_frames:
                file.unlink()
        if not keep_frames:
            frames_dir.unlink()


if __name__ == "__main__":
    description = (
        "Generates a movie (.mp4) of several plots to show how they change over "
        "different time integration periods. The plots are: "
        "a (kx, ky) spectrum of phi2, " # TODO add option for real space
        "a time trace of phi2 and transfer, "
        "and a poloidal structure of phi2 and transfer."
    )
    
    parser = argparse.ArgumentParser(prog="analysis-script.py", 
                                     description=description)
    parser.add_argument("--input-file", type=str, default=None,
        help="GS2 output file (.out.nc).")
    parser.add_argument("--frames-dir", type=str, default=None, 
        help="Directory to save individual frames. If not set, creates " \
            "a directory called 'frames'.")
    parser.add_argument("--from-file", default=False, action="store_true",
        help="If set, creates a movie using the frames in the specified " \
            "frames directory. If '--frames-dir' is not set, looks in 'frames'" \
            "instead. Useful for making cosmetic adjustments to the plot code.")
    parser.add_argument("--movie-name", type=str, default=None,
        help="The name of the output .mp4 file. If not set, defaults to the " \
            "name of the GS2 output file if provided, or 'plots.mp4' if not")
    parser.add_argument("-n", "--nproc", type=int, default=4,
        help="Number of worker processes (default 4).")
    parser.add_argument("--fps", type=int, default=2,
        help="Number of frames per second for the output video (default 2).")
    parser.add_argument("-k", "--keep-frames", default=False, action="store_true",
        help="Keep the frames directory after creating the video (default false).")

    args = parser.parse_args()
    filename = args.input_file
    frames_dir = args.frames_dir
    from_file = args.from_file
    movie_name = args.movie_name
    nproc = args.nproc
    fps = args.fps
    keep_frames = args.keep_frames
    

    if filename is not None:
        filename = Path(filename)

    if frames_dir is not None:
        frames_dir = Path(frames_dir)

    main(filename, frames_dir, from_file, movie_name, nproc, fps, keep_frames)
