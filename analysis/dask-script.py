import json
import numpy as np
import xarray as xr
import dask
from dask.distributed import Client
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import imageio.v3 as iio
from pathlib import Path
from memory_profiler import memory_usage


@dask.delayed
def draw_frame(chunk):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=600)

    t = chunk["t"]
    theta = chunk["theta"]
    Ts = chunk["entropy_transfer_3D"]
    T = t[-1] - t[0]
    Ts_by_theta = np.reciprocal(T) * Ts.integrate(coord="t").sum(
        dim=["kys", "kxs", "kxt"]
    )
    Ts_by_theta /= abs(Ts_by_theta).max() # normalise


    phi2 = chunk["phi_t"].sel(ri=0) ** 2 + chunk["phi_t"].sel(ri=1) ** 2
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
    path = frames_dir / f"t-{int(np.round(t[0].item()))}.png"
    fig.savefig(path)
    plt.close(fig)

    return path


if __name__ == "__main__":
    # TODO create a GUI that enables the user to enter
    # tsat, T, select whether they want a single frame or
    # a movie, if the latter, then set the fps and dt, etc.

    # TODO implement error handling for tsat, T, dt, etc.
    with open("inputs.json", "r") as file:
        data = json.load(file)

    # Required values
    GS2_output_file = data["GS2_output_file"]
    movie = data["movie"]
    tsat = data["tsat"]
    T = data["T"]
    memory_limit = data["memory_limit"]
    n_workers = data["n_workers"]

    # Get optional values 'dt' and 'fps' if making a movie
    if movie:
        dt = data.get("dt", None)
        fps = data.get("fps", None)

    # Create 'frames' directory if it does not already exist
    frames_dir = Path("./frames")
    if not frames_dir.is_dir():
        frames_dir.mkdir()

    ds = xr.open_dataset(GS2_output_file, engine="netcdf4")
    t = ds["t"]
    tlast = t[-1].item()

    # Generate (tstart, tstop) slices for each frame
    if movie:
        start_times = np.arange(tsat, tlast - T, dt)
        time_slices = [slice(t.sel(t=start, method="nearest").item(),
                        t.sel(t=start+T, method="nearest").item()) for start in start_times]
    else:
        time_slices = [slice(t.sel(t=tsat, method="nearest").item(),
                        t.sel(t=tsat + T, method="nearest").item())]


    # Load the first slice into memory to determine whether memory_limit
    # will be exceeded
    sample_slice = ds.sel(t=time_slices[0])
    t = sample_slice["t"]
    theta = sample_slice["theta"]
    Ts = sample_slice["entropy_transfer_3D"]
    phi_t = sample_slice["phi_t"]

    arr_size = (t.nbytes + theta.nbytes + Ts.nbytes + phi_t.nbytes) / (1024 ** 3)

    if (3 * arr_size) > memory_limit:
        fudge_factor = memory_limit / arr_size
        safe_T = np.floor(T * memory_limit / (3 * arr_size)).astype(int)
        raise MemoryError(
            f"Operation would exceed memory limit of {memory_limit} GiB. " \
            f"Estimated memory usage: {3 * arr_size:.2f} GiB. " \
            f"Try reducing 'T' to {safe_T} or smaller."
        )

    client = Client(
        n_workers=n_workers,
        threads_per_worker=1,
        memory_limit=f"{memory_limit}GiB"
    )

    tasks = [draw_frame(ds.sel(t=tslice)) for tslice in time_slices]
    paths = dask.compute(tasks)[0]
    images = [iio.imread(path) for path in paths]
    iio.imwrite("test.mp4", images, fps=fps)
