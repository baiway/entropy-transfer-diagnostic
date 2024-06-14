from gs2ools.utils.ncdf2dict import ncdf2dict
from gs2ools.output.read_restart import read_restart

import imageio.v3 as iio
import xarray as xr
from pathlib import Path

def make_vpar_vperp_plot(dfn_data, nc_data, layout, ik = 0, it = 0, isp = 0, ig = 0,
                         transform = abs, plot_points = False, edgecolor = None,
                         shade = True, return_data = False, outpath=""):

    """Produces a pcolormesh plot showing the distribution function
    as a function of vperp and vpar for the selected theta grid point,
    species and wavenumber.

    Note: this is copied from Ref. [1] (with some small changes).
    
    [1] https://bitbucket.org/gyrokinetics/analysis/src/master/python/gs2ools/analysis/dfn_plots.py
    """

    from matplotlib.pyplot import pcolormesh, colorbar, plot, xlabel, ylabel, title, savefig, close
    from numpy import zeros, sqrt, ones, NaN
    # Sanity check the inputs
    required_entries = ['bmag', 'kx', 'ky', 'theta', 'lambda', 'energy']
    if not all([x in nc_data.keys() for x in required_entries]):
        raise Exception('Not all required keys in nc_data')

    # Work out mapping from dimension to index
    # The dimensions are in reverse order as compared to the layout string
    # so we reverse the layout string here
    layout_list = list(layout)[::-1]
    isp_ind, ie_ind, il_ind, ik_ind, it_ind = [
        layout_list.index(x) for x in ['s','e','l','y','x']]

    bmag = nc_data['bmag']
    has_trapped = bmag.max() > bmag.min()
    pitch = nc_data['lambda']
    energy = nc_data['energy'][isp,:]
    nenergy = energy.shape[0]
    bmag_local = bmag[ig]
    theta = nc_data['theta']
    kx = nc_data['kx']
    ky = nc_data['ky']

    npitch = pitch.shape[0]
    nvalid_pitch = 0
    for ip, pitch_val in enumerate(pitch):
        tmp = bmag_local * pitch_val
        if tmp <= 1:
            nvalid_pitch += 1
        else:
            break

    # Two signs of vpar for each pitch angle
    nvpar = 2*nvalid_pitch
    if has_trapped:
        # Except for vpar = 0
        nvpar -= 1

    from numpy import s_
    def pack_indices(isp,ie,il,ik,it,isig,ig):
        index_pack = list([0,0,0,0,0])
        index_pack[isp_ind] = isp
        index_pack[ie_ind] = ie
        index_pack[il_ind] = il
        index_pack[ik_ind] = ik
        index_pack[it_ind] = it
        return s_[index_pack[0], index_pack[1], index_pack[2],
                  index_pack[3], index_pack[4], isig, ig]

    # We add an extra energy row to provide an energy = 0
    # set of grid points. This helps avoid pcolormesh
    # issues due to not knowing the lower bounds of the
    # cells it draws.
    vpar = zeros([nvpar, nenergy + 1])
    vperp = zeros(vpar.shape)
    data = zeros(vpar.shape, dtype = 'complex')

    for ip, pitch_val in enumerate(pitch):
        tmp = bmag_local * pitch_val
        if tmp <= 1.0:
            for ie, energy_val in enumerate(energy):
                vperp[ip, ie + 1] = sqrt(energy_val*tmp)
                vpar[ip, ie + 1] = -sqrt(energy_val * (1-tmp))
                data[ip, ie + 1] = dfn_data[pack_indices(isp, ie, ip, ik, it, 1, ig)]
                if not has_trapped or (has_trapped and ip < nvalid_pitch):
                    vperp[-ip-1, ie + 1] = vperp[ip, ie + 1]
                    vpar[-ip-1, ie + 1] = -vpar[ip, ie + 1]
                    data[-ip-1, ie + 1] = dfn_data[pack_indices(isp, ie, ip, ik, it, 0, ig)]
    data[:, 0] = NaN

    # plot trapped-passing boundary
    # factors of 2 and the mass appear to have been absorbed into the definitions
    # of energy and lambda, judging from the above
    vpar_bounce = zeros(energy.shape)
    vperp_bounce = zeros(energy.shape)
    for ie, energy_val in enumerate(energy):
        tmp = bmag_local / bmag.max()
        vpar_bounce[ie] = sqrt(energy_val * (1 - tmp))
        vperp_bounce[ie] = sqrt(energy_val * tmp)
    
    def null(a):
        return a

    if transform is None: transform = null

    if shade:
        shading = 'gouraud'
    else:
        shading = 'auto'

    mesh = pcolormesh(vpar, vperp, transform(data), shading=shading, edgecolor=edgecolor)
    colorbar(mesh)

    if plot_points:
        plot(vpar, vperp, 'kx')

    plot(vpar_bounce, vperp_bounce, linestyle="dashed")
    plot(-vpar_bounce, vperp_bounce, linestyle="dashed")

    xlabel(r"$v_\parallel$" + " / " + r"$v_\text{th}$")
    ylabel(r"$v_\perp$" + " / " + r"$v_\text{th}$")
    title(rf"$k_x\rho={kx[it]:.2f}, k_y\rho={ky[ik]:.2f}, \theta={theta[ig]:.2f}$")
    savefig(outpath)
    close()


    if return_data:
        return vpar, vperp, data


if __name__ == "__main__":
    restart_dir = "/users/bc1264/scratch/cbc-distfn-transfer/restart_dir"
    output_filepath = "/users/bc1264/scratch/cbc-distfn-transfer/input.out.nc"

    dfn_data = read_restart(fold=restart_dir, include_fields=False)
    layout = "xyles"
    nc_data = ncdf2dict(filename=output_filepath, 
                        only=["supercell_labels", "kx", "ky", "theta", 
                              "theta0", "bmag", "lambda", "energy"])

    filenames = []
    frames_dir = Path(output_filepath).parent / "frames_dir"
    if not frames_dir.is_dir():
        frames_dir.mkdir()

    ds = xr.open_dataset(output_filepath, engine="netcdf4")
    it = len(ds["kx"].values) // 4

    for ig, theta in enumerate(nc_data["theta"]):
        outpath = str(frames_dir / f"ig-{ig}.png")
        make_vpar_vperp_plot(dfn_data, nc_data, layout, ik=0, it=it, ig=ig, 
                             plot_points=True, outpath=outpath)
        filenames.append(outpath)

    images = [iio.imread(f) for f in filenames]
    iio.imwrite("test.mp4", images, fps=3)
