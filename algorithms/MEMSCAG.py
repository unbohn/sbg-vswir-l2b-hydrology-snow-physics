from itertools import product
import logging
import multiprocessing
import os
import time
from types import SimpleNamespace

import click
import earthaccess
import numpy as np
from spectral.io import envi
import xarray as xr

from isofit import ray
from isofit.core.common import envi_header
from isofit.core.fileio import write_bil_chunk


def generate_filter(wvl_m, wvl, wl_resol):
    num_wvl_m = len(wvl_m)
    num_wvl = len(wvl)

    s_norm_m = np.zeros((num_wvl_m, num_wvl))
    exp_max = 2.
    exp_min = 2.
    exp_arr = exp_max + (exp_min - exp_max) * np.arange(0, 2100, 1) / (num_wvl - 1)
    c_arr = (1 / (2 ** exp_arr * np.log(2))) ** (1 / exp_arr)
    for bd in range(num_wvl):
        li1 = np.logical_and(wvl_m >= (wvl[bd] - 2. * wl_resol[bd]), wvl_m <= (wvl[bd] + 2. * wl_resol[bd]))
        li1 = np.where(li1)
        cnt = len(li1)
        if cnt > 0:
            tmp = np.abs(wvl[bd] - wvl_m[li1]) / (wl_resol[bd] * c_arr[bd])
            s = np.exp(-(tmp ** exp_arr[bd]))
            s_norm_m[li1, bd] = s / np.sum(s)

    return s_norm_m


def check_consecutive_residuals(arr, num_consecutive, threshold):
    n = len(arr)
    for i in range(n - num_consecutive + 1):
        window = arr[i:i + num_consecutive]
        if all(x <= threshold for x in window):
            return True
    return False


def run_memscag(input_data, wl, ndsi, models, Qs, Rs, inds):
    memscag_output = np.zeros(3)
    if ndsi >= 0.0:
        g_fit = None
        x_hat = None
        gs_ind = None
        
        for ii in range(len(models)):
            # compute regression coefficients
            b = Qs[ii].T.dot(input_data)
            # compute endmember fractional cover
            p = np.linalg.inv(Rs[ii]).dot(Qs[ii].T).dot(input_data)
    
            residual = input_data - Qs[ii].dot(b)
            RMSE = np.sqrt(1 / len(wl) * np.sum(residual ** 2))
    
            if ii == 0 or RMSE < g_fit:
                g_fit = RMSE
                x_hat = p
                gs_ind = inds[ii][0]
    
            ccr = check_consecutive_residuals(arr=residual, num_consecutive=7, threshold=0.025)
            
            if RMSE < 0.025 and 0.0 <= np.sum(x_hat) <= 1.0 and -0.01 < np.all(x_hat) < 1.01 and ccr:
                break
    
        snow_fc = x_hat[0] / x_hat.sum()
        grain_size = gs_ind * 10
        uncertainty = RMSE
        
    else:
        snow_fc = np.nan
        grain_size = np.nan
        uncertainty = np.nan

    memscag_output[:] = snow_fc, grain_size, uncertainty

    return memscag_output


def main(args: SimpleNamespace) -> None:
    logging.basicConfig(
        format="%(levelname)s:%(asctime)s ||| %(message)s",
        level=args.loglevel,
        filename=args.logfile,
        datefmt="%Y-%m-%d,%H:%M:%S",
    )

    if os.path.isfile(args.output_memscag_file):
        dat = (
            envi.open(envi_header(args.output_memscag_file))
            .open_memmap(interleave="bip")
            .copy()
        )
        if not np.all(dat == -9999):
            logging.info("Existing MEMSCAG file found, terminating")
            exit()
    
    fid = args.fid

    url = f"https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/EMITL2ARFL.001/EMIT_L2A_RFL_001_{fid}/" \
          f"EMIT_L2A_RFL_001_{fid}.nc "
    
    # Get requests https Session using Earthdata Login Info
    fs = earthaccess.get_requests_https_session()
    # Retrieve granule asset ID from URL (to maintain existing naming convention)
    granule_asset_id = url.split('/')[-1]
    # Define Local Filepath
    fp = f'data/{granule_asset_id}'
    # Download the Granule Asset if it doesn't exist
    if not os.path.isfile(fp):
        logging.info("Download the Granule Asset")
        with fs.get(url,stream=True) as src:
            with open(fp,'wb') as dst:
                for chunk in src.iter_content(chunk_size=64*1024*1024):
                    dst.write(chunk)

    ds = xr.open_dataset(fp)
    rfl = np.array(ds['reflectance'])
    emit_wl = np.array(xr.open_dataset(fp, group='sensor_band_parameters')["wavelengths"])
    emit_fwhm = np.array(xr.open_dataset(fp, group='sensor_band_parameters')["fwhm"])
    input_data = rfl.copy()
    rfls = input_data.shape

    output_metadata = {}
    output_metadata["lines"] = rfls[0]
    output_metadata["samples"] = rfls[1]
    output_metadata["interleave"] = "bil"
    output_metadata["data type"] = "4"
    output_metadata["bands"] = "3"
    output_metadata["description"] = (
        "MEMSCAG Snow fractional cover, snow grain size, and uncertainty (spectral fit RMSE)"
    )

    img = envi.create_image(
        envi_header(args.output_memscag_file), ext="", metadata=output_metadata, force=True
    )
    del img

    logging.info("Load Endmember spectral library and resample to instrument wavelengths")
    spec_lib = envi.open(args.spec_lib)

    emit_norm = generate_filter(wvl_m=np.array(spec_lib.bands.centers), wvl=emit_wl, wl_resol=emit_fwhm)
    
    emit_spec_lib = spec_lib.spectra @ emit_norm
    emit_spec_lib = emit_spec_lib / 10000
    
    snow_spec_lib = emit_spec_lib[1:111, :]
    veg_spec_lib = emit_spec_lib[142:162, :]
    rock_spec_lib = emit_spec_lib[113:142, :]

    logging.info("Store endmember models")
    models = []
    
    for ii in range(snow_spec_lib.shape[0]):
        for jj in range(-1, veg_spec_lib.shape[0]):
            for kk in range(-1, rock_spec_lib.shape[0]):
                snow_spec = snow_spec_lib[ii, :]
                if jj == -1:
                    veg_spec = None
                else:
                    veg_spec = veg_spec_lib[jj, :]
                if kk == -1:
                    rock_spec = None
                else:
                    rock_spec = rock_spec_lib[kk, :]
    
                spec_lib = []
                for lib in [snow_spec, veg_spec, rock_spec]:
                    if type(lib) == np.ndarray:
                        spec_lib.append(lib)

                # take the transpose to calculate the orthogonal instead of the orthonormal matrix
                # => modified Gram-Schmidt
                models.append(np.array(spec_lib).T)

    logging.info("Compute qr factorization (orthogonal and upper-triangular matrices)")
    Qs = []
    Rs = []
    
    for ii in range(len(models)):
        Q, R = np.linalg.qr(models[ii])
        Qs.append(Q)
        Rs.append(R)
    
    snow_ind = np.arange(0, snow_spec_lib.shape[0])
    veg_ind = np.arange(-1, veg_spec_lib.shape[0])
    rock_ind = np.arange(-1, rock_spec_lib.shape[0])
    
    inds = list(product(snow_ind, veg_ind, rock_ind))

    logging.info("Compute normalixed-difference snow index (NDSI) to only include snow-covered pixels")
    ndsi = (input_data[:, :, 25] - input_data[:, :, 170]) / (input_data[:, :, 25] + input_data[:, :, 170])

    logging.info("Set up ray-based parallel processing")
    if args.n_cores == -1:
        n_cores = multiprocessing.cpu_count()
        n_workers = multiprocessing.cpu_count()
    else:
        n_cores = args.n_cores
        n_workers = args.n_cores

    rayargs = {
        "ignore_reinit_error": True,
        "local_mode": args.n_cores == 1,
        "_temp_dir": args.ray_tmp_dir,
        "num_cpus": n_cores,
    }

    ray.init(**rayargs)

    line_breaks = np.linspace(0, rfls[0], n_workers, dtype=int)
    line_breaks = [
        (line_breaks[n], line_breaks[n + 1]) for n in range(len(line_breaks) - 1)
    ]

    params = [
            ray.put(obj)
            for obj in [input_data, emit_wl, ndsi, models, Qs, Rs, inds, args.output_memscag_file, line_breaks,
                        args.loglevel, args.logfile]
        ]
    
    start_time = time.time()
    logging.info("Beginning parallel MEMSCAG runs")
    result_list = [
        run_lines.remote(
            *params, n
        )
        for n in range(len(line_breaks))
    ]
    [ray.get(result) for result in result_list]

    total_time = time.time() - start_time
    logging.info(
        f"MEMSCAG runs complete.  {round(total_time,2)}s total, "
        f"{round(rfls[0]*rfls[1]/total_time,4)} spectra/s, "
        f"{round(rfls[0]*rfls[1]/total_time/n_workers,4)} spectra/s/core"
    )


@ray.remote(num_cpus=1)
def run_lines(rfl, wl, ndsi, models, Qs, Rs, inds, output_memscag_file, startstop, loglevel, logfile, line_id):
    logging.basicConfig(
        format="%(levelname)s:%(asctime)s ||| %(message)s",
        level=loglevel,
        filename=logfile,
        datefmt="%Y-%m-%d,%H:%M:%S",
    )

    start_line, stop_line = startstop[line_id]
    output_memscag = np.zeros((stop_line - start_line, rfl.shape[1], 3)) - 9999

    for r in range(start_line, stop_line):
        for c in range(rfl.shape[1]):
            meas = rfl[r, c, :]
            if np.all(meas < 0):
                continue
            output_memscag[r - start_line, c, :] = run_memscag(meas, wl, ndsi[r, c], models, Qs, Rs, inds)

        logging.info(f"MEMSCAG writing line {r}")

        write_bil_chunk(
            output_memscag[r - start_line, ...].T,
            output_memscag_file,
            r,
            (rfl.shape[0], rfl.shape[1], output_memscag.shape[2]),
        )


@click.command(name="memscag")
@click.argument("fid")
@click.argument("spec_lib")
@click.argument("output_memscag_file", required=False)
@click.option("--loglevel", default="INFO")
@click.option("--logfile")
@click.option("--n_cores", type=int, default=-1)
@click.option("--ray_tmp_dir")
@click.option(
    "--debug-args",
    help="Prints the arguments list without executing the command",
    is_flag=True,
)
def cli_memscag(debug_args, **kwargs):
    """Run MEMSCAG

    Calculate snow fractional cover and snow grain size for
    a set of reflectance data by solving a system of linear equations
    using modified Gram-Schmidt orthogonalization.

    References
    MEMSCAG: Painter et al. (2003) https://www.sciencedirect.com/science/article/abs/pii/S0034425702001876)
    Gram-Schmidt: "https://t-redactyl.io/blog/2020/07/simplifying-the-normal-equation-with-gram-schmidt.html"
    """
    click.echo("Running MEMSCAG from Reflectance")
    if debug_args:
        click.echo("Arguments to be passed:")
        for key, value in kwargs.items():
            click.echo(f"  {key} = {value!r}")
    else:
        # SimpleNamespace converts a dict into dot-notational for backwards compatability with argparse
        main(SimpleNamespace(**kwargs))

    click.echo("Done")


if __name__ == "__main__":
    cli_memscag()
