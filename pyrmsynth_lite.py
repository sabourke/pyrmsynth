#!/usr/bin/env python

"""
pyrmsynth_lite, Stephen Bourke

Based on pyrmsynth (Michael Bell, Henrik Junklewitz and Sarrvesh Sridhar)
"""

from __future__ import print_function
import argparse
from astropy.io import fits
import numpy as np
import csv
import numpy
import rm_tools as R
import math
import sys

VERSION = "0.1"

import astropy.version
if astropy.version.major < 1 or astropy.version.major < 2 and astropy.version.minor < 3:
    overwrite = {'clobber': True}
    no_overwrite = {'clobber': False}
else:
    overwrite = {'overwrite': True}
    no_overwrite = {'overwrite': False}

class Params:
    def __init__(self):
        """
        Params()

        Intitializes the parameter list with default values. Change parameters
        by directly accessing these class variables.
        """
        self.nphi = 200
        self.dphi = 1
        self.phi_min = -100
        temp = numpy.arange(self.nphi)
        self.phi = self.phi_min + temp * self.dphi
        self.niter = 500
        self.gain = 0.1
        self.cutoff = 0.
        self.do_clean = False
        self.isl2 = False
        self.weight = None
        self.outputfn = 'test_100A'
        self.input_dir = './'
        self.ra_lim = []
        self.dec_lim = []
        self.imagemask = ''
        self.threshold = 0.


def parse_input_file(infile):
    """ parse the parameter file.  Does not handle weights yet."""

    params = Params()

    reader = csv.reader(open(infile, 'rb'), delimiter=" ",
                        skipinitialspace=True)
    parset = dict()

    for row in reader:
        if len(row) != 0 and row[0] != '%':
            parset[row[0]] = row[1]

    params.cutoff = float(parset['cutoff'])
    params.dec_lim = [int(parset['dec_min']), int(parset['dec_max'])]
    params.ra_lim = [int(parset['ra_min']), int(parset['ra_max'])]

    params.phi_min = float(parset['phi_min'])
    params.nphi = int(parset['nphi'])
    params.dphi = float(parset['dphi'])
    temp = numpy.arange(params.nphi)
    params.phi = params.phi_min + temp * params.dphi
    
    if 'imagemask' in parset:
        params.imagemask = parset['imagemask']
    else:
        params.imagemask = ''

    if parset['do_clean'].lower() == 'false':
        params.do_clean = False
    else:
        params.do_clean = True

    params.gain = float(parset['gain'])
    params.niter = int(parset['niter'])
    params.outputfn = parset['outputfn']
    params.input_dir = parset['input_dir']

    params.alphafromfile = False
    if 'alpha' in parset:
        try:
            params.alpha = numpy.float(parset['alpha'])
        except ValueError:
            params.alphafromfile = True
            params.alpha = parset['alpha']
    else:
        params.alpha = False

    if 'ref_freq' in parset:
        params.ref_freq = numpy.float(parset['ref_freq'])
    else:
        params.ref_freq = None
    
    if 'do_weight' in parset:
        params.weight = numpy.loadtxt(params.input_dir + parset['do_weight'])
        print('Non-trivial weights enabled! Loaded from ' + parset['do_weight'])

    return params


def progress(width, percent):
    marks = round(width * (percent / 100.0))
    spaces = math.floor(width - marks)
    loader = '[' + ('=' * int(marks)) + (' ' * int(spaces)) + ']'
    sys.stdout.write("%s %d%%\r" % (loader, round(percent)))
    if percent >= 100:
        sys.stdout.write("\n")
    sys.stdout.flush()


def freqs(header, f_axis):
    freq_ref = float(header["CRVAL{}".format(f_axis)])
    freq_ref_i = float(header["CRPIX{}".format(f_axis)]) - 1 # -1 is to do 0 indexed
    freq_delta = float(header["CDELT{}".format(f_axis)])
    freq_num = float(header["NAXIS{}".format(f_axis)])
    return np.array([freq_ref + i * freq_delta for i in np.arange(-freq_ref_i, freq_num+freq_ref_i)])


def channel_flags(cube, clip_lev=5):
    # Algorithm from getnoise_simple_qu
    # I got this from S.O'Sullivan, but I think it's originally by G.Heald.
    n_freq, n_stokes, n_dec, n_ra = cube.shape
    noise_values = np.empty(shape=(n_freq, n_stokes))
    
    for c, channel in enumerate(cube):
        for s, stokes in enumerate(channel):
            inner_third = stokes[n_dec/3:2*n_dec/3,n_ra/3:2*n_ra/3].flatten()
            if np.isnan(inner_third).all():
                noise_values[c,s] = -1
            else:
                rms = inner_third.std()
                mean = inner_third.mean()
                good_positions = np.logical_and(inner_third < mean + 3 * rms, inner_third > mean - 3 * rms)
                good_values = inner_third[good_positions]
                hist = np.histogram(good_values, bins=100)
                if max(hist[0]) == 0:
                    noise_values[c,s] = -1
                
                bin_label = hist[1][:-1] + 0.5*(hist[1][1] - hist[1][0])
                bin_value = hist[0]/float(max(hist[0]))
                noise_values[c,s] = np.sqrt(abs(sum(bin_label**2 * bin_value) / bin_value.sum()))
    
    bad_freqs = np.zeros(shape=(n_freq,), dtype=bool)
    for st in range(n_stokes):
        st_noise = noise_values[:,st]
        avg_noise = np.median(st_noise[abs(st_noise) < 1])
        std_noise = (st_noise[abs(st_noise) < 1]).std()
        st_bad = np.logical_or(st_noise > (avg_noise + clip_lev * std_noise), st_noise == -1)
        bad_freqs = np.logical_or(bad_freqs, st_bad)

    return bad_freqs


def write_cube(empty_cube, data_selection, data_out, header, out_name):
    empty_cube[data_selection] = abs(data_out)
    empty_cube = np.moveaxis(empty_cube, -1, 0)
    hdu = fits.PrimaryHDU(empty_cube)
    hdu.header = header
    hdu.writeto(out_name, **overwrite)
    empty_cube = np.moveaxis(empty_cube, 0, -1)


def phi_range_to_channel_range(phi_range, params):
    ch0_indx = np.abs(params.phi - phi_range[0]).argmin()
    ch1_indx = np.abs(params.phi - phi_range[1]).argmin()
    return (min(ch0_indx, ch1_indx), max(ch0_indx, ch1_indx)+1)


def check_cube_format(header):
    try:
        assert header["CTYPE1"].startswith("RA")
        assert header["CTYPE2"].startswith("DEC")
        assert header["CTYPE3"].startswith("STOKES")
        assert header["CTYPE4"].startswith("FREQ")
    except AssertionError:
        raise ValueError("Input cube must be in order: RA,DEC,STOKES,FREQ")


def main():
    """ Handle all parsing here if started from the command line"""

    parser = argparse.ArgumentParser(description="Rotation Measure Synthesis tool.")
    parser.add_argument("-o", "--outname", type=str, help="output prefix")
    parser.add_argument("-p", "--save-pol-cube", action="store_true",
                        dest="save_pol_cube", help="Save Pol cube", default=False)
    parser.add_argument("-q", "--save-qu-cubes", action="store_true",
                        dest="save_qu_cubes", default=False,
                        help="Save derotated Q and U cubes")
    parser.add_argument("-r", "--save-residual-cubes", action="store_true",
                        dest="save_residual_cubes", default=False,
                        help="Save residual cubes if cleaning")
    parser.add_argument("-d", "--save-dirty-cubes", action="store_true",
                        dest="save_dirty_cubes", default=False,
                        help="Save dirty cubes if cleaning")
    parser.add_argument("-a", "--auto-flag", action="store_true",
                        dest="auto_flag", help="auto flag data", default=False)
    parser.add_argument("-m", "--auto-mask", dest="auto_mask", type=str, nargs=2, metavar=('MASK_FILE', 'CUTOFF'),
                        help="Use stokes I map and cutoff value in Jy. This is in addition to a regular mask. I.E. it is OR'd with the regular mask if provided. Eg. map_i.fits 0.001")
    parser.add_argument("-x", "--exclude-phi", dest="exclude_phi",
                        metavar=('PHI_LOW', 'PHI_HIGH'), nargs=2, type=float, default=(0,0), 
                        help="Exclude this Phi range from 2D maps. Eg: -3 1.5")
    parser.add_argument("-n", "--phi-rms", dest="phi_rms", type=float,
                        metavar=('PHI_LOW', 'PHI_HIGH'), nargs=2, default=(0,0), 
                        help="Make a Phi RMS image from this range. Eg: 100 115")
    parser.add_argument("--single", action="store_true", default=False,
                        help="Use single precision floats internally. Default: False")
    parser.add_argument("--double", action="store_true", default=False,
                        help="Output double precision floats. Default: False. If neither --single or --double is specified, use double precision internally and write out single precision.")
    parser.add_argument("--rmclean-mask", type=str, help="Input mask for RMCLEAN")
    parser.add_argument("--rmclean-sigma", dest="rmclean_sigma",
                        type=float, default=0.0,
                        help="Clean to RMCLEAN_SIGMA times the phi noise. This is used in combination with the cutoff value. Eg. 10")
    parser.add_argument("param_file", type=str, help="Input Parameter file")
    parser.add_argument("fits_cube", type=str, help="QU FITS cube")

    args = parser.parse_args()

    print("Parsing parameter file...")
    params = parse_input_file(args.param_file)

    if args.outname:
        params.outputfn = args.outname

    if args.single:
        in_dtype = np.complex64
    else:
        in_dtype = np.complex128

    if args.double:
        out_dtype = np.complex128
    else:
        out_dtype = np.complex64

    print("Loading data...")
    hdu = fits.open(args.fits_cube)[0]
    data = hdu.data
    header = hdu.header
    check_cube_format(header)

    if params.imagemask:
        mask = fits.open(params.imagemask)[0].data.squeeze().astype(bool)
    else:
        mask = np.ones(shape=data.shape[2:], dtype=bool)
    if args.auto_mask:
        stokesI_map = fits.open(args.auto_mask[0])[0].data.squeeze()
        cutoff = float(args.auto_mask[1])
        stokesI_mask = np.zeros_like(stokesI_map, dtype=np.uint8)
        stokesI_mask[np.where(stokesI_map > cutoff)] = 1
        if params.imagemask:
            mask = np.logical_or(mask, stokesI_mask)
        else:
            mask = stokesI_mask
        del stokesI_map, stokesI_mask

    if args.rmclean_mask:
        rmclean_mask = fits.open(args.rmclean_mask)[0].data.squeeze().astype(bool)
        params.do_clean = True
    else:
        rmclean_mask = np.ones(shape=data.shape[2:], dtype=bool)

    if args.auto_flag:
        bad_chans = channel_flags(data)
        flag_weights = np.logical_not(bad_chans).astype(float)
        if params.weight is None:
            params.weight = flag_weights
        else:
            params.weight = np.minimum(params.weight, flag_weights)
    
    if args.exclude_phi[0] != args.exclude_phi[1]:
        exclude_phi_range = phi_range_to_channel_range(args.exclude_phi, params)
        print("Excluding phi range from 2D maps: {} => {}".format(args.exclude_phi, exclude_phi_range))
    else:
        exclude_phi_range = None

    if args.phi_rms[0] != args.phi_rms[1]:
        phi_rms_range = phi_range_to_channel_range(args.phi_rms, params)
        print("RMS phi range: {} => {}".format(args.phi_rms, phi_rms_range))
    else:
        phi_rms_range = None

    if args.rmclean_sigma != 0:
        if phi_rms_range is None:
            args.error("--phi-rms required for --rmclean-sigma")

    # Survey cubes will be ordered: [Freqs, Stokes, Dec, RA]
    data = np.moveaxis(data, 0, -1) # Move Freq to end
    data = np.moveaxis(data, 0, -1) # Move Stokes to end
    # Now order is [Dec, RA, Freqs, Stokes]
    (DEC,RA,FREQ,STOKES) = range(4)

    data_selection = np.where(mask)
    non_masked_data = data[data_selection]
    # Unflagged NaNs will cause RMSynth to return NaNs
    non_masked_data = np.nan_to_num(non_masked_data.view(dtype=">c8").squeeze()).astype(in_dtype)
    clean_regions = rmclean_mask[data_selection].squeeze()
    
    params.nu = freqs(header, 4)
    params.dnu = params.nu[1] - params.nu[0]
    rms = R.RMSynth(params.nu, params.dnu, params.phi, params.weight)
    if params.do_clean:
        rmc = R.RMClean(rms, params.niter, params.gain, params.cutoff)

    rmsfout = numpy.zeros((len(rms.rmsf), 3))
    rmsfout[:, 0] = rms.rmsf_phi
    rmsfout[:, 1] = rms.rmsf.real
    rmsfout[:, 2] = rms.rmsf.imag
    numpy.savetxt(params.outputfn + "_rmsf.txt", rmsfout)

    data_out = np.empty(shape=(len(non_masked_data), params.nphi), dtype=out_dtype)
    if params.do_clean:
        if args.save_dirty_cubes:
            dirty_out = np.empty(shape=(len(non_masked_data), params.nphi), dtype=out_dtype)
        else:
            dirty_out = None
        if args.save_residual_cubes:
            residual_out = np.empty(shape=(len(non_masked_data), params.nphi), dtype=out_dtype)
        else:
            residual_out = None

    print("Doing RM synthesis")
    n = len(non_masked_data)
    width = len(str(n))
    for i in range(len(non_masked_data)):
        if params.do_clean and clean_regions[i]:
            rmc.reset()
            rmc.residual_map = rms.compute_dirty_image(non_masked_data[i])
            rmc.dirty_image = rmc.residual_map
            if args.rmclean_sigma != 0:
                q_rms = rmc.residual_map[phi_rms_range[0]:phi_rms_range[1]].real.std()
                u_rms = rmc.residual_map[phi_rms_range[0]:phi_rms_range[1]].imag.std()
                cutoff = args.rmclean_sigma * 0.5 * (q_rms + u_rms)
                rmc.cutoff = max(cutoff, params.cutoff)
            rmc.perform_clean()
            rmc.restore_clean_map()
            data_out[i] = rmc.clean_map
            if args.save_residual_cubes:
                residual_out[i] = rmc.residual_map
            if args.save_dirty_cubes:
                dirty_out[i] = rmc.dirty_image
        else:
            data_out[i] = rms.compute_dirty_image(non_masked_data[i])
        if i % 100 == 0:
            progress(20, 100.0 * i / n)
    print("")

    header_axes_attributes = ["NAXIS", "CTYPE",  "CRVAL", "CRPIX", "CDELT", "CUNIT", "CROTA"]
    for attr in header_axes_attributes:
        attr += "4"
        if attr in header:
            del header[attr]
    
    header["CTYPE3"] = "Phi"
    header["CUNIT3"] = "rad/m/m"
    header["CRPIX3"] = 1
    header["CRVAL3"] = params.phi_min
    header["CDELT3"] = params.dphi


    if args.save_pol_cube or args.save_qu_cubes or params.do_clean and (args.save_residual_cubes or args.save_dirty_cubes):
        cube_out = np.full(fill_value=np.NaN, shape=(data.shape[DEC], data.shape[RA], params.nphi), dtype=data_out.real.dtype)
    else:
        cube_out = None

    if args.save_pol_cube:
        print("Saving phi cube")
        if not params.do_clean:
            write_cube(cube_out, data_selection, abs(data_out), header, params.outputfn + "_di_p.fits")
        else:
            write_cube(cube_out, data_selection, abs(data_out), header, params.outputfn + "_clean_p.fits")
            if args.save_dirty_cubes:
                print("Saving phi dirty cube")
                write_cube(cube_out, data_selection, abs(dirty_out), header, params.outputfn + "_di_p.fits")
            if args.save_residual_cubes:
                print("Saving phi residual cube")
                write_cube(cube_out, data_selection, abs(residual_out), header, params.outputfn + "_residual_p.fits")

    if args.save_qu_cubes:
        print("Saving Q & U cubes")
        if not params.do_clean:
            write_cube(cube_out, data_selection, data_out.real, header, params.outputfn + "_di_q.fits")
            write_cube(cube_out, data_selection, data_out.imag, header, params.outputfn + "_di_u.fits")
        else:
            write_cube(cube_out, data_selection, data_out.real, header, params.outputfn + "_clean_q.fits")
            write_cube(cube_out, data_selection, data_out.imag, header, params.outputfn + "_clean_u.fits")
            if args.save_dirty_cubes:
                print("Saving Q & U dirty cubes")
                write_cube(cube_out, data_selection, dirty_out.real, header, params.outputfn + "_di_q.fits")
                write_cube(cube_out, data_selection, dirty_out.imag, header, params.outputfn + "_di_u.fits")
            if args.save_residual_cubes:
                print("Saving Q & U residual cubes")
                write_cube(cube_out, data_selection, residual_out.real, header, params.outputfn + "_residual_q.fits")
                write_cube(cube_out, data_selection, residual_out.imag, header, params.outputfn + "_residual_u.fits")
    
    del cube_out
    if params.do_clean:
        del dirty_out, residual_out
    
    for attr in header_axes_attributes:
        attr += "3"
        if attr in header:
            del header[attr]

    hdu.data = np.full(fill_value=np.NaN, shape=(data.shape[DEC], data.shape[RA]), dtype=data_out.real.dtype)

    # RMS Image
    if phi_rms_range != None:
        q_rms = data_out[:,phi_rms_range[0]:phi_rms_range[1]].real.std(axis=1)
        u_rms = data_out[:,phi_rms_range[0]:phi_rms_range[1]].imag.std(axis=1)
        p_rms = 0.5 * (q_rms + u_rms)
        hdu.data[data_selection] = p_rms
        print("Saving rms map")
        hdu.writeto(params.outputfn + "_rms.fits", **overwrite)

    if exclude_phi_range != None:
        data_out[:,exclude_phi_range[0]:exclude_phi_range[1]] = 0

    # 2D maps
    hdu.data[data_selection] = np.amax(abs(data_out), axis=-1)
    print("Saving int pol map")
    hdu.writeto(params.outputfn + "_polint.fits", **overwrite)

    indx_max = np.argmax(abs(data_out), axis=-1)
    hdu.data[data_selection] = data_out.real[range(len(data_out)), indx_max]
    print("Saving q map")
    hdu.writeto(params.outputfn + "_qmap.fits", **overwrite)
    
    hdu.data[data_selection] = data_out.imag[range(len(data_out)), indx_max]
    print("Saving u map")
    hdu.writeto(params.outputfn + "_umap.fits", **overwrite)
    
    header["BUNIT"] = "rad/m/m"
    hdu.data[data_selection] = params.phi[indx_max]
    print("Saving phi map")
    hdu.writeto(params.outputfn + "_phi.fits", **overwrite)
    

if __name__ == "__main__":
    main()
