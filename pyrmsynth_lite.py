#!/usr/bin/env python

"""
pyrmsynth_lite, Stephen Bourke

Stripped down from pyrmsynth (Michael Bell, Henrik Junklewitz and Sarrvesh Sridhar)
"""

from optparse import OptionParser
from astropy.io import fits
import numpy as np
import csv
import numpy
import rm_tools as R
from timeit import default_timer as tic
import math
import sys

VERSION = "0.1"


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
        print 'Non-trivial weights enabled! Loaded from ' + parset['do_weight']

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


def main():
    """ Handle all parsing here if started from the command line"""

    parser = OptionParser(usage="%prog <input parameter file>",
                          version="%prog " + VERSION)
    parser.add_option("-p", "--save-pol-cube", action="store_true",
                      dest="save_pol_cube",
                      help="Save Pol cube", default=False)
    parser.add_option("-q", "--save-qu-cubes", action="store_true",
                      dest="save_qu_cubes",
                      help="Save derotated Q and U cubes", default=False)
    parser.add_option("-a", "--auto-flag", action="store_true",
                      dest="auto_flag",
                      help="auto flag data", default=False)
    parser.add_option("-x", "--exclude_phi", metavar='phi_range', nargs=2,
                      type=float, default=(0,0), 
                      help="exclude this range from moment maps. Eg: -3 1.5")
    parser.add_option("-n", "--phi_rms", metavar='phi_rms_range', nargs=2,
                      type=float, default=(0,0), 
                      help="Make an RMS image from this range. Eg: 100 115")

    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error("Incorrect number of arguments")

    print "Parsing parameter file..."
    params = parse_input_file(args[0])

    timing = []
    timing.append(tic())
    
    print "Loading data..."
    hdu = fits.open(args[1])[0]
    data = hdu.data
    header = hdu.header
    if params.imagemask:
        mask = fits.open(params.imagemask)[0].data.squeeze().astype(bool)
    else:
        mask = np.ones(shape=data.shape[2:], dtype=bool)

    if options.auto_flag:
        bad_chans = channel_flags(data)
        flag_weights = np.logical_not(bad_chans).astype(float)
        if params.weight is None:
            params.weight = flag_weights
        else:
            params.weight = np.minimum(params.weight, flag_weights)
    
    # Survey cubes will be ordered: [Freqs, Stokes, Dec, RA]
    data = np.moveaxis(data, 0, -1) # Move Freq to end
    data = np.moveaxis(data, 0, -1) # Move Stokes to end
    # Now order is [Dec, RA, Freqs, Stokes]
    (DEC,RA,FREQ,STOKES) = range(4)

    data_selection = np.where(mask)
    data_sel = data[data_selection]
    # Unflagged NaNs will cause RMSynth to return NaNs
    data_sel = np.nan_to_num(data_sel.view(dtype=">c8").squeeze()).astype(np.complex128)
    
    params.nu = freqs(header, 4)
    params.dnu = params.nu[1] - params.nu[0]
    rms = R.RMSynth(params.nu, params.dnu, params.phi, params.weight)

    rmsfout = numpy.zeros((len(rms.rmsf), 3))
    rmsfout[:, 0] = rms.rmsf_phi
    rmsfout[:, 1] = rms.rmsf.real
    rmsfout[:, 2] = rms.rmsf.imag
    numpy.savetxt(params.outputfn + "_rmsf.txt", rmsfout)

    timing.append(tic())

    data_out = np.empty(shape=(len(data_sel), params.nphi), dtype=np.complex128)
    print "Doing RM synthesis"
    n = len(data_sel)
    width = len(str(n))
    for i in range(len(data_sel)):
        data_out[i] = rms.compute_dirty_image(data_sel[i]) 
        if i % 100 == 0:
            progress(20, 100.0 * i / n)
    print ""

    timing.append(tic())
    
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

    if options.save_pol_cube:
        cube_out = np.full(fill_value=np.NaN, shape=(data.shape[DEC], data.shape[RA], params.nphi), dtype=data_out.real.dtype)
        cube_out[data_selection] = abs(data_out)
        cube_out = np.moveaxis(cube_out, -1, 0)
        hdu.data = cube_out
        print "Saving phi cube"
        hdu.writeto(params.outputfn + "_di_p.fits", overwrite=True)

        del cube_out
    
    if options.save_qu_cubes:
        cube_out = np.full(fill_value=np.NaN, shape=(data.shape[DEC], data.shape[RA], params.nphi), dtype=data_out.real.dtype)
        cube_out[data_selection] = data_out.real
        cube_out = np.moveaxis(cube_out, -1, 0)
        hdu.data = cube_out
        print "Saving q cube"
        hdu.writeto(params.outputfn + "_di_q.fits", overwrite=True)
    
        cube_out = np.moveaxis(cube_out, 0, -1)
        cube_out[data_selection] = data_out.imag
        cube_out = np.moveaxis(cube_out, -1, 0)
        hdu.data = cube_out
        print "Saving u cube"
        hdu.writeto(params.outputfn + "_di_u.fits", overwrite=True)
    
        del cube_out
    
    for attr in header_axes_attributes:
        attr += "3"
        if attr in header:
            del header[attr]

    hdu.data = np.full(fill_value=np.NaN, shape=(data.shape[DEC], data.shape[RA]), dtype=data_out.real.dtype)

    # RMS Image
    if options.phi_rms[0] != options.phi_rms[1]:
        phi0_idx = np.abs(params.phi - options.phi_rms[0]).argmin()
        phi1_idx = np.abs(params.phi - options.phi_rms[1]).argmin()
        x_phi = (min(phi0_idx, phi1_idx), max(phi0_idx, phi1_idx)+1)
        print("RMS phi range: {} => {}".format(options.phi_rms, x_phi))
        q_rms = data_out[:,x_phi[0]:x_phi[1]].real.std(axis=1)
        u_rms = data_out[:,x_phi[0]:x_phi[1]].imag.std(axis=1)
        p_rms = 0.5 * (q_rms + u_rms)

        hdu.data[data_selection] = p_rms
        print "Saving rms map"
        hdu.writeto(params.outputfn + "_rms.fits", overwrite=True)

    # Exclude range:
    if options.exclude_phi[0] != options.exclude_phi[1]:
        phi0_idx = np.abs(params.phi - options.exclude_phi[0]).argmin()
        phi1_idx = np.abs(params.phi - options.exclude_phi[1]).argmin()
        x_phi = (min(phi0_idx, phi1_idx), max(phi0_idx, phi1_idx)+1)
        print("Excluding phi range: {} => {}".format(options.exclude_phi, x_phi))
        data_out[:,x_phi[0]:x_phi[1]] = 0

    hdu.data[data_selection] = np.amax(abs(data_out), axis=-1)
    print "Saving int pol map"
    hdu.writeto(params.outputfn + "_polint.fits", overwrite=True)

    indx_max = np.argmax(abs(data_out), axis=-1)
    hdu.data[data_selection] = data_out.real[range(len(data_out)), indx_max]
    print "Saving q map"
    hdu.writeto(params.outputfn + "_qmap.fits", overwrite=True)
    
    hdu.data[data_selection] = data_out.imag[range(len(data_out)), indx_max]
    print "Saving u map"
    hdu.writeto(params.outputfn + "_umap.fits", overwrite=True)
    
    header["BUNIT"] = "rad/m/m"
    hdu.data[data_selection] = params.phi[indx_max]
    print "Saving phi map"
    hdu.writeto(params.outputfn + "_phi.fits", overwrite=True)
    
    timing.append(tic())

if __name__ == "__main__":
    main()
