# pyrmsynth_lite

`pyrmsynth_lite` is a modified version of `pyrmsynth` intended for LoTSS DR2 data. It should be much faster than standard `pyrmsynth` for masked data and slightly faster for full cube synthesis.

Only the main `pyrmsynth` routine is modified. `pyrmsynth_lite` uses the same code to do the actual RM synthesis. The speed-up is achived by optimising the memory access pattern.

Additional features over standard `pyrmsynth`:
- Use LoTSS-DR2 cubes directly
- Significantly faster when using a mask
- Save only wanted output cubes
- Auto flag data after loading
- Exclude a specified Phi range from the integrated maps
- Generate a noise map from a specified Phi range
- Save output in single or double precision format
- Use an input mask to specify areas to do RMCLEAN

## Installing

If you already have the original `pyrmsynth` installed and working, you only need the `pyrmsynth_lite.py` file. Otherwise do the following applies.

Dependancies: `git python-astropy cython libgsl-dev`
```
sudo apt install git python-astropy cython libgsl-dev
```
Build and install the rm_tools component:
```
cd rm_tools
python setup.py install --user
```
Copy the pyrmsynth_lite.py program to a directory in your path:
```
mkdir -p ~/.local/bin
cp pyrmsynth_lite.py ~/.local/bin
```
## Usage
```
pyrmsynth_lite.py <parfile> <QU_cube.fits>
```
By default, `pyrmsynth_lite` outputs only the 2D maps. To save the cubes, use the command line arguments as listed below.

`pyrmsynth_lite` includes an auto-flagging rountine as developed by G. Heald. To enable it see the command line options below.

Default behaviour is to use double precision internally and save output as single precision. This can be modifiled with the `--single` and `--double` command line arguments.
```
usage: pyrmsynth_lite.py [-h] [-o OUTNAME] [-p] [-q] [-r] [-d] [-a]
                         [-x phi_range phi_range]
                         [-n phi_rms_range phi_rms_range] [--single]
                         [--double] [-c CLEAN_MASK]
                         param_file fits_cube

Rotation Measure Synthesis tool.

positional arguments:
  param_file            Input Parameter file
  fits_cube             QU FITS cube

optional arguments:
  -h, --help            show this help message and exit
  -o OUTNAME, --outname OUTNAME
                        output prefix
  -p, --save-pol-cube   Save Pol cube
  -q, --save-qu-cubes   Save derotated Q and U cubes
  -r, --save-residual-cube
                        Save residual cubes if cleaning
  -d, --save-dirty-cube
                        Save dirty cubes if cleaning
  -a, --auto-flag       auto flag data
  -x phi_range phi_range, --exclude_phi phi_range phi_range
                        Exclude this Phi range from 2D maps. Eg: -3 1.5
  -n phi_rms_range phi_rms_range, --phi_rms phi_rms_range phi_rms_range
                        Make a Phi RMS image from this range. Eg: 100 115
  --single              Use single precision floats internally. Default: False
  --double              Output double precision floats. Default: False
  --clean_mask CLEAN_MASK
                        Input mask for RM clean
```
A utility called `mask_map.py` is included to make masks based on a cutoff value in Janskys. This is inteded to be used with a Stokes I map that corresponds to the cube that will be processed. E.g.:
```
mask_map.py --cutoff 0.001 imgI.fits
```
`pyrmsynth_lite` is licensed under the [GPLv3.](http://www.gnu.org/licenses/gpl.html)
