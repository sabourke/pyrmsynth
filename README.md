# `pyrmsynth_lite`

`pyrmsynth_lite` is a stripped-down and modified version of `pyrmsynth`
intended for LoTSS DR2 data.
It should be much faster than standard `pyrmsynth`
for masked data and slightly faster for full cube synthesis.

Only the main `pyrmsynth` routine is modified. `pyrmsynth_lite` uses
the same code to do the actual RM synthesis. The speed up is achived
by optimising the memory access pattern.

For development and testing the RM-Clean capability was removed. This
should be restored at some point.

## Installing

Dependancies: `git python-astropy cython libgsl-dev`
```
sudo apt install git python-astropy cython libgsl-dev
```

The `rm_tools` component needs to be build and installed. If you
already have a working pyrmsynth you may not need to do this,
in which case you can skip this step.
```
cd rm_tools
python setup.py install --user
```

Copy the pyrmsynth_lite.py program to a directory in your path.
```
mkdir -p ~/.local/bin
chmod a+x pyrmsynth_lite.py
cp pyrmsynth_lite.py ~/.local/bin
```


## Usage
```
python rmsynthesis.py <parfile> QU_cube.fits
```

By default, `pyrmsynth_lite` outputs only the 2D maps. To save the qubes,
use the command line options as listed below.

`pyrmsynth_lite` includes an auto-flagging rountine as developed by G. Heald.
To enable it see the command line options below.


Options:
```
  -p, --save-pol-cube  Save Pol cube
  -q, --save-qu-cubes  Save derotated Q and U cubes
  -a, --auto-flag      auto flag data
```

`pyrmsynth_lite` is licensed under the [GPLv3](http://www.gnu.org/licenses/gpl.html).
