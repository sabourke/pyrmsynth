#!/usr/bin/env python

# Create a fits mask based on an input map and a threshold value
# Stephen Bourke
# Onsala Space Observatory / Chalmers
# Nov 2018

import numpy as np
from astropy.io import fits
import os.path
import argparse

DEFAULT_CUTOFF = 0.001 # 1 mJy

def main():
    parser = argparse.ArgumentParser(description='Generate a mask file from an input map')
    parser.add_argument('-c', '--cutoff', metavar='flux', type=float, default=DEFAULT_CUTOFF, help='Flux cutoff in Jy. Eg. 0.001')
    parser.add_argument('fitsfile', type=str, help='Input fits file')
    args = parser.parse_args()

    image = fits.open(args.fitsfile)[0]
    mask = np.zeros_like(image.data, dtype=np.uint8)
    mask[np.where(image.data > args.cutoff)] = 1
    image.data = mask.squeeze()
    image.writeto('{}.mask{}'.format(*os.path.splitext(args.fitsfile)))


if __name__ == "__main__":
    main()
