#!/home/gerold/anaconda3/envs/iraf27/bin/python

# Author: Gerold Busch
# Created:      2018-08-20 G. Busch




# import modules

import os, sys, shutil
import re
from argparse import ArgumentParser

from pyraf import iraf
from iraf import noao
from iraf import mscred

import astropy.io.fits as fits
from astropy import units as u
from astropy.coordinates import SkyCoord

import numpy as np





class Reduction(object):
    """
    Class to reduce LUCI images
    
    Description:
    ------------
    1. 
    
    
    TODO:
    * 
    
    
    """
    
    
    def __init__(self, inpfiles, outpname="./tmp/out.fits", tempdir="./tmp/"):
        """
        Init
        
        Input:
        ------
        * inpfiles: list
            list with input files
                or path to directory with files
        * outpname: string
            filename of the output master dome flat
        * tempdir: string
            directory to store temp files
        
        
        --> 
        """
        
        if type(inpfiles)==type(list()):
            imglist = [os.path.abspath(f) for f in sorted(inpfiles) if re.match(r"luci\d.\d{8}.\d{4}.fits", f)]
        elif os.path.isdir(inpfiles):
            frlist = os.listdir(inpfiles)
            path = os.path.abspath(inpfiles)
            imglist = [os.path.join(path,f) for f in sorted(frlist) if re.match(r"luci\d.\d{8}.\d{4}.fits", f)]
        else:
            raise Exception("Cannot read input list")
        
        self.inpfiles = imglist
        self.outpname = os.path.abspath(outpname)
        self.tempdir = os.path.abspath(tempdir)
        
        
        
        
    def getWCSOffsets(self, offset_file="WCSoffsets.dat"):
        
        offset_pathfile = os.path.join(self.tempdir, offset_file)
        offsets = np.zeros((len(self.inpfiles), 2))
        ref_img = self.inpfiles[0]
        h = fits.open(ref_img)[0].header
        pixsc = h['PIXSCALE']
        c0 = SkyCoord(h['TELRA'], h['TELDEC'], unit=(u.hourangle, u.deg))
        
        with open(offset_pathfile, 'w') as txt:
            for i, img in enumerate(self.inpfiles):
                h = fits.open(img)[0].header
                c = SkyCoord(h['TELRA'], h['TELDEC'], unit=(u.hourangle, u.deg))
                
                offsets[i,0] = (c.ra - c0.ra).arcsec/pixsc
                offsets[i,1] = (c0.dec - c.dec).arcsec/pixsc
                txt.write("{}\t{:6f}\t{:6f}\n".format(img, offsets[i][0], offsets[i][1]))
        
        return offsets
        
        
        
        
        
if __name__ == "__main__":
    
    usage = "reduce.py [options]"
    desc = """This module takes file lists (or directories) and outputs a reduced image."""
    
    parser = ArgumentParser(usage, description=desc)
    
    parser.add_argument("-in", "--files_in", action="store", dest="inpfiles", help="input files")
    parser.add_argument("-out", "--output_filename", action="store", dest="outpname", help="path+filename for output file", default="./tmp/out.fits")
    parser.add_argument("-temp", "--temp_dir", action="store", dest="tempdir", help="directory for temporary files", default="./tmp/")

    
    
    args = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)
       
    if not args.inpfiles:
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    red = Reduction(args.inpfiles, args.outpname, args.tempdir)
    red.getWCSOffsets()
