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
import numpy as np


def tryrem(path, filen):
    try:
        os.remove(os.path.join(path, filen))
    except:
        print("file {} did not exist".format(filen))




class MasterDomeFlat(object):
    """
    Class to build master dome flat.
    Assumes that flat files are already sorted according to on/off, band, and LUCI1/2
    (Sorting script todo, for now this should be done by hand)
    
    Description:
    ------------
    1. Combine lamp-on frames
    2. Combine lamp-off frames
    3. subtract on-off
    4. normalize with median
    
    
    TODO:
    * Check mess with Paths...
    
    
    """
    
    
    def __init__(self, inpfiles_on, inpfiles_off, outpname="./tmp/flat.fits", tempdir="./tmp/", normal=True):
        """
        Init
        
        Input:
        ------
        * inpfiles_on: list
            list with lamp-on input files
                or path to directory with files
        * inpfiles_off: list
            list with lamp-off input files
        * tempdir: string
            directory to store temp files
        * outpname: string
            filename of the output master dome flat
        * normal: bool
            if true, flat will be normalized
        
        --> produces master dome flat, that will be saved in file, as given by outpname
        """
        
        self.__inpfiles_on = inpfiles_on
        self.__inpfiles_off = inpfiles_off
        self.__outpname = os.path.abspath(outpname)
        self.__tempdir = os.path.abspath(tempdir)
        self.__normal = normal
        
    
    def createMaster(self):
        """
        Create master dome flat from dome flat file lists
        """
        
        # remove old file
        list_lampon = []
        list_lampoff = []
        
        if type(self.__inpfiles_on)==type(list()):
            list_lampon = self.__inpfiles_on
        elif os.path.isdir(self.__inpfiles_on):
            frlist = os.listdir(self.__inpfiles_on)
            path = os.path.abspath(self.__inpfiles_on)
            list_lampon = [os.path.join(path,f) for f in sorted(frlist) if re.match(r"luci\d.\d{8}.\d{4}.fits", f)]
        else:
            raise Exception("Cannot read flat-off files")
            
        if type(self.__inpfiles_off)==type(list()):
            list_lampoff = self.__inpfiles_off
        elif os.path.isdir(self.__inpfiles_off):
            frlist = os.listdir(self.__inpfiles_off)
            path = os.path.abspath(self.__inpfiles_off)
            list_lampoff = [os.path.join(path,f) for f in sorted(frlist) if re.match(r"luci\d.\d{8}.\d{4}.fits", f)]
        else:
            raise Exception("Cannot read flat-off files")
            
        
        
        if not os.path.exists(self.__tempdir):
            os.makedirs(self.__tempdir)
        
        ### IRAF will overwrite existing files
        iraf.clobber='yes'
        base, infile = os.path.split(self.__outpname)
        
        if not os.path.exists(base):
            os.makedirs(base)
        
        iraf.chdir(base)
            
        ### combine on-frames
        flat_lampon = os.path.join(self.__tempdir, "flat_lampON.fits")
        with open(os.path.join(self.__tempdir, "files_on.list"), 'w') as fo:
            for f in list_lampon:
                fo.write(str(f)+"\n")
        
        iraf.mscred.flatcombine(input="@"+(os.path.join(self.__tempdir, "files_on.list")).replace('//','/'),
                        output=flat_lampon,
                        combine='median',
                        ccdtype='',
                        process='no',
                        reject='sigclip',
                        subset='no',
                        scale='mode'
                        )
        
        
        ### combine off-frames
        flat_lampoff = os.path.join(self.__tempdir, "flat_lampOFF.fits")
        with open(os.path.join(self.__tempdir, "files_off.list"), 'w') as fo:
            for f in list_lampoff:
                fo.write(str(f)+"\n")
        
        iraf.mscred.flatcombine(input="@"+(os.path.join(self.__tempdir, "files_off.list")).replace('//','/'),
                        output=flat_lampoff,
                        combine='median',
                        ccdtype='',
                        process='no',
                        reject='sigclip',
                        subset='no',
                        scale='mode'
                        )
        
        
        ### subtract files
        flat_diff = os.path.join(self.__tempdir, "flat_lampON_OFF.fits")
        iraf.imarith(operand1 = flat_lampon,
                    operand2 = flat_lampoff,
                    op = '-',
                    result = flat_diff
                    )
                    
                    
        ### Normalize flat
        if self.__normal:
            f = fits.open(flat_diff)
            naxis1 = f[0].header['NAXIS1']
            naxis2 = f[0].header['NAXIS2']
            offset1 = int(naxis1*0.1)
            offset2 = int(naxis2*0.1)
            median = np.median(f[0].data[offset2:(naxis2-offset2), offset1:(naxis1-offset1)])
            
            iraf.imarith(operand1=flat_diff, operand2=median, op='/', result=self.__outpname.replace("//","/"))
        
        else:
            shutil.move(flat_diff, self.__outpname)
            
        
        iraf.chdir()
        
        
        ### Clean up
        tryrem(self.__tempdir, "flat_lampOFF.fits")
        tryrem(self.__tempdir, "flat_lampON.fits")
        tryrem(self.__tempdir, "flat_lampON_OFF.fits")
        tryrem(self.__tempdir, "files_on.list")
        tryrem(self.__tempdir, "files_off.list")
        
        
        
if __name__ == "__main__":
    
    usage = "DomeFlat.py [options]"
    desc = """This module takes file lists (or directories) for Dome Flats (on and off) and creates a Master Dome Flat from this."""
    
    parser = ArgumentParser(usage, description=desc)
    
    parser.add_argument("-on", "--files_on", action="store", dest="inpfiles_on", help="input files for lamp on, can be list or path to directory with lamp-on files")
    parser.add_argument("-off", "--files_off", action="store", dest="inpfiles_off", help="input files for lamp off, can be list or path to directory with lamp-off files")
    parser.add_argument("-out", "--output_filename", action="store", dest="outpname", help="path+filename for output file", default="./tmp/flat.fits")
    parser.add_argument("-temp", "--temp_dir", action="store", dest="tempdir", help="directory for temporary files", default="./tmp/")
    parser.add_argument("-norm", "--normalize", action="store_true", dest="norm", help="normalize flat (default: True)", default=True)
    
    
    
    args = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)
       
    if not args.inpfiles_on or not args.inpfiles_off:
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    mdFlat = MasterDomeFlat(args.inpfiles_on, args.inpfiles_off, args.outpname, args.tempdir, args.norm)
    mdFlat.createMaster()
