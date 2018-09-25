#////////////////////////////////////////////////////////////////////////////////////////////////////#
# * An example code to show fullsky map
# - Last Modified: Mon 25 Jul 2016 06:16:29 PM EDT
#
# Example
#   python showmap_full.py --file=tmp/test.fits
#
#////////////////////////////////////////////////////////////////////////////////////////////////////#

import argparse

# read arguments
parser = argparse.ArgumentParser( description='plot 2D map' )
parser.add_argument( '--file', required=True)
parser.add_argument( '--field', required=False, default=0)
parser.add_argument( '--output', required=False, default='test.png')
parser.add_argument( '--coord_in', required=False, default='G')
parser.add_argument( '--coord_out', required=False, default='C')
args = parser.parse_args()

import healpy as hp
from matplotlib.pyplot import *

map = hp.fitsfunc.read_map(args.file,field=int(args.field))
hp.mollview(map,coord=[args.coord_in,args.coord_out])
hp.graticule()
savefig(args.output)
show()

