#////////////////////////////////////////////////////////////////////////////////////////////////////#
# * An example code to show 2D projected map
# - Last Modified: Sat 15 Jul 2017 12:05:12 PM PDT
#
# Example
#   python showmap_proj.py --file=tmp/map_1.dat --ncol=2 --scale=1e6
#
#////////////////////////////////////////////////////////////////////////////////////////////////////#

import argparse
import numpy as np

# read arguments
parser = argparse.ArgumentParser( description='plot 2D map' )
parser.add_argument( '--file', required=True)
parser.add_argument( '--ncol', required=False, default='')
parser.add_argument( '--nx', required=False, default=236)
parser.add_argument( '--ny', required=False, default=100)
parser.add_argument( '--scale', required=False, default=1.)
parser.add_argument( '--label', required=False, default='')
parser.add_argument( '--xlabel', required=False, default='')
parser.add_argument( '--ylabel', required=False, default='')
parser.add_argument( '--vmin', required=False, default='')
parser.add_argument( '--vmax', required=False, default='')
parser.add_argument( '--vmaxnorm', required=False, default=False)
parser.add_argument( '--ofig', required=False, default='test.png')
parser.add_argument( '--show', required=False, default=False)
parser.add_argument( '--log', required=False, default=False)
parser.add_argument( '--xflip', required=False, default=False)
parser.add_argument( '--yflip', required=False, default=False)
args = parser.parse_args()

# confirm inputs
nx = int(args.nx)
ny = int(args.ny)
print 'map grids: ', nx, 'x', ny

if args.vmaxnorm: vmax = np.max(W)

if args.ncol!='':
  W = np.loadtxt(args.file,unpack=True)[int(args.ncol)]
else:
  W = np.loadtxt(args.file,unpack=True)

# check
print 'input array', W, np.shape(W)

# convert into 2D array
W = np.reshape(W,(nx,ny)).T * np.float(args.scale)

# start plot
import matplotlib as mpl
from matplotlib.pyplot import *
from matplotlib.colors import LogNorm

mpl.rcParams.update({'font.size': 15})
subplot(1,1,1,aspect=nx/ny)

# coord
xlim(0,nx)
ylim(0,ny)
if args.xflip: xlim(nx,0)
if args.yflip: ylim(ny,0)

# label
xlabel(args.xlabel)
ylabel(args.ylabel)

if args.log: W = np.log(W)
if args.vmin!='' and args.vmax=='': map = pcolor(W,vmin=float(args.vmin))
if args.vmin=='' and args.vmax!='': map = pcolor(W,vmax=float(args.vmax))
if args.vmin=='' and args.vmax=='': map = pcolor(W)
if args.vmin!='' and args.vmax!='': map = pcolor(W,vmin=float(args.vmin),vmax=float(args.vmax))

cb  = colorbar(map,orientation='horizontal',fraction=.25)
cb.set_label(r''+args.label,labelpad=5)
savefig(args.ofig)
if args.show: show()

