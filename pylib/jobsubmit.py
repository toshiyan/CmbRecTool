#!/usr/bin/env python

import argparse
import os

def add(text,f):
  os.system('echo "'+text+'" >> '+f)

print 'creating slurm jobsub file'

# read arguments
parser = argparse.ArgumentParser( description='job submit' )
parser.add_argument( '--exe', required=True)
parser.add_argument( '--ini', required=False, default='')
parser.add_argument( '--mem', required=False, default='4')
parser.add_argument( '--time', required=False, default='0-12:00')
args = parser.parse_args()

# set tag
jtag = args.exe+'_'+args.ini

# create job file
f = 'submit_'+jtag+'.sh'
os.system('echo "#!/bin/bash" > '+f)
add('#SBATCH -n 4',f)
add('#SBATCH -N 1',f)
add('#SBATCH -p serial_requeue',f)
add('#SBATCH --mem='+args.mem+'G',f)
add('#SBATCH -t '+args.time,f)
add('#SBATCH --job-name='+jtag,f)
add('#SBATCH --mail-type=ALL',f)
add('#SBATCH --mail-user=toshiyan@stanford.edu',f)
add(args.exe+' '+args.ini,f)
os.system('sbatch '+f)

