#!/usr/bin/env python

import os, argparse

def add(text,f):
  os.system('echo "'+text+'" >> '+f)

print 'creating slurm jobsub file'

# read arguments
parser = argparse.ArgumentParser( description='matlab job submit' )
parser.add_argument( '--func', required=True)
parser.add_argument( '--args', required=False, default='')
parser.add_argument( '--mem', required=False, default='40')
parser.add_argument( '--time', required=False, default='0-12:00')
args = parser.parse_args()

# set tag
jtag = args.func+str(args.args)

# create job file
f = 'submit_'+args.func+'_'+str(args.args)+'.sh'
os.system('echo "#!/bin/bash" > '+f)
add('#SBATCH -n 1',f)
add('#SBATCH -N 1',f)
add('#SBATCH -p serial_requeue',f)
add('#SBATCH --mem='+args.mem+'G',f)
add('#SBATCH -t '+args.time,f)
add('#SBATCH --job-name='+jtag,f)
add('#SBATCH --mail-type=ALL',f)
add('#SBATCH --mail-user=toshiyan@stanford.edu',f)
add('module load matlab',f)
add('echo \''+args.func+'('+args.args+')'+'\' > tmp_'+jtag+'.m',f)
add('cat tmp_'+jtag+'.m',f)
add('matlab -nojvm -nodisplay -nosplash < tmp_'+jtag+'.m',f)
add('rm tmp_'+jtag+'.m',f)
#os.system('sbatch '+f)

