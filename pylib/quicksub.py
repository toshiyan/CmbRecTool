# * quick submit
import os

def add(text,f,ini=False):
  if ini:     os.system('echo "'+text+'" > '+f)
  if not ini: os.system('echo "'+text+'" >> '+f)

def set_sbatch_params(f,tag,n=1,N=1,queue='regular',mem='20G',t='0-12:00',email=False,host='nersc'):
  add('#!/bin/bash',f,True)
  add('#SBATCH -n '+str(n),f)
  add('#SBATCH -N '+str(N),f)
  if queue=='regular' and host=='pdy':  queue='general'
  add('#SBATCH -p '+queue,f)
  add('#SBATCH --mem='+mem,f)
  add('#SBATCH -t '+t,f)
  add('#SBATCH --job-name='+tag,f)
  if host=='nersc':
    add('#SBATCH -C haswell',f)
  if host=='ody':
    add('#SBATCH --open-mode=append',f)
  if email==True:
    add('#SBATCH --mail-type=END,FAIL',f)
    add('#SBATCH --mail-user=toshiyan@stanford.edu',f)

