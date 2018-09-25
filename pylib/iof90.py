
import numpy as np
from scipy.io import FortranFile

def writef90map(fname,d1):
  f = FortranFile(fname,'w')
  #f.write_record(np.int16(dim))
  f.write_record(np.float64(d1))
  f.close()

def readf90map(fname):
  f = FortranFile(fname,'r')
  dat = f.read_ints(dtype=np.float64)
  f.close()
  return dat

def readfile(fname,dt=np.float64,hdrnum=0):
  f = open(fname,'r')
  if hdrnum>0: 
    hdr = np.fromfile(f,dtype=np.int32,count=hdrnum)
    print hdr
  dat = np.fromfile(f,dtype=dt)
  f.close()
  return dat

