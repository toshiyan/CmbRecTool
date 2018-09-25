
import numpy as np, healpy as hp

class analysis_choice: 
  # set variables here

  print 'set analysis parameters'

  # rlz num
  sn = 100

  # map definition
  nside = 2048
  npix  = 12*2048**2

  # mask window
  #wind  = 'apo5'
  #fsky  = 'fsky70'
  wind = 'smica'

  # number of bins
  bmax  = 20

  # data type
  dtype = 'smica'
  #freq  = 143

  # multipoles
  rL    = [100,2048]
  eL    = [0,3000]
  oL    = [0,3000]


class filename:
  # set filename

  # load variables
  ac = analysis_choice()

  # directory
  pldir = '/project/projectdirs/cmb/data/planck2015/'
  mydir = '/global/cscratch1/sd/toshiyan/plk/'

  # tags
  ltag  = 'l'+str(ac.rL[0])+'-'+str(ac.rL[1])

  # mask window
  #wind = mydir+'/map/HFI_Mask_GalPlane-'+ac.wind+'_2048_R2.00.fits'
  wind = mydir+'/map/COM_Mask_CMB-PointSrcGalplane-smica-harmonic-mask_2048_R2.00.fits'
  ptsr = mydir+'/map/HFI_Mask_PointSrc_2048_R2.00.fits'

  # cmb and noise maps
  smap = ['']*1001
  nmap = ['']*1001
  #smap[0] = pldir+'pr2/PR2/frequencymaps/HFI_SkyMap_'+str(ac.freq)+'_2048_R2.02_full.fits'
  smap[0] = pldir+'pr2/PR2/cmbmaps/COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits'
  for i in range(1,1001):
    #smap[i+1] = pldir+'ffp8/mc_cmb/'+str(ac.freq)+'/ffp8_cmb_scl_143_full_map_mc_'+str(i).zfill(5)+'.fits'
    #nmap[i+1] = pldir+'ffp8/mc_noise/'+str(ac.freq)+'/ffp8_noise_143_full_map_mc_'+str(i).zfill(5)+'.fits'
    smap[i] = pldir+'ffp8/compsep/mc_cmb/ffp8_smica_int_cmb_mc_'+str(i-1).zfill(5)+'_005a_2048.fits'
    nmap[i] = pldir+'ffp8/compsep/mc_noise/ffp8_smica_int_noise_mc_'+str(i-1).zfill(5)+'_005a_2048.fits'

  # f90 map file
  f90s = ['']*1001
  f90n = ['']*1001
  f90s[0] = mydir+'map/real_'+ac.dtype+'.bin'
  for i in range(1,1001):
    f90s[i] = mydir+'map/lcdm_'+ac.dtype+'_'+str(i)+'.bin'
    f90n[i] = mydir+'map/nois_'+ac.dtype+'_'+str(i)+'.bin'
  f90m = mydir+'map/mask_'+ac.dtype+'.bin'

  # alms
  salm = ['']*1001
  nalm = ['']*1001
  calm = ['']*1001
  wtag = ac.wind
  calm[0] = mydir+'alm/real_'+wtag+'.fits'
  for i in range(1000):
    salm[i+1] = mydir+'alm/lcdm_'+wtag+'_'+str(i).zfill(5)+'.fits'
    nalm[i+1] = mydir+'alm/nois_'+wtag+'_'+str(i).zfill(5)+'.fits'
    calm[i+1] = mydir+'alm/comb_'+wtag+'_'+str(i).zfill(5)+'.fits'

  # tlm
  tlm = ['']*1001
  tlm[0] = mydir+'tlm/'+wtag+'_real.fits'
  for i in range(1000):
    tlm[i+1] = mydir+'tlm/'+wtag+'_'+str(i).zfill(5)+'.fits'

  # aps
  ilcl = mydir+'aps/lensedfid_P15.dat'
  mcl  = mydir+'aps/mcl_'+wtag+'.dat'
  ocl  = mydir+'aps/ocl_'+wtag+'.dat'
  fl   = mydir+'aps/fl_'+wtag+'.dat'


class recfunc:
  # set arrays

  # set parameters
  ac = analysis_choice()
  f  = filename()

  # galactic mask
  #if ac.fsky=='fsky20': mask = hp.fitsfunc.read_map(f.wind,field=0)
  #if ac.fsky=='fsky40': mask = hp.fitsfunc.read_map(f.wind,field=1)
  #if ac.fsky=='fsky60': mask = hp.fitsfunc.read_map(f.wind,field=2)
  #if ac.fsky=='fsky70': mask = hp.fitsfunc.read_map(f.wind,field=3)
  mask = hp.fitsfunc.read_map(f.smap[0],field=1)

  # ptsr mask 
  #if ac.freq==143: ptsr = hp.fitsfunc.read_map(f.ptsr,field=1)

  # window norm
  w2 = np.average(mask**2)
  print w2
  #w3 = np.average(mask**3)
  #w4 = np.average(mask**4)

  # theoretical cl
  ilcl = np.zeros((2,ac.eL[1]+1))
  L, TT = (np.loadtxt(f.ilcl,usecols=(0,1))).T
  TT = TT/(L**2+L)*2*np.pi
  ilcl[0,2:] = TT[:ac.eL[1]-1]
  ilcl[1,2:] = 1./TT[:ac.eL[1]-1]

  # multipole array
  L = np.linspace(0,ac.eL[1],ac.eL[1]+1)

  # rec multipole filter
  fil = np.zeros(ac.eL[1]+1)
  fil[ac.rL[0]:ac.rL[1]+1] = 1

