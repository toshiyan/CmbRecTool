
import matplotlib as mpl
import numpy as np
from matplotlib.pyplot import *
from matplotlib.colors import LogNorm

nx = 236 
ny = 100
x  = np.linspace(-55,55,nx)
y  = np.linspace(-70,-45,ny)
r  = 110./236.*110./25. # aspect ratio
xt = [-50,0,50]
yt = [-70,-65,-60,-55,-50,-45]
f  = ['grec1_g.dat','grec1_c.dat']
v  = .04
dir0 = './'

mpl.rcParams.update({'font.size': 12})
fig, ax = subplots(nrows=1,ncols=2,figsize=(10,3),dpi=100)
for i in range(2):
  ax[i].set_xlim(55.,-55.)
  ax[i].set_ylim(-70.,-45.)
  ax[i].set_xticks(xt)
  if i == 0: ax[i].set_yticks(yt)
  if i == 1: ax[i].set_yticks([])
  if i == 0: ax[i].set_xlabel('Right Ascension [deg.]')
  if i == 0: ax[i].set_ylabel('Declination [deg.]')
  ax[i].set_aspect(r)
  W = np.reshape(np.loadtxt(dir0+f[i]).T,(nx,ny)).T
  a = ax[i].pcolor(x,y,W,vmin=-v,vmax=v)

cax = fig.add_axes([0.77,0.25,.01,0.5])
cb = fig.colorbar(a,cax=cax,ticks=[-v,.0,v])
cb.set_label(r'',labelpad=5,rotation=270,fontsize=20)
fig.text(0.15,0.7,r'$\kappa$',fontsize=10)
fig.text(0.45,0.7,r'curl',fontsize=10)
subplots_adjust(bottom=0.1,right=0.75,top=0.9,wspace=0)
savefig("test.png")
show()

