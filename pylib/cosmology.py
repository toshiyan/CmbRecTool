# Python tools for cosmology
import numpy as np
import scipy.constants as spc
from scipy.integrate import quad

c0kms = spc.c/1e3


class cosmofunc:


  def __init__(self,Obh2=.0223,Omh2=.142,H0=67.5,OK=0.,As=2.13e-9,ns=.965,k0=.05,tau=.06):
    self.Obh2 = Obh2
    self.Omh2 = Omh2
    self.H0   = H0/c0kms #[1/Mpc]
    self.OK   = OK
    self.As   = As
    self.ns   = ns
    self.k0   = k0       #[Mpc]
    self.tau  = tau

    # derived cp
    self.OmH2 = Omh2*(100./c0kms)**2 #[1/Mpc^2]
    self.h0   = H0/100.
    self.Om   = self.Omh2/self.h0**2
    self.Ob   = self.Obh2/self.h0**2
    self.Ov   = 1. - self.Om - self.OK


  # Transfer function (BBKS)
  def tk_bbks(self,k):
    q  = k/self.Omh2
    return np.log(1.+2.34*q)/(2.34*q) / np.power(1.+3.89*q+(16.2*q)**2+(5.47*q)**3+(6.71*q)**4,.25)


  # H(z)
  def Hz(self,z):
    return self.H0*np.power(self.Om*(1+z)**3+self.OK*(1+z)**2+self.Ov,.5)


  # H(z)/H0
  def expansion(self,z):
    return np.power(self.Om*(1+z)**3+self.OK*(1+z)**2+self.Ov,.5)


  # growth factor
  def growth_D(self,z):
    I = lambda x: (1.+x)/cosmofunc.expansion(self,x)**3
    return 2.5*self.Om*cosmofunc.expansion(self,z) * quad(I,z,100.)[0]


  # comoving distance
  def chi(self,z):
    I = lambda x: 1./(self.H0*cosmofunc.expansion(self,x))
    return quad(I,0.,z)[0]


  # Matter power spectrum
  def Pmk(self,k,z):
    fc = 2*np.pi**2/np.power(self.H0,3+self.ns)
    return (8.*np.pi**2/(25.*self.OmH2**2))*self.As*k*np.power(k/self.k0,self.ns-1.)*cosmofunc.tk_bbks(self,k)**2*cosmofunc.growth_D(self,z)**2



# reionization

def xe_sym(z,zre=8.,Dz=4.):
  y   = np.power(1.+z,1.5)
  yre = np.power(1.+zre,1.5)
  Dy  = 1.5*np.sqrt(1.+zre)*Dz
  return .5*(1.-np.tanh((y-yre)/Dy))


def xe_asym(z,alpha=6.,zend=6.):
  if z<zend:
    return 1.
  if z>=zend and z<20.:
    return np.power((20.-z)/(20.-zend),alpha)
  if z>=20:
    return 1e-30


