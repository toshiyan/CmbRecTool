import numpy as np

# * aps amplitude histogram
class statistics:

  def __init__(self,ocl='',scl=''):
    #   ocl --- observed cl array like [bin]
    #   scl --- simulated cl array like [sim,bin]
    self.ocl = ocl
    self.scl = scl
    self.A  = 0.
    self.mA = 0.
    self.sA = 0.
    self.MA = 0.
    self.p  = 0.
    self.px2 = 0.
    self.ox2 = 0.
    self.sx2 = 0.
    self.px1 = 0.
    self.ox1 = 0.
    self.sx1 = 0.

    self.onlydiag = False


  def x2PTE(self,simamp=1,diag=False):
    # compute chi^2 PTE of ocl using scl
    # simamp /=1 if needed
    n   = len(self.scl[:,0])
    # for real data
    mx  = np.mean(self.scl,axis=0)
    dxi = self.scl - mx
    dx0 = self.ocl - mx*simamp
    cov = np.cov(dxi,rowvar=0)
    if diag: cov = np.diag(np.diag(cov))
    oX2 = np.dot(dx0,np.dot(np.linalg.inv(cov),dx0))
    # for sim (exclude self rlz)
    dxi = np.array([self.scl[i,:]-np.mean(np.delete(self.scl,i,0),axis=0) for i in range(n)])
    sX2 = np.array([np.dot(dxi[i,:],np.dot(np.linalg.inv(np.cov(np.delete(dxi,i,0),rowvar=0)),dxi[i,:])) for i in range(n)])
    # output
    self.px2 = (sX2>oX2).sum()/np.float(n)
    self.ox2 = oX2
    self.sx2 = sX2


  def x1PTE(self,simamp=1,twoside=True):
    # compute chi PTE of ocl using scl
    n   = len(self.scl[:,0])
    # for real data
    mx  = np.mean(self.scl,axis=0)
    sx  = np.std(self.scl,axis=0)
    oX1 = np.sum((self.ocl-mx*simamp)/sx)
    # for sim (exclude self rlz)
    dxi = np.array([self.scl[i,:]-np.mean(np.delete(self.scl,i,0),axis=0) for i in range(n)])
    sX1 = np.array([np.sum(dxi[i,:]/np.std(np.delete(self.scl,i,0),axis=0)) for i in range(n)])
    # output
    px1 = (sX1>oX1).sum()/np.float(n)
    if twoside:
      self.px1 = 1.-2*np.abs(px1-0.5)
    else:
      self.px1 = px1
    self.ox1 = oX1
    self.sx1 = sX1
    #print np.std(sX1), oX1-np.mean(sX1)


  def get_amp(self,fcl='',scale=1.,diag=False):
    # estimate the amplitude of the power spectrum
    #   fcl --- fiducial cl used to define amplitude array like [bin]

    # mean
    if fcl=='': fcl = np.mean(self.scl,axis=0)*scale
    amp = self.scl/fcl

    # covariance
    cov = np.cov(amp,rowvar=0)
    cov[np.isnan(cov)] = 0.
    if diag: cov = np.diag(np.diag(cov))
    wb = np.sum(np.linalg.inv(cov),axis=0)
    wt = np.sum(wb)

    # amplitude estimator
    A  = np.sum(wb*amp,axis=1)/wt
    mA = np.mean(A)
    sA = np.sqrt(np.var(A))

    # observed amplitude
    oA = np.sum(wb*self.ocl/fcl)/wt
    #print 'obs/sim amp', oA, mA, 'sigma amp', sA, 'ratio', oA/sA

    self.A  = A
    self.mA = mA
    self.sA = sA
    self.oA = oA
    self.p  = (A>oA).sum()/np.float(len(A))
    self.MA = np.median(A)


  def get_amp_simfix(self,ocls):
    # estimating bias in the amplitude of the power spectrum with fixed cl and covariance
    # scl is a simulated covariance and fiducial cl
    #   ocls --- mock observed cls array like [sim,bin]

    # fiducial spectrum
    fcl = np.mean(self.scl,axis=0)

    # amplitude parameters for each mock observed and simulated cl
    ampo = ocls/fcl
    amps = self.scl/fcl

    # optimal weighting evaluated from simulation
    wb, wt = opt_weight(amps,self.onlydiag)

    # amplitude estimator
    A  = np.sum(wb*ampo,axis=1)/wt
    mA = np.mean(A)
    sA = np.sqrt(np.var(A))

    print 'mean =', mA, ', sigma =', sA, ', S/N =', np.sqrt(wt)

    return mA, sA


  def plot_hist(self,l=6,histbn=20,showfig=False,f=''):
    import matplotlib.pyplot as plt
    #plt.xlim(xran)
    #plt.ylim(yran)
    plt.xlabel(r'Power spectrum amplitude')
    plt.ylabel(r'Probability distribution')
    plt.hist(self.A,bins=histbn,normed=1,weights=np.ones_like(self.A)/float(len(self.A)),alpha=.5,lw=0)
    plt.figtext(0.60,0.80,r'mean ='+str(self.mA)[:l])
    plt.figtext(0.60,0.75,r'median ='+str(self.MA)[:l])
    plt.figtext(0.60,0.70,r'$\sigma(A)$ ='+str(self.sA)[:l])
    plt.figtext(0.60,0.65,r'$A^{obs}$ ='+str(self.oA)[:l])
    plt.figtext(0.60,0.60,r'$P(A>A^{obs})$ ='+str(self.p)[:l])
    plt.axvline(self.oA,color='r',lw=2)
    #plot(x,np.exp(-(x-w.mA)**2/(2.*w.sA**2))/np.sqrt(2.*np.pi)/w.sA,label=r'normal dist.')
    #axvline(w.mA-2*w.sA,ls='--',color='k')
    #axvline(w.mA-1*w.sA,ls='--',color='k')
    #axvline(w.mA,ls='--',color='k')
    if f!='': plt.savefig(f,bbox_inches='tight')
    if showfig:  plt.show()
    plt.clf()



#////////// data analysis functions //////////#
def quickplot_hist(dat,obs='',bins=20,label='',histtype='step',density=True,gpdf=True,err_type=''):
  import matplotlib.pyplot as plt
  import inspect
  from scipy.stats import norm

  plt.xlim(min(dat),max(dat))
  h0, b, p = plt.hist(dat,bins,histtype=histtype,label=label,density=density)

  h1, b = np.histogram(dat,density=False)
  h2, b = np.histogram(dat,density=True)
  # scaling
  if density:
    R = 1.
  else:
    R = h1[0]/h2[0]

  X = (b[1:]+b[:-1])/2
  x   = np.arange(min(dat),max(dat),(max(dat)-min(dat))/100.)
  mu  = np.mean(dat)
  sig = np.std(dat)
  y   = norm.pdf(x,loc=mu,scale=sig)*R
  Y   = norm.pdf(X,loc=mu,scale=sig)*R

  # plot Gaussian PDF
  if gpdf:
    plt.plot(x,y)

  # plot error bars of histogram
  if err_type!='': 

    if err_type=='gaussian': 
      yerr = np.sqrt(Y)

    if err_type=='hist':
      yerr = np.sqrt(h1)

    if density: yerr = yerr/R
    plt.errorbar(X,h0,yerr=yerr,fmt='o')

    # chi-square to gaussian pdf
    chi2 = np.sqrt(np.sum((h0-Y)**2/yerr**2))
    plt.figtext(.8,.8,'$\chi^2=$'+str(chi2)[:4])

  # observed value
  if obs!='': 
    plt.axvline(obs)

  # add legend
  plt.legend(loc=0,frameon=False)




def hist_errorbars( data, xerrs=True, *args, **kwargs) :
    """Plot a histogram with error bars. Accepts any kwarg accepted by either numpy.histogram or pyplot.errorbar"""
    import matplotlib.pyplot as plt
    import inspect

    # pop off normed kwarg, since we want to handle it specially
    norm = False
    if 'normed' in kwargs.keys() :
        norm = kwargs.pop('normed')

    # retrieve the kwargs for numpy.histogram
    histkwargs = {}
    for key, value in kwargs.iteritems() :
        if key in inspect.getargspec(np.histogram).args :
            histkwargs[key] = value

    histvals, binedges = np.histogram( data, **histkwargs )
    a, binedges = np.histogram( data)
    yerrs = np.sqrt(a)*histvals[0]/a[0]
    print yerrs, histvals

    if norm :
        nevents = float(sum(histvals))
        binwidth = (binedges[1]-binedges[0])
        histvals = histvals/nevents/binwidth
        yerrs = yerrs/nevents/binwidth

    bincenters = (binedges[1:]+binedges[:-1])/2

    if xerrs :
        xerrs = (binedges[1]-binedges[0])/2
    else :
        xerrs = None

    # retrieve the kwargs for errorbar
    ebkwargs = {}
    for key, value in kwargs.iteritems() :
        if key in inspect.getargspec(plt.errorbar).args :
            histkwargs[key] = value
    out = plt.errorbar(bincenters, histvals, yerrs, xerrs, fmt=".", **ebkwargs)

    if 'log' in kwargs.keys() :
        if kwargs['log'] :
            plt.yscale('log')

    if 'range' in kwargs.keys() :
        plt.xlim(*kwargs['range'])

    return out


def apofunc(distance,aposcale):
  # apodization window
  #   distance --- unit in rad
  #   aposcale --- scale of apodization in deg
  x = (1.-np.cos(distance))/(1.-np.cos(aposcale*np.pi/180.))
  x[x>=1.] = 1.
  y = np.sqrt(x)
  return  y - np.sin(2*np.pi*y)/(2*np.pi)


def opt_weight(x,diag=False):
  # optimal weighting
  #   x --- data like [sim,bin]
  if x.shape[1]>=2: 
    cov  = np.cov(x,rowvar=0)
    if diag: cov = np.diag(np.diag(cov)) # set off-diag to zero
    cov[np.isnan(cov)] = 0.
    cinv = np.linalg.inv(cov)   # inverse covariance
  else:
    cinv = 1./np.var(x)
  wb   = np.sum(cinv,axis=0)  # optimal weight for each x
  wt   = np.sum(wb)           # normalization
  return wb, wt


#////////// Cl binning //////////#

def cl_binning(OL,WL,bp):
  # binning of power spectrum
  bn = np.size(bp) - 1
  bc = (bp[1:]+bp[:-1])*.5
  cb = np.zeros(bn)
  for i in range(bn):
    b0 = int(bp[i])
    b1 = int(bp[i+1])
    cl = OL[b0:b1]
    wl = 1./WL[b0:b1]**2
    N0 = np.count_nonzero(wl)
    if N0==b1-b0:
      norm = 1./np.sum(wl)
      cb[i] = norm*np.sum(wl*cl)
    elif N0<b1-b0 and N0>0:
      norm = 1./np.sum(wl[wl!=0])
      cb[i] = norm*np.sum(wl[wl!=0]*cl[wl!=0])
    else:
      cb[i] = 0
  return bc, cb

