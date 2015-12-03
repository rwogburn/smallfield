#!/opt/local/bin/python2.7

import numpy as np
import healpy as hp
from matplotlib import pyplot as plt
from scipy import optimize
import fieldmask
import healmap

# default_bins = np.concatenate([[0],np.arange(20,(581+35*2),35)])
default_bins = np.array([40,70,110,160,220,290,370])

def bin_cl(clall,bins):
  try:
    n = len(clall)
  except TypeError:
    clall = [clall]
    n = 1 
  ell = np.arange(clall[0].size)
  ellfac = ell * (ell+1) / (2.0 * np.pi)
  dl = np.zeros([np.array(bins).size-1,6])
  for i in xrange(6):
    if n>i:
      tmpcl = clall[i]
    else:
      tmpcl = 0.0 * clall[0]
    tmpcl = tmpcl * ellfac
    for j in xrange(len(bins)-1):
      idx = np.where((ell>=bins[j]) & (ell<bins[j+1]))
      dl[j,i] = np.mean(tmpcl[idx])
  lc = 0.5 * (bins[:-1] + bins[1:])
  nell = (bins[1:] - bins[:-1])
  return dl,lc,nell

def patch_dl(lon,lat,map1,map2=None,bins=default_bins,nside=512,dlat=None):
  if dlat is None:
    msk = fieldmask.GaussMask(lon,lat,11.3,fwhm=2.0,nside=nside)
  else:
    msk = fieldmask.GaussRacetrackMask(lon,lat,dlat=dlat,rad_equiv=11.3,fwhm=2.0,nside=nside)
  apmap1 = msk.apply_ap(map1)
  if map2 is not None:
    apmap2 = msk.apply_ap(map2)
  else:
    apmap2 = None
  lmax = bins[-1]
  clall = hp.anafast(apmap1,map2=apmap2,lmax=lmax,pol=True)
  dl,ell,nell = bin_cl(clall,bins)
  dl = dl * msk.apfac()
  return dl,ell,nell

def calc_lcdm_dl(clfile='camb_66469116_scalcls.fits',cldir='.',bins=default_bins):
  cl = hp.read_cl(cldir+'/'+clfile)
  dl,ell,nell = bin_cl(cl,bins)
  return dl,ell

def patch_dl_xpol(lon,lat,map1,map2=None,bins=default_bins,nside=512,dlat=None):
  if dlat is None:
    msk = fieldmask.GaussMask(lon,lat,11.3,fwhm=2.0,nside=nside)
  else:
    msk = fieldmask.GaussRacetrackMask(lon,lat,dlat=dlat,rad_equiv=11.3,fwhm=2.0,nside=nside)
  w1,w2,w4 = msk.weights()
  apmap1 = msk.apply_ap(map1)
  if map2 is None:
    raise ValueError('Xpol requires cross spectra, not auto!')
  apmap2 = msk.apply_ap(map2)
  lmax = bins[-1]
  #
  # Need auto and cross spectra
  cl11 = hp.anafast(apmap1,lmax=lmax,pol=True)
  dl11 = bin_cl(cl11,bins)[0]
  cl22 = hp.anafast(apmap2,lmax=lmax,pol=True)
  dl22 = bin_cl(cl22,bins)[0]
  cl12 = hp.anafast(apmap1,map2=apmap2,lmax=lmax,pol=True)
  dl12,ell,nell = bin_cl(cl12,bins)
  #
  # Sky area correction factors
  dl11 = dl11 / w2
  dl22 = dl22 / w2
  dl12 = dl12 / w2
  #
  # Xpol recipe: use cross spectrum for angular power spectra
  dl = dl12
  # Xpol recipe: analytic expression for error bars
  # Tristram et al. Eq. 29
  nu_l = (2*ell+1.0)*nell*np.square(w2)/w4
  # Tristram et al. Eq. 32 (OK to work in D_l instead of C_l because all terms quadratic in C_l)
  var = (1.0/nu_l) * (np.square(dl12) + dl11*dl22)
  # Paper XXX Eq. 2 for sample variance (Try f_eff = w2^2/w4)
  # sfac = (2*ell+1.0)*nell*w1
  sfac = nu_l
  svar = (2.0/sfac) * np.square(dl12)
  #
  # Keep variances only for TT, EE, BB, not for TE, TB, EB where
  # these analytic formulae don't make sense.
  var[:,3:] = 0.0
  svar[:,3:] = 0.0
  # Paper XXX Sec. 3.1.1 says sample variance subtracted from Xpol
  # analytic expression for variance, so do that.
  dvar = np.square(np.sqrt(var) - np.sqrt(svar))
  if np.any(var<0) or np.any(svar<0):
    print "Negative value found in variance."
    for i in xrange(var.size):
      print "bin="+str(i)+", var="+str(var[i])+", svar="+str(svar[i])+", dvar="+str(dvar[i])
  #
  return dl,dvar,ell

def fit_powerlaw(lc,dl,var,alpha=-2.42,ell0=80.0):
  fitfunc = lambda x, a: a * (x/ell0) ** (alpha+2)
  if np.any(var<=0):
    print "At least one D_l variance is zero or negative."
  sigma = np.sqrt(var)
  out,cov = optimize.curve_fit(fitfunc,lc,dl,p0=[1.0],sigma=sigma)
  return out[0],np.sqrt(cov[0,0])

def calc_r_equiv(lc,dl,var,sub_lcdm=None,alpha=-2.42):
  if sub_lcdm is not None:
    if np.all(dl.shape == sub_lcdm.shape):
      dl_lcdm = sub_lcdm
    else:
      dl_lcdm = calc_lcdm_dl()
    dl = dl - dl_lcdm
  x,sig = fit_powerlaw(lc,dl,var,alpha,ell0=80)
  f = np.square(0.0395)        # Scaling 353 to 150 GHz uK_CMB, so square it in D_l
  f /= 6.71e-2                 # D_l at ell=80 for r=1, 6.71e-2 uK_CMB^2
  u = (0.004 / 0.0395)         # Frac uncertainty on 353 to 150 scaling from Paper XXX Sec. 5.3)
  sig = np.sqrt(np.square(sig) + np.square(2.0*x*u))
  return x*f,sig*f

def gen_field_centers(nside=8,minlat=35):
  npix = 12 * nside * nside
  th,ph=hp.pix2ang(nside,range(0,npix),nest=False)
  lat = np.array(90.0 - th*180.0/np.pi)
  lon = np.array(ph*180.0/np.pi)
  ipix = np.arange(0,npix)
  idx = np.where(np.abs(lat)>=minlat)
  lat = lat[idx]
  lon = lon[idx]
  ipix = ipix[idx]
  return lat,lon,ipix
