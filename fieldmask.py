import numpy as np
import healpy as hp

class FieldMask(object):

  def apfac(self):
    npix = self.msk.size
    wsum = np.sum(np.square(self.msk))
    return npix/wsum

  # Functions relating to Xspect / Xpol
  # See Tristram et al.
  # MNRAS 358 (2005) 833
  # http://arxiv.org/abs/astro-ph/0405575

  # w2 is the inverse of apfac from mask_apfac
  def weights(self):
    npix = self.msk.size
    w1 = np.sum(self.msk) / npix
    w2 = np.sum(np.square(self.msk)) / npix
    w4 = np.sum(self.msk ** 4) / npix
    return w1,w2,w4

  def apply_ap(self,map):
    n = map.map.shape
    apmap = 0.0 * np.array(map.map)
    if len(n)==1:
      mu = np.mean(map.map[np.where(self.msk>0)])
      apmap = (map.map - mu) * self.msk
    else:
      for i in xrange(n[0]):
        mu = np.mean(map.map[i,np.where(self.msk>0)])
        apmap[i,:] = (map.map[i,:] - mu) * self.msk
    return apmap

  def __init__(self,lon,lat,nside=512,nest=False):
    self.lon = lon
    self.lat = lat
    self.nside = nside
    self.ap = ''
    self.nest = nest
    npix = 12 * self.nside * self.nside
    self.msk = np.ones(npix)
    return
 
class CircMask(FieldMask):
  def _make_circle_r(self):
    # Calculate radius of each pixel center in field-centered coordinates
    npix = 12 * self.nside * self.nside
    thph=hp.pix2ang(self.nside,range(0,npix),nest=self.nest)
    ddeg=180.0/np.pi*hp.rotator.angdist(thph,[(90-self.lat)*np.pi/180,self.lon*np.pi/180],lonlat=False)
    return ddeg

  def __init__(self,lon,lat,rad,nside=512,nest=False):
    self.lon = lon
    self.lat = lat
    self.rad = rad
    self.nside = nside
    self.nest = nest
    self.ap = 'hard'
    r = self._make_circle_r()
    self.msk = np.array (0.0 + r <= rad)
    return

class CosMask(CircMask):
  def __init__(self,lon,lat,rad,nside=512,nest=False,n=2,rtuk=10.0):
    self.lon = lon
    self.lat = lat
    self.rad = rad
    self.nside = nside
    self.nest = nest
    self.ap = 'cos'
    self.n = n
    self.rtuk = rtuk
    r = self._make_circle_r()
    self.msk = np.array (0.0 + r <= rad)
    x = (rad - r) / (0.0 + rad*rtuk)
    self.msk = np.select([(x>=0) & (x<=1)],[0.5*(1-np.cos(np.pi*x))],[self.msk]).flatten()
    return 

class GaussMask(CircMask):
  # Not 100% sure what they mean by Gaussian apodization...
  def __init__(self,lon,lat,rad,nside=512,nest=False,fwhm=2.0):
    self.lon = lon
    self.lat = lat
    self.rad = rad
    self.nside = nside
    self.nest = nest
    self.ap = 'gauss'
    self.fwhm = 2.0
    r = self._make_circle_r()
    x = (rad - r) / (0.0 + fwhm) *  2.355
    self.msk = np.select([x>=0],[1.0-np.exp(-0.5 * np.square(x))],0.0).flatten()
    return

