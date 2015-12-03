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

class RacetrackMask(FieldMask):
  def _make_racetrack_r(self):
    # Need to build a pseudo-"radius" for the racetrack shape
    # Take the minimum of these three values:
    #   - distance from left center point
    #   - distance from right center point
    #   - distance in declination from central declination
    npix = 12 * self.nside * self.nside
    thph = hp.pix2ang(self.nside,range(0,npix),nest=self.nest)
    if self.galmap:
      # Convert to equatorial coordinates to do the calc
      r = hp.rotator.Rotator(coord='GC')
      thph = r(thph)
      lon0,lat0 = r(self.lon,self.lat,lonlat=True)
    else:
      lon0 = self.lon
      lat0 = self.lat

    # ra = thph(1)
    # de = np.pi/2 - thph(0)

    # The logic here has to get a little complicated.  We want to make a racetrack with fixed width
    # (that is, declination range) and with circular ends, adjusting the length (az range between
    # the centers of the two circles) to keep the area constant.  This may become impossible near the
    # poles.
    a_circ = 2.*np.pi * (1. - np.cos(self.dlat/2. * np.pi/180.))
    a_annulus = 2.*np.pi * (np.sin((lat0 + self.dlat/2.) * np.pi/180.) - np.sin((lat0 - self.dlat/2.) * np.pi/180.)) 
    frac = (self.a_equiv - a_circ) / a_annulus
    if (frac<0) or (frac>1.):
      raise ValueError('Impossible racetrack shape requested.')
    dra = 360. * frac
    ddeg1 = 180./np.pi * hp.rotator.angdist(thph,[(90.-lat0)*np.pi/180.,(lon0+dra/2.)*np.pi/180.],lonlat=False)
    ddeg2 = 180./np.pi * hp.rotator.angdist(thph,[(90.-lat0)*np.pi/180.,(lon0-dra/2.)*np.pi/180.],lonlat=False)
    ddeg3 = np.abs(lat0 - 180./np.pi*(np.pi/2. - thph[0])) + 180. * (np.abs(np.mod((lon0 - 180/np.pi*thph[1]) + 180.,360.)-180.) >= dra/2.)
    ddeg = np.minimum(np.minimum(ddeg1,ddeg2),ddeg3)

    # If we're close to a pole and the racetrack isn't big enough, fatten it out with a circle
    a_actual = np.mean(ddeg <= self.dlat/2.) * 4.*np.pi
    if (a_actual < self.a_equiv): 
      print a_actual, self.a_equiv
      ddeg4 = 180./np.pi * hp.rotator.angdist(thph,[(90.-lat0)*np.pi/180.,lon0*np.pi/180.],lonlat=False)
      rpad = dra/2.
      while (a_actual < self.a_equiv):
        ddeg = np.minimum(ddeg, ddeg4 - rpad + dra/2.)
        a_actual = np.mean(ddeg <= self.dlat/2.) * 4.*np.pi
        rpad = rpad + 0.01

    return ddeg

  def __init__(self,lon,lat,dlat,rad_equiv=11.3,a_equiv=None,nside=512,nest=False,galmap=True):
    self.lon = lon
    self.lat = lat
    self.dlat = dlat
    self.nside = nside
    self.nest = nest
    self.galmap = galmap
    self.ap = 'hard'
    if a_equiv is None:
      self.a_equiv = 2*np.pi * (1 - np.cos(rad_equiv * np.pi/180.))
    else:
      self.a_equiv = a_equiv
    r = self._make_racetrack_r()
    self.msk = np.array (0.0 + r <= self.dlat/2.)
    return

class CosRacetrackMask(RacetrackMask):
  def __init__(self,lon,lat,dlat,rad_equiv=11.3,a_equiv=None,nside=512,nest=False,n=2,rtuk=10.0,galmap=True):
    self.lon = lon
    self.lat = lat
    self.dlat = dlat
    self.nside = nside
    self.galmap = galmap
    self.nest = nest
    self.ap = 'cos'
    self.n = n
    self.rtuk = rtuk
    if a_equiv is None:
      self.a_equiv = 2*np.pi * (1 - np.cos(rad_equiv * np.pi/180.))
    else:
      self.a_equiv = a_equiv
    r = self._make_racetrack_r()
    x = (self.dlat/2. - r) / (0.0 + self.dlat/2.*rtuk)
    self.msk = np.select([(x>=0) & (x<=1)],[0.5*(1-np.cos(np.pi*x))],[self.msk]).flatten()
    return

class GaussRacetrackMask(RacetrackMask):
  # Not 100% sure what they mean by Gaussian apodization...
  def __init__(self,lon,lat,dlat,rad_equiv=11.3,a_equiv=None,nside=512,nest=False,fwhm=2.0,galmap=True):
    self.lon = lon
    self.lat = lat
    self.dlat = dlat
    self.nside = nside
    self.galmap = galmap
    self.nest = nest
    self.ap = 'gauss'
    self.fwhm = 2.0
    if a_equiv is None:
      self.a_equiv = 2*np.pi * (1 - np.cos(rad_equiv * np.pi/180.))
    else:
      self.a_equiv = a_equiv
    r = self._make_racetrack_r()
    x = (self.dlat/2. - r) / (0.0 + fwhm) * 2.355
    self.msk = np.select([x>=0],[1.0-np.exp(-0.5 * np.square(x))],0.0).flatten()
    return

