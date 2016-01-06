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

  def tqu_mean(self,map):
    n = map.map.shape
    if len(n)==1:
      mu = [np.mean(map.map[np.where(self.msk>0)])]
    else:
      mu = np.zeros(n[0])
      for i in xrange(n[0]):
        mu[i] = np.mean(map.map[i,np.where(self.msk>0)])
    return mu

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

  def boundary(self,reso=0.1):
    lon = []
    lat = []
    return lon,lat

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

  def boundary(self,reso=0.1):
    r = hp.Rotator(rot=[self.lon+180.,90.-self.lat,0],eulertype='XYZ',inv=True)
    ph0 = np.arange(0,361,reso/max([np.sin(1.*np.pi/180.),np.sin(self.rad*np.pi/180.)]))*np.pi/180.
    th,ph = r(self.rad*np.pi/180. + 0*ph0,ph0)
    lon = ph * 180./np.pi
    lat = 90. - th * 180./np.pi
    return lon,lat

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
    if (frac<0) or (frac>1):
      print "Unusual frac value outside 0-1 range."
      print "a_circ="+str(a_circ)+", a_annulus="+str(a_annulus)+", a_equiv="+str(self.a_equiv)+", frac="+str(frac)
    if (frac<0):
      raise ValueError('Impossible racetrack shape requested.')
    elif (frac<=1):
      if (frac<=0.5):
        dra = 360. * frac
      else:
        dra = 180.
      ddeg1 = 180./np.pi * hp.rotator.angdist(thph,[(90.-lat0)*np.pi/180.,(lon0+dra/2.)*np.pi/180.],lonlat=False)
      ddeg2 = 180./np.pi * hp.rotator.angdist(thph,[(90.-lat0)*np.pi/180.,(lon0-dra/2.)*np.pi/180.],lonlat=False)
      ddeg3 = np.abs(lat0 - 180./np.pi*(np.pi/2. - thph[0])) + 180. * (np.abs(np.mod((lon0 - 180/np.pi*thph[1]) + 180.,360.)-180.) >= dra/2.)
      ddeg = np.minimum(np.minimum(ddeg1,ddeg2),ddeg3)
    else:
      print "No racetrack, using a spot"
      dra = 0.
      ddeg = 180. + 0. * thph[0]

    # If we're close to a pole and the racetrack isn't big enough, fatten it out with a circle
    a_actual = np.mean(ddeg <= self.dlat/2.) * 4.*np.pi
    rpad = 0.
    if (a_actual < self.a_equiv): 
      print a_actual, self.a_equiv
      ddeg4 = 180./np.pi * hp.rotator.angdist(thph,[(90.-lat0)*np.pi/180.,lon0*np.pi/180.],lonlat=False)
      rpad = self.dlat/2.
      while (a_actual < self.a_equiv):
        ddeg = np.minimum(ddeg, ddeg4 - rpad + self.dlat/2.)
        a_actual = np.mean(ddeg <= self.dlat/2.) * 4.*np.pi
        rpad = rpad + 0.01

    # Store shape info
    self.dra = dra
    self.rpad = rpad
    return ddeg

  def boundary(self,reso=0.1):
    lon = np.array([])
    lat = np.array([])

    if self.galmap:
      # Convert to equatorial coordinates to do the calc
      r = hp.rotator.Rotator(coord='GC')
      lon0,lat0 = r(self.lon,self.lat,lonlat=True)
    else:
      lon0 = self.lon
      lat0 = self.lat

    if (self.dra>0):
      tmplon = np.arange(-self.dra/2.,self.dra/2.,reso/max([np.sin(1.*np.pi/180.),np.cos(lat0*np.pi/180.)]))
      if (lat0+self.dlat/2. > -90.) and (lat0+self.dlat/2. < 90.):
        lon = np.concatenate((lon,lon0+tmplon))
        lat = np.concatenate((lat,lat0+self.dlat/2. + 0.*tmplon))

      r = hp.Rotator(rot=[180.+lon0+self.dra/2.,90.-lat0,0],eulertype='XYZ',inv=True)
      ph0 = np.arange(360.,180.,-reso/max([np.sin(1.*np.pi/180.),np.sin(self.dlat/2.*np.pi/180.)]))*np.pi/180.
      th,ph = r(self.dlat/2.*np.pi/180. + 0*ph0,ph0)
      w = np.where((th>=0) & (th<=np.pi))
      lon = np.concatenate((lon,np.asarray(ph[w])*180./np.pi))
      lat = np.concatenate((lat,90.-np.asarray(th[w])*180./np.pi))

      if (lat0-self.dlat/2. > -90.) and (lat0-self.dlat/2. < 90.):
        lon = np.concatenate((lon,lon0-tmplon))
        lat = np.concatenate((lat,lat0-self.dlat/2. + 0.*tmplon))

      r = hp.Rotator(rot=[180.+lon0-self.dra/2.,90.-lat0,0],eulertype='XYZ',inv=True)
      ph0 = np.arange(180.,0.,-reso/max([np.sin(1.*np.pi/180.),np.sin(self.dlat/2.*np.pi/180.)]))*np.pi/180.
      th,ph = r(self.dlat/2.*np.pi/180. + 0*ph0,ph0)
      w = np.where((th>=0) & (th<=np.pi))
      lon = np.concatenate((lon,ph[w]*180./np.pi))
      lat = np.concatenate((lat,90.-th[w]*180./np.pi))

      lon = np.concatenate((lon,[lon[0]]))
      lat = np.concatenate((lat,[lat[0]]))
    else:
      r = hp.Rotator(rot=[lon0,lat0,0],eulertype='XYZ',inv=True)
      ph0 = np.arange(0.,360.,reso/max([np.sin(1.*np.pi/180.),np.sin(self.rpad*np.pi/180.)]))*np.pi/180.
      th,ph = r(self.rpad*np.pi/180. + 0*ph0,ph0)
      w = np.where((th>=0) and (th<=np.pi))
      lon = np.concatenate(lon,ph[w]*180./np.pi)
      lat = np.concatenate(lat,90.-th[w]*180./np.pi)

    if self.galmap:
      # Convert results from cel to gal if needed
      r = hp.Rotator(coord=['C','G'])
      lon,lat = r(lon,lat,lonlat=True)

    return lon,lat

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

