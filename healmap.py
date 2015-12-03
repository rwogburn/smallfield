import healpy as hp
import numpy as np
from os import path

class HealMap(object):
    
  def _get_a_map(self,field=0):
    print "Loading healpix map field="+str(field)+"."
    h = hp.read_map(path.join(self.dir,self.fname),field)
    if self.nside != []:
      if (hp.get_map_size(h) > hp.nside2npix(self.nside)):
        print "Downgrading."
        h = hp.ud_grade(h,nside_out=self.nside)
    print "Scaling from K_CMB to uK_CMB."
    h = h * self.cal
    return h

  def _get_a_map_tqu(self):
    h = self._get_a_map(field=0)
    try:
      hq = self._get_a_map(field=1)
      hu = self._get_a_map(field=2)
      # h = np.concatenate((h,hq,hu),axis=1)
      h = np.vstack((h,hq,hu))
    except ValueError:
      print "No polarization information in fits map."
    print h.shape
    return h

  def __init__(self,fname,dir='.',nside=[],cal=1.0):
    self.dir=dir
    self.fname=fname
    self.nside=nside
    self.cal=cal
    self.map = self._get_a_map_tqu()
    if self.nside == []:
      self.nside = hp.get_nside(self.map)
    return

def get_input_maps(split='D',use_nside=512,use_cal=1.0,dir='.'):
  if split == 'D':
    print "Using det set split"
    h1 = HealMap('HFI_SkyMap_353-ds1_2048_R2.02_full.fits',nside=use_nside,cal=1.0e6,dir=dir)
    h2 = HealMap('HFI_SkyMap_353-ds2_2048_R2.02_full.fits',nside=use_nside,cal=1.0e6,dir=dir)
  elif split == 'Y':
    print "Using year split"
    h1 = HealMap('HFI_SkyMap_353_2048_R2.02_year-1.fits',nside=use_nside,cal=1.0e6,dir=dir)
    h2 = HealMap('HFI_SkyMap_353_2048_R2.02_year-2.fits',nside=use_nside,cal=1.0e6,dir=dir)
  elif split == 'R':
    print "Using half-ring split"
    h1 = HealMap('HFI_SkyMap_353_2048_R2.02_full-ringhalf-1.fits',nside=use_nside,cal=1.0e6,dir=dir)
    h2 = HealMap('HFI_SkyMap_353_2048_R2.02_full-ringhalf-2.fits',nside=use_nside,cal=1.0e6,dir=dir)
  elif split == 'M':
    print "Using half-mission split"
    h1 = HealMap('HFI_SkyMap_353_2048_R2.02_halfmission-1.fits',nside=use_nside,cal=1.0e6,dir=dir)
    h2 = HealMap('HFI_SkyMap_353_2048_R2.02_halfmission-2.fits',nside=use_nside,cal=1.0e6,dir=dir)
  elif split == 'F':
    print "Using 353/30 frequency cross"
    h1 = HealMap('HFI_SkyMap_353_2048_R2.02_full.fits',nside=use_nside,cal=1.0e6,dir=dir)
    h2 = HealMap('HFI_SkyMap_30-1024_R2.01_full.fits',nside=use_nside,cal=1.0e6,dir=dir)
  else:
    raise ValueError("Unrecognized split code "+split+".")
  print "Maps are ready."
  return h1,h2
