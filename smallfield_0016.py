#!/bin/env python

import smallfield
import healmap
import healpy as hp
import numpy as np

if __name__ == "__main__":
  use_split = 'D'
  use_nside = 512
  out_nside = 32
  outfname = 'smallfield_0016.fits'
  dl_lcdm,lc = smallfield.calc_lcdm_dl()
  print dl_lcdm

  h1,h2 = healmap.get_input_maps(use_split,use_nside=use_nside,dir='/n/bicepfs3/general/input_maps/planck_pr2/')
  lat,lon,ipix = smallfield.gen_field_centers(nside=out_nside,minlat=-1.0)
  i0 = 0
  if i0==0:
    rb = np.zeros([12*out_nside*out_nside])
    sb = 0.0 * rb
    re = 0.0 * rb
    se = 0.0 * rb
  else:
    rb = hp.read_map(outfname,field=0)
    sb = hp.read_map(outfname,field=1)
    re = hp.read_map(outfname,field=2)
    se = hp.read_map(outfname,field=3)

  for i in xrange(i0,lat.size):
    if (rb[ipix[i]]!=0.) or (sb[ipix[i]]!=0.):
      continue
    print str(i)+": ipix="+str(ipix[i])+", lat="+str(lat[i])+", lon="+str(lon[i])
    dl,var,lc = smallfield.patch_dl_xpol(lon[i],lat[i],h1,map2=h2,nside=use_nside)
    rb[ipix[i]],sb[ipix[i]] = smallfield.calc_r_equiv(lc,dl[:,2],var[:,2])
    re[ipix[i]],se[ipix[i]] = smallfield.calc_r_equiv(lc,dl[:,1],var[:,1],sub_lcdm=dl_lcdm[:,1])
    print "   r_B  "+str(rb[ipix[i]])+" +- "+str(sb[ipix[i]])+", r_E  "+str(re[ipix[i]])+" +- "+str(se[ipix[i]])
    if ((i % 10)==0):
      hp.write_map(outfname,(rb,sb,re,se),nest=False,coord='G')
  hp.write_map(outfname,(rb,sb,re,se),nest=False,coord='G')
