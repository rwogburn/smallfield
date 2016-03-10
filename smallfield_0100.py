#!/bin/env python

import smallfield
import healmap
import healpy as hp
import numpy as np

if __name__ == "__main__":
  use_split = 'M'
  use_nside = 512
  out_nside = 32
  outfname = 'smallfield_hm_353_circ_dl_0100.fits'

  h1,h2 = healmap.get_input_maps(use_split,use_nside=use_nside,dir='/n/bicepfs3/general/input_maps/planck_pr2/')
  lat,lon,ipix = smallfield.gen_field_centers(nside=out_nside,minlat=-1.)
  i0 = 0
  nbins = len(smallfield.get_default_bins())-1
  dl = ()
  for j in xrange(2*3*nbins):
    dl = dl + (np.zeros([12*out_nside*out_nside]),)
  if i0>0:
    for j in xrange(3*nbins):
      dl[j] = hp.read_map(outfname,field=j)

  for i in xrange(i0,lat.size):
    print str(i)+": ipix="+str(ipix[i])+", lat="+str(lat[i])+", lon="+str(lon[i])
    # set dlat for racetrack mask
    tmpdl,var,lc = smallfield.patch_dl_xpol(lon[i],lat[i],h1,map2=h2,nside=use_nside)
    for j in xrange(3):
      for k in xrange(nbins):
        dl[j*nbins+k][ipix[i]] = tmpdl[k,j]
        dl[(j+3)*nbins+k][ipix[i]] = var[k,j]
    if ((i % 10)==0):
      hp.write_map(outfname,dl,nest=False,coord='G')
  hp.write_map(outfname,dl,nest=False,coord='G')
