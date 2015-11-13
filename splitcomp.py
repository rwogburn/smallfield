#!/opt/local/bin/python2.7

import smallfield
import healmap
import healpy as hp
import numpy as np

if __name__ == "__main__":
  use_split = 'D'
  use_nside = 512

  h1,h2 = healmap.get_input_maps(use_split,use_nside=use_nside,dir='..')
  sfac = 0.0395
  squ = 10
  st = 200
  rot=[0,90,0]
  #
  hp.orthview(h1.map[0,:]*sfac,rot=rot,nest=False,title="Split 1 T",min=0,max=st)
  hp.orthview(h1.map[1,:]*sfac,rot=rot,nest=False,title="Split 1 Q",min=-squ,max=squ)
  hp.orthview(h1.map[2,:]*sfac,rot=rot,nest=False,title="Split 1 U",min=-squ,max=squ)
  #
  hp.orthview(h2.map[0,:]*sfac,rot=rot,nest=False,title="Split 2 T",min=0,max=st)
  hp.orthview(h2.map[1,:]*sfac,rot=rot,nest=False,title="Split 2 Q",min=-squ,max=squ)
  hp.orthview(h2.map[2,:]*sfac,rot=rot,nest=False,title="Split 2 U",min=-squ,max=squ)
  #
  plt.show()

