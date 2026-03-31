#!/home/charmm-gui/local/bin/python
import sys, os
import cPickle
import numpy as np
from sklearn.neighbors import KernelDensity

segid = sys.argv[1]
#segid = segid.lower()
#saved_pickle = 'glycan_rotlib.dat'
saved_pickle = sys.argv[2]

kde = cPickle.load( open(saved_pickle, 'rb') )
sample = kde[segid].sample(n_samples=5000)

np.savetxt('./phi_'+segid.lower()+'.dat', sample[:,0], fmt='%1.4f')
np.savetxt('./psi_'+segid.lower()+'.dat', sample[:,1], fmt='%1.4f')
