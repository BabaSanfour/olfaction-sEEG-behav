# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../behavior') 
from connection import *

import numpy as np



subject = session.query(Subject).filter_by(name = 'VACJ').first()

anasig = subject.blocks[0].segments[0].analogsignals[0]
sig = anasig.signal.magnitude
print anasig.signal.units

print np.max(sig)
print np.min(sig)

bit_depth1 =  np.min(np.diff(np.sort(np.unique(sig))))
print bit_depth1
all = np.diff(np.sort(sig))
bit_depth2 = np.min(all[all!=0])
print bit_depth2

print 2**16 * bit_depth2