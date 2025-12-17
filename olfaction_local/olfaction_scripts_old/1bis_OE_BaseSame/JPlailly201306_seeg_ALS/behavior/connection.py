# -*- coding: utf-8 -*-
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import sys,os
import getpass

if getpass.getuser() in ('karim', 'samuel') and  sys.platform.startswith('linux'):
    p1 = r'/home/%s/Documents/projet/OpenElectrophy3/OpenElectrophy/'% getpass.getuser()
    p2 = '/home/karim/Documents/Backup_Dropbox/al_perso/Git_Geek/timefreqtools/'
elif sys.platform=='win32' :
    p1 = 'C:/Users/Anne-Lise/Anaconda3/envs/OE/Lib/site-packages/OpenElectrophy3'
    p2 = 'C:/Users/Anne-Lise/Anaconda3/envs/OE/Lib/site-packages/timefreqtools'
elif getpass.getuser() == 'david' and  sys.platform.startswith('linux'):
    p1 = r'/home/david/smb4k/CRNLDATA.UNIV-LYON1.FR/crnldata/cmo/scripts/OpenElectrophy3'
    p2 = r'/home/david/smb4k/CRNLDATA.UNIV-LYON1.FR/crnldata/cmo/scripts/timefreqtools_v0.1'
else :
    p1 = r'/home/%s/mnt/CRNLDATA/crnldata/cmo/scripts/OpenElectrophy3' % getpass.getuser()
    p2 = r'/home/%s/mnt/CRNLDATA/crnldata/cmo/scripts/timefreqtools_v0.1' % getpass.getuser()


sys.path = [p1,p2,] + sys.path

from distutils import version
import OpenElectrophy
assert version.LooseVersion(OpenElectrophy.__version__)>'0.3', 'Bad version of OpenElectrophy {}'.format( OpenEletrophy.__version__ )
from OpenElectrophy import *

from distutils import version
import timefreqtools
print timefreqtools.__file__
#print timefreqtools.__version__
assert version.LooseVersion(timefreqtools.__version__)<'0.2', 'Bad version of timefreqtools {}'.format(timefreqtools.__version__ )




sys.path.append('../seeg_analysis') 


import quantities as pq
#~ from datetime import datetime
import datetime
import numpy as np

from scipy import *
from scipy.signal import medfilt

medfilt_kernel_size = 5
max_percent_for_respiration_cycle = 0.6


from sqlalchemy.sql.expression import asc, desc, text


#This is OE0.3
from OpenElectrophy.core.sqlmapper import open_db, sql
from classes import classes

if getpass.getuser() == 'sgarcia' and  sys.platform.startswith('linux'):
    #~ url = 'sqlite:///../database JPlailly201306_seeg_ALS.sqlite'
    url = 'mysql+mysqldb://episodic:quiquandou@neuro001.univ-lyon1.fr/JPlailly201306_seeg_ALS'
elif getpass.getuser() == 'samuel' and  sys.platform.startswith('linux'):
    #~ url = 'sqlite:///../database JPlailly201306_seeg_ALS.sqlite'
    url = 'mysql+mysqldb://episodic:quiquandou@neuro001.univ-lyon1.fr/JPlailly201306_seeg_ALS'
elif getpass.getuser() == 'Sam' and  sys.platform.startswith('win'):
    #~ url = 'sqlite:///../database JPlailly201306_seeg_ALS.sqlite'
    url = 'mysql+mysqldb://episodic:quiquandou@neuro001.univ-lyon1.fr/JPlailly201306_seeg_ALS'
    
elif getpass.getuser() == 'karim' and  sys.platform.startswith('linux'):
    url = 'sqlite:///../database JPlailly201306_seeg_ALS.sqlite'
    #url = 'mysql+mysqldb://episodic:quiquandou@neuro001.univ-lyon1.fr/JPlailly201306_seeg_ALS'
    
elif getpass.getuser() == 'alsaive' and  sys.platform.startswith('win'):
    url = 'sqlite:///../database JPlailly201306_seeg_ALS.sqlite'
    #~ url = 'mysql+mysqldb://episodic:quiquandou@neuro001.univ-lyon1.fr/JPlailly201306_seeg_ALS'
    
elif getpass.getuser() == 'Anne-Lise' and  sys.platform.startswith('win'):
    #~ url = 'sqlite:///../database JPlailly201306_seeg_ALS.sqlite'
    url = 'mysql+mysqldb://episodic:quiquandou@neuro001.univ-lyon1.fr/JPlailly201306_seeg_ALS'

elif getpass.getuser() == 'david' and  sys.platform.startswith('linux'):
    #~ url = 'sqlite:///../database JPlailly201306_seeg_ALS.sqlite'
    url = 'mysql+mysqldb://episodic:quiquandou@neuro001.univ-lyon1.fr/JPlailly201306_seeg_ALS'
    
    #~ url = 'sqlite:///../database JPlailly201306_seeg_ALS.sqlite'

# print "#"*30
# print "Attention ça tappe sur la base", url[:5]
# print url
# print "#"*30

if sys.platform.startswith('linux') and url.startswith('sqlite'):
    #~ compress = 'lz4'
    compress = 'blosc'
else:
    compress = 'blosc'

dbinfo = open_db(url, myglobals = globals(), use_global_session = False, 
                        object_number_in_cache = 10000,  numpy_storage_engine = 'sqltable', compress = compress,
                        hdf5_filename = None,
                        relationship_lazy = 'dynamic', predefined_classes = classes)
session = dbinfo.Session()


#~ print session.query(AnalogSignal)[5].sampling_rate


