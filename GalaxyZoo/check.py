import seaborn as sns
import pandas as pd
import termcolor
import numpy as np
from astropy.coordinates import SkyCoord as sky
import astropy.units as unit


def feature_farction(feature,sample):
    sep = 0.01
    for z in np.arange(0, 0.2, sep):
        whole = sample[(sample.z <= z + sep) & (sample.z > z)]
        whole = whole.sort_values(by='Mr')
        whole.index = range(len(whole))
        pass
        # for m in np.arange(-19,-21,-0.5)


zoo = pd.read_csv('zoo2MainSpecz.csv', header=0)
uni = pd.read_csv('gz2sample.csv', header=0, index_col='OBJID')
zoo['Mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_MR.values
zoo['mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_R.values
zoo['z'] = uni.ix[zoo.dr7objid.values].REDSHIFT.values
# z, Mr, mr
gal = pd.DataFrame(columns=['z','Mr','fraction','feature'])
gal['z'] = np.repeat(np.arange(0,0.2,0.01),4)
gal['Mr'] = np.tile(np.arange(-19,-21,-0.5),20)
gal['featrue'] = np.tile(['bar','bulge','odd','spiral'],20)
print(gal)