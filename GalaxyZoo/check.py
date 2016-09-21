import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import termcolor
import numpy as np
from astropy.coordinates import SkyCoord as sky
import astropy.units as unit


def selection(feature, sample):
    if feature == 'bar':
        return sample[(sample.t01_smooth_or_features_a02_features_or_disk_fraction > 0.227) &
                      (sample.t02_edgeon_a05_no_fraction > 0.519) &
                      (sample.t03_bar_a06_bar_count > 10) &
                      (sample.t03_bar_a06_bar_fraction > 0.3)]


def calculate_fraction():
    for z in z_range:
        whole = zoo[(zoo.z > z) & (zoo.z < z+0.01) & (zoo.Mr < -20.89) & (zoo.mr < 17)]
        bar = whole[(whole.t01_smooth_or_features_a02_features_or_disk_fraction > 0.430) &
                    (whole.t02_edgeon_a05_no_fraction > 0.715) &
                    (whole.t02_edgeon_a05_no_count >= 20)]
        result = bar[bar.t03_bar_a06_bar_fraction >= 0.8]
        rid = gal[(gal.Mr == -20.89) & (gal.z == z) & (gal.root == 'bar')].index[0]
        gal.at[rid,'fraction'] = len(result)/len(bar)
    gal.to_csv('gal.csv', index=False)


zoo = pd.read_csv('zoo2MainSpecz.csv', header=0)
uni = pd.read_csv('gz2sample.csv', header=0, index_col='OBJID')
zoo['Mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_MR.values
zoo['mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_R.values
zoo['z'] = uni.ix[zoo.dr7objid.values].REDSHIFT.values
# z, Mr, mr
gal = pd.DataFrame(columns=['z', 'Mr', 'fraction', 'feature', 'root'])

features = ['bar', 'odd', 'spiral']
z_range = np.arange(0.01, 0.2, 0.01)
nz = len(z_range)
mag_range = np.arange(-19.89, -21.89, -0.5)
nm = len(mag_range)
root_range = ['whole', 'bar']
nr = len(root_range)
gal['z'] = np.repeat(z_range, nm * nr)
gal['Mr'] = np.tile(np.repeat(mag_range, nr), np.array(nz))
gal['root'] = np.tile(root_range, np.array(nz * nm))

calculate_fraction()
plt.plot(z_range,gal[(gal.Mr == -20.89) & (gal.root == 'bar')]['fraction'].values)
plt.show()
