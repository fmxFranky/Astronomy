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
    elif feature == 'spiral':
        return sample[(sample.t01_smooth_or_features_a02_features_or_disk_fraction > 0.227) &
                      (sample.t02_edgeon_a05_no_fraction > 0.519) &
                      (sample.t04_spiral_a08_spiral_count > 10) &
                      (sample.t04_spiral_a08_spiral_fraction > 0.3)]
    elif feature == 'odd':
        return sample[(sample.t06_odd_a14_yes_count > 10) &
                      (sample.t06_odd_a14_yes_fraction > 0.3)]


def feature_farction(sample):
    sep = 0.01
    for z in np.arange(0, 0.2, sep):
        whole = sample[(sample.z <= z + sep) & (sample.z > z)]
        whole = whole.sort_values(by='Mr')
        whole.index = range(len(whole))
        complete_index = whole[whole.mr > 17].index[0] - 1
        print(complete_index)
        com_set = whole.ix[:complete_index]
        for fea in features:
            reset_index = gal[(gal.z == z) & (gal.feature == fea)].index[0]
            gal.at[reset_index, 'Mr'] = whole.ix[complete_index, 'Mr']
            gal.at[reset_index, 'fraction'] = len(selection(fea, com_set)) / len(com_set)
    print(gal)
    gal.to_csv('gal.csv', index=False)


zoo = pd.read_csv('zoo2MainSpecz.csv', header=0)
uni = pd.read_csv('gz2sample.csv', header=0, index_col='OBJID')
zoo['Mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_MR.values
zoo['mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_R.values
zoo['z'] = uni.ix[zoo.dr7objid.values].REDSHIFT.values
# z, Mr, mr
gal = pd.DataFrame(columns=['z', 'Mr', 'fraction', 'feature'])

features = ['bar', 'odd', 'spiral']

gal['z'] = np.repeat(np.arange(0, 0.2, 0.01), 3)
gal['feature'] = np.tile(features, 20)
# gal.index = gal.z

# print(gal)
feature_farction(zoo)
f = sns.FacetGrid(gal,col='feature',sharey=False,xlim=[0,0.2])
f.map(sns.tsplot, 'z','fraction')
plt.show()
