import pandas as pd
import termcolor
import numpy as np
from astropy.coordinates import SkyCoord as sky
import astropy.units as unit
import seaborn as sns


def match(sample_ra, sample_dec, catalog_ra, catalog_dec):
    a = sky(ra=sample_ra, dec=sample_dec, unit=unit.deg, frame='icrs')
    b = sky(ra=catalog_ra, dec=catalog_dec, unit=unit.deg, frame='icrs')
    idx, d2d, d3d = a.match_to_catalog_sky(b, nthneighbor=1)
    print('%d/%d  %f' % (len(np.where(d2d.arcsec < 1)[0]),len(a),len(np.where(d2d.arcsec < 5)[0])/len(a)))
    # return np.where(d2d.arcsec < 5)[0]
    # return len(np.where(d2d.arcsecond < 1)[0]) / len(a)


def select_bar(sample):
    return sample[(sample.t01_smooth_or_features_a02_features_or_disk_fraction >= 0.227) &
                  (sample.t02_edgeon_a05_no_fraction >= 0.519) &
                  (sample.t03_bar_a06_bar_count >= 10) &
                  (sample.t03_bar_a06_bar_fraction > 0.3)]


def check_mag(sample):
    sep = 0.01
    for z in np.arange(0,0.2,sep):
        # 从 z 到 z+0.01的区间
        whole = sample[(sample.z <= z+sep) & (sample.z > z)]
        whole = whole.sort_values(by='Mr')
        whole.index = range(len(whole))
        bar_index = whole[whole.mr > 17].index[0]
        barred_sample = select_bar(whole.ix[:bar_index-1])
        print(len(barred_sample)/bar_index, whole.ix[bar_index-1, 'Mr'])

print(sns.load_dataset('titanic'))
# one = pd.read_table('sdssDR4_type1_galaxy.txt', header=0, sep='\s+', engine='python')
# # one = one[one.Z<0.05]
# two = pd.read_table('sdssDR4_type2_galaxy.txt', header=0, sep='\s+', engine='python')
# # two = two[two.Z<0.05]
# zoo = pd.read_csv('zoo2MainSpecz.csv', header=0)
#
# # 预处理,在 Main Sample 加上红移和绝对星等
# uni = pd.read_csv('gz2sample.csv', header=0, index_col='OBJID')
# zoo['Mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_MR.values
# zoo['mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_R.values
# zoo['z'] = uni.ix[zoo.dr7objid.values].REDSHIFT.values
# check_mag(zoo)

