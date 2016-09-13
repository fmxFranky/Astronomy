import pandas as pd
import termcolor
import numpy as np
from astropy.coordinates import SkyCoord as sky
import astropy.units as unit


def match(sample_ra, sample_dec, catalog_ra, catalog_dec):
    a = sky(ra=sample_ra, dec=sample_dec, unit=unit.deg, frame='icrs')
    b = sky(ra=catalog_ra, dec=catalog_dec, unit=unit.deg, frame='icrs')
    idx, d2d, d3d = a.match_to_catalog_sky(b, nthneighbor=1)
    print('%d/%d' % (len(np.where(d2d.arcsecond < 1)[0]),len(a)))
    return len(np.where(d2d.arcsecond < 1)[0]) / len(a)


one = pd.read_table('sdssDR4_type1_galaxy.txt', header=0, sep='\s+', engine='python')
# one = one[one.Z<0.05]
two = pd.read_table('sdssDR4_type2_galaxy.txt', header=0, sep='\s+', engine='python')
# two = two[two.Z<0.05]
zoo = pd.read_csv('zoo2bar.csv', header=0)
x = match(one.RA, one.DEC, zoo.ra, zoo.dec)
y = match(two.RA, two.DEC, zoo.ra, zoo.dec)
print(x*100)
print(y*100)
print(x/y)
# type12 = pd.read_table('type12_match.asc', header=0, sep='\s+', engine='python')
# type12 = type12[type12.Z1<0.05]
# print(match(type12.RA1, type12.DEC1, zoo.ra, zoo.dec))
# print(match(type12.RA2, type12.DEC2, zoo.ra, zoo.dec))

