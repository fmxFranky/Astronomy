import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.cosmology as COS
from astropy.coordinates import SkyCoord
from astropy.cosmology import LambdaCDM
import warnings
from scipy import interpolate
import pyprind
import sys
from termcolor import colored
import warnings

# Mr 的bin:np.arange(-22, -17.5, 0.5)
warnings.filterwarnings('ignore')
flatui = ["#95a5a6", "#3498db", "#9b59b6", "#e74c3c", "#34495e", "#2ecc71", "#1de9b6","#827717"]

# zooniverse = pd.read_csv('zoo.csv')
bar = pd.read_csv('bar.csv')
kind = {
    -1: 'Unclassifiable',
    0: 'Not used',
    1: 'SF',
    2: 'low S/N SF',
    3: 'Composite',
    4: 'AGN non-Liner',
    5: 'Low S/N Liner'
}
# zoo = pd.read_table('GalaxyType.dat', engine='python', sep='\s+', header=None)
# zoo.columns = ['kind']
# zoo.columns = ['col' + str(i) for i in range(1, 10)]
# zoo.index += 1
# bar.at[bar[bar.kind != 'unmatched'].index, 'kind'] = [kind[i] for i in zoo.ix[bar[bar.kind != 'unmatched']['galaxyid'].values, 'kind'].values]
# typ = pd.read_table('GalaxyType.dat', header=None)
# typ.columns = ['t']
# zoo['kind'] = typ.t.values
# one = pd.read_csv('initial_type1.csv')
# two = pd.read_csv('initial_type2.csv')
# gal = pd.read_csv('full.csv')
# idx, sep2d, dist3d = SkyCoord(ra=bar.ra, dec=bar.dec, unit=u.degree).match_to_catalog_sky(
#     SkyCoord(ra=zoo.col3, dec=zoo.col4, unit=u.degree), nthneighbor=1)
# mid = np.where(sep2d.arcsec < 1)[0]
# bar.at[mid, 'galaxyid'] = zoo.ix[idx[mid], 'col1'].values
# zoo = pd.read_csv('SDSS7_ST', engine='python', sep='\s+', header=None)
#
# bar.at[bar[bar.kind != 'unmatched'].index, 'mass'] = zoo.ix[bar[bar.kind != 'unmatched']['galaxyid'].values, 'col4'].values
# bar.to_csv('bar.csv', index=None)

# scale-Mr(带中值)
# med = []
# for i in np.arange(-22, -17.5, 0.5):
#     arr = bar[(bar.abs_mag <= i + 0.5) & (bar.abs_mag > i)]['length_scaled'].values
#     med.append(np.median(arr))
#
# g = sns.FacetGrid(bar, hue='kind', size=7, palette=flatui)
# g.map(plt.scatter, 'abs_mag', 'length_scaled', s=7)
# g.set(ylim=(0, 1.), title='Length_scaled vs Absolutely_magnitude', xlim=(-17, -23))
# plt.plot(np.arange(-22, -17.5, 0.5), med, 'd', color=flatui[2])
# sns.plt.legend(['Median of Scaled Lengths','Normal Galaxies','Strong Agns'])
# plt.savefig('/Users/franky/Desktop/Length_scaled vs Absolutely_magnitude', dpi=100)
# bar = bar[bar.kind != 'unmatched']


# frac-scale/Mr/Mass
# bar['fraction'] = 0
# bar.at[bar[bar.kind == 'agn'].index, 'fraction'] = 1
# bar1 = bar[bar.abs_mag < -20]
# bar1.length_scaled.plot.hist()
# bar2 = bar[bar.abs_mag > -20]
# bar2.length_scaled.plot.hist()

# for mag in np.arange(-22, -17.5, 0.5):
#     bar.at[bar[(bar.abs_mag <= mag + 0.5) & (bar.abs_mag > mag)].index, 'Mr'] = mag
# sns.plt.subplot(2, 2, 1)
# sns.pointplot(x='Mr', y='fraction', data=bar, size=9, estimator=np.mean, ci=68).set(ylim=(-.1, .6), ylabel='fraction')
# for mass in np.arange(9, 11.2, 0.2):
#     bar.at[bar[(bar.mass <= mass + 0.2) & (bar.mass > mass)].index, 'stellar_mass'] = mass
# sns.plt.subplot(2, 2, 2)
#
# sns.pointplot(x='stellar_mass', y='fraction', data=bar, size=9, estimator=np.mean, ci=68).set(ylim=(-.1, .6), ylabel='fraction')
# for ls in np.arange(0.1, .9, 0.2):
#     bar1.at[bar1[(bar1.length_scaled <= ls + 0.2) & (bar1.length_scaled > ls)].index, 'scaled_length'] = ls
#     x = bar1[(bar1.length_scaled <= ls + 0.2) & (bar1.length_scaled > ls)]
#     print(len(x),len(x[x.kind == 'agn']))
#     bar2.at[bar2[(bar2.length_scaled <= ls + 0.2) & (bar2.length_scaled > ls)].index, 'scaled_length'] = ls
# sns.plt.subplot(2, 2, 3)
# sns.pointplot(x='scaled_length', y='fraction', data=bar1,color=flatui[1], size=9, estimator=np.mean, ci=68).set(ylim=(-.1, .6), ylabel='fraction')
# sns.pointplot(x='scaled_length', y='fraction', data=bar2,color=flatui[3], size=9, estimator=np.mean, ci=68).set(ylim=(-.1, .6), ylabel='fraction')
# sns.plt.legend([plt.plot(np.NaN,np.NaN,color=flatui[1])[0],plt.plot(np.NaN,np.NaN,color=flatui[3])[0]],['Mr < -20', 'Mr > -20'])
# plt.savefig('/Users/franky/Desktop/fraction vs scaled_length', dpi=150)
# bar = bar[(bar.kind == 'low S/N SF') ]
bar_copy = bar[bar.kind != 'unmatched']
plt.figure('frac-ls in mag bins')
for t in np.arange(-19,-22, -1):
    bar = bar_copy[(bar_copy.abs_mag > t-1) & (bar_copy.abs_mag < t)]
    for ls in np.arange(0.1, .9, 0.2):
        bar.at[bar[(bar.length_scaled <= ls + 0.2) & (bar.length_scaled > ls)].index, 'scaled_length'] = ls
    # print(bar)
    plt.subplot(2,2,22+t)
    for i in [-1,0,1,2,3,4,5,]:
        bar['fraction'] = 0
        bar.at[bar[bar.kind == kind[i]].index, 'fraction'] = 1
        sns.pointplot(x='scaled_length', y='fraction', data=bar, color=flatui[i], size=9, estimator=np.mean, ci=68).set(ylim=(-.1,0.7))
    plt.legend([plt.plot(np.NaN, np.NaN, color=flatui[j])[0] for j in [-1,0,1,2,3,4,5]], [kind[j] for j in [-1,0,1,2,3,4,5,]])
# plt.savefig('/Users/franky/Desktop/fraction vs scaled_length(all kinds,Mr>-20)', dpi=150)
# TODO 1.
# zooniverse = pd.read_csv('gz2sample.csv')
# zooniverse.index = zooniverse.OBJID.values
# bar['Mr'] = zooniverse.ix[bar.dr7objid.values, 'PETROMAG_MR'].values
# bar['Redshift'] = zooniverse.ix[bar.dr7objid.values, 'REDSHIFT'].values
# bar.to_csv('bar.csv', index=None)

# bar = bar[bar.kind=='unmatched']
# bar.z.plot.hist()
# print(bar.z.describe())
# print(bar[bar.Mr < -19.38])
# bar = bar[bar.kind != 'unmatched']
# bar.abs_mag += 5*np.log10(0.7)
# plt.figure('mag_diff distribution')
# bar['mag_diff'] = bar.abs_mag-bar.Mr
# print(bar.mag_diff.describe())
# bar.mag_diff.plot.hist(bins=np.linspace(-0.2,0.2,10))
# sns.regplot(x='abs_mag',y='Mr',data=bar,fit_reg=False,marker='o',color=flatui[1]).set(ylim=(-17, -23), xlim=(-17,-23))
# sns.regplot(x='Mr',y='length_scaled',data=bar,fit_reg=False,marker='*',color=flatui[1]).set(ylim=(0, 1), xlim=(-17,-23))
# plt.legend([plt.plot(np.NaN,np.NaN,'*',color=flatui[1])[0],plt.plot(np.NaN,np.NaN,'o',color=flatui[0])[0]],['gz2','group'])
# t.mag_diff.plot.hist(bins=np.linspace(-1.5,-0.5,20))
# sns.lmplot(x='Mr',y='length_scaled',data=a.append(b),hue='mag_type',fit_reg=False,size=9,palette=flatui).set(ylim=(0, 1), xlim=(-17,-23))
# sns.lmplot(x='abs_mag',y='length_scaled',data=bar,fit_reg=False,size=9,markers='*').set(ylim=(0, 1), xlim=(-17,-23))
# sns.plt.plot([-19.38, -19.38], [0, 1], 'k--')
# for mass in np.arange(10,11.8,0.3):
#     bar.at[bar[(bar.mass <= mass + 0.3) & (bar.mass > mass)].index, 'Stellar Mass'] = mass
#     print(len(bar[(bar.mass <= mass + 0.3) & (bar.mass > mass)]))
# sns.factorplot(x='Stellar Mass',y='length_scaled',data=bar,estimator=np.mean,ci=68)
sns.plt.show()
