import pandas as pd
from scipy import interpolate
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import LambdaCDM
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import pyprind
import sys
from termcolor import cprint, colored


sns.set_style("darkgrid")
pd.set_option('display.width', 400)
warnings.filterwarnings('ignore')
# tips = sns.load_dataset("tips")
# df = pd.DataFrame(columns=['cnt','agn_type','mor'])
# df.cnt = [10,40,20,30]
# df.agn_type = np.repeat(['type1','type2'],2)
# df.mor = ['E','S','E','S']
# tips.total_bill = np.random.randint(0,2,244)
# print(tips)

# ax = sns.barplot(y="cnt",x='agn_type',hue='mor', data=df)


# ax = sns.countplot(x="day", data=tips)
# arrays = [np.array(['a','bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux']),
#           np.array(['a','one', 'two', 'one', 'two', 'one', 'two', 'one', 'two'])]
# print(pd.DataFrame(np.random.randn(9, 9), columns=arrays))
#
# zoo['Mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_MR.values-np.log10(0.7)*5
# zoo['mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_R.values

# a = pd.read_csv('matched_type2_galaxy.csv')
# a = a[a.gz2class == 'A']
# x = zoo.ix[a.gz2id,['ra','dec','Mr','z']]
# x.index = range(len(x))
# x = x.rename(columns={'Mr':'petro_abs_mag'})
# y=a[['ra','dec','petro_abs_mag','z']]
# y.index=range(len(y))
# print(x-y)
# sns.set(style="darkgrid")
# titanic = sns.load_dataset("titanic")
# print(titanic)
# ax = sns.countplot(x="class", data=titanic)


# zoo = pd.read_csv('zoo2MainSpecz.csv', header=0)
# uni = pd.read_csv('gz2sample.csv', header=0, index_col='OBJID')
# zoo['z'] = uni.ix[zoo.dr7objid.values].REDSHIFT.values
# zoo['Mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_MR.values
# zoo['mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_R.values
#
# bar = zoo[(zoo.gz2class.str.contains('SB')) & (zoo.z < 0.06)]
# bar.to_csv('bar.csv',index=None,columns=['ra','dec','gz2class','z'])


b = pd.read_csv('old_matched_type12.csv').drop_duplicates('name1')
b = b[b.z1 < 0.05]
c = pd.read_csv('type12_controlled_initial.csv')
c = c[(c.z1 < 0.05) & (c.central1 == 1)]
d = pd.read_csv('matched_type2_galaxy.csv')
d = d[d.z < 0.2]
# print(c)
# print(len(a[a.z < 0.05]),len(b[b.Z1 < 0.05]),len(c[(c.z1 < 0.05) & (c.central1 == 1)]))
# x1 = b[b.name2.isin(d.name.values)]
# x2 = c[c.name2.isin(d.name.values)]
# print(len(x1),len(b),len(x2),len(c))
x = b[b.name2.isin(d.name.values)]
# y = zoo.ix[x.gz2id.values,['Mr','mr','z']]
# print(y)
# y['new'] = y.mr-LambdaCDM(H0=70, Om0=0.3, Ode0=0.7).distmod(y.z).value-5*np.log10(0.7)


y = b[b.name2.isin(d.name.values) == False]
y['mr'] = y.mrpetro2+LambdaCDM(H0=70, Om0=0.3, Ode0=0.7).distmod(y.z1).value+5*np.log10(0.7)
x['mr'] = x.mrpetro2+LambdaCDM(H0=70, Om0=0.3, Ode0=0.7).distmod(x.z1).value+5*np.log10(0.7)
a = y[['z1','mr']].append(x[['z1','mr']])
a['t'] = np.concatenate((np.repeat(['not_in'],len(y)),np.repeat(['in'],len(x))))
sns.lmplot('z1','mr',fit_reg=False,data=a,hue='t')
# sns.lmplot('z2','new',fit_reg=False,data=x)
# b.at[x.index,'belong'] = 'in_GZ2'
# b.at[y.index,'belong'] = 'not_in_GZ2'
# sns.lmplot('z2','mrpetro2',data=b,hue='belong',fit_reg=False)
# idx,sep2d,d3d = SkyCoord(ra=x.ra2,dec=x.dec2,unit=u.deg,frame='icrs').match_to_catalog_sky(SkyCoord(ra=bar.ra,dec=bar.dec,unit=u.deg,frame='icrs'),nthneighbor=1)
# print(len(np.where(sep2d.arcsec < 1)[0]),len(x),len(b))
# x = (17-0*np.log10(0.7))*u.mag-LambdaCDM(H0=70, Om0=0.3, Ode0=0.7).distmod(0.050002)
# print(x)
# a = pd.read_csv('old_matched_type12.csv')
# print(a[(a.z2 < 0.05+0.0001) & (a.z2 > 0.05-0.0001)].mrpetro)
# b['app_mag2'] = b.mrpetro2.values+LambdaCDM(H0=70, Om0=0.3, Ode0=0.7).distmod(y.z2.values).value
y.index = range(len(y))
# print(y.app_mag2)
plt.show()
