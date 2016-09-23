import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import termcolor
import numpy as np
from astropy.coordinates import SkyCoord as sky
import astropy.units as unit
import astropy.cosmology as COS
from astropy.cosmology import LambdaCDM
import warnings


def calculate_fraction():
    for m in mag_range:
        for z in z_range:
            whole = zoo[(zoo.z > z) & (zoo.z < z + 0.01) & (zoo.Mr < m) & (zoo.mr < 17)]
            bar = whole[(whole.t01_smooth_or_features_a02_features_or_disk_fraction > 0.430) &
                        (whole.t02_edgeon_a05_no_fraction > 0.715) &
                        (whole.t02_edgeon_a05_no_count >= 20)]
            result = bar[bar['gz2class'].str.contains('SB')]
            rid = gal[(gal.Mr == m) & (gal.z == z) & (gal.root == 'bar')].index[0]
            wid = gal[(gal.Mr == m) & (gal.z == z) & (gal.root == 'whole')].index[0]
            gal.at[rid, 'fraction'] = len(result) / len(bar)
            gal.at[wid, 'fraction'] = len(result) / len(whole)
    gal.to_csv('gal2.csv', index=False)


warnings.filterwarnings('ignore')
flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]

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
mag_range = np.arange(-18.89, -21.89, -0.5)
nm = len(mag_range)
root_range = ['whole', 'bar']
nr = len(root_range)
gal['z'] = np.repeat(z_range, nm * nr)
gal['Mr'] = np.tile(np.repeat(mag_range, nr), np.array(nz))
gal['root'] = np.tile(root_range, np.array(nz * nm))

calculate_fraction()
redshift_cut = [COS.z_at_value(LambdaCDM(H0=70, Om0=0.3, Ode0=0.7).distmod, (17 - mag) * unit.mag)
                for mag in mag_range]
gal = pd.read_csv('gal.csv')
g = sns.FacetGrid(gal, hue='root', col='Mr', col_wrap=3, col_order=mag_range, palette=flatui)
g.map(plt.plot, 'z', 'fraction')
g.add_legend()
g.set(ylim=(0, .5))
i = 0
for ax in g.axes.ravel():
    ax.plot([redshift_cut[i],redshift_cut[i]],[0,1],'k--')
    ax.text(redshift_cut[i]+0.01,0.45,'z=%f'% redshift_cut[i],
            fontsize=10,bbox=dict(facecolor='#95a5a6', alpha=0.4))
    i += 1

plt.show()
