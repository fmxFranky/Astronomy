import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings

warnings.filterwarnings('ignore')
flatui = ["#95a5a6", "#3498db", "#9b59b6", "#e74c3c", "#34495e", "#2ecc71", "#1de9b6", "#827717"]
bar = pd.read_csv('bar.csv')
bar = bar[bar.kind != 'unmatched']
kind = {
    -1: 'Unclassifiable',
    0: 'Not used',
    1: 'SF',
    2: 'low S/N SF',
    3: 'Composite',
    4: 'AGN non-Liner',
    5: 'Low S/N Liner'
}
# zoo = pd.read_csv('zoo.csv', engine='c', header=0)
# gz = pd.read_csv('gz.csv', engine='c', index_col='OBJID', header=0)
# zoo['z'] = gz.ix[zoo.dr7objid.values].REDSHIFT.values
# zoo['Mr'] = gz.ix[zoo.dr7objid.values].PETROMAG_MR.values


# mag-length_scaled with galaxy_kind
def mag_vs_length():
    # bar = bar[(bar.kind == 'Composite')]
    # sns.lmplot(x='Mr', y='length_scaled', data=bar, hue='kind', palette=flatui, scatter_kws={'s': 9}, fit_reg=False, size=10).set(ylim=(0,1), xlim=(-18, -23))
    for ax in range(1, 6):
        plt.subplot(2, 3, ax)
        sample = bar[bar.kind == kind[ax]]
        sns.kdeplot(sample.length_scaled, sample.Mr, cmap=sns.light_palette(color=flatui[1], as_cmap=True), shade=True, shade_lowest=True).set(xlim=(0, 1.1), ylim=(-18, -23), title=kind[ax])

        # sns.kdeplot(bar.z, bar.Mr, cmap="Blues", shade=True, shade_lowest=False).set(xlim=(0, 0.07), ylim=(-18, -23))


# calculate cross fraction
def cross_fraction(a, b):
    print()


# print(gz[(gz.REDSHIFT < 0.06) & (gz.REDSHIFT > 0.01)])
# print(len(gz), len(gz[~np.isnan(gz.REDSHIFT)]), len(gz[~np.isnan(gz.REDSHIFTERR)]))
print(len(bar[(bar.z > 0.01) & (bar.z < 0.06) & (bar.abs_mag < -19.38)]))
print(len(bar[(bar.z > 0.01) & (bar.z < 0.06) & (bar.abs_mag < -20.15)]))
print(len(bar[(bar.z > 0.01) & (bar.z < 0.06) & (bar.Mr < -19.38)]))
print(len(bar[(bar.z > 0.01) & (bar.z < 0.06) & (bar.Mr < -20.15)]))
# sns.regplot(x='z', y='Mr', data=zoo, scatter_kws={'s': 4, 'c': flatui[0]}, fit_reg=False).set(xlim=(0, 0.07), ylim=(-17, -24))
# sns.regplot(x='z', y='Mr', data=zoo.ix[bar.gz2id], scatter_kws={'s': 11, 'c': flatui[1]}, fit_reg=False, scatter='*').set(xlim=(0, 0.07), ylim=(-17, -24))
# selection = zoo[(zoo.z > 0.01) & (zoo.z < 0.06) & (zoo.t03_bar_a06_bar_fraction >= 0.8)]
# sns.regplot(x='z', y='Mr', data=bar, scatter_kws={'s': 4, 'c': flatui[1]}, fit_reg=False).set(xlim=(0, 0.07), ylim=(-17, -24))

sns.plt.show()
