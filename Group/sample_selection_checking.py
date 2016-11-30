import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings

warnings.filterwarnings('ignore')
flatui = ["#95a5a6", "#3498db", "#9b59b6", "#e74c3c", "#34495e", "#2ecc71", "#1de9b6", "#827717"]
bar = pd.read_csv('bar_details.csv')
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


# mag-length_scaled with galaxy_kind
def mag_vs_length():
    # bar = bar[(bar.kind == 'Composite')]
    # sns.lmplot(x='Mr', y='length_scaled', data=bar, hue='kind', palette=flatui, scatter_kws={'s': 9}, fit_reg=False, size=10).set(ylim=(0,1), xlim=(-18, -23))
    for ax in range(1, 6):
        plt.subplot(2, 3, ax)
        sample = bar[bar.kind == kind[ax]]
        sns.kdeplot(sample.length_scaled, sample.Mr, cmap=sns.light_palette(color=flatui[1], as_cmap=True), shade=True, shade_lowest=True).set(xlim=(0, 1.1), ylim=(-18, -23), title=kind[ax])

    # sns.kdeplot(bar.z, bar.Mr, cmap="Blues", shade=True, shade_lowest=False).set(xlim=(0, 0.07), ylim=(-18, -23))


# calculate completeness fraction

mag_vs_length()
# sns.kdeplot(bar.z, bar.Mr, cmap=sns.light_palette(color=flatui[1], as_cmap=True), shade=True, shade_lowest=False).set(xlim=(0, 0.07), ylim=(-18, -23))

sns.plt.show()

