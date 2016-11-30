import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings


# mag-length_scaled with galaxy_kind
warnings.filterwarnings('ignore')
flatui = ["#95a5a6", "#3498db", "#9b59b6", "#e74c3c", "#34495e", "#2ecc71", "#1de9b6", "#827717"]

bar = pd.read_csv('bar_details.csv')
# bar = bar[(bar.kind == 'Composite')]
# sns.lmplot(x='Mr', y='length_scaled', data=bar, hue='kind', palette=flatui, scatter_kws={'s': 9}, fit_reg=False, size=10).set(ylim=(0,1), xlim=(-18, -23))
sns.kdeplot(bar.z, bar.Mr, cmap="Reds", shade=True, shade_lowest=False).set(xlim=(0, 0.07), ylim=(-18, -23))
sns.kdeplot(bar.z, bar.Mr, cmap="Reds", shade=True, shade_lowest=False).set(xlim=(0, 0.07), ylim=(-18, -23))
sns.plt.show()
