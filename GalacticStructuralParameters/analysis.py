import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import PIL.Image as pi
import subprocess


# sns.set_style('white')
# type1 = pd.read_csv('type1.csv').ix[:784]
# type2 = pd.read_csv('type2.csv').ix[:784]
# sample1 = type1.drop_duplicates(['NAME1'])
# print(sample1)
# sample1 = type1[(type1.I > 0.2) & (type1.I <= 0.999)]
# sample2 = type2.drop_duplicates(['NAME2'])
# bins = np.linspace(0.0, 1, 51)
# print(len(sample1), len(sample2))
# sample2.to_csv('t2.csv', index=None, columns=['NAME2', 'I'])
# print(len(sample1[sample1.INIT > 0.25]), len(sample2[sample2.INIT > 0.25]))

# sns.distplot(type1.A, hist=True, color='r', kde=False, bins=bins, hist_kws={"histtype": "step","linewidth": 2})
# sns.distplot(type1.I, hist=True, color='k', kde=False, bins=bins, hist_kws={"histtype": "step","linewidth": 2})
# sns.distplot(type2.I, hist=True, kde=False, bins=bins, color='b', hist_kws={"histtype": "step","linewidth": 2})
# sns.distplot(type2.A, hist=True, kde=False, bins=bins, color='y', hist_kws={"histtype": "step","linewidth": 2})
# plt.show()
#
ds9 = '/Applications/SAOImage\ DS9.app/Contents/MacOS/ds9'
cmd = '-log -zoom to fit -minmax -invert -cmap value 1.75 0.6'
ct = pd.read_csv('t2.csv')
ct = ct[ct.type == '4']
n = len(ct)
ct.index = range(n)
print(ct)
for i in range(1):
    objs = ''
    for j in range(n):
        objs += '/Users/franky/Desktop/type2/image/%s_r.fits ' % ct.ix[j]['NAME2']
    subprocess.Popen('%s %s %s' % (ds9, cmd, objs), executable='/bin/zsh', shell=True)
