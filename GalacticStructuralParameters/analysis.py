import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import PIL.Image as pi
import subprocess


sns.set_style('white')
init1 = pd.read_csv('init1.csv').ix[:784]
init2 = pd.read_csv('init2.csv').ix[:784]
log1 = pd.read_csv('log1.csv').ix[:784]
log2 = pd.read_csv('log2.csv').ix[:784]
sample1 = log1[(log1.LOG > 0.2) & (log1.LOG < 1)]
sample2 = log2[(log2.LOG > 0.2) & (log2.LOG < 1)]
bins = np.linspace(0.0, 1, 51)
print(len(sample1), len(sample2))
# print(len(sample1[sample1.INIT > 0.25]), len(sample2[sample2.INIT > 0.25]))

# sns.distplot(log1['LOG'], hist=True, color='r', kde=False, bins=bins)
# sns.distplot(init1['INIT'], hist=True, color='b', kde=False, bins=bins)
# sns.distplot(init2['INIT'], hist=True, kde=False, bins=bins, color='r')
# sns.distplot(log2['LOG'], hist=True, kde=False, bins=bins)
# plt.show()

ds9 = '/Applications/SAOImage\ DS9.app/Contents/MacOS/ds9'
cmd = '-log -zoom to fit -minmax -invert -cmap value 1.75 0.5'
ct = sample1
n = len(ct)
ct.index = range(n)
print(ct)
for i in range(1):
    objs = ''
    for j in range(0, n, 5):
        objs += '/Users/franky/Desktop/type1/image/%s_r.fits ' % ct.ix[j]['NAME1']
    subprocess.Popen('%s %s %s' % (ds9, cmd, objs), executable='/bin/zsh', shell=True)
