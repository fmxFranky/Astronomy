import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import PIL.Image as pi
import subprocess
import pygal

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

# sns.distplot(type1.A, hist=True, color='r', kde=False, bins=bins, hist_kws={"histtype": "step","barwidth": 2})
# sns.distplot(type1.I, hist=True, color='k', kde=False, bins=bins, hist_kws={"histtype": "step","barwidth": 2})
# sns.distplot(type2.I, hist=True, kde=False, bins=bins, color='b', hist_kws={"histtype": "step","barwidth": 2})
# sns.distplot(type2.A, hist=True, kde=False, bins=bins, color='y', hist_kws={"histtype": "step","barwidth": 2})
# plt.show()
#
#
#
# ds9 = '/Applications/SAOImage\ DS9.app/Contents/MacOS/ds9'
# cmd = '-log -zoom to fit -minmax -invert -cmap value 1.45 0.4 -frame first '
# ct = pd.read_csv('t1.csv')
# # ct = pd.read_csv('t1.csv/Sheet 1-Table 1.csv')
# # ct = ct.drop_duplicates(['NAME2'])
# ct = ct[(ct.bar == 1)]
# ct = ct.sort_index(by=['NAME1'])
# n = len(ct)
# ct.index = range(n)
# print(ct)
# for i in range(1):
#     objs = ''
#     for j in range(n):
#         objs += '/Users/franky/Desktop/type1/image/%s_r.fits ' % ct.ix[j]['NAME1']
#     subprocess.Popen('%s %s %s' % (ds9, cmd, objs), executable='/bin/zsh', shell=True)


def get_ratio(arr):
    return float('{:.2f}'.format(len(arr.values) / n * 100.0))


items = ['hubble', 'sat', 'bar', 'core', 'asa', 'ring', 'gud', 'tl', 'mg', 'onsa']
n = 785
t1 = pd.read_csv('t1.csv', na_filter=0, index_col='NAME1')
i1 = pd.read_csv('type1.csv', na_filter=0).ix[:n-1]
t2 = pd.read_csv('t2.csv', na_filter=0, index_col='NAME2')
i2 = pd.read_csv('type2.csv', na_filter=0).ix[:n-1]
for i in range(n):
    for j in items:
        i1.at[i, j] = t1.at[i1.NAME1[i], j]
        i2.at[i, j] = t2.at[i2.NAME2[i], j]
r1 = [0 for x in range(6)]
r2 = [0 for y in range(6)]
# # HUBBLE
# r1[0] = get_ratio(i1[(i1.hubble == 'e')])
# r1[1] = get_ratio(i1[(i1.hubble == 's') | (i1.sat == '1') | (i1.bar == '1') | (i1.asa == '1') | (i1.onsa == '1')])
# r1[3] = get_ratio(i1[(i1.bar == '1')])
# r1[2] = r1[1] - r1[3]
# r1[4] = get_ratio(i1[(i1.hubble == 'l')])
# r1[5] = 100-(r1[0]+r1[1]+r1[4])
#
# r2[0] = get_ratio(i2[(i2.hubble == 'e')])
# r2[1] = get_ratio(i2[(i2.hubble == 's') | (i2.sat == '1') | (i2.bar == '1') | (i2.asa == '1') | (i2.onsa == '1')])
# r2[3] = get_ratio(i2[(i2.bar == '1')])
# r2[2] = r2[1] - r2[3]
# r2[4] = get_ratio(i2[(i2.hubble == 'l')])
# r2[5] = 100-(r2[0]+r2[1]+r2[4])

# # ASYMMETRY
# r1[0] = get_ratio(i1[(i1.bar == '1')])
# r1[1] = get_ratio(i1[(i1.asa == '1') | (i1.onsa == '1')])
# r1[2] = get_ratio(i1[(i1.gud == '1')])
# r1[3] = get_ratio(i1[(i1.tl == 'l') | ((i1.bar != '1') & (i1.ring == '1'))])
# r1[4] = get_ratio(i1[(i1.mg == '1')])
#
# r2[0] = get_ratio(i2[(i2.bar == '1')])
# r2[1] = get_ratio(i2[(i2.asa == '1') | (i2.onsa == '1')])
# r2[2] = get_ratio(i2[(i2.gud == '1')])
# r2[3] = get_ratio(i2[(i2.tl == 'l') | ((i2.bar != '1') & (i2.ring == '1'))])
# r2[4] = get_ratio(i2[(i2.mg == '1')])

# FEATURE
r1[0] = get_ratio(i1[(i1.bar == '1') & (i1.core == '1')])
r1[1] = get_ratio(i1[(i1.bar == '1') & (i1.ring == '1')])
r1[2] = get_ratio(i1[(i1.onsa == '1')])
r1[3] = get_ratio(i1[(i1.asa == '1')])
r1[4] = get_ratio(i1[(i1.tl == '1')])
r1[5] = get_ratio(i1[(i1.bar != '1') & (i1.ring == '1')])

r2[0] = get_ratio(i2[(i2.bar == '1') & (i2.core == '1')])
r2[1] = get_ratio(i2[(i2.bar == '1') & (i2.ring == '1')])
r2[2] = get_ratio(i2[(i2.onsa == '1')])
r2[3] = get_ratio(i2[(i2.asa == '1')])
r2[4] = get_ratio(i2[(i2.tl == '1')])
r2[5] = get_ratio(i2[(i2.bar != '1') & (i2.ring == '1')])

bar_chart = pygal.Bar(legend_at_bottom=True, legend_at_bottom_columns=2)
bar_chart._series_margin = .2
bar_chart._serie_margin = .1
bar_chart.title = 'The Asymmetric Features Distribution(in %)'
bar_chart.x_labels = map(str, ['Bar&Core', 'Bar&Ring', 'Odd Arms', 'Asymmetric Arms', 'Tail', 'Peripheral Circle'])
bar_chart.add('Type-I', r1)
bar_chart.add('Type-II', r2)
bar_chart.render_in_browser()
