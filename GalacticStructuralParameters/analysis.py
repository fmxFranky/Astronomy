import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import PIL.Image as pi

sns.set_style('white')
init2 = pd.read_csv('init2.csv').ix[:784]
sample = init2[init2.INIT > 0.3]
# sample['index'] = range(len(sample))
print(sample.drop_duplicates(['NAME2']))
# for fit in sample.drop_duplicates(['NAME2']).NAME2.values:
#     pi.open('/Users/franky/Desktop/templates/init/%s.jpg' % fit).show()

