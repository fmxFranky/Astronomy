import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import termcolor
import numpy as np
from astropy.coordinates import SkyCoord as sky
import astropy.units as unit
from astropy.cosmology import LambdaCDM
import astropy.cosmology as cosmo
import pyprind
import time
import sys
from termcolor import cprint, colored
import warnings
import re


def split_features(agn_type):
    """
    对 type1/2 的 gz2class 进行拆分, 并返回文件
    :parameter: agn_type
    :return:
    """
    agn_type_morph = pd.DataFrame(np.full((len(control_sample), len(feature_columns)), np.nan), columns=feature_columns)
    agn_type_morph[['name', 'gz2class']] = control_sample[['name1', 'gz2class1']] if agn_type == 'type1' else control_sample[['name2', 'gz2class2']]
    progress_bar = pyprind.ProgBar(len(agn_type_morph), stream=sys.stdout, bar_char='█', width=47, title=colored('controlling sample', color='blue', attrs=['bold']))

    for i in agn_type_morph.index:
        gz2class = agn_type_morph.ix[i, 'gz2class']
        if 'E' in gz2class:
            agn_type_morph.at[i, 'smooth'] = 1
            if 'Er' in gz2class:
                agn_type_morph.at[i, 'completely_round'] = 1
            elif 'Ei' in gz2class:
                agn_type_morph.at[i, 'in_between'] = 1
            elif 'Ec' in gz2class:
                agn_type_morph.at[i, 'cigar_shaped'] = 1
        elif 'S' in gz2class:
            agn_type_morph.at[i, 'features_or_disks'] = 1
            if 'Se' in gz2class:
                agn_type_morph.at[i, 'edge_on'] = 1
                if 'Ser' in gz2class:
                    agn_type_morph.at[i, 'round_bulge_shape'] = 1
                if 'Seb' in gz2class:
                    agn_type_morph.at[i, 'boxy_bulge_shape'] = 1
                if 'Sen' in gz2class:
                    agn_type_morph.at[i, 'none_bulge_shape'] = 1
            elif 'SB' in gz2class:
                agn_type_morph.at[i, 'bar'] = 1
            if 'Sd' in gz2class.replace('B', ''):
                agn_type_morph.at[i, 'none_bulge_prominence'] = 1
            elif 'Sc' in gz2class.replace('B', ''):
                agn_type_morph.at[i, 'just_noticeable_bulge_prominence'] = 1
            elif 'Sb' in gz2class.replace('B', ''):
                agn_type_morph.at[i, 'obvious_bulge_prominence'] = 1
            elif 'Sa' in gz2class.replace('B', ''):
                agn_type_morph.at[i, 'dominant_bulge_prominence'] = 1

            unodd_gz2class = gz2class
            for ch in 'rldiomu':
                unodd_gz2class = gz2class.replace('(%s)' % ch, '')
            if re.findall(r'[0-9?+lmt]', unodd_gz2class):
                agn_type_morph.at[i, 'spiral_structure'] = 1
        elif 'A' in gz2class:
            agn_type_morph.at[i, 'star'] = 1

        if '(' in gz2class and ')' in gz2class:
            agn_type_morph.at[i, 'odd_feature'] = 1
            if '(r)' in gz2class:
                agn_type_morph.at[i, 'ring'] = 1
            if '(l)' in gz2class:
                agn_type_morph.at[i, 'lens_arc'] = 1
            if '(d)' in gz2class:
                agn_type_morph.at[i, 'disturbed'] = 1
            if '(i)' in gz2class:
                agn_type_morph.at[i, 'irregular'] = 1
            if '(o)' in gz2class:
                agn_type_morph.at[i, 'other'] = 1
            if '(m)' in gz2class:
                agn_type_morph.at[i, 'merge'] = 1
            if '(u)' in gz2class:
                agn_type_morph.at[i, 'dust'] = 1

        progress_bar.update()

    agn_type_morph.to_csv('matched_morph_detail_%s_galaxy.csv' % agn_type, index=None)


def get_counts(z_cut,color):
    """
    统计 type1/2 的各个特征的数量
    :return:
    """
    # z_cut = 0.05
    idm = control_sample[(control_sample.z1 < z_cut) & (control_sample.central1 == 1)].index
    morph_type1 = pd.read_csv('matched_morph_detail_type1_galaxy.csv').ix[idm]
    cnt_type1 = morph_type1.count(axis=0)
    morph_type2 = pd.read_csv('matched_morph_detail_type2_galaxy.csv').ix[idm]
    cnt_type2 = morph_type2.count(axis=0)
    cnt = pd.DataFrame(columns=['agn_type', 'feature','count'])
    cnt['feature'] = np.tile(feature_columns[2:], [2])
    cnt['agn_type'] = np.repeat(['type1', 'type2'], len(feature_columns[2:]))
    cnt['count'] = np.concatenate((cnt_type1.values[2:],cnt_type2.values[2:]))/len(idm)
    # print(cnt)
    print(z_cut,len(idm))
    s = sns.barplot('feature','count',hue='agn_type',ci=False,palette=color,
                data=cnt[cnt['feature'].isin(['features_or_disks','edge_on','bar','odd_feature','spiral_structure','merge'])])
    s.set(ylabel='percentage',title='z=%f' % z_cut)

if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    pd.set_option('display.width', 400)
    flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
    feature_columns = ['name', 'gz2class',
                       'smooth', 'completely_round', 'in_between', 'cigar_shaped',
                       'features_or_disks',
                       'edge_on', 'round_bulge_shape', 'boxy_bulge_shape', 'none_bulge_shape',
                       'bar',
                       'none_bulge_prominence', 'just_noticeable_bulge_prominence', 'obvious_bulge_prominence', 'dominant_bulge_prominence',
                       'spiral_structure',
                       'odd_feature', 'ring', 'lens_arc', 'disturbed', 'irregular', 'other', 'merge', 'dust',
                       'star']
    control_sample = pd.read_csv('type12_controlled_initial.csv')
    # for t in ['type1', 'type2']:
    #     split_features(t)
    plt.subplot(221)
    get_counts(0.05,sns.color_palette('deep'))
    plt.subplot(222)
    get_counts(0.1, sns.color_palette('hls'))
    plt.subplot(223)
    get_counts(0.15, sns.color_palette('husl'))
    plt.subplot(224)
    get_counts(0.2, sns.color_palette(flatui))
    plt.show()
