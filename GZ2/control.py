import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as unit
import astropy.cosmology as COS
from astropy.cosmology import LambdaCDM
import warnings
from scipy import interpolate


def merge_init_and_mat(agn_type):
    # Merge type1/2 initial data and matched data
    x = pd.read_csv('initial_%s_galaxy.csv' % agn_type, index_col='name')
    y = pd.read_csv('matched_%s_galaxy.csv' % agn_type, index_col='name')
    for col in ['lum_oiii', 'gz2id', 'gz2class']:
        x.at[y.index, col] = y[col]
    zfl_type1 = pd.read_csv('zfl_type1_information.csv')
    zfl_type2 = pd.read_csv('zfl_type2_information.csv')
    zfl = zfl_type1.append(zfl_type2, ignore_index=True).drop_duplicates('z').sort_values(by='z')
    itp = interpolate.UnivariateSpline(zfl.z, np.sqrt(np.power(10, zfl.luminous - zfl.flux)), bbox=[0, 0.25], k=1)
    x['petro_app_mag'] = x.petro_abs_mag + LambdaCDM(H0=70, Om0=0.3, Ode0=0.7).distmod(x.z).value + 5 * np.log10(0.7)
    x.lum_oiii = np.log10(np.square(itp(x.z))) + x.flux_oiii
    x.reset_index(inplace=True)
    x.to_csv('initial_%s.csv' % agn_type, index=None, columns=['name', 'ra', 'dec', 'z', 'lum_oiii', 'petro_abs_mag', 'petro_app_mag', 'central', 'gz2id', 'gz2class'])


def control(z_cut, pre_match=False):
    mag_cut = 17 - LambdaCDM(H0=70, Om0=.3, Ode0=.7).distmod(z_cut).value - 5 * np.log10(0.7)

    x = pd.read_csv('initial_type1.csv')
    y = pd.read_csv('initial_type2.csv')
    if pre_match:
        x = x[~x.gz2id.isnull()]
        y = y[~y.gz2id.isnull()]
    x = x[(x.z < z_cut) & (x.petro_abs_mag < mag_cut)]
    y = y[(y.z < z_cut) & (y.petro_abs_mag < mag_cut)]

    import pyprind
    import sys
    from termcolor import colored
    progress_bar = pyprind.ProgBar(len(x), stream=sys.stdout, bar_char='█', width=47, title=colored('Controlling sample   ', color='blue', attrs=['bold']) +
                                                                                            colored('z_cut=%f mag_cut=%f pre_match=%s' % (z_cut, mag_cut, pre_match), color='cyan', attrs=['bold']))
    ctrl = pd.DataFrame(columns=['name_1', 'ra_1', 'dec_1', 'z_1', 'lum_oiii_1', 'petro_abs_mag_1', 'petro_app_mag_1', 'central_1', 'gz2id_1', 'gz2class_1',
                                 'name_2', 'ra_2', 'dec_2', 'z_2', 'lum_oiii_2', 'petro_abs_mag_2', 'petro_app_mag_2', 'central_2', 'gz2id_2', 'gz2class_2'])
    for i in x.index:
        agn = x.ix[i]
        cs = y[(abs(y.petro_abs_mag - agn.petro_abs_mag) < 0.1) &
               (abs(y.z - agn.z) < 0.01) &
               (y.central == agn.central)]
        if len(cs) >= 5:
            cs['diff_lum'] = abs(cs.lum_oiii - agn.lum_oiii)
            cs.sort_values(inplace=True, by='diff_lum')
            cs.drop('diff_lum', 1, inplace=True)
            for j in cs.index[:5]:
                ctrl = ctrl.append(pd.DataFrame([np.concatenate((agn.values, cs.ix[j].values), axis=0)],
                                                columns=ctrl.columns), ignore_index=True)
        progress_bar.update()

    ctrl.to_csv('match_type12_%.2f_%s.csv' % (z_cut, pre_match), index=False)


def transform_format(z_cut, pre_match):
    x = pd.read_csv('match_type12_%.2f_%s.csv' % (z_cut, pre_match))
    y = pd.DataFrame(data=np.full((len(x) * 2, 11), np.nan), columns=['name', 'ra', 'dec', 'z', 'lum_oiii', 'petro_abs_mag', 'petro_app_mag', 'central', 'gz2id', 'gz2class', 'agn'])
    y.agn = np.tile(['type1', 'type2'], [len(x)])

    for col in cls[:len(cls) - 1]:
        y.at[0:len(x) * 2:2, col] = x.ix[:len(x), col + '_1'].values
        y.at[1:len(x) * 2:2, col] = x.ix[:len(x), col + '_2'].values
    return y


def distribution_check():
    od = 1
    z_cut_range = np.linspace(0.05, 0.09, 3)
    plt.figure(1).suptitle('UnMatched Sample')
    z_bins = np.linspace(0., 0.09, 20)
    l_bins = np.linspace(38.5, 42.5, 20)
    m_bins = np.linspace(-22.5, -18.5, 20)
    for z in z_cut_range:
        # for pm in [False, True]:
        a = transform_format(z, False)
        a = a[a.central == 1]

        plt.subplot(3, 3, od)
        sns.distplot(a[a.agn == 'type1'].z,
                     color='b', hist=True, kde=False, bins=z_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1}).set(ylabel='z_cut=%.2f\n\ncount' % z)
        sns.distplot(a[a.agn == 'type2'].z,
                     color='r', hist=True, kde=False, bins=z_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1})  # .legend(['type1','type2'])
        if od == 1:
            plt.legend(['type1', 'type2'], loc='upper left')

        plt.subplot(3, 3, od + 1)
        sns.distplot(a[a.agn == 'type1'].lum_oiii,
                     color='b', hist=True, kde=False, bins=l_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1})
        sns.distplot(a[a.agn == 'type2'].lum_oiii,
                     color='r', hist=True, kde=False, bins=l_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1})  # .legend(['type1', 'type2'])

        plt.subplot(3, 3, od + 2)
        sns.distplot(a[a.agn == 'type1'].petro_abs_mag,
                     color='b', hist=True, kde=False, bins=m_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1})
        sns.distplot(a[a.agn == 'type2'].petro_abs_mag,
                     color='r', hist=True, kde=False, bins=m_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1})  # .legend(['type1', 'type2'])
        od += 3


def show_diff():
    z_cut_range = np.linspace(0.05, 0.09, 3)
    xy = pd.DataFrame(columns=['agn', 'fraction', 'z', 'pre_matched','total'])
    for z in z_cut_range:
        for pm in [False, True]:
            x = transform_format(z, pm)
            x = x[x.central == 1]
            y = pd.DataFrame([['type1', len(x[(x.agn == 'type1') & (x.gz2class.str.contains('SB'))]) / len(x[x.agn=='type1']), z, pm,len(x)/2],
                              ['type2', len(x[(x.agn == 'type2') & (x.gz2class.str.contains('SB'))]) / len(x[x.agn=='type2']), z, pm,len(x)/2]],
                             columns=['agn', 'fraction', 'z', 'pre_matched','total'])
            xy = xy.append(y)
    xy['num'] = xy.total*xy.fraction
    print(xy)
    sns.set(style="ticks")
    sns.factorplot('z', 'fraction', hue='agn', col='pre_matched',ci=None, data=xy, size=7, alpha=.7,
                   palette=sns.color_palette(flatui))


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    pd.set_option('display.width', 200)
    flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
    cls = ['name', 'ra', 'dec', 'z', 'lum_oiii', 'petro_abs_mag', 'petro_app_mag', 'central', 'gz2id', 'gz2class', 'agn']

    # w = {}
    # for i in ['name', 'ra', 'dec', 'z', 'lum_oiii', 'petro_abs_mag', 'petro_app_mag', 'central', 'gz2id', 'gz2class']:
    #     w[i+'1'] = i+'_1'
    #     w[i+'2'] = i+'_2'
    # distribution_check()
    show_diff()
    plt.show()
    # print(a)
