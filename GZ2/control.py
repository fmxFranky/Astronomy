import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.cosmology as COS
from astropy.coordinates import SkyCoord
from astropy.cosmology import LambdaCDM
import warnings
from scipy import interpolate
import pyprind
import sys
from termcolor import colored
from scipy.stats import poisson


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


def mag_at_value(z):
    return 17 - LambdaCDM(H0=70, Om0=.3, Ode0=.7).distmod(z).value - 5 * np.log10(0.7)


def control(z_cut, mag_cut):
    x = pd.read_csv('initial_type1.csv')
    y = pd.read_csv('initial_type2.csv')
    z = pd.read_csv('full.csv')

    x = x[(x.z < z_cut) & (x.petro_abs_mag < mag_cut) & (~x.gz2id.isnull())]
    y = y[(y.z < z_cut) & (y.petro_abs_mag < mag_cut) & (~y.gz2id.isnull())]
    z = z[(z.z < z_cut) & (z.petro_abs_mag < mag_cut) & (~z.gz2id.isnull())]

    title = colored('Controlling sample   ', color='blue', attrs=['bold'])
    info = colored('z_cut=%f mag_cut=%f' % (z_cut, mag_cut), color='cyan', attrs=['bold'])
    progress_bar = pyprind.ProgBar(len(x), stream=sys.stdout, bar_char='█', width=47, title=title + info)

    ctrl = pd.DataFrame()
    for i in x.index:
        one = x.ix[i]
        two = y[(abs(y.petro_abs_mag - one.petro_abs_mag) < 0.1) &
                (abs(y.z - one.z) < 0.01) &
                (y.central == one.central)]

        if len(two) >= 5:
            two['diff_lum'] = abs(two.lum_oiii - one.lum_oiii)
            two.sort_values(inplace=True, by='diff_lum')
            two.drop('diff_lum', 1, inplace=True)
            two = two[:5]

            normtwo = pd.DataFrame()
            for j in two.index:
                tw = two.ix[j]
                ntw = z[(abs(z.z - tw.z) < 0.01) & (z.central == tw.central) & (z.galaxyid != tw.galaxyid)]
                ntw['diff_mag'] = abs(ntw.petro_abs_mag - tw.petro_abs_mag)
                ntw.sort_values(inplace=True, by='diff_mag')
                ntw.drop('diff_mag', 1, inplace=True)
                normtwo = normtwo.append(ntw[:1])

            normone = z[(abs(z.z - one.z) < 0.01) & (z.central == one.central) & (z.galaxyid != one.galaxyid)]
            normone['diff_mag'] = abs(normone.petro_abs_mag - one.petro_abs_mag)
            normone.sort_values(inplace=True, by='diff_mag')
            normone.drop('diff_mag', 1, inplace=True)
            normone = normone[:5]

            ctrl = ctrl.append(one.to_frame().T.append(two.append(normtwo.append(normone))))

        progress_bar.update()

    ctrl.to_csv('match_%.2f_%.2f.csv' % (z_cut, mag_cut), index=False)


def distribution_check():
    od = 1
    z_cut_range = np.linspace(0.05, 0.09, 3)
    plt.figure(1).suptitle('UnMatched Sample')
    z_bins = np.linspace(0., 0.09, 20)
    l_bins = np.linspace(38.5, 42.5, 20)
    m_bins = np.linspace(-22.5, -18.5, 20)
    for z in z_cut_range:
        # for pm in [False, True]:
        a = pd.read_csv('match.csv')
        a = a[a.central == 1]

        plt.subplot(3, 3, od)
        sns.distplot(a[a.kind == 'type1'].z,
                     color='b', hist=True, kde=False, bins=z_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1}).set(ylabel='z_cut=%.2f\n\ncount' % z)
        sns.distplot(a[a.kind == 'type2'].z,
                     color='r', hist=True, kde=False, bins=z_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1})  # .legend(['type1','type2'])
        if od == 1:
            plt.legend(['type1', 'type2'], loc='upper left')

        plt.subplot(3, 3, od + 1)
        sns.distplot(a[a.kind == 'type1'].lum_oiii,
                     color='b', hist=True, kde=False, bins=l_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1})
        sns.distplot(a[a.kind == 'type2'].lum_oiii,
                     color='r', hist=True, kde=False, bins=l_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1})  # .legend(['type1', 'type2'])

        plt.subplot(3, 3, od + 2)
        sns.distplot(a[a.kind == 'type1'].petro_abs_mag,
                     color='b', hist=True, kde=False, bins=m_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1})
        sns.distplot(a[a.kind == 'type2'].petro_abs_mag,
                     color='r', hist=True, kde=False, bins=m_bins,
                     hist_kws={"histtype": "step", 'alpha': 1, 'linewidth': 1})  # .legend(['type1', 'type2'])
        od += 3


def show_diff():
    ctrl = pd.DataFrame()
    for z_cut in np.linspace(0.05, 0.09, 3):
        for mag_cut in [0, mag_at_value(z_cut)]:
            x = pd.read_csv('match_%.2f_%.2f.csv' % (z_cut, mag_cut))
            x = x[x.central == 1]
            x['bar'] = 0
            for k in ['type1', 'type2', 'normal']:
                spiral_fraction = len(x[x.kind == k]) / len(x[(x.kind == k) & (x.gz2class.str.contains('S')) & (~x.gz2class.str.contains('Se'))])
                x['bar'][(~x.gz2id.isnull()) & (x.gz2class.str.contains('SB'))] = 1
            x['cut_z'] = z_cut
            x['cut_mag'] = 'mag(z)' if mag_cut == mag_at_value(z_cut) else ('mag(z+0.02)' if mag_cut == mag_at_value(z_cut + 0.02) else 'none')
            ctrl = ctrl.append(x)
    ctrl.to_csv('hh.csv', index=None)
    sns.factorplot(x='cut_z', y='bar', hue='kind', col='cut_mag', data=ctrl, size=7,
                   estimator=np.mean, ci=68, capsize=.1, n_boot=1000).fig.suptitle('Bar/Spiral')
    return


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    pd.set_option('display.width', 200)
    flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
    cls = ['galaxyid', 'kind', 'ra', 'dec', 'z', 'app_mag', 'petro_abs_mag', 'central', 'lum_oiii', 'gz2id', 'gz2class']
    sns.set(style="ticks")
    # a = pd.read_csv('match_0.05_0.00.csv')
    # print(len(a[a.kind == 'type1']), len(a[(a.kind=='type1') & (a.gz2class.str.contains('S'))]), len(a[(a.kind=='type1') & (a.gz2class.str.contains('SB'))]))
    # a = pd.read_csv('hh.csv')
    # a = a[a.cut_mag != 'mag(z+0.02)']
    # sns.factorplot(x='cut_z', y='bar', hue='kind', col='cut_mag', data=a, size=7,
    #                estimator=np.mean, ci=68, capsize=.1, n_boot=1000).fig.suptitle('Bar/Spiral')
    show_diff()
    # print(a[(a.kind == 'type1') & (a.gz2class.str.contains('S')) & (~a.gz2class.str.contains('Se'))])
    plt.show()
    # print(a)
