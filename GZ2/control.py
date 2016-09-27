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
    mag_cut = 17 - LambdaCDM(H0=70, Om0=.3, Ode0=.7).distmod(z_cut).value

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
    ctrl = pd.DataFrame(columns=['name1', 'ra1', 'dec1', 'z1', 'lum_oiii1', 'petro_abs_mag1', 'petro_app_mag1', 'central1', 'gz2id1', 'gz2class1',
                                 'name2', 'ra2', 'dec2', 'z2', 'lum_oiii2', 'petro_abs_mag2', 'petro_app_mag2', 'central2', 'gz2id2', 'gz2class2'])
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


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    pd.set_option('display.width', 200)
    # for z in np.linspace(0.05, 0.09, 3):
    #     for pm in [False, True]:
    a = pd.read_csv('match_type12_%.2f_%s.csv' % (0.05, False))
    a[['z1','z2']].plot.hist()
    plt.show()
