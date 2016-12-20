import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings


def creat_cross_sample():
    group_galaxies = pd.read_csv('group_galaxies.csv', engine='c')
    group_galaxy_type = pd.read_table('group_galaxy_type.dat', engine='c', header=None)
    group_galaxies['typ'] = [galaxy_types[i[0]] for i in group_galaxy_type.values]
    cross = group_galaxies[~np.isnan(group_galaxies.gz2id)]
    gz2 = pd.read_csv('GZ2.csv', engine='c')
    cross['dr7id'] = gz2.ix[cross.gz2id, 'dr7objid'].values
    cross['bar_debiased_fraction'] = gz2.ix[cross.gz2id, 't03_bar_a06_bar_debiased'].values
    n4_bdd = pd.read_csv('n4_bdd.csv', engine='c', index_col='col0')
    cross['B_T_n4'] = n4_bdd.ix[cross.dr7id, 'col14'].values
    cross['Re_n4'] = n4_bdd.ix[cross.dr7id, 'col22'].values
    cross['be_n4'] = n4_bdd.ix[cross.dr7id, 'col24'].values
    cross['Rd_n4'] = n4_bdd.ix[cross.dr7id, 'col28'].values
    cross['dia_n4'] = n4_bdd.ix[cross.dr7id, 'col30'].values
    fn_bdd = pd.read_csv('fn_bdd.csv', engine='c', index_col='col0')
    cross['B_T_fn'] = fn_bdd.ix[cross.dr7id, 'col14'].values
    cross['Re_fn'] = fn_bdd.ix[cross.dr7id, 'col22'].values
    cross['be_fn'] = fn_bdd.ix[cross.dr7id, 'col24'].values
    cross['Rd_fn'] = fn_bdd.ix[cross.dr7id, 'col28'].values
    cross['dia_fn'] = fn_bdd.ix[cross.dr7id, 'col30'].values
    cross_gz2_group = cross[['galaxyid', 'gz2id', 'z', 'petro_abs_mag', 'gz2class', 'bar_debiased_fraction', 'typ',
                             'B_T_n4', 'Re_n4', 'be_n4', 'Rd_n4', 'dia_n4',
                             'B_T_fn', 'Re_fn', 'be_fn', 'Rd_fn', 'dia_fn']]
    cross_gz2_group.to_csv('cross_gz2_group_bdd.csv', index=None)


def control_sample():
    cross = pd.read_csv('cross_gz2_group_bdd.csv', engine='c')
    vls = cross[(~np.isnan(cross.B_T_n4)) & (cross.B_T_n4 > -1) &
                (cross.z < 0.06) & (cross.petro_abs_mag < -19.38)]
    vls['kind'] = 'unbar'
    bar_sample = vls[(vls.bar_debiased_fraction > 0.5) & (vls.gz2class.str.contains('B'))]
    bar_sample['kind'] = 'bar'

    # g = sns.regplot(x='z', y='petro_abs_mag', data=bar_sample, scatter_kws={'s': 7}, color=flatui[1], fit_reg=False)
    # g.set(ylim=(-19, -22), title='Scatter of Mr-z and mean(B/T)')
    # h = g.twinx()
    # h.plot(np.arange(0.02, 0.06, 0.005)+0.0025,
    #                [np.mean(bar_sample[(bar_sample.z < z+0.005) & (bar_sample.z > z)]['B_T_n4']) for z in np.arange(0.02, 0.06, 0.005)],
    #                'kv--')
    # h.set(ylim=(0, 0.6), ylabel='mean(B/T)')

    ctrl = pd.DataFrame()
    for i in bar_sample.index[:]:
        bs = bar_sample.ix[i]
        cubs = vls[(abs(vls.z - bs.z) < 0.01) &
                   (abs(vls.petro_abs_mag - bs.petro_abs_mag) < 0.1) &
                   (abs(vls.B_T_n4 - bs.B_T_n4) < 0.01) &
                   (vls.gz2id != bs.gz2id)]
        if len(cubs) >= 3:
            cubs['diff_B_T_n4'] = abs(cubs.B_T_n4 - bs.B_T_n4)
            cubs['diff_Mr'] = abs(cubs.petro_abs_mag - bs.petro_abs_mag)
            cubs.sort_values(inplace=True, by=['diff_B_T_n4', 'diff_Mr'])
            cubs.drop(['diff_B_T_n4', 'diff_Mr'], 1, inplace=True)
            ctrl = ctrl.append(bs.to_frame().T.append(cubs[:3]))

    print(len(ctrl[(ctrl.kind == 'unbar') & (ctrl.gz2clss.str.contain('B'))]), len(ctrl))
    # print(control_unbar_sample)


def verify():
    gz2 = pd.read_csv('GZ2.csv', engine='c')
    n4_bdd = pd.read_csv('n4_bdd.csv', engine='c')
    # meta = pd.read_csv('gz2sample.csv', engine='c')
    # from astropy.coordinates import SkyCoord
    # from astropy import units as u
    # idx,sep2d,d3d = SkyCoord(ra=gz2.ra2, dec=gz2.dec2, unit=u.deg, frame='icrs').match_to_catalog_sky(SkyCoord(ra=n4_bdd.ra,dec=bar.dec,unit=u.deg,frame='icrs'),nthneighbor=1)


def plot_distribution():
    cross = pd.read_csv('cross_gz2_group_bdd.csv', engine='c')
    vls = cross[(cross.z < 0.06) & (cross.petro_abs_mag < -19.38) & (~np.isnan(cross.B_T_n4)) & (cross.B_T_n4 > -1)]
    # B_T distribution
    vls = vls[~vls.gz2class.str.contains('B')]
    sns.distplot(vls[vls.petro_abs_mag > -20].B_T_n4, bins=np.linspace(0, 1, 11), rug=False, kde=False, norm_hist=False,
                 hist=True, hist_kws={'histtype': 'step', 'linewidth': 3}, color='blue').set(xlim=(-0.04, 1.04))
    sns.distplot(vls[(vls.petro_abs_mag > -21) & (vls.petro_abs_mag < -20)].B_T_n4, bins=np.linspace(0, 1, 11), rug=False, kde=False, norm_hist=False,
                 hist=True, hist_kws={'histtype': 'step', 'linewidth': 3}, color='red').set(xlim=(-0.04, 1.04))
    sns.distplot(vls[vls.petro_abs_mag < -21].B_T_n4, bins=np.linspace(0, 1, 11), rug=False, kde=False, norm_hist=False,
                 hist=True, hist_kws={'histtype': 'step', 'linewidth': 3}, color='green').set(xlim=(-0.04, 1.04))
    plt.legend(('-20 < Mr < -19.38', '-21 < Mr < -20', 'Mr < -21'))


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    pd.set_option('display.width', 200)
    flatui = ["#95a5a6", "#3498db", "#9b59b6", "#e74c3c", "#34495e", "#2ecc71", "#1de9b6", "#827717"]
    galaxy_types = {
        -1: 'Unclassifiable',
        0: 'Not used',
        1: 'SF',
        2: 'low S/N SF',
        3: 'Composite',
        4: 'AGN non-Liner',
        5: 'Low S/N Liner'
    }
    # creat_cross_sample()
    control_sample()
    # plot_distribution()
    sns.plt.show()
