import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings
import matplotlib_venn as mv


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
    cross_gz2_group_n4bdd = cross[['galaxyid', 'gz2id', 'z', 'petro_abs_mag', 'gz2class', 'bar_debiased_fraction',
                                   'B_T_n4', 'typ', 'Re_n4', 'be_n4', 'Rd_n4', 'dia_n4']]
    cross_gz2_group_n4bdd.to_csv('cross_gz2_group_n4bdd.csv', index=None)


def control_sample():
    cross = pd.read_csv('cross_gz2_group_n4bdd.csv', engine='c')
    cross = cross[(~np.isnan(cross.B_T_n4)) & (cross.B_T_n4 > -1) &
                  (cross.z < 0.06) & (cross.petro_abs_mag < -19.38)]
    bar_sample = cross[(cross.bar_debiased_fraction > 0.5) & (cross.gz2class.str.contains('B'))]

    # plt.plot(np.arange(0.02, 0.06, 0.005), [np.mean(bar_sample[(bar_sample.z < z+0.005) & (bar_sample.z > z)]['B_T_n4']) for z in np.arange(0.02, 0.06, 0.005)], 'kv--')
    # unbar_sample = cross[(cross.bar_debiased_fraction < 0.3) & (~cross.gz2class.str.contains('B'))]
    # print(len(bar_sample), len(unbar_sample), len(cross))
    # create controll sample
    # k = 0
    # for i in bar_sample.index[:1]:
    #     bs = bar_sample.ix[i]
    #     cbs = unbar_sample[(abs(unbar_sample.z-bs.z) < 0.01) &
    #                        (abs(unbar_sample.petro_abs_mag - bs.petro_abs_mag) < 0.1) &
    #                        (abs(unbar_sample.B_T_n4-bs.B_T_n4) < 0.01)]
    #     if len(cbs) >= 5:
    #         k += 1
    # print(k/len(bar_sample))
    # print(control_unbar_sample)


def verify():
    gz2 = pd.read_csv('GZ2.csv', engine='c')
    n4_bdd = pd.read_csv('n4_bdd.csv', engine='c')
    # meta = pd.read_csv('gz2sample.csv', engine='c')
    # from astropy.coordinates import SkyCoord
    # from astropy import units as u
    # idx,sep2d,d3d = SkyCoord(ra=gz2.ra2, dec=gz2.dec2, unit=u.deg, frame='icrs').match_to_catalog_sky(SkyCoord(ra=n4_bdd.ra,dec=bar.dec,unit=u.deg,frame='icrs'),nthneighbor=1)


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
    # gns_ids = pd.read_table('group_galaxy_id', header=None, engine='python', sep='\s+')
    # gns_ids.to_csv('group_vagc_dr7_ids.csv', engine='c', index=None)
    # creat_bar_sample()
    # cross = pd.read_csv('cross_gz2_group_n4bdd.csv', engine='c')
    # print(cross)
    control_sample()
    sns.plt.show()
