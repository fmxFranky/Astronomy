import pandas as pd
from scipy import interpolate
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import seaborn as sns
import matplotlib.pyplot as plt
import warnings


def add_central_flag_from_group(agn_type):
    """
    在 sample 中加入是否是中心星系的判断标示并且精简数据,重新生成的文件的列表里面只剩下match需要用到的列,
    这个时候文件里面只有 flux_oiii并没有 lum_oiii
    :param agn_type
    """
    initial_agn_type = pd.read_table('sdssDR4_%s_galaxy.txt' % agn_type, engine='python', sep='\s+')
    initial_agn_type.columns = initial_agn_type.columns.str.lower()
    initial_agn_type.rename(columns={'logf_5007': 'flux_oiii', 'col1': 'galaxy_id', 'col9': 'petro_abs_mag'}, inplace=True)

    group = pd.read_table('imodelC_1', engine='python', sep='\s+', header=None,
                          names=['galaxy_id', 'other_id', 'group_id', 'brightest', 'most_massive'],
                          index_col='galaxy_id')
    # 确定每个星系是中心星系还是伴星系
    initial_agn_type['brightest'] = group.ix[initial_agn_type.galaxy_id, 'brightest'].values
    initial_agn_type['most_massive'] = group.ix[initial_agn_type.galaxy_id, 'brightest'].values
    initial_agn_type['central'] = 0
    initial_agn_type['central'][(initial_agn_type.brightest == 1) & (initial_agn_type.most_massive == 1)] = 1
    # 在我们的样本中选取最亮的当做中心星系和选取质量最大的当做中心星系是没有区别的
    # print(len(initial_agn_type[initial_agn_type.most_massive == 1]))
    # print(len(initial_agn_type[initial_agn_type.brightest == 1]))
    # print(len(initial_agn_type[(initial_agn_type.brightest == 1) & (initial_agn_type.most_massive == 1)]))

    initial_agn_type.to_csv('initial_%s_galaxy.csv' % agn_type, index=None,
                            columns=['name', 'ra', 'dec', 'z', 'flux_oiii', 'petro_abs_mag', 'central'])


def match_initial_data(agn_type):
    """
    将原始数据和 GZ2的 Full-Catalog 进行匹配,在重新生成另一个 catalog,
    同时将 flux 替换为 lum_oiii,并且加入这个星系在 GZ2 Full-Catalog 里面的 id
    type1匹配率: 43%
    type2匹配率: 53%
    :param agn_type:
    :return:
    """
    initial_agn_type = pd.read_csv('initial_%s_galaxy.csv' % agn_type)
    idx, sep2d, dist3d = SkyCoord(ra=initial_agn_type.ra, dec=initial_agn_type.dec, unit=u.degree, frame='icrs').match_to_catalog_sky(
        SkyCoord(ra=zoo.ra, dec=zoo.dec, unit=u.degree, frame='icrs'), nthneighbor=1)
    matched_index = np.where(sep2d.arcsec < 1)[0]
    matched_agn_type = initial_agn_type.ix[matched_index, :]
    matched_agn_type['gz2id'] = idx[matched_index]
    flux2luminous(matched_agn_type).to_csv('matched_%s_galaxy.csv' % agn_type, index=None,
                                           columns=['name', 'ra', 'dec', 'z', 'lum_oiii', 'petro_abs_mag', 'central', 'gz2id'])


def control_sample(match_num=5):
    """
    在原始的匹配到 GZ2 的两个 sample 里面进行条件 control,最后把结果存储到文件type12_controlled_initial. csv
    控制的条件(Jnac paper):
        |∆Mr| < 0.1
        |∆z| < 0.01
        5 closest Lum[OIII]
    :param match_num:
    :return:
    """
    matched_type1 = pd.read_csv('matched_type1_galaxy.csv')
    matched_type2 = pd.read_csv('matched_type2_galaxy.csv')

    ctrl = pd.DataFrame(columns=['name1', 'ra1', 'dec1', 'z1', 'lum_oiii1', 'petro_abs_mag1', 'central1', 'gz2id1',
                                 'name2', 'ra2', 'dec2', 'z2', 'lum_oiii2', 'petro_abs_mag2', 'central2', 'gz2id2'])

    for i in matched_type1.index:
        agn = matched_type1.ix[i]
        cs = matched_type2[(abs(matched_type2.petro_abs_mag - agn.petro_abs_mag) < 0.1) &
                           (abs(matched_type2.z - agn.z) < 0.01) &
                           (matched_type2.central == agn.central)]
        if len(cs) >= match_num:
            cs['diff_lum'] = abs(cs.lum_oiii - agn.lum_oiii)
            cs.sort_values(inplace=True, by='diff_lum')
            cs.drop('diff_lum', 1, inplace=True)
            for j in cs.index[:match_num]:
                ctrl = ctrl.append(pd.DataFrame([np.concatenate((agn.values, cs.ix[j].values), axis=0)],
                                                columns=ctrl.columns), ignore_index=True)

    ctrl.to_csv('type12_controlled_initial.csv', index=None)


def create_zfl(agn_type):
    """
    将 type1/2原始数据里面的  redshift,luminous,flux 信息单独提取出来
    后面会把所有文件里面的 flux 替换为 luminous,所以这个函数不会用了
    :param agn_type:
    :return:
    """
    name = 'NAME1' if agn_type == 'type1' else 'NAME2'
    z = 'Z1' if agn_type == 'type1' else 'Z2'
    lum = 'LOGL1_5007' if agn_type == 'type1' else 'LOGL2_5007'
    old_matched_type12 = pd.read_csv('old_matched_type12.csv').drop_duplicates(name).sort_values(by=z)
    initial_agn_type = pd.read_csv('initial_%s_galaxy.csv' % agn_type, index_col='name')
    mat = initial_agn_type.ix[old_matched_type12[name]]
    zfl = pd.DataFrame(columns=['z', 'flux', 'luminous'])
    zfl.z = mat.z
    zfl.flux = mat.flux_oiii.astype(float)
    zfl.luminous = old_matched_type12[lum].values.astype(float)
    zfl.to_csv('zfl_%s_information.csv' % agn_type)


def flux2luminous(sample):
    """
    通过插值将 flux 转换成 luminous
    :param: sample
    :return:
    """
    zfl_type1 = pd.read_csv('zfl_type1_information.csv')
    zfl_type2 = pd.read_csv('zfl_type2_information.csv')
    # # 这几行代码是检查一下用 type1/2 做出的差值有无差异,事实证明没有差异
    # plt.plot(zfl_type1.z,np.sqrt(np.power(10,zfl_type1.luminous-zfl_type1.flux)),'r')
    # plt.plot(zfl_type2.z,np.sqrt(np.power(10,zfl_type2.luminous-zfl_type2.flux)),'b')
    zfl = zfl_type1.append(zfl_type2, ignore_index=True).drop_duplicates('z').sort_values(by='z')
    itp = interpolate.UnivariateSpline(zfl.z, np.sqrt(np.power(10, zfl.luminous - zfl.flux)), bbox=[0, 0.25], k=1)
    df = sample.sort_values(by='z')
    df['lum_oiii'] = np.log10(np.square(itp(df.z))) + df.flux_oiii
    return df


def initialize_gz2():
    """
    预处理,在 Main Sample 加上红移和绝对星等
    :return:
    """
    uni = pd.read_csv('gz2sample.csv', header=0, index_col='OBJID')
    zoo['Mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_MR.values
    zoo['mr'] = uni.ix[zoo.dr7objid.values].PETROMAG_R.values
    zoo['z'] = uni.ix[zoo.dr7objid.values].REDSHIFT.values


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    pd.set_option('display.width', 400)

    # zoo 读取不要修改 index,用 range(len(zoo))
    zoo = pd.read_csv('zoo2MainSpecz.csv', header=0)
    initialize_gz2()
    for t in ['type1', 'type2']:
        match_initial_data(t)
    control_sample(5)

    plt.show()
