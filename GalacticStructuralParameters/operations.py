import subprocess

import astropy.io.fits as ft
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skimage.exposure as eps
import skimage.morphology as mph
from skimage.restoration import inpaint
from termcolor import cprint as tc


def load(tp, obj):
    # load necessary data
    init_data = np.copy(ft.open('/Users/franky/Desktop/type%s/image/%s_r.fits' % (tp, obj))[0].data)
    seg_data = np.copy(ft.open('/Users/franky/Desktop/type%s/seg/%s_r.fits' % (tp, obj))[0].data)
    cat_data = pd.read_table('/Users/franky/Desktop/type%s/catalog/%s_r.txt' % (tp, obj), sep='\s+', header=None, names=['mag', 'x', 'y', 'fi', 'fp'])
    # find center
    y, x = init_data.shape[0] / 2, init_data.shape[1] / 2
    offset = np.argmax(init_data[y - 5:y + 6, x - 5:x + 6])
    return y + (offset // 11 - 5), x + (offset % 11 - 5), init_data, seg_data, cat_data


def norm(arr):
    return (arr - arr.min()) / (arr.max() - arr.min())


def tr(arr, y, x, r):
    return np.copy(arr[y - r:y + r + 1, x - r:x + r + 1])


def log(gal, exponent=1000):
    return np.log10(norm(gal) * exponent + 1) / np.log10(exponent)


def adapt_hist(gal):
    return eps.equalize_adapthist(norm(gal), nbins=2 ** 14, kernel_size=2 ** 7, clip_limit=0.01)


def truncate(init_data, seg_data, cat_data, py, px):
    gal_flag = seg_data[py][px]
    one = np.ones_like(seg_data, dtype=bool)
    zero = np.zeros_like(seg_data, dtype=bool)
    pollutoin = np.zeros_like(seg_data)
    for i in np.arange(py-150, py+151):
        for j in np.arange(px-150, px+151):
            s = seg_data[i][j]
            if s and s != gal_flag:
                pollutoin[i][j] = 1 if cat_data.ix[s - 1].fi / cat_data.ix[s - 1].fp < 2.5 else 0
    belong_to_galaxy = mph.remove_small_objects(np.where((abs(seg_data - gal_flag) < 2) & (pollutoin == 0) & (seg_data != 0), one, zero))
    y, x = np.where(belong_to_galaxy)
    dist = np.sqrt((y - py) ** 2 + (x - px) ** 2)
    radius = min(150, np.max(dist))
    return tr(init_data, py, px, radius), tr(seg_data, py, px, radius), tr(belong_to_galaxy, py, px, radius), tr(pollutoin, py, px, radius)


def display():
    ds9 = '/Applications/SAOImage\ DS9.app/Contents/MacOS/ds9'
    cmd = '-log -zoom to fit -minmax -invert -cmap value 1.75 0.5'
    ct = pd.read_csv('list.csv')
    for i in range(16):
        objs = ''
        for j in range(i * 50, min((i + 1) * 50, 784)):
            objs += '/Users/franky/Desktop/type2/image/%s_r.fits ' % ct.ix[j]['NAME2']
        subprocess.Popen('%s %s %s' % (ds9, cmd, objs), executable='/bin/zsh', shell=True)
    return


def seg2mask(segmentation, pollution):
    sy, sx = segmentation.shape

    flag = segmentation[sy / 2][sx / 2]
    one = np.ones_like(segmentation)
    zero = np.zeros_like(segmentation)
    return np.where(((pollution == 0) & (abs(segmentation - flag) < 2)) | (segmentation == 0), zero, one)


def inpainting(image, mask_map):
    ex = norm(image)
    ex[np.where(mask_map)] = 0
    return inpaint.inpaint_biharmonic(ex, mask_map)


def cal_ay(image, belong_to_galaxy):
    zero = np.zeros_like(image)
    I = np.where(belong_to_galaxy, image, zero)
    I180 = np.rot90(I, 2)
    return np.sum(np.abs(I - I180)) * 0.5 / np.sum(np.abs(I))


if __name__ == '__main__':
    fits = pd.read_csv('list.csv')
    pro = 0
    for fit in fits[fits.Z1 < 0.05]['NAME2'].values:
        print(fit, end='    ')
        cy, cx, init, seg, cat = load(2, fit)
        img, sg, btg, poll = truncate(init, seg, cat, cy, cx)
        inpaint_img = inpainting(img, seg2mask(sg, poll))
        plt.imsave('/Users/franky/Desktop/templates/init/%s.jpg' % fit, log(inpaint_img), cmap='gray', dpi=200)
        fits.at[pro, 'INIT'] = cal_ay(inpaint_img, btg)
        print(pro)
        pro += 1
    fits.to_csv('init2.csv', index=None, columns=['NAME2', 'INIT'])
    # plt.imshow(-inpainting(norm(img), seg2mask(sg)), cmap='gray')
    # plt.show()
