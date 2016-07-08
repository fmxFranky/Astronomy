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
    return eps.equalize_adapthist(norm(gal), nbins=512, kernel_size=32, clip_limit=0.1)


def truncate(init_data, seg_data, cat_data, py, px):
    gal_flag = seg_data[py][px]
    real_galaxy = mph.remove_small_objects(np.where(seg_data == gal_flag, np.ones_like(seg_data, dtype=bool), np.zeros_like(seg_data, dtype=bool)),
                                           connectivity=2, min_size=100)
    y, x = np.where(real_galaxy)
    radius = np.max(np.sqrt((y - py) ** 2 + (x - px) ** 2))
    init_tr = tr(init_data, py, px, radius)
    seg_tr = tr(seg_data, py, px, radius)
    poll_tr = np.zeros_like(seg_tr)

    def is_pollution(p):
        return True if cat_data.ix[p - 1].fi / cat_data.ix[p - 1].fp < 2.5 or abs(p-gal_flag) > 1 else False
    for i in np.arange(poll_tr.shape[0]):
        for j in np.arange(poll_tr.shape[1]):
            s = seg_tr[i][j]
            if s and s != gal_flag:
                poll_tr[i][j] = poll_tr[-i-1][-j-1] = 1 if is_pollution(s) else 0
    belong_to_galaxy = mph.remove_small_objects(np.where((abs(seg_tr - gal_flag) < 2) & (poll_tr == 0) & (seg_tr != 0), np.ones_like(seg_tr, dtype=bool), np.zeros_like(seg_tr, dtype=bool)),
                                                connectivity=2, min_size=100)
    # plt.imshow(belong_to_galaxy)
    # plt.show()
    return init_tr, seg_tr, belong_to_galaxy, poll_tr


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
    plt.imsave('/Users/franky/Desktop/templates/init/1/%s.jpg' % fit, -log(I), cmap='gray', dpi=400)
    I180 = np.rot90(I, 2)
    return np.sum(np.abs(I - I180)) * 0.5 / np.sum(np.abs(I))


if __name__ == '__main__':
    fits = pd.read_csv('list.csv')
    pro = 0
    for fit in fits[fits.Z1 < 0.05]['NAME1'].values[:1]:
        print(fit, end='    ')
        # fit = 'J114433.07+613200.7'
        cy, cx, init, seg, cat = load(1, fit)
        img, sg, btg, poll = truncate(init, seg, cat, cy, cx)
        # inpaint_img = inpainting(img, seg2mask(sg, poll))
        # plt.imsave('/Users/franky/Desktop/templates/init/%s.jpg' % fit, -adapt_hist(img), cmap='gray', dpi=200)
        fits.at[pro, 'A'] = cal_ay(img, btg)
        print(cal_ay(img, btg), end='   ')
        print(pro)
        pro += 1
        fits.to_csv('init1.csv', index=None, columns=['NAME1', 'A'])
