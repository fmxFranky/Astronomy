import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.io.fits as ft
import skimage.exposure as eps
import skimage.data as skd
import pandas as pd
import numpy as np
import skimage.morphology as mph
import termcolor


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


def log(gal, exponent=1000):
    return norm(np.log10(norm(gal)*exponent+1)/np.log10(exponent))


def adapt_hist(gal):
    return eps.equalize_adapthist(norm(gal), nbins=2**14, kernel_size=2**7, clip_limit=0.01)


# TODO: select correct region including the galaxy
def truncate():
    return

if __name__ == '__main__':
    cy, cx, init, seg, cat = load(1, 'J083732.70+284218.7')
    plt.imshow(-adapt_hist(log(init)), cmap='gray')
    plt.show()


