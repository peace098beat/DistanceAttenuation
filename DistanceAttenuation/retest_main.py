#! coding:utf-8
"""
retest_main.py
幾何音響計算プログラムの検証.
参考論文を追試．

"""
__version__ = '0.0'

import os
import sys

from distance_attenuation import *
import matplotlib.pyplot as plt


def paper_retest():
    """参考文献を追試"""

    # パラメータ
    # -------------------------
    h1 = 0.31
    h2 = 1.22
    L = 15.2
    T = 20.
    sigmas = [10., 100., 1000., 10000]
    # パラメータ
    # -------------------------
    f_min = 1
    f_max = 10000.
    df = 1.

    for sigma in sigmas:
        print '* Caluculate Transfer Respose'
        print '** Paramater : h1=%0.2f, h2=%0.2f, L=%0.2f, T=%0.1f' % (h1, h2, L, T)

        # 伝達関数の計算
        # -------------------------
        Pr_s, freq = distance_TF(h1, h2, L, f_min, f_max, df, temp=T)
        Pr_s_dump, freq = distance_TF_ref_dumping(h1, h2, L, f_min, f_max, df, temp=T, omega=sigma)
        P0_s, freq = distance_TF_non_refrect(h1, h2, L, f_min, f_max, df, temp=T)

        # グラフ
        # -------------------------
        fig = plt.figure(1)
        gdata = -1 * (20 * np.log10(np.abs(Pr_s_dump)) - 20 * np.log10(np.abs(P0_s)))
        plt.semilogx(freq, gdata)
        plt.xlim([100, 10000])

        # Anotation
        _x = 0.05
        font_size = 13
        txt = 'Hs=%0.2f[m]\nHr=%0.2f[m]\nL=%0.2f[m]\nT=%0.1f[m]\nsigma=%0.1f[Pa*s/m]' % (h1, h2, L, T, sigma)
        plt.annotate(txt, xy=(_x, 0.95), fontsize=font_size, xycoords='axes fraction',
                     horizontalalignment='left', verticalalignment='top')
        plt.ylabel('Excess attenuation P/P0 [dB]')
        plt.xlabel('Frequency [Hz]')
        plt.xlim([100, 10000])

        # 画像保存
        fig_name = './retest/TF_h1=%0.2f[m]_h2=%0.2f[m]_L=%0.2f[m]_T=%0.1f[T]_sigma=%d.jpg' % (h1, h2, L, T, sigma)
        plt.savefig(fig_name)
        # plt.show()

        # データ保存
        txt_data = np.asarray((freq, gdata)).transpose()
        txt_name = './retest/TF_h1=%0.2f[m]_h2=%0.2f[m]_L=%0.2f[m]_T=%0.1f[T]_sigma=%d.csv' % (h1, h2, L, T, sigma)
        np.savetxt(txt_name, txt_data, delimiter=',', header='Frequency, Excess attenuation [db]', comments='')


if __name__ == '__main__':
    paper_retest()
