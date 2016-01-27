#! coding:utf-8
"""
distance_attenuation.py
鏡像音源による音の干渉をシミュレート

Created by 0160929 on 2016/01/27 17:01
"""
__version__ = '0.0'

import os
import sys
import numpy as np


def sound_speed(temp=20., stric=True):
    """
    :param temp: 摂氏{-273.15 < tmp < ∞}
    :return: 音速[ms]
    """
    if stric:
        # 厳密
        c = 20.055 * (float(temp) + 273.15) ** (1. / 2.)
    else:
        # テイラー展開一時近時
        c = 331.5 + 0.6 * float(temp)
    return c


def distance_TF(h1=0.1, h2=1.8, L=2., f_min=1, f_max=24000, df=10, temp=20.):
    """
    点音源と壁面における鏡像音源との干渉による伝達関数を測定
    :param h1: 音源S高さ[m]
    :param h2: 受音点R高さ[m]
    :param L: 音源Sと受音点Rとの水平距離[m]
    :param f_min: 解析する最低周波数[Hz]
    :param f_max: 解析する最高周波数[Hz]
    :param df: 解析の周波数幅[Hz]
    :param temp: 温度[℃]
    :return: Pr{Conmplex}, Freq{周波数列}
    """
    # 音速
    c = sound_speed(temp, True)
    print c

    # 音源の振幅
    A = 1.
    # 音源のベクトル
    S = np.array([0, h1])
    # 虚像音源のベクトル
    Si = np.array([0, -h1])
    # 受音点位置へのベクトル
    R = np.array([L, h2])
    # 音源から受音点へのベクトル
    vec_r1 = R - S
    # 虚像音源から受音点へのベクトル
    vec_r2 = R - Si
    # 音源から受音点までの距離
    r1 = np.linalg.norm(vec_r1)
    # 虚像音源から受音点までの距離
    r2 = np.linalg.norm(vec_r2)
    # 解析する周波数幅
    f_num = f_max / df
    freq = np.linspace(f_min, f_max, f_num)

    # 結果格納用バッファ
    Pr_s = list()
    for f in freq:
        # TODO:速度改善のためnumpy型での計算に変更
        # 波数
        k = 2 * np.pi * f / c
        Pr1 = A * np.exp(1j * k * r1) / r1
        Pr2 = A * np.exp(1j * k * r2) / r2
        # 複素数
        Pr = Pr1 + Pr2
        Pr_s.append(Pr)

    return np.asarray(Pr_s), freq


def main(h1,h2,L):
    # 音源の高さ[m]
    # h1 = 0.1
    # 受音点の高さ[m]
    # h2 = 1.7
    # 音源と受音点の水平距離[m]
    # L = 0.25

    # 伝達関数の計算
    Pr_s, freq = distance_TF(h1, h2, L, f_min=1., f_max=24000., df=1., temp=10)

    # 解析結果
    abs_Pr = np.abs(Pr_s)
    ang_Pr = np.angle(Pr_s)
    real_Pr = np.real(Pr_s)
    imag_Pr = np.imag(Pr_s)

    txt_data = np.asarray((freq, abs_Pr, ang_Pr, real_Pr, imag_Pr)).transpose()
    txt_name =  './TF_h1=%0.2f[m]_h2=%0.2f[m]_L=%0.2f[m].csv' % (h1, h2, L)
    np.savetxt(txt_name, txt_data, delimiter=',',header='Frequency,Amplitude,Angle,Real,Image',comments='')



    # グラフプロット

    fig = plt.figure(1)
    mm = 4;
    nn = 1
    ax1 = fig.add_subplot(mm, nn, 1)
    ax2 = fig.add_subplot(mm, nn, 2)
    ax3 = fig.add_subplot(mm, nn, 3)
    ax4 = fig.add_subplot(mm, nn, 4)

    plt.axes(ax1)
    plt.plot(freq, abs_Pr)
    plt.ylabel('Amplitude')

    plt.axes(ax2)
    plt.plot(freq, ang_Pr)
    plt.ylabel('Angle')

    plt.axes(ax3)
    plt.plot(freq, real_Pr)
    plt.ylabel('Real')

    plt.axes(ax4)
    plt.plot(freq, imag_Pr)
    plt.ylabel('Image')
    plt.xlabel('Frequency [Hz]')

    plt.suptitle('h1=%0.2f[m], h2=%0.2f[m], L=%0.2f[m]' % (h1, h2, L))

    fig_name = './TF_h1=%0.2f[m]_h2=%0.2f[m]_L=%0.2f[m].jpg' % (h1, h2, L)
    plt.savefig(fig_name)

    # plt.show()
    plt.clf()

    pass


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    for L in [0,0.25,0.5,0.75,1.0,2.0,3.0]:
        for h1 in [0.01, 0.1, 0.2, 0.3]:
            for h2  in [0.3, 1.0, 1.8, 2.0]:
                main(h1,h2,L)
