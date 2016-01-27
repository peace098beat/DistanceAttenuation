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

import time
import math

from scipy import special

erfc = special.erfc
print erfc(1)




# 時間計測用でコレータ
def time_func(func):
    """時間計測用でこれーた"""

    def decorator(*args, **kwargs):
        start = time.time()
        # funcの実行
        ret = func(*args, **kwargs)
        print '%15s was executed, it took %-010.6f sec' % (func.func_name, (time.time() - start) * 1.)
        return ret

    return decorator


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


@time_func
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
    # print 'Sound Speed',c

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
        k = 2. * np.pi * f / c
        Pr1 = A * np.exp(1j * k * r1) / r1
        Pr2 = A * np.exp(1j * k * r2) / r2
        # 複素数
        Pr = (Pr1 + Pr2)
        Pr_s.append(Pr)

    return np.asarray(Pr_s), freq


def distance_TF_ref_dumping(h1=0.1, h2=1.8, L=2., f_min=1, f_max=24000, df=10, temp=20., omega=30.):
    """
    点音源と壁面における鏡像音源との干渉による伝達関数を測定
    床面のインピーダンス差による減衰を考慮
    http://apmr.matelys.com/PropagationModels/MotionlessSkeleton/DelanyBazleyModel.html
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

    # 虚像から音源までのベクトルと、地表との角度
    # phy = np.pi/4.
    x, y = vec_r2
    phy = np.arctan2(y, x)

    # -------------------------------------
    # f = float(freq[5])
    # omega = float(30)
    # f = np.asarray(freq)
    # def cal_(f):
    # 空気中の音響インピーダンス
    rho1 = 1.2
    c1 = c
    k1 = 2. * np.pi * freq / c1
    Z1 = rho1 * c1

    # -------------------------------------
    # 地表面の音響インピーダンス(多孔質インピーダンス)
    # R2 = (rho1 * c1) * (1. + 9.08 * pow(float(f) / float(omega), -0.75))
    # X2 = (rho1 * c1) * 11.9 * pow(float(f) / omega, -0.73)
    R2 = (rho1 * c1) * (1. + 9.08 * np.power(freq / float(omega), -0.75))
    X2 = (rho1 * c1) * 11.9 * np.power(freq / float(omega), -0.73)
    Z2d = R2 + 1j * X2

    # -------------------------------------
    # 地表面の波長定数(多孔質インピーダンス)
    # alpha2 = k1 * (1 + 10.8 * pow(float(f) / omega, -0.73))
    # beta2 = k1 * (10.3 * pow(float(f) / omega, -0.59))
    alpha2 = k1 * (1 + 10.8 * np.power(freq / float(omega), -0.73))
    beta2 = k1 * (10.3 * np.power(freq / float(omega), -0.59))
    k2d = alpha2 + 1j * beta2


    # -------------------------------------
    # 平面波の反射係数
    def Rd(k1_, Z1_, k2_, Z2_, phy):
        Rd_A = np.sin(phy)
        Rd_B = Z1_ * np.sqrt(1. - ((k1_ / k2_) ** 2) * (np.cos(phy) ** 2)) / Z2_
        Rd_ = (Rd_A - Rd_B) / (Rd_A + Rd_B)
        return Rd_

    _Rd = Rd(k1, Z1, k2d, Z2d, phy)
    # print '_Rd', _Rd

    # -------------------------------------
    # wd = 1j * (2 * k1 * r2) * pow(Z1 / Z2d, 2) * (1. - pow(k1 / k2d, 2) * pow(np.cos(phy), 2)) / pow((1 - _Rd), 2)
    wd = 1j * (2. * k1 * r2) * np.power(Z1 / Z2d, 2) * (
        1. - np.power(k1 / k2d, 2) * np.power(np.cos(phy), 2)) / np.power((1 - _Rd), 2)
    # print 'wd', wd

    # -------------------------------------
    # 境界損失係数
    def F(w_):
        ret = 1 + 1j * np.sqrt(np.pi) * np.sqrt(w_) * np.exp(-1. * w_) * erfc(-1j * np.sqrt(w_))
        return ret

    # -------------------------------------
    Qd = _Rd + (1 - _Rd) * F(wd)
    # print 'Qd', Qd

    # -------------------------------------
    # 合成波の音圧(複素数)
    Pr = (1. / r1 * np.exp(1j * k1 * r1)) + Qd * (1. / r2 * np.exp(1j * k1 * r2))
    # -------------------------------------

    return Pr, freq


# @time_func
def main(h1, h2, L):
    # 音源の高さ[m]
    # h1 = 0.1
    # 受音点の高さ[m]
    # h2 = 1.7
    # 音源と受音点の水平距離[m]
    # L = 0.25

    # 伝達関数の計算
    Pr_s, freq = distance_TF(h1, h2, L, f_min=1., f_max=24000., df=1., temp=10)
    Pr_s_str, freq = distance_TF_ref_dumping(h1, h2, L, f_min=1., f_max=24000., df=1., temp=10, omega=300)
    a = Pr_s - Pr_s_str
    b = np.abs(Pr_s) - np.abs(Pr_s_str)
    print 'a', a, np.sum(a)
    print 'b', b, np.sum(b)

    # グラフ
    fig = plt.figure(1)
    plt.plot(freq, np.abs(Pr_s))
    plt.plot(freq, np.abs(Pr_s_str))
    plt.show()
    return

    # 解析結果
    abs_Pr = np.abs(Pr_s)
    ang_Pr = np.angle(Pr_s)
    real_Pr = np.real(Pr_s)
    imag_Pr = np.imag(Pr_s)

    txt_data = np.asarray((freq, abs_Pr, ang_Pr, real_Pr, imag_Pr)).transpose()
    txt_name = './TF_h1=%0.2f[m]_h2=%0.2f[m]_L=%0.2f[m].csv' % (h1, h2, L)
    np.savetxt(txt_name, txt_data, delimiter=',', header='Frequency,Amplitude,Angle,Real,Image', comments='')



    # グラフプロット
    @time_func
    def graph_plot():
        fig = plt.figure(1)
        mm = 4
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

    # graph_plot()
    pass


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    main(h1=0.1, h2=1.8, L=2)

    # for L in [0,0.25,0.5,0.75,1.0,2.0,3.0]:
    #     for h1 in [0.01, 0.1, 0.2, 0.3]:
    #         for h2  in [0.3, 1.0, 1.8, 2.0]:
    #             main(h1,h2,L)
    # for L in [0.5, 1]:
    #     for h1 in [0.1]:
    #         for h2 in [1.0, 1.8]:
    #             main(h1, h2, L)
