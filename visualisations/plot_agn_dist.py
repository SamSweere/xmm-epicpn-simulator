import matplotlib.pyplot as plt
import numpy as np

# S [erg / cm ** 2 / s]
# N( > S) [deg ** -2]
from simput.agn_distribution import get_S_N_from_file

S, N = get_S_N_from_file()

plt.plot(S, N, 'g', label="data")


plt.xlabel("S [erg / cm ** 2 / s]")
plt.ylabel("N( > S) [deg ** -2]")
plt.xscale('log')
plt.yscale('log')
# plt.legend()
# plt.savefig("xray_agn_number_count.pdf")
plt.show()


# C = N/np.sum(N) # Chance
#
# np.interp(1.1e-12, S, C)
#
# points = []
#
# while len(points) < 100:
#     x = np.random.uniform(np.min(S), np.max(S))
#     chance = np.interp(x, S, C)
#     if np.random.rand() < chance:
#         points.append((x, np.interp(x, S, N)))
#
# print(points)

#
# def func(x, a, b, c):
#      res = a * (x**b) + c
#      return res
#
#
# def inv_func(x, a, b, c):
#     res = ((x-c)/a)**(1/b)
#     return res
#
#
# def round_to_n(x, n):
#     " Round x to n significant figures "
#     return round(x, -int(np.floor(np.sign(x) * np.log10(abs(x)))) + n)
#
# def str_fmt(x, n=2):
#     " Format x into nice Latex rounding to n"
#     power = int(np.log10(round_to_n(x, 0)))
#     f_SF = round_to_n(x, n) * pow(10, -power)
#     return r"${}\cdot 10^{}$".format(f_SF, power)
#
# # init_vals = [1484.3392378755964, 977035624376729.8, 18.05549012693649]
# init_vals = [1, -1, 1]
# init_vals = [-0.0001, 10, -0.0001]
# # popt, pcov = curve_fit(f=func, xdata=N, ydata=S, p0=init_vals, maxfev=100000)
# popt, pcov = curve_fit(f=inv_func, xdata=S, ydata=N, p0=init_vals)
# print(popt)
# a, b, c = popt
#
#
#
#
# print(f"Fitted formula: S = {a} * N^({b}) + {c}")
# print(f"Inverse of formula: N = ((S-{c})/{a})**(1/{b})")
#
#
#
# # plt.plot(N, func(N, *popt), 'g--', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
# # plt.plot(N, func(N, *popt), 'g--', label=f'fit: a={popt[0]:.3e}, b={popt[1]:.3e}, c={popt[2]:.3e}')
# # plt.plot(N, func(N, *popt), 'r--', label=f'fit: S = {popt[0]:.3e} * N^({popt[1]:.3e}) + {popt[2]:.3e}')
# plt.plot(S, inv_func(S, *popt), 'r--', label=f'fit: N = ((S-{popt[2]:.3e})/{popt[0]:.3e})**(1/{popt[1]:.3e})')
# # plt.plot(N, func(N, *popt), 'g--', label=f'fit: a={str_fmt(, n=2)}, b={popt[1]}, c={popt[2]}')
# # from random import random
# # x_min = 1.0000e-15
# alpha = 1/b
# r = random()
# # x_smp = x_min * (1 - r) ** (-1 / (alpha - 1))
#
# # def fitted_func(x):
# #     return 9.213261619713252e-14 * (x ** (-0.661765007798531)) + 1.844972267368158e-15
#
# points = []
#
# for i in range(100):
#     r = random()
#     x_smp = x_min * (1 - r) ** (-1 / (alpha - 1))
#     point = (x_smp, inv_func(x_smp, *popt))
#     print(point)
#     points.append(point)
# #
# plt.scatter(x=list((x[0] for x in points)), y=list((x[1] for x in points)), label="points")


# plt.ylabel("S [erg / cm ** 2 / s]")
# plt.xlabel("N( > S) [deg ** -2]")
# plt.xscale('log')
# plt.yscale('log')
# plt.legend()
# # plt.savefig("xray_agn_number_count.pdf")
# plt.show()