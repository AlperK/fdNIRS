import matplotlib.pyplot as plt
import numpy as np

v_in = np.linspace(start=50, stop=840, num=80)
amp_1 = np.array([219, 251, 284, 318, 350,
                  383, 415, 448, 479, 511,
                  545, 579, 612, 644, 675,
                  708, 749, 770, 802, 831,
                  865, 896, 925, 954, 985,
                  1015, 1048, 1080, 1107, 1136,
                  1166, 1192, 1219, 1244, 1273,
                  1298, 1323, 1345, 1373, 1395,
                  1422, 1446, 1467, 1488, 1513,
                  1532, 1553, 1572, 1587, 1609,
                  1624, 1642, 1662, 1678, 1691,
                  1704, 1717, 1731, 1738, 1747,
                  1791, 1808, 1819, 1825, 1836,
                  1844, 1852, 1861, 1866, 1870,
                  1879, 1889, 1892, 1901, 1908,
                  1915, 1918, 1924, 1921, 1923])
amp_1_linear = v_in * 3.25 + 56.73
amp_2 = np.array([254, 295, 337, 381, 425,
                 454, 495, 536, 577, 618,
                 651, 691, 733, 775, 812,
                 848, 878, 917, 956, 990,
                 1019, 1055, 1089, 1127,1162,
                 1198, 1242, 1277, 1309, 1344,
                 1377, 1411, 1443, 1475, 1505,
                 1522, 1554, 1579, 1612, 1641,
                 1662, 1689, 1710, 1745, 1760,
                 1788, 1814, 1835, 1864, 1883,
                 1905, 1928, 1954, 1961, 1992,
                 2006, 2030, 2046, 2062, 2099,
                 2120, 2127, 2143, 2161, 2161,
                 2166, 2175, 2188, 2199, 2214,
                 2232, 2239, 2259, 2263, 2277,
                 2299, 2313, 2327, 2336, 2342])

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True)
fig.suptitle('Input vs Output Amplitudes with Linear Fit\nwith 60dB attenuation on the Input')

axes[0].scatter(x=v_in, y=amp_1, s=40, label='Channel-1', color='r')
m, b = np.polyfit(v_in[:10], amp_1[:10], 1)
axes[0].plot(v_in, m * v_in + b, '--r',
             linewidth=2, alpha=0.5,
             label=f'Channel-1 Fit\nSlope:{round(m, 4)}')
axes[0].set_ylabel('Output (mV)')

axes[1].scatter(x=v_in, y=amp_2, s=40, label='Channel-2', color='b')
m, b = np.polyfit(v_in[:10], amp_2[:10], 1)
print(m, b)
axes[1].plot(v_in, m * v_in + b, '--b',
             linewidth=2, alpha=0.5,
             label=f'Channel-2 Fit\nSlope:{round(m, 4)}')
axes[1].set_xlabel('Input (mV)')
axes[1].set_ylabel('Output (mV)')
axes[0].legend()
axes[1].legend()
axes[0].grid()
axes[1].grid()
fig.tight_layout()


closest = min(amp_1, key=lambda x: abs(x-1552))
index_closest = np.where(amp_1 == closest)[0]
print(amp_1_linear[index_closest])
plt.show()
