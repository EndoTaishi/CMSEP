import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

site = 'KWG'
year = 2000
path = f'./../out/{site}/{year}_wind11/'

# CSVファイルの読み込み
df_clim = pd.read_csv(path + 'meanclim/' + f'meanclim_{year}.csv', header=None)
df_lai = pd.read_csv(path + 'daily/' + f'allocation_{year}.csv', header=None)
df_gpp = pd.read_csv(path + 'daily/' + f'daily_{year}.csv', header=None)

# グラフの作成
fig, axs = plt.subplots(4, 1, figsize=(14, 14), sharex=True) # 4行1列のレイアウトと共有X軸
fig.subplots_adjust(hspace = 0)  # マージンをなくす

# 第1グラフ: DOY vs T_a_C_mean, R_s_sum
ax1 = axs[0].twinx()
axs[0].plot(df_clim[0], df_clim[1], color='black')
axs[0].set_ylabel('Air Temperature (℃)')
ax1.plot(df_clim[0], df_clim[2], color='gray', linestyle='dashed')
ax1.set_ylabel('Sunshine (MJm$^{-2}$dy$^{-1}$)')

# 第2グラフ: DOY vs rainfall_sum, W_mean
ax2 = axs[1].twinx()
axs[1].bar(df_clim[0], df_clim[4], color='black')
axs[1].set_ylabel('Rainfall (mmdy$^{-1}$)')
ax2.plot(df_clim[0], df_clim[7], color='black')
ax2.set_ylabel('Soil Water (mm)')
ax2.set_yticks(np.arange(300, 450, 50))

# 第3グラフ: DOY vs LAI
axs[2].plot(df_lai[0], df_lai[1], color='black')
axs[2].set_ylabel('LAI (m$^2$m$^{-2}$)')

# 第4グラフ: DOY vs GPP
axs[3].plot(df_gpp[0], df_gpp[1], color='black')
axs[3].set_ylabel('GPP (gCm$^{-2}$dy$^{-1}$))')

# 最後のグラフのX軸に名前を付ける
axs[3].set_xlabel('Day of year')  # X軸の名前

# X軸の範囲を指定
xmin = df_gpp[0].min()  # X軸の最小値
xmax = df_gpp[0].max()  # X軸の最大値
for ax in axs:
    ax.set_xlim(xmin, xmax)

# グラフを表示
plt.show()
