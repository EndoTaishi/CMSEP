import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

site = 'TYO'
start_year, end_year = 2023, 2023

for year in range(start_year, end_year+1):
    #year = 2000
    path = f'./../../out/{site}/{year}/'

    # read csv
    df_clim = pd.read_csv(path + 'meanclim/' + f'meanclim_{year}.csv', header=None)
    df_lai = pd.read_csv(path + 'daily/' + f'allocation_{year}.csv', header=None)
    df_gpp = pd.read_csv(path + 'daily/' + f'daily_{year}.csv', header=None)

    # graph setting
    fig, axs = plt.subplots(4, 1, figsize=(14, 14), sharex=True) # 4行1列のレイアウトと共有X軸
    fig.subplots_adjust(hspace = 0)  # マージンをなくす

    # DOY vs Air Temperature, Sunshine
    ax1 = axs[0].twinx()
    axs[0].plot(df_clim[0], df_clim[1], color='black')
    axs[0].set_ylabel('Air Temperature (℃)')
    ax1.plot(df_clim[0], df_clim[2], color='gray', linestyle='dashed')
    ax1.set_ylabel('Sunshine (MJm$^{-2}$dy$^{-1}$)')

    # DOY vs rainfall_sum, W_mean
    ax2 = axs[1].twinx()
    axs[1].bar(df_clim[0], df_clim[4], color='black')
    axs[1].set_ylabel('Rainfall (mmdy$^{-1}$)')
    ax2.plot(df_clim[0], df_clim[7], color='black')
    ax2.set_ylabel('Soil Water (mm)')
    ax2.set_yticks(np.arange(300, 450, 50))

    # DOY vs LAI
    axs[2].plot(df_lai[0], df_lai[1], color='black')
    axs[2].set_ylabel('LAI (m$^2$m$^{-2}$)')

    # DOY vs GPP
    ax3 = axs[3].twinx()
    axs[3].plot(df_gpp[0], df_gpp[2], color='black')
    axs[3].set_ylabel('GPP (gCm$^{-2}$dy$^{-1}$))')
    ax3.plot(df_gpp[0], df_gpp[1], color='pink')
    ax3.set_ylabel('GPP obs (gCm$^{-2}$dy$^{-1}$))')

    # axis setting
    axs[3].set_xlabel('Day of year')  # X軸の名前

    p = np.arange(0, 365, 365/12)
    l = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nob', 'Dec']
    plt.xticks(p, l)
    plt.xlim(0, 365)

    # plot
    #plt.show()

    # save as pdf
    fig.savefig(path + f'PlantGrowthPlot_{year}.pdf', bbox_inches='tight')
    plt.close(fig)