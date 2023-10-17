import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sum_gpp = [0]*365

site = 'TYO'
start_year = 2001
end_year = 2022

path = f'./../../out/{site}'

# Create the figure and axes outside the loop
fig, axs = plt.subplots(2, 1, figsize=(14, 14), sharex=True) # 2行1列のレイアウトと共有X軸
fig.subplots_adjust(hspace = 0)  # マージンをなくす

for year in range(start_year, end_year+1):
    path = f'./../../out/{site}/{year}/'
    
    # read csv
    df_gpp = pd.read_csv(path + 'daily/' + f'daily_{year}.csv', header=None)
    df_lai = pd.read_csv(path + 'daily/' + f'allocation_{year}.csv', header=None)
    
    axs[0].plot(df_lai[0], df_lai[1], color='yellowgreen', linewidth=0.5, label=f'{year}')
    axs[0].set_ylabel('LAI (m$^2$m$^{-2}$)')
    
    axs[1].plot(df_gpp[0], df_gpp[2], color='yellowgreen', linewidth=0.5, label=f'{year}')
    axs[1].set_ylabel('GPP (gCm$^{-2}$dy$^{-1}$))')

# axis setting
p = np.arange(0, 365, 365/12)
l = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nob', 'Dec']
plt.xticks(p, l)
plt.xlim(0, 365)

# Labeling and presentation
plt.xlabel('Day of year')
plt.ylabel('GPP (gCm$^{-2}$dy$^{-1}$))')
plt.grid(True)
plt.tight_layout()
plt.legend() # Add a legend to distinguish between years
plt.show()
