import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

site = 'TYO'
start_year = 2001
end_year = 2022

path = f'./../../out/{site}'

# read csv
df_npp = pd.read_csv(path + '/yearly/yearly.csv', header=None)

prec_data = []
npp_data = []

plt.boxplot((df_npp[1]-df_npp[2])/10, showmeans=True, sym='')
v = np.var((df_npp[1]-df_npp[2])/10)
print(v)
print(max((df_npp[1]-df_npp[2])/10))
print(min((df_npp[1]-df_npp[2])/10))
print(np.mean((df_npp[1]-df_npp[2])/10))
# graph setting
fig, ax = plt.subplots(figsize=(14, 7))
for year in range(start_year, end_year+1):
    df_clim = pd.read_csv(f'{path}/{year}/meanclim/meanclim_{year}.csv', header=None)
    # sum of df_clim[4] is yearly sum of rainfall
    yearly_prec = df_clim[4].sum()
    yearly_npp = (df_npp[1][year-start_year]-df_npp[2][year-start_year])/10
    ax.plot(yearly_prec, yearly_npp, 'o', color='black')

    # Store data for regression
    prec_data.append(yearly_prec)
    npp_data.append(yearly_npp)

# Compute linear regression
coefficients = np.polyfit(prec_data, npp_data, 1)
regression_line = np.poly1d(coefficients)
x_reg = np.linspace(min(prec_data), max(prec_data), 100)
ax.plot(x_reg, regression_line(x_reg), color='black', linestyle='dashed')

# axis setting
ax.set_xlabel('Precipitation (mm)')
ax.set_ylabel('NPP (kgCm$^{-2}$yr$^{-1}$)')
# ax.set_xlim(0, 4500)
# ax.set_ylim(0, 3.5)

# plot
plt.show()
