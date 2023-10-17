import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

site = 'TYO'
start_year = 1981
end_year = 2022

path = f'./../../out/{site}'

# read csv
files = path + '/{year}/daily/daily_{year}.csv'
dfs = []

for year in range(start_year, end_year+1):
    df = pd.read_csv(files.format(year=year), header=None)
    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

plt.figure(figsize=(14, 7))
plt.plot(df[2], linewidth=0.5)

# axis setting
plt.xlabel('Day of year')
plt.ylabel('GPP (gCm$^{-2}$dy$^{-1}$))')


plt.ylim(bottom=0)
plt.xlim(left=0)

plt.show()