import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

site = 'TYO'
start_year = 2001
end_year = 2023

path = f'./../../out/{site}'

# read csv
files = path + '/{year}/daily/allocation_{year}.csv'
dfs = []

for year in range(start_year, end_year+1):
    df = pd.read_csv(files.format(year=year), header=None)
    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

plt.figure(figsize=(14, 7))
plt.plot(df[1])

# axis setting
plt.xlabel('Day of year')
plt.ylabel('LAI (m$^2$m$^{-2}$)')


plt.ylim(bottom=0)
plt.xlim(left=0)


plt.show()