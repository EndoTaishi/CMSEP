import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

site = 'TYO'
start_year = 1981
end_year = 2022

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

def is_leap_year(year):
    if (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0):
        return True
    return False

# データの日数を格納するリスト
days = []

# 年ごとの日数を追加（閏年の場合は366日、それ以外は365日）
for year in range(start_year, end_year+1):
    if is_leap_year(year):
        days.append(366)
    else:
        days.append(365)

# np.cumsumで累積和を取ることで、年ごとの日数の累積を取得
# ただし、最初の値は0にするために、先頭に0を追加
p = np.insert(np.cumsum(days), 0, 0)

# 年のリストを作成
l = [str(year) for year in range(start_year, end_year+1)]

p = p[:-1]
plt.xticks(p, l)

plt.show()