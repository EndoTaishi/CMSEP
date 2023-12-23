import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

site = 'TYO'
start_year = 2000
end_year = 2023

path = f'./../../out/{site}'

# read yearly NPP data
df_npp = pd.read_csv(path + '/yearly/yearly.csv', header=None)

# initialize lists to store data
temp_data = []
sunlight_data = []
humidity_data = []
prec_data = []
npp_data = []

# prepare the figure for 2x2 subplots
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
axs = axs.flatten()  # flatten the 2D array of axes

for year in range(start_year, end_year+1):
    df_clim = pd.read_csv(f'{path}/{year}/meanclim/meanclim_{year}.csv', header=None)

    # Extract data
    yearly_temp = df_clim[1].mean()  # average temperature
    yearly_sunlight = df_clim[2].mean()  # average total solar radiation
    yearly_humidity = df_clim[5].mean()  # average relative humidity
    yearly_prec = df_clim[4].sum()  # total precipitation
    yearly_npp = (df_npp[1][year-start_year]-df_npp[2][year-start_year])/10  # calculated NPP

    # Store data for plotting and regression
    temp_data.append(yearly_temp)
    sunlight_data.append(yearly_sunlight)
    humidity_data.append(yearly_humidity)
    prec_data.append(yearly_prec)
    npp_data.append(yearly_npp)

    # Plot each year's data with text
    axs[0].plot(yearly_temp, yearly_npp, 'o', color='red')
    axs[0].text(yearly_temp, yearly_npp, str(year), fontsize=8)

    axs[1].plot(yearly_sunlight, yearly_npp, 'o', color='green')
    axs[1].text(yearly_sunlight, yearly_npp, str(year), fontsize=8)

    axs[2].plot(yearly_humidity, yearly_npp, 'o', color='blue')
    axs[2].text(yearly_humidity, yearly_npp, str(year), fontsize=8)

    axs[3].plot(yearly_prec, yearly_npp, 'o', color='purple')
    axs[3].text(yearly_prec, yearly_npp, str(year), fontsize=8)

# Compute linear regression for each climate factor and plot
def plot_regression(ax, x_data, y_data, color, xlabel):
    coefficients = np.polyfit(x_data, y_data, 1)
    regression_line = np.poly1d(coefficients)
    x_reg = np.linspace(min(x_data), max(x_data), 100)
    ax.plot(x_reg, regression_line(x_reg), color=color, linestyle='dashed')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('NPP (kgCm$^{-2}$yr$^{-1}$)')

plot_regression(axs[0], temp_data, npp_data, 'red', 'Average Temperature (°C)')
plot_regression(axs[1], sunlight_data, npp_data, 'green', 'Average Total Solar Radiation (MJ/m²)')
plot_regression(axs[2], humidity_data, npp_data, 'blue', 'Average Relative Humidity (%)')
plot_regression(axs[3], prec_data, npp_data, 'purple', 'Total Precipitation (mm)')

# Set titles for each subplot
axs[0].set_title('Temperature vs NPP')
axs[1].set_title('Sunlight vs NPP')
axs[2].set_title('Humidity vs NPP')
axs[3].set_title('Precipitation vs NPP')

plt.tight_layout()
plt.show()
