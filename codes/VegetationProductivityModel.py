#!/usr/bin/env python3
import math
import os
import csv
import time
import numpy as np
import mypkg.functions as func

start_time = time.time()

site = 'TYO'

if site == 'TYO':
    lat, lon = 35.415, 139.450
    t_offset_GMT = 9  # Site local time offset from UTC (GMT) (hours)
    elv = 25  # elv: elevation (m)
    z = 30  # z: reference height (m)
    PFT = 'DBF'
    vegetation_type = 0
    leaf_type = 0
    t_step = 60  # (min)
    year_start, year_end = 2000, 2000

if site == 'KWG':
    lat, lon = 35.8725, 139.4869
    t_offset_GMT = 9  # Site local time offset from UTC (GMT) (hours)
    elv = 26  # elv: elevation (m)
    z = 30  # z: reference height (m)
    PFT = 'DBF'
    vegetation_type = 0
    leaf_type = 0
    t_step = 30  # (min)
    year_start, year_end = 1995, 2004

t_start = int(t_step / 2)
lon_LST = 15 * t_offset_GMT

## vegetatiom type ######################
# 0: forest
# 1: grass
# vegetation_type = 0
#########################################

## leaf type ############################
# 0: deciduous
# 1: conifer-evergreen
# leaf_type = 1
#########################################

# input parameters ####################################
emissivity_c = 0.97  # emissivity of canopy
emissivity_s = 0.97  # emissivity of soil
extinction_coef = 0.5  # extinction coeffcient
u_attenuation_coef = 2.5 # in-canopy extinction coefficient ~ (Shuttleworth & Wallace, 1985)
d_leaf = np.array([0.05, 0.05, 0.03, 0.03])  # characteristic dimension of leaf (m) ~ (Jones, 1992)
albedo_c = 0.15
reflectance = np.array([0.09, 0.3])  # PAR, SR
transmissivity = np.array([0.06, 0.2])  # PAR, SR
albedo_soil = np.array([0.06, 0.26])  # PAR, SR
scattering_coefficient = np.array([0.20, 0.85])  # scattering coefficient of the leaf (PAR, SR)
reflectance_b = 0.09  # effective canopy-soil reflectance for direct beam radiaiton
reflectance_d = 0.09  # effective canopy-soil reflectance for diffusive radiaiton
k_b = 0.5  # extinction_coeffcient of a canopy for direct beam radiation
k_d = 0.4  # extinction_coeffcient of a canopy for direct diffuse radiation
k_b_black = 0.5  # extinction_coeffcient of a canopy of black leaves for direct beam radiation
k_d_black = 0.4 # extinction_coeffcient of a canopy of black leaves for direct diffuse radiation
k_n = 0.5  # extinction_coeffcient of a canopy for leaf nitrogen
k_u = 0.5  # extinction_coeffcient of a canopy for wind speed
absorptance = 1 - reflectance[0] - transmissivity[0]  # PAR
alpha = np.array([5, 2.5, 2.5, 2.5])
beta = np.array([2, 4.5, 4.5, 4.5])
C_a_out = 380  # atmospheric CO2 (ppm = μmol/mol)
O_a_out = 209490  # oxygenc partial pressure in chloroplast (ppm) (205–209 μbar)
g_s_max = np.array([[7.5, 8.3, 8.3, 8.3], [12, 10]], dtype=object)  # maximum stomatal conductance (mm/s) ~ (IGBP)
R_g_parameter = np.array([0.25, 0.25])  # growth respiration parameter ~ (Knorr, 2000)
allocation_parameter = np.array([[95, 52, 30, 35], [0.8, 3]], dtype=object)
allocation_parameter_index = np.array([[1.4, 1.3, 1.4, 1.6], [1.2, 1.2]], dtype=object)
r_0 = 0.3  # fractional carbon allocation to root for non-limiting conditions
s_0 = 0.3  # fractional carbon allocation to stem for non-limiting conditions
r_00 = np.array([[0.60, 0.60, 0.55, 0.6], [0.8, 0.8]], dtype=object) # fractional carbon allocation to root for non-limiting conditions
s_00 = np.array([0.1, 0.05, 0.2, 0.05]) # fractional carbon allocation to stem for non-limiting conditions
sensitivity_allocation = np.array([[0.8, 0.8, 0.8, 0.5], [1, 1]], dtype=object) # sensitivity of allocation to changes in W & L
specified_R_rate = 0.218 # specified respiration rate based on N content (kgC/kgN/dy) ~ (Keyser et al., 2000)
CN_ratio_stem = 50  # C:N ratio of stem ~ (Arora, 2003)
CN_ratio_root = 50  # C:N ratio of root ~ (Arora, 2003)
V_c_max_org = np.array([[80, 96, 50, 94], [80, 80]], dtype=object)
Ha_V = np.array([[62000, 71000, 71000, 71000], [65000, 70000]], dtype=object) # activation_energy (J/mol)
J_m_25 = np.array([[80, 84, 60, 82], [100, 60]], dtype=object)  # μmol/m^2/s
Ha_J = np.array([[48000, 59000, 52000, 52000], [48000, 52000]], dtype=object) # activation_energy (J/mol)
Hd_J = np.array([[219400, 219400, 219400, 219400], [219400, 219400]], dtype=object)  # deactivation_energy (J/mol)
Delta_S = np.array([635.2, 650, 650, 650])  # entropy factor (J/mol/K)
Gamma_25 = np.array([[37, 44.7, 40, 44.7], [45, 42]], dtype=object) # CO2 compensation point without dark respiration (μbar)
quantum_efficiency = 0.385 # quantum efficiency of RuBP regeneration (mol e mol-1 quanta)
a_1 = 10
g_0_sunlit = 0.1
g_0_shaded = 0.1
VPD_0 = 15  # Leuning, 1995
R_m_base_stem = np.array([0.00005, 0.0008, 0.00002, 0.0008])
R_m_base_root = np.array([[0.0012, 0.0008, 0.0012, 0.001], [0.0002, 0.0015]], dtype=object)
death_rate_leaf = 0.0045  # death rate of leaf (/dy)
push_down_rate = 0.05  # rate of leaf pushed down (/dy)
# push_down_rate = 0.15  # rate of leaf pushed down (/dy)
loss_rate_W_max = np.array([[0.005, 0.005, 0.005, 0.015], [0.025, 0.015]], dtype=object)  # maximum drought leaf loss rate (/dy)
loss_rate_T_max = 3 # maximum cold leaf loss rate (/dy)
b_W = np.array([[1.1, 3, 3, 2], [2, 3]], dtype=object)  # shape parameter for leaf loss (drought)
b_T = 3  # shape parameter for leaf loss (cold)
T_cold = np.array([[7, 0, 5, -5], [5, 8]], dtype=object)  # Temperature threshold for leaf loss because of cold stress (°C)
loss_rate_stem = np.array([[0.000047, 0.000061, 0.000061, 0.000055], [0.00005, 0.00004]], dtype=object) # stem turnover rate (/dy)
loss_rate_root = np.array([[0.00023, 0.00034, 0.00034, 0.00031], [0.00025, 0.00025]], dtype=object) # root turnover rate (/dy)
loss_rate_standby_leaf = np.array([[0, 0.00052, 0.00052, 0.00055], [0, 0.0001]], dtype=object)
loss_rate_standby_stem = np.array([0.00012, 0.00012, 0.00022, 0.000055])
loss_rate_standby_root = np.array([[0.00012, 0.00012, 0.00022, 0.000055], [0.01, 0.0006]], dtype=object)
virtual_LAI_succeed_day = 7 # successive days necessary for virtual leaf to survive at onset

dormancy_terminate_day = 90
porosity = 0.5  # porosity (i.e. saturated moisture content)

########################################################

## constants ###########################################
pi = math.pi
R = 8.31447  # $R: gas constant (J/K/mol)
rho = 1.20  # $rho: air density (kg/m^3)
M_d = 0.02897  # $M_d: molecular mass of dry air (kg/mol)
lapse_rate = 0.0065  # $lapse_rate: lapse-rate of air temperature (K/m)
refraction = 0.01  # $refraction: horizontal refraction (rad)
# $SteBol_const: Stefan-Bolzman constatnt (W/m^2/K^4)
SteBol_const = 5.67*10**(-8)
C_pd = 1005  # $C_pd: specific heat of dry air at constant pressure (J/kg/K)
gamma = 0.68  # $gamma: psychrometric constant (hPa=mbar)
# empirical parameters ~ (Tetens's equation, 1930)
a0, b0, c0 = 17.67, 243.5, 6.112
Karman_const = 0.41  # $Karman_const: the von-Karman constant
f = 0.15  # to correct for spectral quality of the light (Evans, 1987)
# empirical curvature factor (0.7 is a good average value, Evans, 1989)
a = 0.7
f_s = 0.7  # Direct solar fraction of solar radiation incident on the canopy
########################################################

if vegetation_type == 0:
    N_layer = 20
elif vegetation_type == 1:
    N_layer = 10

N_angle = 10
N2_angle = 40
diff_level = 0.01

latitude = lat / 180 * pi

if not os.path.isdir(f'../out/{site}'):
    os.mkdir(f'../out/{site}')
if not os.path.isdir(f'../out/{site}/yearly'):
    os.mkdir(f'../out/{site}/yearly')

with open(f'./../out/{site}/yearly/yearly.csv', 'a') as YEARLY:
    '''
    with open(f'./../data/{site}/soil/fieldcap.out', 'r') as FCAP:
        FCAP_line = FCAP.readline().strip()
        FCAP_data = FCAP_line.split('\t')
        W_capacity = FCAP_data[2]

    with open(f'./../data/{site}/soil/wiltpoint.out', 'r') as WILT:
        WILT_line = WILT.readline().strip()
        WILT_data = WILT_line.split('\t')
        W_wilting = WILT_data[2]
        W_critical = W_capacity * 0.75  # W_critical: critical point (mm)
    '''
    W_retention = 0
    W_capacity = 400
    W_wilting = 200
    W_critical = W_capacity * 0.75  # W_critical: critical point (mm)

    leaf_onset = 0
    leaf_fall = 0
    leaf_dormant = 0
    leaf_normal = 0
    C_increase_dy = 0
    C_decrease_dy = 0

    if leaf_type == 0:
        phenophase = 0
        C_leaf = 0.0
        Cg_leaf = 0.0
        Cd_leaf = 0.0
        C_stem = 8.0
        C_root = 1.8
        C_all = C_stem+C_root
    elif leaf_type == 1:
        phenophase = 2
        C_leaf = 0.0
        Cg_leaf = 0.0
        Cd_leaf = 0.0
        C_stem = 0.5
        C_root = 0.2
        C_all = C_stem+C_root
    elif leaf_type == 2:
        phenophase = 2
        C_leaf = 0.15
        Cg_leaf = 0.13
        Cd_leaf = 0.02
        C_stem = 6.5
        C_root = 1.8
        C_all = C_leaf+C_stem+C_root
    elif leaf_type == 3:
        phenophase = 2
        C_leaf = 0.18
        Cg_leaf = 0.16
        Cd_leaf = 0.02
        C_stem = (0.6-0.03*(lat-40))*10
        C_all = C_leaf+C_stem #+C_root

    if vegetation_type == 1:
        C_leaf = 0.0
        Cg_leaf = 0.0
        Cd_leaf = 0.0
        C_stem = 0
        if leaf_type == 0:
            C_root = 0
        elif leaf_type == 1:
            C_root = 0.55
        C_all = C_leaf+C_root

    phase2to3_dy = 0
    phase3to2_dy = 0
    phase3to2_dy2 = 0
    phaseto4_dy = 0

    W_bare = 0
    A_sum_yearly_bare = 0
    R_m_leaf_sum_yearly_bare = 0
    ET_yearly_bare = 0
    ET_c_yearly_bare = 0
    ET_eq_yearly_bare = 0
    C_leaf_bare = 0
    Cg_leaf_bare = 0
    Cd_leaf_bare = 0
    C_stem_bare = 0
    C_root_bare = 0
    C_all_bare = 0

    sun_duration_sum = 0
    sun_duration_count = -9
    solar_elevation = 0
    solar_elevation_pre = -30
    T_a_C_pre = 0
    rainfall_pre = 0
    rh_pre = 0
    u_z_pre = 0
    R_s_total_pre = 0

    if vegetation_type == 0:
        W = W_capacity * 0.9
    elif vegetation_type == 1:
        W = W_capacity * 0.8
    W_pre = W

    LAI_g = 0
    virtual_LAI = 0
    virtual_LAI_day = 0
    dormancy_dys = 0
    below_minus5 = 0
    above_minus5 = 0
    below_0 = 0
    above_0 = 0

    Q_n_sunlit_pre = 0
    Q_n_shaded_pre = 0
    Gamma_respiration_sunlit = 0
    Gamma_respiration_shaded = 0
    Gamma_respiration_sunlit_pre = 0
    Gamma_respiration_shaded_pre = 0
    C_i_sunlit = 0
    C_i_shaded = 0
    A_c_sunlit = 0
    A_c_shaded = 0
    J_sunlit_pre = 0
    J_shaded_pre = 0
    V_j_sunlit_pre = 0
    V_j_shaded_pre = 0

    g_b_free_sunlit_pre = 0
    g_b_free_shaded_pre = 0
    g_b_free_soil_pre = 0
    G_0_sunlit_pre = 0
    G_0_shaded_pre = 0
    g_b_forced_sunlit_pre = 0
    g_b_forced_shaded_pre = 0
    g_b_sunlit_pre = 0
    g_b_shaded_pre = 0
    g_h_sunlit_pre = 0
    g_h_shaded_pre = 0
    g_r_sunlit_pre = 0
    g_r_shaded_pre = 0
    r_a_c_pre = 0
    r_s_c_pre = 0

    itrn = 0
    pass200 = 0

    diff_count_yr = 0
    sum_C_stem = 0
    sum_C_root = 0
    ave_C_stem_pre = 0
    ave_C_root_pre = 0
    A_sum_yearly_pre = 0
    R_m_leaf_sum_yearly_pre = 0
    R_a_c_yearly_pre = 0
    NPP_yearly_pre = 0
    ET_yearly_pre = 0
    ET_c_yearly_pre = 0
    T_a_C_mean = 0
    pressure_mean = 0
    rh_mean = 0
    u_z_mean = 0
    W_mean = 0
    T_a_C_mean_7 = [0] * 7
    virtual_DOY = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_T_a_C_mean = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_R_s_sum = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_pressure_mean = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_rainfall_sum = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_rh_mean = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_u_z_mean = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_W_mean = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_A_sum_daily = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_R_a_c = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_R_m_leaf_sum_daily = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_ET_daily = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_ET_c_daily = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_ET_eq_daily = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_LAI_list = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_LAI_g = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_C_leaf = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_Cg_leaf = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_Cd_leaf = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_C_stem = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_C_root = np.zeros(virtual_LAI_succeed_day + 1)
    virtual_C_all = np.zeros(virtual_LAI_succeed_day + 1)

    # iterate calculations until steady conditions
    BreakLOOP = False
    while True:
        if BreakLOOP:
            break
        for year in range(year_start, year_end + 1):
            if BreakLOOP:
                break
            month_days = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
            dy_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            DOY_max = 365
            if year % 4 == 0:
                month_days = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
                dy_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
                DOY_max = 366

            if leaf_type == 0:
                growing_dys = 365  # $growing_dys: leaf life span (dy)
                loss_rate = 1 / growing_dys  # loss_rate: normal turnover (/dy)
            elif leaf_type in [1, 2, 3]:
                if vegetation_type == 0:
                    loss_rate = 1 / DOY_max / 1.75
                elif vegetation_type == 1:
                    loss_rate = 1 / DOY_max

            if virtual_LAI != 1:
                diff_count_yr = 0
                A_sum_yearly = 0
                R_a_c_yearly = 0
                R_m_leaf_sum_yearly = 0
                NPP_yearly = 0
                ET_yearly = 0
                ET_c_yearly = 0
                ET_eq_yearly = 0

                if vegetation_type == 0:
                    if leaf_type == 0:
                        SLA = 27  # $SLA: specidic leaf area (m^2/kgC)
                    elif leaf_type == 1:
                        SLA = 20
                    elif leaf_type == 2 or leaf_type == 3:
                        SLA = 12
                elif vegetation_type == 1:
                    if leaf_type == 0:
                        SLA = 30
                    elif leaf_type == 1:
                        SLA = 25

            func.MakeDirectory(site, year)
            with open(f'../out/{site}/{year}/daily/allocation_{year}.csv', 'w') as ALLOCATION, \
                    open(f'../out/{site}/{year}/daily/daily_{year}.csv', 'w') as DAILY, \
                    open(f'../out/{site}/{year}/meanclim/meanclim_{year}.csv', 'w') as MEANCLIM:

                y2 = -1
                for month in range(1, 13):
                    if BreakLOOP:
                        break
                    day = 0

                    while day <= dy_month[month - 1]:
                        day += 1

                        if day > dy_month[month - 1]:
                            break

                        DOY = month_days[month - 1] + day
                        if DOY == 200:
                            pass200 = 1
                        
                        # initial conditions
                        if virtual_LAI == -9:
                            LAI = 0
                            virtual_LAI = 0
                            W = W_bare
                            A_sum_yearly = A_sum_yearly_bare
                            R_m_leaf_sum_yearly = R_m_leaf_sum_yearly_bare
                            ET_yearly = ET_yearly_bare
                            ET_c_yearly = ET_c_yearly_bare
                            ET_eq_yearly = ET_eq_yearly_bare
                            C_leaf = C_leaf_bare
                            Cg_leaf = Cg_leaf_bare
                            Cd_leaf = Cd_leaf_bare
                            C_stem = C_stem_bare
                            C_root = C_root_bare
                            C_all = C_all_bare
                        elif leaf_dormant == 0 and virtual_LAI != 1 and (
                            (vegetation_type == 0 and LAI_g < 0.3) or 
                            (vegetation_type == 1 and LAI_g < 0.2) or 
                            Cg_leaf < ((C_stem + C_root) / allocation_parameter[vegetation_type][leaf_type]) ** (1 / allocation_parameter_index[vegetation_type][leaf_type]) * 0.075
                            ):
                            if leaf_type == 0:
                                leaf_fall = 0
                                W_bare = W
                                A_sum_yearly_bare = A_sum_yearly
                                R_m_leaf_sum_yearly_bare = R_m_leaf_sum_yearly
                                ET_yearly_bare = ET_yearly
                                ET_c_yearly_bare = ET_c_yearly
                                ET_eq_yearly_bare = ET_eq_yearly
                                C_leaf_bare = C_leaf
                                Cg_leaf_bare = Cg_leaf
                                Cd_leaf_bare = Cd_leaf
                                C_stem_bare = C_stem
                                C_root_bare = C_root
                                C_all_bare = C_all
                                C_leaf_0 = ((C_stem + C_root) / allocation_parameter[vegetation_type][leaf_type]) ** (1 / allocation_parameter_index[vegetation_type][leaf_type]) * 0.075
                                if vegetation_type == 0:
                                    if C_leaf_0 < 0.3 / SLA:
                                        C_leaf_0 = 0.3 / SLA
                                        Cg_leaf = C_leaf_0
                                        Cd_leaf = 0
                                elif vegetation_type == 1:
                                    if C_leaf_0 < 0.2 / SLA:
                                        C_leaf_0 = 0.2 / SLA
                                        Cg_leaf = C_leaf_0
                                        Cd_leaf = 0
                                        C_root = C_leaf_0
                                C_leaf = C_leaf_0
                                Cg_leaf = C_leaf_0
                                Cd_leaf = 0
                                C_all += C_leaf + C_root
                                virtual_LAI = 1
                                virtual_LAI_day = 0
                                virtual_A_n = 0
                                if leaf_type == 0 and growing_dys > virtual_LAI_succeed_day:
                                    growing_dys_pre = growing_dys
                                    growing_dys = 0
                                print('>', growing_dys, '\t', Cg_leaf, '\t', Cd_leaf, '\t', C_leaf)
                            elif leaf_type in [1, 2, 3]:
                                if W >= W_wilting:
                                    leaf_normal = 0
                                    phenophase = 1
                            
                        LAI = SLA * C_leaf
                        LAI_g = SLA * Cg_leaf
                        LAI_d = SLA * Cd_leaf
                        N_layer = int(LAI*4)
                        if N_layer <= 4:
                            N_layer = 4
                        
                        if vegetation_type == 0:
                            a_ch = 6.5 # for forest ~ (Lawrence & Slingo, 2004)
                        elif vegetation_type == 1:
                            a_ch = 0.5 # for grass ~ (Lawrence & Slingo, 2004)
                        
                        # height_c: canopy height (m)
                        if LAI == 0:
                            if vegetation_type == 0:
                                height_c = 10.0 * C_stem ** 0.385
                            elif vegetation_type == 1:
                                height_c = 0.1
                        else:
                            if vegetation_type == 0:
                                height_c = 10.0 * C_stem ** 0.385
                            elif vegetation_type == 1:
                                height_c = (a_ch * LAI) ** (2/3)
                            if height_c > z:
                                z = height_c + 1

                        d_z_h = 1 / N_layer
                        d_z = d_z_h * height_c

                        A_sum_daily = 0
                        A_n_sum_daily = 0
                        R_m_leaf_sum_daily = 0
                        ET_daily = 0
                        ET_c_daily = 0
                        ET_eq_daily = 0
                        T_a_C_sum = 0
                        R_s_sum = 0
                        pressure_sum = 0
                        rainfall_sum = 0
                        rh_sum = 0
                        u_z_sum = 0
                        W_sum = 0
                        A_n_obs_sum = 0
                        S_b_down = np.zeros(N_layer+1)
                        S_d_down = np.zeros(N_layer+1)
                        S_d_up = np.zeros(N_layer+1)

                        t = -10 # time of the day (min)
                        with open(f'./../data/{site}/forcing/{year}/{site}_{year}_{str(month).zfill(2)}_{str(day).zfill(2)}.csv', 'r') as IN, \
                            open(f'./../out/{site}/{year}/canopy/canopy_{year}_{str(month).zfill(2)}_{str(day).zfill(2)}.csv', 'w') as CANOPY:

                            reader = csv.reader(IN)
                            for t, flux_data in zip(range(t_start, 1441, t_step), reader):

                                hr = int(t / 60) # (hour)
                                min = int(t % 60) # (min)
                                flux_data[1:] = [float(data) for data in flux_data[1:]]

                                T_a_C, T_a_K, rainfall, u_z, pressure, R_s_total, A_n_obs, rh, VPD_a, rho_a, C_p, c_p, gamma, Delta, sun_duration_pos, COSTHETA, R_s_d_total, sun_duration, S_max, R_s_b_total, solar_elevation, k_b, k_d, k_b_black, k_d_black, L_down, R_s_b, R_s_d, T_a_C_pre, rainfall_pre, rh_pre, u_z_pre, R_s_total_pre = func.inputClimateData(t, t_start, t_step, elv, flux_data, solar_elevation, a0, b0, c0, R, M_d, C_pd, k_b_black, DOY_max, DOY, latitude, refraction, lon, lon_LST, SteBol_const, sun_duration_sum, sun_duration_count, T_a_C_pre, rainfall_pre, rh_pre, u_z_pre, R_s_total_pre)
                                ### Main ##########

                                T_a_C_sum, R_s_sum, pressure_sum, rainfall_sum, rh_sum, u_z_sum, W_sum, A_n_obs_sum, T_a_C_mean, pressure_mean, rh_mean, u_z_mean, W_mean = func.DailyMeanClimate(t, t_start, t_step, T_a_C, rainfall, pressure, R_s_total, A_n_obs, rh, u_z, W, T_a_C_sum, R_s_sum, pressure_sum, rainfall_sum, rh_sum, u_z_sum, W_sum, A_n_obs_sum, T_a_C_mean, pressure_mean, rh_mean, u_z_mean, W_mean)
                                
                                if LAI > 0:
                                    f_list, z_c, LA_z, d_LA, d_LA_g, LA_z_max, z_c = func.LAdistribution(N_layer, d_z_h, alpha, beta, leaf_type, height_c, d_z, Cg_leaf, C_leaf, LAI)
                                    if S_max <= 0:
                                        for i_l in range(N_layer, -1, -1):
                                            S_b_down[i_l] = 0
                                            S_d_down[i_l] = 0
                                            S_d_up[i_l] = 0
                                    I_b, I_d, G_layer = func.CalI(N_layer, solar_elevation, N_angle, N2_angle, d_LA)
                                
                                if S_max <= 0:
                                    for i_l in range(N_layer, -1, -1):
                                        S_b_down[i_l] = 0
                                        S_d_down[i_l] = 0
                                        S_d_up[i_l] = 0
                                    LAI_sunlit = 0
                                    LAI_shaded = LAI
                                    Q_sunlit = [0]*2
                                    Q_shaded = [0]*2
                                else:
                                    if LAI > 0:
                                        d_LA_su = np.zeros(N_layer+1)
                                        d_LA_sh = np.zeros(N_layer+1)
                                        S_d_down_pre = np.zeros(N_layer+1)
                                        S_d_up_pre = np.zeros(N_layer+1)
                                        diff_S_d_down = np.zeros(N_layer+1)
                                        diff_S_d_up = np.zeros(N_layer+1)
                                        Q_sunlit = [0, 0]
                                        Q_shaded = [0, 0]
                                        L = 1 # near infra red radiation
                                        d_LA_su , d_LA_sh , S_d_down_pre , S_d_up_pre , diff_S_d_down , diff_S_d_up , Q_shaded , Q_sunlit = func.S_balance(reflectance, transmissivity, albedo_soil, N_layer, diff_level, S_max, R_s_b, R_s_d, solar_elevation, d_LA, S_b_down, S_d_down, S_d_up, L, I_b, I_d, G_layer, d_LA_su , d_LA_sh , S_d_down_pre , S_d_up_pre , diff_S_d_down , diff_S_d_up , Q_shaded , Q_sunlit)
                                        LAI_sunlit = 0
                                        if R_s_b_total > 0:
                                            for i_l in range(N_layer, 0, -1):
                                                LAI_sunlit += d_LA_su[i_l]
                                        LAI_shaded = LAI - LAI_sunlit
                                        L = 0
                                        d_LA_su , d_LA_sh , S_d_down_pre , S_d_up_pre , diff_S_d_down , diff_S_d_up , Q_shaded , Q_sunlit = func.S_balance(reflectance, transmissivity, albedo_soil, N_layer, diff_level, S_max, R_s_b, R_s_d, solar_elevation, d_LA, S_b_down, S_d_down, S_d_up, L, I_b, I_d, G_layer, d_LA_su , d_LA_sh , S_d_down_pre , S_d_up_pre , diff_S_d_down , diff_S_d_up , Q_shaded , Q_sunlit)
                                if LAI > 0:
                                    if R_s_total > 0 and S_max > 0:
                                        g_s_sunlit = 1 * W_retention
                                        g_s_shaded = 0.5 * W_retention
                                        C_s_sunlit = C_a_out * 0.8
                                        C_s_shaded = C_a_out * 0.8
                                        C_i_sunlit = C_a_out * 0.7
                                        C_i_shaded = C_a_out * 0.7
                                        diff_T_sunlit = 0
                                        diff_T_shaded = 0
                                        T_c_C_soil = T_a_C
                                        T_c_K_soil = T_c_C_soil + 273.15
                                        no_converged = 0
                                        count_i = 0

                                        while True:
                                            diff_T_sunlit_pre = diff_T_sunlit
                                            diff_T_shaded_pre = diff_T_shaded
                                            C_s_sunlit_pre = C_s_sunlit
                                            C_s_shaded_pre = C_s_shaded
                                            C_i_sunlit_pre = C_i_sunlit
                                            C_i_shaded_pre = C_i_shaded

                                            ### SimultaneousEquations ##########
                                            T_c_C_sunlit = T_a_C + diff_T_sunlit
                                            T_c_K_sunlit = T_c_C_sunlit + 273.15
                                            T_c_C_shaded = T_a_C + diff_T_shaded
                                            T_c_K_shaded = T_c_C_shaded + 273.15
                                            Q_long_isothermal_sunlit, Q_long_isothermal_shaded = func.L_balance(emissivity_c , k_b_black , k_d_black , SteBol_const , LAI , R_s_total , T_a_K , L_down , LAI_sunlit , LAI_shaded)
                                            Q_n_isothermal_sunlit, Q_n_isothermal_shaded = func.CalRn(Q_sunlit , Q_shaded , Q_long_isothermal_sunlit , Q_long_isothermal_shaded)
                                            r_a, r_a_soil, G_0_sunlit, G_0_shaded, g_a, g_b_free_sunlit, g_b_forced_sunlit, g_b_sunlit, g_b_forced_shaded, g_b_free_shaded, g_b_shaded, g_h_sunlit, g_h_shaded, g_b_free_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, h_beta, g_b_free_sunlit_pre, g_b_free_shaded_pre, g_b_free_soil_pre, G_0_sunlit_pre, G_0_shaded_pre, g_b_forced_sunlit_pre, g_b_forced_shaded_pre, g_b_sunlit_pre, g_b_shaded_pre, g_h_sunlit_pre, g_h_shaded_pre, g_r_sunlit_pre, g_r_shaded_pre, r_a_c_pre, r_s_c_pre = func.CalConductance(vegetation_type, leaf_type, d_leaf, R, Karman_const, SteBol_const, emissivity_c, k_b_black, k_d_black, k_n, k_u, g_0_sunlit, g_0_shaded, porosity, W, W_capacity, z, c_p, u_z, u_attenuation_coef, LAI, height_c, pressure, T_a_K, LAI_sunlit, LAI_shaded, T_c_K_soil, T_c_K_sunlit, T_c_K_shaded, g_s_sunlit, g_s_shaded, g_b_free_sunlit_pre , g_b_free_shaded_pre , g_b_free_soil_pre , G_0_sunlit_pre , G_0_shaded_pre , g_b_forced_sunlit_pre , g_b_forced_shaded_pre , g_b_sunlit_pre , g_b_shaded_pre , g_h_sunlit_pre , g_h_shaded_pre , g_r_sunlit_pre , g_r_shaded_pre , r_a_c_pre , r_s_c_pre)
                                            gamma_modified_sunlit, gamma_modified_shaded, C_s_sunlit, C_s_shaded, VPD_s_sunlit, VPD_s_shaded = func.ResConst(C_a_out , VPD_a , gamma , Delta , g_s_sunlit , g_s_shaded , C_i_sunlit , C_i_shaded , diff_T_sunlit , diff_T_shaded , g_a , g_b_forced_sunlit , g_b_free_sunlit , g_b_sunlit , g_b_forced_shaded , g_b_free_shaded , g_b_shaded , g_h_sunlit , g_h_shaded , g_r_sunlit , g_r_shaded)
                                            k_n, Gamma_respiration_sunlit, Gamma_respiration_shaded, C_i_sunlit, C_i_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, f_w, Gamma_respiration_sunlit_pre, Gamma_respiration_shaded_pre, W_retention = func.SimultaneousPhotosynthesis(vegetation_type, leaf_type, reflectance, transmissivity, scattering_coefficient, reflectance_b, reflectance_d, k_b, k_d, k_b_black, absorptance, C_a_out, V_c_max_org, Ha_V, J_m_25, Ha_J, Hd_J, Delta_S, Gamma_25, a_1, VPD_0, T_cold, R, f, a, W_capacity, W_wilting, W, R_s_b, R_s_d, R_s_b_total, R_s_d_total, R_s_total, pressure, S_max, T_a_C_mean, LAI, LAI_g, LAI_sunlit, LAI_shaded, Q_sunlit, C_s_shaded, C_i_sunlit, C_i_shaded, T_c_C_sunlit, T_c_K_sunlit, T_c_C_shaded, T_c_K_shaded, G_0_sunlit, G_0_shaded, h_beta, C_s_sunlit, VPD_s_sunlit, VPD_s_shaded, J_sunlit_pre, J_shaded_pre, V_j_sunlit_pre, V_j_shaded_pre, Gamma_respiration_sunlit_pre, Gamma_respiration_shaded_pre)
                                            r_a, r_a_soil, G_0_sunlit, G_0_shaded, g_a, g_b_free_sunlit, g_b_forced_sunlit, g_b_sunlit, g_b_forced_shaded, g_b_free_shaded, g_b_shaded, g_h_sunlit, g_h_shaded, g_b_free_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, h_beta, g_b_free_sunlit_pre, g_b_free_shaded_pre, g_b_free_soil_pre, G_0_sunlit_pre, G_0_shaded_pre, g_b_forced_sunlit_pre, g_b_forced_shaded_pre, g_b_sunlit_pre, g_b_shaded_pre, g_h_sunlit_pre, g_h_shaded_pre, g_r_sunlit_pre, g_r_shaded_pre, r_a_c_pre, r_s_c_pre = func.CalConductance(vegetation_type, leaf_type, d_leaf, R, Karman_const, SteBol_const, emissivity_c, k_b_black, k_d_black, k_n, k_u, g_0_sunlit, g_0_shaded, porosity, W, W_capacity, z, c_p, u_z, u_attenuation_coef, LAI, height_c, pressure, T_a_K, LAI_sunlit, LAI_shaded, T_c_K_soil, T_c_K_sunlit, T_c_K_shaded, g_s_sunlit, g_s_shaded, g_b_free_sunlit_pre , g_b_free_shaded_pre , g_b_free_soil_pre , G_0_sunlit_pre , G_0_shaded_pre , g_b_forced_sunlit_pre , g_b_forced_shaded_pre , g_b_sunlit_pre , g_b_shaded_pre , g_h_sunlit_pre , g_h_shaded_pre , g_r_sunlit_pre , g_r_shaded_pre , r_a_c_pre , r_s_c_pre)
                                            Q_n_sunlit, Q_n_shaded, Q_n_sunlit_pre, Q_n_shaded_pre, R_n_soil, G, R_n_sum, lE_soil, lE_sum, ET_c_sum, ET_soil, ET_sum, H_sum = func.CalET_SW(emissivity_c, emissivity_s, extinction_coef, albedo_soil, SteBol_const, T_a_K, R_s_total, VPD_a, rho_a, C_p, c_p, gamma, Delta, S_max, L_down, COSTHETA, LAI, LAI_g, T_c_K_sunlit, T_c_K_shaded, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, r_a, r_a_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, Q_n_sunlit_pre, Q_n_shaded_pre)
                                            diff_T_sunlit, diff_T_shaded = func.CalDiffT(VPD_a, c_p, Delta, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, g_h_sunlit, g_h_shaded, g_r_sunlit, g_r_shaded, gamma_modified_sunlit, gamma_modified_shaded)
                                            if abs(diff_T_sunlit) > 10:
                                                diff_T_sunlit = 0
                                                no_converged = 1
                                                break
                                            if abs(diff_T_shaded) > 10:
                                                diff_T_shaded = 0
                                                no_converged = 1
                                                break
                                            #####################################

                                            count_i += 1
                                            if S_max > 0 and R_s_total > 0:
                                                if (abs(diff_T_sunlit - diff_T_sunlit_pre) < 0.01 and 
                                                    abs(diff_T_shaded - diff_T_shaded_pre) < 0.01 and 
                                                    abs(C_s_sunlit - C_s_sunlit_pre) < 0.01 and 
                                                    abs(C_s_shaded - C_s_shaded_pre) < 0.01 and 
                                                    abs(C_i_sunlit - C_i_sunlit_pre) < 0.01 and 
                                                    abs(C_i_shaded - C_i_shaded_pre) < 0.01):
                                                    break
                                            else:
                                                if (abs(diff_T_sunlit - diff_T_sunlit_pre) < 0.01 and
                                                    abs(diff_T_shaded - diff_T_shaded_pre) < 0.01):
                                                    break

                                            if count_i > 100:
                                                break

                                        if no_converged == 1:
                                            ### NoConverged ###################
                                            g_s_sunlit = 1
                                            g_s_shaded = 0.5
                                            C_s_sunlit = C_a_out * 0.8
                                            C_s_shaded = C_a_out * 0.8
                                            C_i_sunlit = C_a_out * 0.7
                                            C_i_shaded = C_a_out * 0.7
                                            diff_T_sunlit = 0
                                            diff_T_shaded = 0
                                            T_c_C_soil = T_a_C
                                            T_c_K_soil = T_c_C_soil + 273.15

                                            Q_long_isothermal_sunlit, Q_long_isothermal_shaded = func.L_balance(emissivity_c , k_b_black , k_d_black , SteBol_const , LAI , R_s_total , T_a_K , L_down , LAI_sunlit , LAI_shaded)
                                            Q_n_isothermal_sunlit, Q_n_isothermal_shaded = func.CalRn(Q_sunlit , Q_shaded , Q_long_isothermal_sunlit , Q_long_isothermal_shaded)
                                            Q_n_sunlit, Q_n_shaded, R_n_soil, G, R_n_sum, Q_n_sunlit_pre, Q_n_shaded_pre = func.CalRnSoil(emissivity_c, emissivity_s, extinction_coef, albedo_soil, SteBol_const, T_a_K, R_s_total, c_p, S_max, L_down, COSTHETA, LAI, T_c_K_sunlit, T_c_K_shaded, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, g_r_sunlit, g_r_shaded, Q_n_sunlit_pre, Q_n_shaded_pre)

                                            W_retention = (W - W_wilting) / (W_capacity - W_wilting)
                                            beta_water = max(0, W_retention)
                                            G_water = max(0.001, 1 - (1-beta_water) ** 2)
                                            k_n, K_c_sunlit, K_c_shaded, K_o_sunlit, K_o_shaded, Gamma_sunlit, Gamma_shaded, O_i, V_c_max_sunlit, V_c_max_shaded, V_c_sunlit, V_c_shaded, J_sunlit, J_shaded, V_j_sunlit, V_j_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, J_sunlit_pre, J_shaded_pre, V_j_sunlit_pre, V_j_shaded_pre = func.Photosynthesis(vegetation_type, leaf_type, transmissivity, reflectance, scattering_coefficient, reflectance_b, reflectance_d, k_b, k_d, k_b_black, absorptance, V_c_max_org, Ha_V, J_m_25, Ha_J, Hd_J, Delta_S, Gamma_25, T_cold, R, f, a, R_s_b, R_s_d, R_s_total, pressure, S_max, T_a_C_mean, LAI, LAI_g, LAI_sunlit, LAI_shaded, Q_sunlit, C_i_sunlit, C_i_shaded, T_c_C_sunlit, T_c_K_sunlit, T_c_C_shaded, T_c_K_shaded, J_sunlit_pre, J_shaded_pre, V_j_sunlit_pre, V_j_shaded_pre, called_from='NoConverged')
                                            extinction_coefficient_n = func.ExtinctionCoeffcientFunction(LAI,k_n)
                                            extinction_coeffcient_b_black_plus_n = func.ExtinctionCoeffcientFunction(LAI,k_b_black+k_n)
                                            G_0_sunlit = g_0_sunlit * extinction_coeffcient_b_black_plus_n
                                            G_0_shaded = g_0_shaded * (extinction_coefficient_n - extinction_coeffcient_b_black_plus_n)
                                            g_s_sunlit = G_0_sunlit + a_1 * f_w * A_c_sunlit / (C_s_sunlit - Gamma_respiration_sunlit) / (1 + VPD_s_sunlit / VPD_0)  # mol m^2/s
                                            g_s_shaded = G_0_shaded + a_1 * f_w * A_c_shaded / (C_s_shaded - Gamma_respiration_shaded) / (1 + VPD_s_shaded / VPD_0)  # mol m^2/s
                                            r_s_c = 2 / (g_s_sunlit + g_s_shaded) * pressure * 100 / R / T_a_K  # (s/m)
                                            lE_c = (Delta * (R_n_sum - R_n_soil) + (rho_a * C_p * VPD_a) / (r_a + r_a_c)) / (Delta + gamma * (1 + r_s_c / (r_a + r_a_c))) * LAI_g / LAI
                                            lE_soil = (Delta * (R_n_soil - G) + (rho_a * C_p * VPD_a) / (r_a + r_a_soil)) / (Delta + gamma * (1 + r_s / (r_a + r_a_soil))) * f_w
                                            if abs(lE_sum) > abs(R_n_sum * 1.5):
                                                lE_sum = R_n_sum
                                                lE_c = lE_sum
                                                lE_soil = 0
                                            lambda_kg = (56780.3 - 42.84 * T_a_K) / 0.018 # latent heat of vaporization (J/kg)
                                            ET_c_sum = lE_c / lambda_kg
                                            ET_soil = lE_soil / lambda_kg
                                            ET_sum = lE_sum / lambda_kg
                                            H_sum = R_n_sum - lE_sum - G
                                            #####################################

                                    else: # R_s_total <= 0 or S_max <= 0
                                        g_s_sunlit = 1
                                        g_s_shaded = 0.5
                                        C_s_sunlit = C_a_out * 0.8
                                        C_s_shaded = C_a_out * 0.8
                                        C_i_sunlit = C_a_out * 0.7
                                        C_i_shaded = C_a_out * 0.7
                                        diff_T_sunlit = 0
                                        diff_T_shaded = 0
                                        T_c_C_soil = T_a_C
                                        T_c_K_soil = T_c_C_soil + 273.15
                                        no_converged = 0
                                        count = 0
                                        while True:
                                            ### SimultaneousEquationsNight ##########
                                            T_c_C_sunlit = T_a_C + diff_T_sunlit
                                            T_c_K_sunlit = T_c_C_sunlit + 273.15
                                            T_c_C_shaded = T_a_C + diff_T_shaded
                                            T_c_K_shaded = T_c_C_shaded + 273.15
                                            Q_long_isothermal_sunlit, Q_long_isothermal_shaded = func.L_balance(emissivity_c , k_b_black , k_d_black , SteBol_const , LAI , R_s_total , T_a_K , L_down , LAI_sunlit , LAI_shaded)
                                            Q_n_isothermal_sunlit = Q_long_isothermal_sunlit
                                            Q_n_isothermal_shaded = Q_long_isothermal_shaded
                                            r_a, r_a_soil, G_0_sunlit, G_0_shaded, g_a, g_b_free_sunlit, g_b_forced_sunlit, g_b_sunlit, g_b_forced_shaded, g_b_free_shaded, g_b_shaded, g_h_sunlit, g_h_shaded, g_b_free_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, h_beta, g_b_free_sunlit_pre, g_b_free_shaded_pre, g_b_free_soil_pre, G_0_sunlit_pre, G_0_shaded_pre, g_b_forced_sunlit_pre, g_b_forced_shaded_pre, g_b_sunlit_pre, g_b_shaded_pre, g_h_sunlit_pre, g_h_shaded_pre, g_r_sunlit_pre, g_r_shaded_pre, r_a_c_pre, r_s_c_pre = func.CalConductance(vegetation_type, leaf_type, d_leaf, R, Karman_const, SteBol_const, emissivity_c, k_b_black, k_d_black, k_n, k_u, g_0_sunlit, g_0_shaded, porosity, W, W_capacity, z, c_p, u_z, u_attenuation_coef, LAI, height_c, pressure, T_a_K, LAI_sunlit, LAI_shaded, T_c_K_soil, T_c_K_sunlit, T_c_K_shaded, g_s_sunlit, g_s_shaded, g_b_free_sunlit_pre , g_b_free_shaded_pre , g_b_free_soil_pre , G_0_sunlit_pre , G_0_shaded_pre , g_b_forced_sunlit_pre , g_b_forced_shaded_pre , g_b_sunlit_pre , g_b_shaded_pre , g_h_sunlit_pre , g_h_shaded_pre , g_r_sunlit_pre , g_r_shaded_pre , r_a_c_pre , r_s_c_pre)
                                            gamma_modified_sunlit, gamma_modified_shaded, C_s_sunlit, C_s_shaded, VPD_s_sunlit, VPD_s_shaded = func.ResConst(C_a_out , VPD_a , gamma , Delta , g_s_sunlit , g_s_shaded , C_i_sunlit , C_i_shaded , diff_T_sunlit , diff_T_shaded , g_a , g_b_forced_sunlit , g_b_free_sunlit , g_b_sunlit , g_b_forced_shaded , g_b_free_shaded , g_b_shaded , g_h_sunlit , g_h_shaded , g_r_sunlit , g_r_shaded)
                                            k_n, Gamma_respiration_sunlit, Gamma_respiration_shaded, C_i_sunlit, C_i_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, f_w, Gamma_respiration_sunlit_pre, Gamma_respiration_shaded_pre, W_retention = func.SimultaneousPhotosynthesis(vegetation_type, leaf_type, reflectance, transmissivity, scattering_coefficient, reflectance_b, reflectance_d, k_b, k_d, k_b_black, absorptance, C_a_out, V_c_max_org, Ha_V, J_m_25, Ha_J, Hd_J, Delta_S, Gamma_25, a_1, VPD_0, T_cold, R, f, a, W_capacity, W_wilting, W, R_s_b, R_s_d, R_s_b_total, R_s_d_total, R_s_total, pressure, S_max, T_a_C_mean, LAI, LAI_g, LAI_sunlit, LAI_shaded, Q_sunlit, C_s_shaded, C_i_sunlit, C_i_shaded, T_c_C_sunlit, T_c_K_sunlit, T_c_C_shaded, T_c_K_shaded, G_0_sunlit, G_0_shaded, h_beta, C_s_sunlit, VPD_s_sunlit, VPD_s_shaded, J_sunlit_pre, J_shaded_pre, V_j_sunlit_pre, V_j_shaded_pre, Gamma_respiration_sunlit_pre, Gamma_respiration_shaded_pre)
                                            r_a, r_a_soil, G_0_sunlit, G_0_shaded, g_a, g_b_free_sunlit, g_b_forced_sunlit, g_b_sunlit, g_b_forced_shaded, g_b_free_shaded, g_b_shaded, g_h_sunlit, g_h_shaded, g_b_free_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, h_beta, g_b_free_sunlit_pre, g_b_free_shaded_pre, g_b_free_soil_pre, G_0_sunlit_pre, G_0_shaded_pre, g_b_forced_sunlit_pre, g_b_forced_shaded_pre, g_b_sunlit_pre, g_b_shaded_pre, g_h_sunlit_pre, g_h_shaded_pre, g_r_sunlit_pre, g_r_shaded_pre, r_a_c_pre, r_s_c_pre = func.CalConductance(vegetation_type, leaf_type, d_leaf, R, Karman_const, SteBol_const, emissivity_c, k_b_black, k_d_black, k_n, k_u, g_0_sunlit, g_0_shaded, porosity, W, W_capacity, z, c_p, u_z, u_attenuation_coef, LAI, height_c, pressure, T_a_K, LAI_sunlit, LAI_shaded, T_c_K_soil, T_c_K_sunlit, T_c_K_shaded, g_s_sunlit, g_s_shaded, g_b_free_sunlit_pre , g_b_free_shaded_pre , g_b_free_soil_pre , G_0_sunlit_pre , G_0_shaded_pre , g_b_forced_sunlit_pre , g_b_forced_shaded_pre , g_b_sunlit_pre , g_b_shaded_pre , g_h_sunlit_pre , g_h_shaded_pre , g_r_sunlit_pre , g_r_shaded_pre , r_a_c_pre , r_s_c_pre)
                                            Q_n_sunlit, Q_n_shaded, Q_n_sunlit_pre, Q_n_shaded_pre, R_n_soil, G, R_n_sum, lE_soil, lE_sum, ET_c_sum, ET_soil, ET_sum, H_sum = func.CalET_SW(emissivity_c, emissivity_s, extinction_coef, albedo_soil, SteBol_const, T_a_K, R_s_total, VPD_a, rho_a, C_p, c_p, gamma, Delta, S_max, L_down, COSTHETA, LAI, LAI_g, T_c_K_sunlit, T_c_K_shaded, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, r_a, r_a_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, Q_n_sunlit_pre, Q_n_shaded_pre)
                                            diff_T_sunlit, diff_T_shaded = func.CalDiffT(VPD_a, c_p, Delta, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, g_h_sunlit, g_h_shaded, g_r_sunlit, g_r_shaded, gamma_modified_sunlit, gamma_modified_shaded)
                                            #####################################
    
                                            diff_T_sunlit_pre = diff_T_sunlit
                                            diff_T_shaded_pre = diff_T_shaded
                                            C_s_sunlit_pre = C_s_sunlit
                                            C_s_shaded_pre = C_s_shaded
                                            if abs(diff_T_shaded - diff_T_shaded_pre) < 0.01:
                                                break
                                            if count_i > 10000:
                                                break

                                else: # LAI <= 0
                                    T_c_C_soil = T_a_C
                                    T_c_K_soil = T_c_C_soil + 273.15
                                    Q_n_sunlit = 0
                                    Q_n_shaded = 0
                                    r_a, r_a_soil, G_0_sunlit, G_0_shaded, g_a, g_b_free_sunlit, g_b_forced_sunlit, g_b_sunlit, g_b_forced_shaded, g_b_free_shaded, g_b_shaded, g_h_sunlit, g_h_shaded, g_b_free_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, h_beta, g_b_free_sunlit_pre, g_b_free_shaded_pre, g_b_free_soil_pre, G_0_sunlit_pre, G_0_shaded_pre, g_b_forced_sunlit_pre, g_b_forced_shaded_pre, g_b_sunlit_pre, g_b_shaded_pre, g_h_sunlit_pre, g_h_shaded_pre, g_r_sunlit_pre, g_r_shaded_pre, r_a_c_pre, r_s_c_pre = func.CalConductance(vegetation_type, leaf_type, d_leaf, R, Karman_const, SteBol_const, emissivity_c, k_b_black, k_d_black, k_n, k_u, g_0_sunlit, g_0_shaded, porosity, W, W_capacity, z, c_p, u_z, u_attenuation_coef, LAI, height_c, pressure, T_a_K, LAI_sunlit, LAI_shaded, T_c_K_soil, T_c_K_sunlit, T_c_K_shaded, g_s_sunlit, g_s_shaded, g_b_free_sunlit_pre , g_b_free_shaded_pre , g_b_free_soil_pre , G_0_sunlit_pre , G_0_shaded_pre , g_b_forced_sunlit_pre , g_b_forced_shaded_pre , g_b_sunlit_pre , g_b_shaded_pre , g_h_sunlit_pre , g_h_shaded_pre , g_r_sunlit_pre , g_r_shaded_pre , r_a_c_pre , r_s_c_pre)
                                    Q_n_sunlit, Q_n_shaded, Q_n_sunlit_pre, Q_n_shaded_pre, R_n_soil, G, R_n_sum, lE_soil, lE_sum, ET_c_sum, ET_soil, ET_sum, H_sum = func.CalET_SW(emissivity_c, emissivity_s, extinction_coef, albedo_soil, SteBol_const, T_a_K, R_s_total, VPD_a, rho_a, C_p, c_p, gamma, Delta, S_max, L_down, COSTHETA, LAI, LAI_g, T_c_K_sunlit, T_c_K_shaded, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, r_a, r_a_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, Q_n_sunlit_pre, Q_n_shaded_pre)
                                
                                R_n_sum = Q_n_sunlit + Q_n_shaded + R_n_soil
                                A_sum = V_n_sunlit + V_n_shaded
                                A_n_sum = A_c_sunlit + A_c_shaded
                                R_m_leaf_sum = R_d_sunlit + R_d_shaded

                                W, W_retention, G_water, SW_retention, W_pre = func.WaterBalance(t_step, W, W_pre, W_capacity, W_wilting, rainfall, ET_sum)
                                
                                #####################
                                # InitialTemp()
                                T_c_C_list = [0]* (N_layer + 1); g_s_W_top = 0; g_s_W_middle = 0; ET_eq = 0; g_s_W_sum = 0
                                print(f'{hr}:{min},{T_a_C},{R_s_total},{rainfall},{rh},{u_z},{ET_soil},{lE_soil},{ET_sum},{lE_sum},{A_sum},{R_m_leaf_sum},{T_c_C_list[N_layer]},{T_c_C_list[N_layer//2]},{g_s_W_top},{g_s_W_middle},{H_sum},{G},{R_n_sum},{ET_eq},{g_s_W_sum},{sun_duration}', file = CANOPY)

                                A_n_sum_daily += A_n_sum * 60 * t_step # (μmol/m^2/dy)
                                R_m_leaf_sum_daily += R_m_leaf_sum * 60 * t_step # (μmol/m^2/dy)
                                A_sum_daily += A_sum * 60 * t_step # (μmol/m^2/dy)
                                ET_daily += ET_sum * 60 * t_step
                                ET_c_daily += ET_c_sum * 60 * t_step
                                ET_eq_daily += ET_eq * 60 * t_step
                        
                        result = func.CalLAI(vegetation_type, leaf_type, k_n, R_g_parameter, allocation_parameter, allocation_parameter_index, r_00, s_00, sensitivity_allocation, R_m_base_stem, R_m_base_root, death_rate_leaf, push_down_rate, loss_rate_W_max, loss_rate_T_max, b_W, b_T, T_cold, loss_rate_stem, loss_rate_root, loss_rate_standby_leaf, loss_rate_standby_stem, loss_rate_standby_root, virtual_LAI_succeed_day, dormancy_terminate_day, leaf_onset, leaf_normal, leaf_dormant, C_increase_dy, phenophase, C_leaf, Cg_leaf, Cd_leaf, C_stem, C_root, C_all, phase2to3_dy, phase3to2_dy, phase3to2_dy2, phaseto4_dy, W, W_wilting, W_retention, virtual_LAI, virtual_LAI_day, dormancy_dys, pass200, dy_month, DOY_max, growing_dys, loss_rate, y2, month, day, DOY, virtual_DOY, virtual_T_a_C_mean, virtual_R_s_sum, virtual_pressure_mean, virtual_rainfall_sum, virtual_rh_mean, virtual_u_z_mean, virtual_W_mean, virtual_A_sum_daily, virtual_R_a_c, virtual_R_m_leaf_sum_daily, virtual_ET_daily, virtual_ET_c_daily, virtual_ET_eq_daily, virtual_LAI_list, virtual_LAI_g, virtual_C_leaf, virtual_Cg_leaf, virtual_Cd_leaf, virtual_C_stem, virtual_C_root, virtual_C_all, virtual_A_n, LAI, LAI_g, A_sum_daily, A_n_sum_daily, R_m_leaf_sum_daily, ET_daily, ET_c_daily, ET_eq_daily, R_s_sum, rainfall_sum, sun_duration_pos, T_a_C_mean, pressure_mean, rh_mean, u_z_mean, W_mean)
                        if result[0] == True:
                            continueLoop, dormancy_dys, growing_dys, y2, month, day, LAI, virtual_LAI, virtual_LAI_day, A_sum_daily, A_n_sum_daily, R_m_leaf_sum_daily, virtual_A_n = result
                            continue
                        if result[0] == False:
                            continueLoop, dormancy_dys, growing_dys, y2, month, day, LAI, virtual_LAI, virtual_LAI_day, A_sum_daily, A_n_sum_daily, A_stem_daily, A_root_daily, R_m_leaf_sum_daily, R_g_leaf_daily, virtual_A_n, leaf_onset, leaf_normal, C_increase_dy, R_a_c, L_leaf, L_leaf_d, L_all, Cg_leaf, Cd_leaf, C_leaf, C_stem, C_root, C_all, phenophase, phase2to3_dy, phase3to2_dy, phase3to2_dy2, phaseto4_dy, leaf_dormant, virtual_DOY, virtual_T_a_C_mean, virtual_R_s_sum, virtual_pressure_mean, virtual_rainfall_sum, virtual_rh_mean, virtual_u_z_mean, virtual_W_mean, virtual_A_sum_daily, virtual_R_a_c, virtual_R_m_leaf_sum_daily, virtual_ET_daily, virtual_ET_c_daily, virtual_ET_eq_daily, virtual_LAI_list, virtual_LAI_g, virtual_C_leaf, virtual_Cg_leaf, virtual_Cd_leaf, virtual_C_stem, virtual_C_root, virtual_C_all = result
                        
                        if virtual_LAI == 0:
                            print(f'{year}-{itrn}-{DOY}\t{LAI}[{phenophase}]', end = '\t') 
                            if leaf_onset == 1:
                                print('*', end='')
                            print(f'\n\t{A_sum_daily*1000}\t{R_a_c*1000}\t{ET_daily}\t{ET_c_daily}\t{W_mean}')
                            print(f'{DOY},{T_a_C_mean},{R_s_sum},{pressure_mean},{rainfall_sum},{rh_mean},{u_z_mean},{W_mean}', file = MEANCLIM)
                            print(f'{DOY},{A_n_obs_sum},{A_sum_daily*1000},{R_a_c*1000},{R_m_leaf_sum_daily*1000},{ET_daily},{ET_c_daily},{ET_eq_daily}', file = DAILY)
                            print(f'{DOY},{LAI},{C_leaf},{C_stem},{C_root},{C_all},{LAI_g},{phenophase}', file = ALLOCATION)

                            A_sum_yearly += A_sum_daily
                            R_m_leaf_sum_yearly += R_m_leaf_sum_daily
                            R_a_c_yearly += R_a_c
                            ET_yearly += ET_daily
                            ET_c_yearly += ET_c_daily
                            ET_eq_yearly += ET_eq_daily
                        elif virtual_LAI == 10:
                            virtual_LAI_day = 1
                            while virtual_LAI_day <= virtual_LAI_succeed_day:
                                print(f'{virtual_DOY[virtual_LAI_day-1]},{virtual_A_sum_daily[virtual_LAI_day-1]*1000},{virtual_R_a_c[virtual_LAI_day-1]*1000},{virtual_ET_daily[virtual_LAI_day-1]},{virtual_ET_c_daily[virtual_LAI_day-1]},{virtual_ET_eq_daily[virtual_LAI_day-1]}')
                                print(f'{virtual_DOY[virtual_LAI_day-1]},{virtual_T_a_C_mean[virtual_LAI_day-1]},{virtual_R_s_sum[virtual_LAI_day-1]},{virtual_pressure_mean[virtual_LAI_day-1]},{virtual_rainfall_sum[virtual_LAI_day-1]},{virtual_rh_mean[virtual_LAI_day-1]},{virtual_u_z_mean[virtual_LAI_day-1]},{virtual_W_mean[virtual_LAI_day-1]}', file=MEANCLIM)
                                print(f'{virtual_DOY[virtual_LAI_day-1]},{A_n_obs_sum},{virtual_A_sum_daily[virtual_LAI_day-1]*1000},{virtual_R_a_c[virtual_LAI_day-1]*1000},{virtual_R_m_leaf_sum_daily[virtual_LAI_day-1]*1000},{virtual_ET_daily[virtual_LAI_day-1]},{virtual_ET_c_daily[virtual_LAI_day-1]},{virtual_ET_eq_daily[virtual_LAI_day-1]}', file=DAILY)
                                print(f'{virtual_DOY[virtual_LAI_day-1]},{virtual_LAI_list[virtual_LAI_day-1]},{virtual_C_leaf[virtual_LAI_day-1]},{virtual_C_stem[virtual_LAI_day-1]},{virtual_C_root[virtual_LAI_day-1]},{virtual_C_all[virtual_LAI_day-1]},{virtual_LAI_g[virtual_LAI_day-1]},1', file=ALLOCATION)
                                
                                A_sum_yearly += virtual_A_sum_daily[virtual_LAI_day-1]
                                R_m_leaf_sum_yearly += virtual_R_m_leaf_sum_daily[virtual_LAI_day-1]
                                R_a_c_yearly += virtual_R_a_c[virtual_LAI_day-1]
                                ET_yearly += virtual_ET_daily[virtual_LAI_day-1]
                                ET_c_yearly += virtual_ET_c_daily[virtual_LAI_day-1]
                                ET_eq_yearly += virtual_ET_eq_daily[virtual_LAI_day-1]
                                if virtual_DOY[virtual_LAI_day-1] == DOY_max:
                                    print(f'{DOY},{A_sum_yearly},{R_m_leaf_sum_yearly},{R_a_c_yearly},{ET_yearly},{ET_c_yearly},{W_mean}')
                                    NPP_yearly = (A_sum_yearly - R_a_c_yearly) * 10  # (tC/ha)
                                    itrn += 1
                                    if itrn >= 3 and DOY == DOY_max:
                                        BreakLOOP = True
                                        break
                                    A_sum_yearly = 0
                                    R_a_c_yearly = 0
                                    R_m_leaf_sum_yearly = 0
                                    NPP_yearly = 0
                                    ET_yearly = 0
                                    ET_c_yearly = 0
                                    ET_eq_yearly = 0
                                virtual_LAI_day += 1 #Counting virtual_LAI_day
                            virtual_LAI = 0
                            leaf_onset = 1
                            phenophase = 1
                            growing_dys = virtual_LAI_succeed_day
                        
                        print(f'>{growing_dys}\t{Cg_leaf}({LAI_g})\t{Cd_leaf}({LAI_d})\t{C_leaf}({LAI})')

                #print(f'#SUM, {A_sum_yearly}, {R_m_leaf_sum_yearly}, {R_a_c_yearly}, {ET_yearly}, {ET_c_yearly}, {ET_eq_yearly}', file=DAILY)
                print(f'{year}, {A_sum_yearly*10}, {R_a_c_yearly*10}, {NPP_yearly}, {ET_yearly}, {ET_c_yearly}, {ET_eq_yearly}', file=YEARLY)

            sum_C_stem += C_stem
            sum_C_root += C_root
        
        ave_C_stem = sum_C_stem / (year_end - year_start + 1)
        ave_C_root = sum_C_root / (year_end - year_start + 1)

        diff_C_stem = ave_C_stem - ave_C_stem_pre
        diff_C_root = ave_C_root - ave_C_root_pre

        if diff_C_stem > 0.5:
            diff_count_yr += 1
        if diff_C_root > 0.5:
            diff_count_yr += 1

        print(f'{itrn}\t{diff_C_stem}\t{diff_C_root}')

        if diff_count_yr == 0 and DOY == DOY_max:
            BreakLOOP = True

        ave_C_stem_pre = ave_C_stem
        ave_C_root_pre = ave_C_root

end_time = time.time()
min = (end_time - start_time) // 60
sec = (end_time - start_time) % 60
print(f'Time: {min} min {sec} sec')