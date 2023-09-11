import math
import numpy as np
import os
from typing import Tuple
from numba import jit

# functions ############################################
def MakeDirectory(site: str, year: int):
        directories = [
            f'./../out/{site}/{year}',
            f'./../out/{site}/{year}/canopy',
            f'./../out/{site}/{year}/daily',
            f'./../out/{site}/{year}/meanclim',
        ]

        for directory in directories:
            os.makedirs(directory, exist_ok=True)

def inputClimateData(t: int, t_start: float, t_step: int, elv: int, flux_data: list, solar_elevation: float, a0: float, b0: float, c0: float, R: float, M_d: float, C_pd: int, k_b_black: float, DOY_max: int, DOY: int, latitude: float, refraction: float, lon: float, lon_LST: float, SteBol_const: float, sun_duration_sum: int, sun_duration_count: int, T_a_C_pre: float, rainfall_pre: float, rh_pre: float, u_z_pre: float, R_s_total_pre: float):
    if flux_data[3] == -9999:
        T_a_C = T_a_C_pre
    else:
        T_a_C = flux_data[3] # air temperature (℃)
    T_a_K = T_a_C + 273.15 # air temperature (K)
    
    if flux_data[2] == -9999:
        rainfall = rainfall_pre
    else:
        rainfall = flux_data[2] # * 60 * t_step # (mm/30min)

    if flux_data[4] == -9999:
        rh = rh_pre
    else:
        rh = max(flux_data[4], 0.0001)

    if flux_data[5] == -9999:
        u_z = u_z_pre
    else: 
        u_z = max(flux_data[5], 0.1)

    # pressure = flux_data[9] / 100 # (hPa)
    pressure = 1013.25*(1-0.0065/T_a_K*elv)**5.2552

    if flux_data[6] == -9999:
        R_s_total = R_s_total_pre
    else:
        R_s_total = flux_data[6] * 277.78 # (W/m^2)
    # if flux_data[19] == -9999:
    #     A_n_obs = 0
    # else:
    #     A_n_obs = flux_data[19]
    A_n_obs = 0
    sun_duration = 0
    solar_elevation_pre = solar_elevation

    e_a_Ta, q, VPD_a, rho_m, rho_a, C_p, c_p, gamma, Delta = CalPressVPD(a0, b0, c0, R, M_d, C_pd, T_a_C, T_a_K, pressure, rh)
    sun_duration_pos, COSTHETA, R_s_d_total, sun_duration, S_max, R_l_f, solar_elevation, k_b, k_d, k_b_black, k_d_black = CalSunDuration(e_a_Ta, DOY_max, DOY, latitude, refraction, lon, lon_LST, t, pressure, sun_duration, T_a_K, SteBol_const, R_s_total, k_b_black)

    R_s_b_total = R_s_total - R_s_d_total
    if solar_elevation > 0:
        if solar_elevation_pre < 0:
            sun_duration_sum = 0
            sun_duration_count = 0
        sun_duration_sum += sun_duration
        sun_duration_count += 1
    else:
        if sun_duration_count == -9:
            sun_duration = 0.8
        else:
            sun_duration = sun_duration_sum / sun_duration_count
    if S_max > 0:
        sun_duration = R_s_b_total / (S_max*(1-0.156))
    if sun_duration < 0:
        sun_duration = 0
    elif sun_duration > 1:
        sun_duration = 1
    
    C = 0.826*sun_duration**3 - 1.234*sun_duration**2 + 1.135*sun_duration + 0.298
    if C > 1:
        C = 1
    L_down = SteBol_const*T_a_K**4*(1-(1-R_l_f/(SteBol_const*T_a_K**4))*C) # downward longwave radiation ~ (Kondo & Xu, 1997)

    R_s_b, R_s_d = [0, 0], [0, 0]

    R_s_b[0] = R_s_b_total*0.43 # PAR content: 43%; 4.6 μmol/J * 43 % (~ Larcher, 1995)
    R_s_d[0] = R_s_d_total*0.57 # PAR content: 57%; 4.2 μmol/J * 57 %
    R_s_b[1] = R_s_b_total*0.57 # 4.6 μmol/J * 43 % (~ Larcher, 1995)
    R_s_d[1] = R_s_d_total*0.43 # 4.2 μmol/J * 57 %

    T_a_C_pre = T_a_C
    rainfall_pre = rainfall
    rh_pre = rh
    u_z_pre = u_z
    R_s_total_pre = R_s_total

    return T_a_C, T_a_K, rainfall, u_z, pressure, R_s_total, A_n_obs, rh, VPD_a, rho_a, C_p, c_p, gamma, Delta, sun_duration_pos, COSTHETA, R_s_d_total, sun_duration, S_max, R_s_b_total, solar_elevation, k_b, k_d, k_b_black, k_d_black, L_down, R_s_b, R_s_d, T_a_C_pre, rainfall_pre, rh_pre, u_z_pre, R_s_total_pre

def CalPressVPD(a0: float, b0: float, c0: float, R: float, M_d: float, C_pd: int, T_a_C: float, T_a_K: float, pressure: float, rh: float):
    e_s_Ta = c0 * np.exp(a0 * T_a_C / (T_a_C + b0)) # saturation vapor pressure (hPa)
    e_a_Ta = e_s_Ta * rh / 100 # vapor pressure (hPa)
    q = 0.622 * e_a_Ta / pressure - 0.379*e_a_Ta # specific humidity (kg/kg)
    # vapor pressure deficit (hPa)
    VPD_a = e_s_Ta - e_a_Ta
    # mole fraction of vapor
    e_frac_s = e_s_Ta / pressure
    e_frac_a = e_a_Ta / pressure
    VPD_frac = e_frac_s - e_frac_a

    # rho_m: molar density (mol/m^3)
    rho_m = pressure * 100 / R / T_a_K
    # rho_a: air density (kg/m^3)
    rho_a = rho_m * M_d * (1 - 0.378 * e_a_Ta / pressure)
    # M_a: molecular mass of air (kg/mol)
    M_a = rho_a / rho_m
    # C_p: specific heat of air at constant pressure (J/kg/K)
    C_p = C_pd * (1 + 0.84 * q)
    # c_p (J/mol/K)
    c_p = C_pd * (1 + 0.84 * q) * M_a
    # latent heat of vaporization (J/mol)
    lambda_mol = 56780.3 - 42.84 * T_a_K
    lambda_kg = lambda_mol / 0.018
    # gamma: psychrometric constant (hPa/K)
    gamma = C_p * pressure / 0.622 / lambda_kg

    Delta = a0 * b0 * c0 / (T_a_C + b0) ** 2 * np.exp(a0 * T_a_C / (T_a_C + b0)) # slope of the saturated vapor pressure (hPa/°C)

    return e_a_Ta, q, VPD_a, rho_m, rho_a, C_p, c_p, gamma, Delta

def CalSunDuration(e_a_Ta: float, DOY_max: int, DOY: int, latitude: float, refraction: float, lon: float, lon_LST: float, t: int, pressure: float, sun_duration_pre: float, T_a_K: float, SteBol_const: float, R_s_total: float, k_b_black: float): # to calculate possible sunshine duration

    albedo, dust, LW10 = vapor(e_a_Ta)

    eta = (2 * math.pi / DOY_max) * DOY
    SA2 = 4.871 + eta + 0.033 * math.sin(eta)
    A22 = 0.398 * math.sin(SA2)
    declination = math.atan(A22 / math.sqrt(1 - A22**2))  # declination (rad)
    DD2 = 1.00011 + 0.034221*math.cos(eta) + 0.00128*math.sin(eta) + 0.000719*math.cos(2*eta) + 0.000077*math.sin(2*eta)
    j1 = (0.066 + 0.34*math.sqrt(dust))*(albedo - 0.15)
    C1 = 0.21 - 0.2*dust
    C2 = 0.15 - 0.2*dust
    if dust > 0.3: 
        C1 = 0.15
        C2 = 0.09
    Io = 1365*DD2
    F1 = 0.056 + 0.16*math.sqrt(dust)
    F2 = 0.075 + 0.65*dust

    hour_angle = math.asin(math.sqrt(math.sin(math.pi/4+(latitude-declination+refraction)/2)*math.sin(math.pi/4-(latitude-declination-refraction)/2)/(math.cos(latitude)*math.cos(declination))))
    sun_duration_pos = 4*hour_angle/0.2618

    t_local = t + 4*(lon-lon_LST)
    H = math.pi*(t_local/60/12-1)
    COSTHETA = math.sin(latitude)*math.sin(declination) + math.cos(latitude)*math.cos(declination)*math.cos(H)
    solar_elevation = math.asin(COSTHETA)

    So = Io*COSTHETA
    M = pressure/1013.25/COSTHETA
    if COSTHETA <= 0:
        I = 0
        So = 0
        S = 0
        S_max = 0
        sun_duration = 0
        R_s_total = 0
        R_s_d_total = 0
    else:
        I1 = 0.014*(M+7+2*LW10)*LW10
        I2 = 0.02*(M+5.5+1.5*LW10)*LW10
        I = Io*(C2+0.75*10**(-M*F2))*(1-I2)
        S_max = So*(C1+0.7*10**(-M*F1))*(1-I1)*(1+j1) # 全天日射量 (W/m^2)
        sun_duration = R_s_total/S_max
        K_Tt = R_s_total/So
        if K_Tt < 0:
            K_Tt = 0
        if K_Tt <= 0.15:
            R_s_d_total = R_s_total
        elif K_Tt <= 0.22:
            R_s_d_total = R_s_total*(1-0.09*K_Tt)
        elif K_Tt <= 0.80:
            R_s_d_total = R_s_total*(0.9511-0.1604*K_Tt+4.388*K_Tt**2-16.638*K_Tt**3+12.336*K_Tt**4)
        else:
            R_s_d_total = R_s_total*0.156

    if R_s_d_total > R_s_total:
        R_s_d_total = R_s_total
    elif R_s_total > 0 and R_s_d_total <= 0:
        R_s_d_total = 0.01

    # $R_l_f: downward longwave radiation under clear sky (W/m^2) ~ (Zuo et al., 1991)
    xi = 46.5 * e_a_Ta / T_a_K
    emissivity_a = 1 - (1 + xi) * np.exp(-math.sqrt(1.2 + 3 * xi))
    R_l_f = SteBol_const * T_a_K ** 4 * emissivity_a

    if solar_elevation > 0:
        k_b_black = 0.5 / math.sin(solar_elevation)

    if k_b_black > 1 or solar_elevation <= 0:
        k_b_black = 1

    k_b = math.sqrt(0.5) * k_b_black
    k_d_black = 0.8
    k_d = math.sqrt(0.5) * k_d_black

    return sun_duration_pos, COSTHETA, R_s_d_total, sun_duration, S_max, R_l_f, solar_elevation, k_b, k_d, k_b_black, k_d_black

def vapor(e_a_Ta: float) -> Tuple[float, float, float]:

    vapor_pressure = e_a_Ta # (hPa)
    albedo = 0.15
    dust = 0.05 # 混濁係数

    dew_point = (27 * math.log(vapor_pressure / 6.11)) / (17.27 - math.log(vapor_pressure / 6.11)) # (℃)
    Bo = 0
    if dew_point < -5:
        LNW = 0.0622 * dew_point + 1.958 - Bo
    elif dew_point < 23:
        LNW = 0.0714 * dew_point + 2.003 - Bo
    else:
        LNW = 0.0345 * dew_point + 2.851 - Bo
    WSTOP = np.exp(LNW)
    LW10 = 0.4343 * math.log((1.234 * WSTOP - 0.21) / 10)

    return albedo, dust, LW10

def LAdistribution(N_layer: int, d_z_h: float, alpha: list, beta: list, leaf_type: int, height_c: float, d_z: float, Cg_leaf: float, C_leaf: float, LAI: float) -> list: # describe leaf area distribution for the multi-layer canopy model
    sum = 0
    for i_l in range(N_layer, 0, -1):
        x = d_z_h * (2* i_l - 1) / 2
        sum += x ** (alpha[leaf_type] - 1) * (1 - x) ** (beta[leaf_type] - 1) * d_z_h

    f_sum = 0
    LA_z_max = 0
    f_list = [0]*(N_layer+1)
    z_c = [0]*(N_layer+2)
    LA_z = [0]*(N_layer+1)
    d_LA = [0]*(N_layer+1)
    d_LA_g = [0]*(N_layer+1)
    for i_l in range(N_layer, 0, -1):
        x = d_z_h * (2* i_l - 1) / 2
        f_list[i_l] = x ** (alpha[leaf_type] - 1) * (1 - x) ** (beta[leaf_type] - 1) / sum * LAI # (m^2/m^2)
        z_c[i_l] = x * height_c # (m)
        LA_z[i_l] = f_list[i_l] / height_c # LA_z: leaf area density (m^2/m^3)
        d_LA[i_l] = LA_z[i_l] * d_z # (m^2/m^2)
        d_LA_g[i_l] = d_LA[i_l] * Cg_leaf / C_leaf # (m^2/m^2)
        LA_z_max = max(LA_z_max, LA_z[i_l])

    z_c[N_layer+1]  = d_z_h * (2* (N_layer + 1) - 1) / 2 * height_c

    return f_list, z_c, LA_z, d_LA, d_LA_g, LA_z_max, z_c

def ExtinctionCoeffcientFunction(var1, var2): #一般の葉群における群落光合成
    return (1 - np.exp(-var2 * var1)) / var2

def IB(I_b: list, i_l: int, d_LA: list, x: float, N_angle: int, N2_angle: int) -> float:
    G_layer = GLayer(N_angle, N2_angle, x)
    I_b[i_l] = 1 - G_layer * d_LA[i_l] / math.sin(x)
    I_b[i_l] = max(0, min(I_b[i_l], 1))

    return I_b, G_layer

def ID(I_d: list, I_b: list, d_LA: list, i_l: int, N_angle: int, N2_angle: int) -> list:
    I_d[i_l] = 0
    d_angle = math.pi / 2 / N_angle
    for i_a3 in range(1, N_angle+1):
        x = d_angle * (2 * i_a3 - 1) / 2
        I_b, G_layer = IB(I_b, i_l, d_LA, x, N_angle, N2_angle)
        I_d[i_l] += I_b[i_l] * math.sin(x) * math.cos(x) * d_angle
    I_d[i_l] *= 2
    I_d[i_l] = max(0, min(I_d[i_l], 1))
    return I_d

def GLayer(N_angle: int, N2_angle: int, x: float) -> float:
    angle_ave = 50 * math.pi / 180 # angle_ave: the mean of leaf inclination (degree)
    angle_SD = 18 * math.pi / 180 # angle_SD: the standard deviation of leaf inclination (degree)

    G_layer = 0
    d_angle = math.pi / 2 / N_angle
    d2_angle = 2 * math.pi / N2_angle
    for i_a2 in range(1, N2_angle+1):
        G_layer0 = 0
        leaf_orientation_angle = d2_angle * (2 * i_a2 - 1) / 2
        for i_a in range(1, N_angle+1):
            leaf_inclination = d_angle * (2 * i_a - 1) / 2 # (rad)
            g_alpha = np.exp(-((leaf_inclination - angle_ave) / angle_SD) ** 2 / 2) / (angle_SD * math.sqrt(2 * math.pi))
            G_leaf = abs(math.cos(leaf_inclination) * math.sin(x) + math.sin(leaf_inclination) * math.cos(x) * math.cos(leaf_orientation_angle)) # G_leaf: leaf projection area
            G_layer0 += G_leaf * g_alpha * d_angle # G_layer0: the ratio of the area of leaves, projected into a plane normal to the solar elevation, to the leaf area index within the layer 
        G_layer += G_layer0 / (2 * math.pi) * d2_angle
    
    return G_layer

def CalI(N_layer: int, solar_elevation: float, N_angle: int, N2_angle: int, d_LA: list) -> Tuple[list, list]:
    I_b = [0]*(N_layer+1)
    I_d = [0]*(N_layer+1)
    for i_l in range(N_layer, 0, -1):
        I_d = ID(I_d, I_b, d_LA, i_l, N_angle, N2_angle)
        x = solar_elevation 
        I_b, G_layer = IB(I_b, i_l, d_LA, x, N_angle, N2_angle)
    
    return I_b, I_d, G_layer

def S_b_balance(N_layer: int, S_b_down: list, I_b: list): # calculate direct radiation balance within each layer
    for i_l in range(N_layer, 0, -1):
        S_b_down[i_l-1] = S_b_down[i_l] * I_b[i_l]

    return S_b_down

def S_d_balance(reflectance: list, transmissivity: list, albedo_soil: list, N_layer: int, S_b_down: list, S_d_down: list, S_d_up: list, L: int, I_b: list, I_d: list): # calculate diffuse radiation balance within each layer
    S_d_up[0] = (S_b_down[0] + S_d_down[0]) * albedo_soil[L]
    for i_l in range(N_layer, 0, -1):
        S_d_down[i_l-1] = S_d_down[i_l] * (transmissivity[L] * (1 - I_d[i_l]) + I_d[i_l]) + S_d_up[i_l - 1] * reflectance[L] * (1 - I_d[i_l]) + S_b_down[i_l] * transmissivity[L] * (1 - I_b[i_l])
        S_d_up[N_layer-(i_l-1)] = S_d_up[N_layer-1-(i_l-1)] * (transmissivity[L] * (1 - I_d[N_layer-(i_l-1)]) + I_d[N_layer-(i_l-1)]) + S_d_down[N_layer-(i_l-1)] * reflectance[L] * (1 - I_d[(N_layer-(i_l-1))]) + S_b_down[N_layer-(i_l-1)] * reflectance[L] * (1 - I_b[N_layer-(i_l-1)])
    return S_d_down, S_d_up

def S_balance(reflectance: list, transmissivity: list, albedo_soil: list, N_layer: int, diff_level: float, S_max: float, R_s_b: list, R_s_d: list, solar_elevation: float, d_LA: list, S_b_down: list, S_d_down: list, S_d_up: list, L: int, I_b: list, I_d: list, G_layer: float, d_LA_su: list, d_LA_sh: list, S_d_down_pre: list, S_d_up_pre: list, diff_S_d_down: list, diff_S_d_up: list, Q_shaded: list, Q_sunlit: list): # calculate radiation balance within each layer
    S_b_down[N_layer] = R_s_b[L]
    S_d_down[N_layer] = R_s_d[L]
    S_b_down = S_b_balance(N_layer, S_b_down, I_b)

    if S_max > 0 and G_layer > 0:
        for i_l in range(N_layer, 0, -1):
            if S_b_down[N_layer] <= 0:
                d_LA_su[i_l] = 0
            else:
                d_LA_su[i_l] = (S_b_down[i_l] - S_d_down[i_l]) / S_b_down[N_layer] * math.sin(solar_elevation) / G_layer
                if d_LA_su[i_l] < d_LA[i_l]:
                    d_LA_su[i_l] = d_LA[i_l]
            d_LA_sh[i_l] = d_LA[i_l] - d_LA_su[i_l]

    for i_l in range(N_layer, 0, -1):
        S_d_down[i_l-1] = S_d_down[N_layer] * i_l / N_layer
        S_d_up[i_l-1] = S_b_down[N_layer] * reflectance[L] * reflectance[L] * i_l / N_layer
    S_d_up[N_layer] = R_s_b[L] * reflectance[L] * reflectance[L]
    
    i = 0
    while True: # iterate until calculations are converged
        diff_count = 0
        for i_l in range(N_layer, 0, -1):
            S_d_down_pre[i_l-1] = S_d_down[i_l-1]
            S_d_up_pre[i_l-1] = S_d_up[i_l-1]
        S_d_up_pre[N_layer] = S_d_up[N_layer]
        S_d_down, S_d_up = S_d_balance(reflectance, transmissivity, albedo_soil, N_layer, S_b_down, S_d_down, S_d_up, L, I_b, I_d)
        for i_l in range(N_layer, 0, -1):
            diff_S_d_down[i_l-1] = abs(S_d_down[i_l-1] - S_d_down_pre[i_l-1])
            diff_S_d_up[i_l-1] = abs(S_d_up[i_l-1] - S_d_up_pre[i_l-1])
        diff_S_d_up[N_layer] = abs(S_d_up[N_layer] - S_d_up_pre[N_layer])

        for i_l in range(N_layer, 0, -1):
            if diff_S_d_down[i_l-1] > diff_level or diff_S_d_up[i_l-1] > diff_level:
                diff_count += 1
            if diff_S_d_up[N_layer] > diff_level:
                diff_count += 1
        
        if diff_count == 0:
            break

        i += 1

    Q_shaded[L] = 0
    Q_sunlit[L] = 0
    for i_l in range(N_layer, 0, -1):
        if S_b_down[N_layer] <= 0:
            d_LA_su[i_l] = 0
        else:
            d_LA_su[i_l] = (S_b_down[i_l] - S_b_down[i_l - 1]) / S_b_down[N_layer] * math.sin(solar_elevation) / G_layer
            if d_LA_su[i_l] > d_LA[i_l]:
                d_LA_su[i_l] = d_LA[i_l]
        d_LA_sh[i_l] = d_LA[i_l] - d_LA_su[i_l]
        Q_shaded[L] += (1 - transmissivity[L]) * (1 - transmissivity[L] - reflectance[L]) * (S_d_down[i_l] + S_d_up[i_l - 1]) * d_LA_sh[i_l]
        Q_sunlit[L] += (1 - transmissivity[L]) * (1 - transmissivity[L] - reflectance[L]) * (S_d_down[i_l] + S_d_up[i_l - 1] + S_b_down[N_layer] * G_layer / math.sin(solar_elevation)) * d_LA_sh[i_l]
    return d_LA_su , d_LA_sh , S_d_down_pre , S_d_up_pre , diff_S_d_down , diff_S_d_up , Q_shaded , Q_sunlit

def L_balance(emissivity_c: float, k_b_black: float, k_d_black: float, SteBol_const: float, LAI: float, R_s_total: float, T_a_K: float, L_down: float, LAI_sunlit: float, LAI_shaded: float):
    extinction_coeffcient_d_black = ExtinctionCoeffcientFunction(LAI, k_d_black)
    extinction_coeffcient_2d_black = ExtinctionCoeffcientFunction(LAI, 2*k_d_black)
    extinction_coeffcient_b_black_plus_d_black = ExtinctionCoeffcientFunction(LAI, k_b_black + k_d_black)
    extinction_coeffcient_b_black_minus_d_black = ExtinctionCoeffcientFunction(LAI, k_b_black - k_d_black)
    Q_long_isothermal_sunlit = (L_down - emissivity_c * SteBol_const * T_a_K ** 4) * LAI_sunlit / LAI
    Q_long_isothermal_shaded = (L_down - emissivity_c * SteBol_const * T_a_K ** 4) * LAI_shaded / LAI
    if LAI < 1:
        Q_long_isothermal_sunlit *= LAI
        Q_long_isothermal_shaded *= LAI
    if abs(Q_long_isothermal_sunlit) > R_s_total:
        Q_long_isothermal_sunlit = 0
    if abs(Q_long_isothermal_shaded) > R_s_total:
        Q_long_isothermal_shaded = 0

    return Q_long_isothermal_sunlit, Q_long_isothermal_shaded

def CalRn(Q_sunlit: list, Q_shaded: list, Q_long_isothermal_sunlit: float, Q_long_isothermal_shaded: float):
    Q_n_isothermal_sunlit = Q_sunlit[0] + Q_sunlit[1] + Q_long_isothermal_sunlit
    Q_n_isothermal_shaded = Q_shaded[0] + Q_shaded[1] + Q_long_isothermal_shaded

    return Q_n_isothermal_sunlit, Q_n_isothermal_shaded

def ResConst(C_a_out: int, VPD_a: float, gamma: float, Delta: float, g_s_sunlit: float, g_s_shaded: float, C_i_sunlit: float, C_i_shaded: float, diff_T_sunlit: float, diff_T_shaded: float, g_a: float, g_b_forced_sunlit: float, g_b_free_sunlit: float, g_b_sunlit: float, g_b_forced_shaded: float, g_b_free_shaded: float, g_b_shaded: float, g_h_sunlit: float, g_h_shaded: float, g_r_sunlit: float, g_r_shaded: float):
    g_w_sunlit = 1/(1/g_a + 1/g_b_sunlit + 1/g_s_sunlit)
    g_w_shaded = 1/(1/g_a + 1/g_b_shaded + 1/g_s_shaded)

    gamma_modified_sunlit = gamma*(g_h_sunlit + g_r_sunlit)/g_w_sunlit
    gamma_modified_shaded = gamma*(g_h_shaded + g_r_shaded)/g_w_shaded

    g_c_sunlit = 1/(1/g_a + 1/(0.11/0.147*g_b_forced_sunlit + 0.038/0.055*g_b_free_sunlit) + 1/0.64/g_s_sunlit)
    g_c_shaded = 1/(1/g_a + 1/(0.11/0.147*g_b_forced_shaded + 0.038/0.055*g_b_free_shaded) + 1/0.64/g_s_shaded)
    C_s_sunlit = C_i_sunlit + (C_a_out - C_i_sunlit)*g_c_sunlit/0.64/g_s_sunlit
    C_s_shaded = C_i_shaded + (C_a_out - C_i_shaded)*g_c_shaded/0.64/g_s_shaded

    VPD_s_sunlit = (VPD_a + Delta*diff_T_sunlit)*g_w_sunlit/g_s_sunlit
    VPD_s_shaded = (VPD_a + Delta*diff_T_shaded)*g_w_shaded/g_s_shaded

    return gamma_modified_sunlit, gamma_modified_shaded, C_s_sunlit, C_s_shaded, VPD_s_sunlit, VPD_s_shaded

def CalDiffT(VPD_a: float, c_p: float, Delta: float, Q_n_isothermal_sunlit: float, Q_n_isothermal_shaded: float, g_h_sunlit: float, g_h_shaded: float, g_r_sunlit: float, g_r_shaded: float, gamma_modified_sunlit: float, gamma_modified_shaded: float):
    diff_T_sunlit = gamma_modified_sunlit / (Delta + gamma_modified_sunlit) * Q_n_isothermal_sunlit / c_p / (g_h_sunlit + g_r_sunlit) - VPD_a / (Delta + gamma_modified_sunlit)
    diff_T_shaded = gamma_modified_shaded / (Delta + gamma_modified_shaded) * Q_n_isothermal_shaded / c_p / (g_h_shaded + g_r_shaded) - VPD_a / (Delta + gamma_modified_shaded)

    return diff_T_sunlit, diff_T_shaded

def CalConductance(vegetation_type: int, leaf_type: int, d_leaf: list, R: float, Karman_const: float, SteBol_const: float, emissivity_c: float, k_b_black: float, k_d_black: float, k_n: float, k_u: float, g_0_sunlit: float, g_0_shaded: float, porosity: float, W: float, W_capacity: int, z: float, c_p: float, u_z: float, u_attenuation_coef: float, LAI: float, height_c: float, pressure: float, T_a_K: float, LAI_sunlit: float, LAI_shaded: float, T_c_K_soil: float, T_c_K_sunlit: float, T_c_K_shaded: float, g_s_sunlit: float, g_s_shaded: float, g_b_free_sunlit_pre : float, g_b_free_shaded_pre : float, g_b_free_soil_pre : float, G_0_sunlit_pre : float, G_0_shaded_pre : float, g_b_forced_sunlit_pre : float, g_b_forced_shaded_pre : float, g_b_sunlit_pre : float, g_b_shaded_pre : float, g_h_sunlit_pre : float, g_h_shaded_pre : float, g_r_sunlit_pre : float, g_r_shaded_pre : float, r_a_c_pre : float, r_s_c_pre: float):
    if vegetation_type == 1 and leaf_type == 0 and LAI == 0:
        d = 0.01
        z_0_m = 0.01
    else:
        if vegetation_type == 0:
            z_0_m = 0.10 * height_c # for forest ~ (Verseghy et al., 1993)
            z_0_h = z_0_m / 2.0 # for forest
            d = 0.70 * height_c # for forest ~ (Verseghy et al., 1993)
        elif vegetation_type == 1:
            z_0_m = 0.123 * height_c # for crop and grass ~ (Monteith, 1981)
            z_0_h = z_0_m / 12.0 # for grass
            d = 0.67 * height_c # for crop and grass ~ (Monteith, 1981)
        u_h = u_z * math.log((height_c - d) / z_0_m) * math.log((height_c - d) / z_0_h) / math.log((z - d) / z_0_m) / math.log((z - d) / z_0_h)
        u = u_h * np.exp(-u_attenuation_coef * 0.5) # within-canopy wind speed ~ (Cionco, 1972)
        u = max(u, 0.0001)
        r_a_alpha = math.log((z - d) / z_0_m) / (Karman_const ** 2 * u_z) * (math.log((z - d) / (height_c - d)) + height_c / (2.5 * (height_c - d)) * (np.exp(2.5 * (1 - (d + z_0_m) / height_c)) - 1))
        r_s_a_alpha = math.log((z - d) / z_0_m) / (Karman_const ** 2 * u_z) * height_c / (2.5 * (height_c - d)) * (np.exp(2.5) - np.exp(2.5 * (1 - (d + z_0_m) / height_c)))

    z_0_e = 0.01

    u_soil = u_h * np.exp(-u_attenuation_coef)
    u_soil = max(u_soil, 0.0001)

    r_s_a_0 = math.log(z / z_0_e) * math.log((d + z_0_m) / z_0_e) / (Karman_const ** 2 * u_z)
    r_s_a_0 = max(r_s_a_0, 0)

    r_a_0 = (math.log(z / z_0_e)) ** 2 / (Karman_const ** 2 * u_z) - r_s_a_0
    if r_a_0 <= 0:
        r_a_0 = 0.1

    if LAI <= 4:
        r_a = r_a_alpha * LAI / 4 + r_a_0 * (1 - LAI / 4)
        r_a_soil = r_s_a_alpha * LAI / 4 + r_s_a_0 * (1 - LAI / 4)
    else:
        r_a = r_a_alpha
        r_a_soil = r_s_a_alpha
    
    g_a = 1 / r_a * pressure * 100 / R / T_a_K

    if LAI > 0:
        g_b_free_sunlit = 0.055 * ((T_c_K_sunlit - T_a_K) / d_leaf[leaf_type]) ** (1 / 4) * LAI_sunlit # g_b_free: boundary conductance due to free convection for H2O (mol/m2/s)
        g_b_free_shaded = 0.055 * ((T_c_K_shaded - T_a_K) / d_leaf[leaf_type]) ** (1 / 4) * LAI_shaded # g_b_free: boundary conductance due to free convection for H2O (mol/m2/s)
        g_b_free_soil = 0.055 * ((T_c_K_soil - T_a_K) / d_leaf[leaf_type]) ** (1 / 4) # g_b_free: boundary conductance due to free convection for H2O (mol/m2/s)
        extinction_coeffcient_b_black_plus_d_black = ExtinctionCoeffcientFunction(LAI,k_b_black+k_d_black)
        extinction_coeffcient_d_black = ExtinctionCoeffcientFunction(LAI, k_d_black)
        extinction_coeffcient_u_plus_b_black = ExtinctionCoeffcientFunction(LAI, 0.5*k_u+k_b_black)
        extinction_coeffcient_u = ExtinctionCoeffcientFunction(LAI, 0.5*k_u)
        extinction_coeffcient_b_black_plus_n = ExtinctionCoeffcientFunction(LAI, k_b_black + k_n)
        extinction_coefficient_n = ExtinctionCoeffcientFunction(LAI, k_n)
        G_0_sunlit = g_0_sunlit * extinction_coeffcient_b_black_plus_n
        G_0_shaded = g_0_shaded * (extinction_coefficient_n - extinction_coeffcient_b_black_plus_n)
        g_b_forced = 0.147 * math.sqrt(u_h / d_leaf[leaf_type]) # g_b_forced: boundary conductance due to forced convection for H2O (mol/m2/s)
        g_b_forced_sunlit = g_b_forced * extinction_coeffcient_u_plus_b_black
        if extinction_coeffcient_u - extinction_coeffcient_u_plus_b_black < 0.0001:
            g_b_forced_shaded = g_b_forced * 0.0001
        else:
            g_b_forced_shaded = g_b_forced * (extinction_coeffcient_u - extinction_coeffcient_u_plus_b_black)
        if math.isnan(g_b_free_sunlit):
            g_b_free_sunlit = 0
        if math.isnan(g_b_free_shaded):
            g_b_free_shaded = 0
        g_b_sunlit = g_b_free_sunlit + g_b_forced_sunlit
        g_b_shaded = g_b_free_shaded + g_b_forced_shaded
        g_h_sunlit = 1 / (1 / g_a + 1 / (0.135 / 0.147 * g_b_forced_sunlit + 0.05 / 0.055 * g_b_free_sunlit))
        g_h_shaded = 1 / (1 / g_a + 1 / (0.135 / 0.147 * g_b_forced_shaded + 0.05 / 0.055 * g_b_free_shaded))
        g_r_sunlit = 4 * emissivity_c * SteBol_const * T_a_K ** 3 * k_d_black / c_p * (extinction_coeffcient_b_black_plus_d_black + (np.exp(-k_d_black * LAI) - np.exp(-k_b_black * LAI)) / (k_b_black - k_d_black))
        g_r_shaded = 4 * emissivity_c * SteBol_const * T_a_K ** 3 * k_d_black / c_p * (2 * extinction_coeffcient_d_black - extinction_coeffcient_b_black_plus_d_black - (np.exp(-k_d_black * LAI) + np.exp(-k_b_black * LAI)) / (k_b_black - k_d_black))
        r_a_c = 2 / (g_b_sunlit + g_b_shaded) * pressure * 100 / R / T_a_K # (s/m)
        r_s_c = 2 / (g_s_sunlit + g_s_shaded) * pressure * 100 / R / T_a_K # (s/m)
    else:
        g_b_free_sunlit = g_b_free_sunlit_pre
        g_b_free_shaded = g_b_free_shaded_pre
        g_b_free_soil = g_b_free_soil_pre
        G_0_sunlit = G_0_sunlit_pre
        G_0_shaded = G_0_shaded_pre
        g_b_forced_sunlit = g_b_forced_sunlit_pre
        g_b_forced_shaded = g_b_forced_shaded_pre
        g_b_sunlit = g_b_sunlit_pre
        g_b_shaded = g_b_shaded_pre
        g_h_sunlit = g_h_sunlit_pre
        g_h_shaded = g_h_shaded_pre
        g_r_sunlit = g_r_sunlit_pre
        g_r_shaded = g_r_shaded_pre
        r_a_c = r_a_c_pre
        r_s_c = r_s_c_pre

    if math.isnan(g_b_free_soil):
        g_b_free_soil = 0

    W_content_soil_surface = (W - (W_capacity - 100 * porosity)) / 100
    if W_content_soil_surface < 0: 
        W_content_soil_surface = 0
    # wetness_soil_surface: degree of saturation of the surface layer
    wetness_soil_surface = W_content_soil_surface / porosity
    # r_s: soil resistance (s/m) ~ (Sellers et al., 1992)
    r_s = np.exp(8.206 - 4.255 * wetness_soil_surface)
    
    h_beta = 1 / 4 * (1 - math.cos(W / W_capacity * math.pi)) ** 2

    g_b_free_sunlit_pre = g_b_free_sunlit
    g_b_free_shaded_pre = g_b_free_shaded
    g_b_free_soil_pre = g_b_free_soil
    G_0_sunlit_pre = G_0_sunlit
    G_0_shaded_pre = G_0_shaded
    g_b_forced_sunlit_pre = g_b_forced_sunlit
    g_b_forced_shaded_pre = g_b_forced_shaded
    g_b_sunlit_pre = g_b_sunlit
    g_b_shaded_pre = g_b_shaded
    g_h_sunlit_pre = g_h_sunlit
    g_h_shaded_pre = g_h_shaded
    g_r_sunlit_pre = g_r_sunlit
    g_r_shaded_pre = g_r_shaded
    r_a_c_pre = r_a_c
    r_s_c_pre = r_s_c

    return r_a, r_a_soil, G_0_sunlit, G_0_shaded, g_a, g_b_free_sunlit, g_b_forced_sunlit, g_b_sunlit, g_b_forced_shaded, g_b_free_shaded, g_b_shaded, g_h_sunlit, g_h_shaded, g_b_free_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, h_beta, g_b_free_sunlit_pre, g_b_free_shaded_pre, g_b_free_soil_pre, G_0_sunlit_pre, G_0_shaded_pre, g_b_forced_sunlit_pre, g_b_forced_shaded_pre, g_b_sunlit_pre, g_b_shaded_pre, g_h_sunlit_pre, g_h_shaded_pre, g_r_sunlit_pre, g_r_shaded_pre, r_a_c_pre, r_s_c_pre

"""
def CalConductance(k_n: float, g_0_sunlit: float, g_0_shaded: float, LAI: float, k_b_black: float, k_d_black: float, SteBol_const: float, c_p: float, emissivity_c: float, T_a_K: float, T_c_K_soil: float):
    if LAI > 0:
        extinction_coeffcient_b_black_plus_d_black = ExtinctionCoeffcientFunction(LAI,k_b_black+k_d_black)
        g_r_sunlit = 4 * emissivity_c * SteBol_const * T_a_K ** 3 * k_d_black / c_p * (extinction_coeffcient_b_black_plus_d_black + (np.exp(-k_d_black * LAI) - np.exp(-k_b_black * LAI)) / (k_b_black - k_d_black))

        extinction_coeffcient_d_black = ExtinctionCoeffcientFunction(LAI, k_d_black)
        g_r_shaded = 4 * emissivity_c * SteBol_const * T_a_K ** 3 * k_d_black / c_p * (2 * extinction_coeffcient_d_black - extinction_coeffcient_b_black_plus_d_black - (np.exp(-k_d_black * LAI) + np.exp(-k_b_black * LAI)) / (k_b_black - k_d_black))

        extinction_coeffcient_b_black_plus_n = ExtinctionCoeffcientFunction(LAI, k_b_black + k_n)
        G_0_sunlit = g_0_sunlit * extinction_coeffcient_b_black_plus_n
        extinction_coefficient_n = ExtinctionCoeffcientFunction(LAI, k_n)
        G_0_shaded = g_0_shaded * (extinction_coefficient_n - extinction_coeffcient_b_black_plus_n)

        g_r_soil = 4 * emissivity_c * SteBol_const * T_c_K_soil ** 3

    else:
        g_r_soil = 4 * emissivity_c * SteBol_const * T_c_K_soil ** 3
    
    if A_c_sunlit > 0:
        g_s_sunlit = G_0_sunlit + a_1 * f_w * A_c_sunlit / (C_s_sunlit - Gamma_respiration_sunlit) / (1 + VPD_s_sunlit / VPD_0) # mol m^2/s
    elif R_s_total <= 0 or S_max <= 0:
        g_s_sunlit = G_0_sunlit * 0.1
    else:
        g_s_sunlit = G_0_sunlit

    if g_s_sunlit < 0:
        g_s_sunlit = G_0_sunlit

    if A_c_shaded > 0:
        g_s_shaded = G_0_shaded + a_1 * f_w * A_c_shaded / (C_s_shaded - Gamma_respiration_shaded) / (1 + VPD_s_shaded / VPD_0) # mol m^2/s
    elif R_s_total <= 0 or S_max <= 0:
        g_s_shaded = G_0_shaded * 0.1
    else:
        g_s_shaded = G_0_shaded

    if g_s_shaded < 0:
        g_s_shaded = G_0_shaded

    W_content_soil_surface = (W - (W_capacity - 100 * porosity)) / 100
    if W_content_soil_surface < 0: 
        W_content_soil_surface = 0
    # wetness_soil_surface: degree of saturation of the surface layer
    wetness_soil_surface = W_content_soil_surface / porosity
    # r_s: soil resistance (s/m) ~ (Sellers et al., 1992)
    r_s = np.exp(8.206 - 4.255 * wetness_soil_surface)

    if LAI > 0:
        u_soil = u_h * np.exp(-u_attenuation_coef) # within-canopy wind speed ~ (Cionco, 1972)
        r_a_soil = (math.log(z / 0.01)) ** 2 / (Karman_const ** 2 * u_soil)
        g_b_forced_soil = 0.147 * math.sqrt(u_soil) # r_b: boundary layer resistance (s/m)
        g_b_free_soil = 0.055 * (T_c_K_soil - T_a_K) ** (1 / 4) # g_b_free: boundary conductance due to free convection for H2O (mol/m2/s)
    else:
        r_a_soil = (math.log(z / 0.01)) ** 2 / (Karman_const ** 2 * u_z / 2)
        g_b_forced = 0.147 * math.sqrt(u_z / 2) # r_b: boundary layer resistance (s/m)
        g_b_forced_soil = g_b_forced
        g_b_free_soil = 0.055 * (T_c_K_soil - T_a_K) ** (1 / 4) # g_b_free: boundary conductance due to free convection for H2O (mol/m2/s)

    g_b_soil = g_b_forced_soil + g_b_free_soil
    g_w_mol = 1 / ((r_s + r_a_soil) / pressure / 100 * R * T_c_K_soil + 1 / g_b_soil)
"""
"""
def SimultaneousEquations(vegetation_type: int, leaf_type: int, z: int, emissivity_c: float, emissivity_s: float, extinction_coef: float, u_attenuation_coef: float, d_leaf: list, reflectance: list, transmissivity: list, albedo_soil: list, scattering_coefficient: list, reflectance_b: float, reflectance_d: float, k_b: float, k_d: float, k_b_black: float, k_d_black: float, k_n: float, k_u: float, absorptance: float, C_a_out: int, V_c_max_org: list, Ha_V: list, J_m_25: list, Ha_J: list, Hd_J: list, Delta_S: list, Gamma_25: list, a_1: int, g_0_sunlit: float, g_0_shaded: float, VPD_0: int, T_cold: list, porosity: float, R: float, a0: float, b0: float, c0: float, Karman_const: float, SteBol_const: float, f: float, a: float, W: float, W_wilting: float, W_capacity: float, height_c: float, LAI: float, LAI_g: float, T_a_C: float, T_a_K: float, pressure: float, R_s_total: float, VPD_a: float, rho_a: float, C_p: float, c_p: float, u_z: float, R_s_b: list, R_s_d: list, R_s_b_total: float, R_s_d_total: float, S_max: float, L_down: float, gamma: float, Delta: float, COSTHETA: float, T_a_C_mean: float, LAI_sunlit: float, LAI_shaded: float, g_s_sunlit: float, g_s_shaded: float, C_i_sunlit: float, C_i_shaded: float, diff_T_sunlit: float, diff_T_shaded: float, T_c_C_soil: float, T_c_K_soil: float, Q_sunlit: list, Q_shaded: list):
    T_c_C_sunlit = T_a_C + diff_T_sunlit
    T_c_K_sunlit = T_c_C_sunlit + 273.15
    T_c_C_shaded = T_a_C + diff_T_shaded
    T_c_K_shaded = T_c_C_shaded + 273.15
    
    Q_long_isothermal_sunlit, Q_long_isothermal_shaded = L_balance(emissivity_c , k_b_black , k_d_black , SteBol_const , LAI , R_s_total , T_a_K , L_down , LAI_sunlit , LAI_shaded)
    Q_n_isothermal_sunlit, Q_n_isothermal_shaded = CalRn(Q_sunlit , Q_shaded , Q_long_isothermal_sunlit , Q_long_isothermal_shaded)
    r_a, r_a_soil, G_0_sunlit, G_0_shaded, g_a, g_b_free_sunlit, g_b_forced_sunlit, g_b_sunlit, g_b_forced_shaded, g_b_free_shaded, g_b_shaded, g_h_sunlit, g_h_shaded, g_b_free_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, h_beta = CalConductance(vegetation_type, leaf_type, d_leaf, R, Karman_const, SteBol_const, emissivity_c, k_b_black, k_d_black, k_n, k_u, g_0_sunlit, g_0_shaded, porosity, W, W_capacity, z, c_p, u_z, u_attenuation_coef, LAI, height_c, pressure, T_a_K, LAI_sunlit, LAI_shaded, T_c_K_soil, T_c_K_sunlit, T_c_K_shaded, g_s_sunlit, g_s_shaded)
    gamma_modified_sunlit, gamma_modified_shaded, C_s_sunlit, C_s_shaded, VPD_s_sunlit, VPD_s_shaded = ResConst(C_a_out , VPD_a , gamma , Delta , g_s_sunlit , g_s_shaded , C_i_sunlit , C_i_shaded , diff_T_sunlit , diff_T_shaded , g_a , g_b_forced_sunlit , g_b_free_sunlit , g_b_sunlit , g_b_forced_shaded , g_b_free_shaded , g_b_shaded , g_h_sunlit , g_h_shaded , g_r_sunlit , g_r_shaded)
    Gamma_respiration_sunlit, Gamma_respiration_shaded, C_i_sunlit, C_i_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, f_w = SimultaneousPhotosynthesis(vegetation_type, leaf_type, reflectance, transmissivity, scattering_coefficient, reflectance_b, reflectance_d, k_b, k_d, k_b_black, absorptance, C_a_out, V_c_max_org, Ha_V, J_m_25, Ha_J, Hd_J, Delta_S, Gamma_25, a_1, VPD_0, T_cold, R, f, a, W_capacity, W_wilting, W, R_s_b, R_s_d, R_s_b_total, R_s_d_total, R_s_total, pressure, S_max, T_a_C_mean, LAI, LAI_g, LAI_sunlit, LAI_shaded, Q_sunlit, C_s_shaded, C_i_sunlit, C_i_shaded, T_c_C_sunlit, T_c_K_sunlit, T_c_C_shaded, T_c_K_shaded, G_0_sunlit, G_0_shaded, h_beta, C_s_sunlit, VPD_s_sunlit, VPD_s_shaded)
    r_a, r_a_soil, G_0_sunlit, G_0_shaded, g_a, g_b_free_sunlit, g_b_forced_sunlit, g_b_sunlit, g_b_forced_shaded, g_b_free_shaded, g_b_shaded, g_h_sunlit, g_h_shaded, g_b_free_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, h_beta = CalConductance(vegetation_type, leaf_type, d_leaf, R, Karman_const, SteBol_const, emissivity_c, k_b_black, k_d_black, k_n, k_u, g_0_sunlit, g_0_shaded, porosity, W, W_capacity, z, c_p, u_z, u_attenuation_coef, LAI, height_c, pressure, T_a_K, LAI_sunlit, LAI_shaded, T_c_K_soil, T_c_K_sunlit, T_c_K_shaded, g_s_sunlit, g_s_shaded)
    Q_n_sunlit, Q_n_shaded, R_n_soil, lE_sum, ET_c_sum, ET_soil, ET_sum, H_sum = CalET_SW(emissivity_c, emissivity_s, extinction_coef, albedo_soil, SteBol_const, a0, b0, c0, T_a_K, R_s_total, VPD_a, rho_a, C_p, c_p, gamma, Delta, S_max, L_down, COSTHETA, LAI, LAI_g, T_c_C_soil, T_c_K_sunlit, T_c_K_shaded, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, r_a, r_a_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s)
    diff_T_sunlit, diff_T_shaded = CalDiffT(VPD_a, c_p, Delta, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, g_h_sunlit, g_h_shaded, g_r_sunlit, g_r_shaded, gamma_modified_sunlit, gamma_modified_shaded)
    if abs(diff_T_sunlit) > 10:
        diff_T_sunlit = 0
        no_converged = 1
        return T_c_K_sunlit, T_c_K_shaded, r_a, r_a_soil, r_a_c, r_s, Gamma_respiration_sunlit, Gamma_respiration_shaded, C_i_sunlit, C_i_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, f_w, C_s_sunlit, C_s_shaded, Q_n_sunlit, Q_n_shaded, R_n_soil, lE_sum, ET_c_sum, ET_soil, ET_sum, H_sum, diff_T_sunlit, diff_T_shaded, no_converged, True
    if abs(diff_T_shaded) > 10:
        diff_T_shaded = 0
        no_converged = 1
        return T_c_K_sunlit, T_c_K_shaded, r_a, r_a_soil, r_a_c, r_s, Gamma_respiration_sunlit, Gamma_respiration_shaded, C_i_sunlit, C_i_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, f_w, C_s_sunlit, C_s_shaded, Q_n_sunlit, Q_n_shaded, R_n_soil, lE_sum, ET_c_sum, ET_soil, ET_sum, H_sum, diff_T_sunlit, diff_T_shaded, no_converged, True
    return T_c_K_sunlit, T_c_K_shaded, r_a, r_a_soil, r_a_c, r_s, Gamma_respiration_sunlit, Gamma_respiration_shaded, C_i_sunlit, C_i_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, f_w, C_s_sunlit, C_s_shaded, Q_n_sunlit, Q_n_shaded, R_n_soil, lE_sum, ET_c_sum, ET_soil, ET_sum, H_sum, diff_T_sunlit, diff_T_shaded, no_converged, False
"""
"""
def SimultaneousEquationsNight(vegetation_type: int, leaf_type: int, z: int, emissivity_c: float, emissivity_s: float, extinction_coef: float, u_attenuation_coef: float, d_leaf: list, reflectance: list, transmissivity: list, albedo_soil: list, scattering_coefficient: list, reflectance_b: float, reflectance_d: float, k_b: float, k_d: float, k_b_black: float, k_d_black: float, k_n: float, k_u: float, absorptance: float, C_a_out: int, V_c_max_org: list, Ha_V: list, J_m_25: list, Ha_J: list, Hd_J: list, Delta_S: list, Gamma_25: list, a_1: int, g_0_sunlit: float, g_0_shaded: float, VPD_0: int, T_cold: list, porosity: float, R: float, a0: float, b0: float, c0: float, Karman_const: float, SteBol_const: float, f: float, a: float, W: float, W_wilting: float, W_capacity: float, height_c: float, LAI: float, LAI_g: float, T_a_C: float, T_a_K: float, pressure: float, R_s_total: float, VPD_a: float, rho_a: float, C_p: float, c_p: float, u_z: float, R_s_b: list, R_s_d: list, R_s_b_total: float, R_s_d_total: float, S_max: float, L_down: float, gamma: float, Delta: float, COSTHETA: float, T_a_C_mean: float, LAI_sunlit: float, LAI_shaded: float, g_s_sunlit: float, g_s_shaded: float, C_i_sunlit: float, C_i_shaded: float, diff_T_sunlit: float, diff_T_shaded: float, T_c_C_soil: float, T_c_K_soil: float, Q_sunlit: list, Q_shaded: list):
    
    T_c_C_sunlit = T_a_C + diff_T_sunlit
    T_c_K_sunlit = T_c_C_sunlit + 273.15
    T_c_C_shaded = T_a_C + diff_T_shaded
    T_c_K_shaded = T_c_C_shaded + 273.15
    
    Q_long_isothermal_sunlit, Q_long_isothermal_shaded = L_balance(emissivity_c , k_b_black , k_d_black , SteBol_const , LAI , R_s_total , T_a_K , L_down , LAI_sunlit , LAI_shaded)
    Q_n_isothermal_sunlit = Q_long_isothermal_sunlit
    Q_n_isothermal_shaded = Q_long_isothermal_shaded
    r_a, r_a_soil, G_0_sunlit, G_0_shaded, g_a, g_b_free_sunlit, g_b_forced_sunlit, g_b_sunlit, g_b_forced_shaded, g_b_free_shaded, g_b_shaded, g_h_sunlit, g_h_shaded, g_b_free_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, h_beta = CalConductance(vegetation_type, leaf_type, d_leaf, R, Karman_const, SteBol_const, emissivity_c, k_b_black, k_d_black, k_n, k_u, g_0_sunlit, g_0_shaded, porosity, W, W_capacity, z, c_p, u_z, u_attenuation_coef, LAI, height_c, pressure, T_a_K, LAI_sunlit, LAI_shaded, T_c_K_soil, T_c_K_sunlit, T_c_K_shaded, g_s_sunlit, g_s_shaded)
    gamma_modified_sunlit, gamma_modified_shaded, C_s_sunlit, C_s_shaded, VPD_s_sunlit, VPD_s_shaded = ResConst(C_a_out , VPD_a , gamma , Delta , g_s_sunlit , g_s_shaded , C_i_sunlit , C_i_shaded , diff_T_sunlit , diff_T_shaded , g_a , g_b_forced_sunlit , g_b_free_sunlit , g_b_sunlit , g_b_forced_shaded , g_b_free_shaded , g_b_shaded , g_h_sunlit , g_h_shaded , g_r_sunlit , g_r_shaded)
    Gamma_respiration_sunlit, Gamma_respiration_shaded, C_i_sunlit, C_i_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, f_w = SimultaneousPhotosynthesis(vegetation_type, leaf_type, reflectance, transmissivity, scattering_coefficient, reflectance_b, reflectance_d, k_b, k_d, k_b_black, absorptance, C_a_out, V_c_max_org, Ha_V, J_m_25, Ha_J, Hd_J, Delta_S, Gamma_25, a_1, VPD_0, T_cold, R, f, a, W_capacity, W_wilting, W, R_s_b, R_s_d, R_s_b_total, R_s_d_total, R_s_total, pressure, S_max, T_a_C_mean, LAI, LAI_g, LAI_sunlit, LAI_shaded, Q_sunlit, C_s_shaded, C_i_sunlit, C_i_shaded, T_c_C_sunlit, T_c_K_sunlit, T_c_C_shaded, T_c_K_shaded, G_0_sunlit, G_0_shaded, h_beta, C_s_sunlit, VPD_s_sunlit, VPD_s_shaded)
    r_a, r_a_soil, G_0_sunlit, G_0_shaded, g_a, g_b_free_sunlit, g_b_forced_sunlit, g_b_sunlit, g_b_forced_shaded, g_b_free_shaded, g_b_shaded, g_h_sunlit, g_h_shaded, g_b_free_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s, h_beta = CalConductance(vegetation_type, leaf_type, d_leaf, R, Karman_const, SteBol_const, emissivity_c, k_b_black, k_d_black, k_n, k_u, g_0_sunlit, g_0_shaded, porosity, W, W_capacity, z, c_p, u_z, u_attenuation_coef, LAI, height_c, pressure, T_a_K, LAI_sunlit, LAI_shaded, T_c_K_soil, T_c_K_sunlit, T_c_K_shaded, g_s_sunlit, g_s_shaded)
    Q_n_sunlit, Q_n_shaded, R_n_soil, lE_sum, ET_c_sum, ET_soil, ET_sum, H_sum = CalET_SW(emissivity_c, emissivity_s, extinction_coef, albedo_soil, SteBol_const, a0, b0, c0, T_a_K, R_s_total, VPD_a, rho_a, C_p, c_p, gamma, Delta, S_max, L_down, COSTHETA, LAI, LAI_g, T_c_C_soil, T_c_K_sunlit, T_c_K_shaded, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, r_a, r_a_soil, g_r_sunlit, g_r_shaded, r_a_c, r_s_c, r_s)
    diff_T_sunlit, diff_T_shaded = CalDiffT(VPD_a, c_p, Delta, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, g_h_sunlit, g_h_shaded, g_r_sunlit, g_r_shaded, gamma_modified_sunlit, gamma_modified_shaded)
    
    return T_c_K_sunlit, T_c_K_shaded, r_a, r_a_soil, r_a_c, r_s, C_s_sunlit, C_s_shaded, Gamma_respiration_sunlit, Gamma_respiration_shaded, C_i_sunlit, C_i_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, f_w, Q_n_sunlit, Q_n_shaded, R_n_soil, lE_sum, ET_c_sum, ET_soil, ET_sum, H_sum, diff_T_sunlit, diff_T_shaded
"""

def SimultaneousPhotosynthesis(vegetation_type: int, leaf_type: int, reflectance: list, transmissivity: list, scattering_coefficient: list, reflectance_b: float, reflectance_d: float, k_b: float, k_d: float, k_b_black: float, absorptance: float, C_a_out: int, V_c_max_org: list, Ha_V: list, J_m_25: list, Ha_J: list, Hd_J: list, Delta_S: list, Gamma_25: list, a_1: int, VPD_0: int, T_cold: list, R: float, f: float, a: float, W_capacity: float, W_wilting: float, W: float, R_s_b: list, R_s_d: list, R_s_b_total: float, R_s_d_total: float, R_s_total: float, pressure: float, S_max: float, T_a_C_mean: float, LAI: float, LAI_g: float, LAI_sunlit: float, LAI_shaded: float, Q_sunlit: list, C_s_shaded: float, C_i_sunlit: float, C_i_shaded: float, T_c_C_sunlit: float, T_c_K_sunlit: float, T_c_C_shaded: float, T_c_K_shaded: float, G_0_sunlit: float, G_0_shaded: float, h_beta: float, C_s_sunlit: float, VPD_s_sunlit: float, VPD_s_shaded: float, J_sunlit_pre: float, J_shaded_pre: float, V_j_sunlit_pre: float, V_j_shaded_pre: float, Gamma_respiration_sunlit_pre: float, Gamma_respiration_shaded_pre: float, called_from=None):
    k_n, K_c_sunlit, K_c_shaded, K_o_sunlit, K_o_shaded, Gamma_sunlit, Gamma_shaded, O_i, V_c_max_sunlit, V_c_max_shaded, V_c_sunlit, V_c_shaded, J_sunlit, J_shaded, V_j_sunlit, V_j_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, J_sunlit_pre, J_shaded_pre, V_j_sunlit_pre, V_j_shaded_pre = Photosynthesis(vegetation_type, leaf_type, transmissivity, reflectance, scattering_coefficient, reflectance_b, reflectance_d, k_b, k_d, k_b_black, absorptance, V_c_max_org, Ha_V, J_m_25, Ha_J, Hd_J, Delta_S, Gamma_25, T_cold, R, f, a, R_s_b, R_s_d, R_s_total, pressure, S_max, T_a_C_mean, LAI, LAI_g, LAI_sunlit, LAI_shaded, Q_sunlit, C_i_sunlit, C_i_shaded, T_c_C_sunlit, T_c_K_sunlit, T_c_C_shaded, T_c_K_shaded, J_sunlit_pre, J_shaded_pre, V_j_sunlit_pre, V_j_shaded_pre, called_from)
    W_retention = (W - W_wilting) / (W_capacity - W_wilting)

    if vegetation_type == 0:
        f_w = 10 * (W - W_wilting) / 3 / (W_capacity - W_wilting)
    elif vegetation_type == 1:
        f_w = (W - W_wilting) / (W_capacity - W_wilting)

    if f_w < 0:
        f_w = 0

    f_w = h_beta

    if R_s_total > 0 and S_max > 0: # analytical solution of photosynthesis rate and incellular CO2 concenteration
        if R_s_b_total > 0 and A_c_sunlit >= 0:
            Gamma_respiration_sunlit = (R_d_sunlit * (K_c_sunlit * pressure * 1000) * (1 + O_i / K_o_sunlit) + V_c_max_sunlit * Gamma_sunlit * (pressure / 1000)) / (V_c_max_sunlit - R_d_sunlit) / (pressure / 1000)
            d_1 = a_1 * f_w / (C_s_sunlit - Gamma_respiration_sunlit) / (1 + VPD_s_sunlit / VPD_0)

            if V_c_sunlit < V_j_sunlit:
                d_2 = V_c_max_sunlit
                d_3 = K_c_sunlit * (1 + O_i / K_o_sunlit)
            else:
                d_2 = 0.25 * J_sunlit
                d_3 = 2 * Gamma_sunlit

            b_0 = -(1 - d_1 * C_s_sunlit) * (d_2 * Gamma_sunlit + d_3 * R_d_sunlit) - G_0_sunlit * d_3 * C_s_sunlit
            b_1 = (1 - d_1 * C_s_sunlit) * (d_2 - R_d_sunlit) + G_0_sunlit * (d_3 - C_s_sunlit) - d_1 * (d_2 * Gamma_sunlit + d_3 * R_d_sunlit)
            b_2 = G_0_sunlit + d_1 * (d_2 - R_d_sunlit)

            if b_1 ** 2 > 4 * b_0 * b_2 and b_2 != 0:
                C_i_sunlit = (-b_1 + math.sqrt(b_1 ** 2 - 4 * b_0 * b_2)) / (2 * b_2)
            else:
                C_i_sunlit = C_a_out * 0.7

            if math.isnan(C_i_sunlit) or C_i_sunlit + d_3 == 0:
                C_i_sunlit = C_a_out * 0.7
            else:
                A_c_sunlit = d_2 * (C_i_sunlit - Gamma_sunlit) / (C_i_sunlit + d_3) - R_d_sunlit
        else:
            Gamma_respiration_sunlit = Gamma_respiration_sunlit_pre

        if R_s_d_total > 0 and A_c_shaded >= 0:
            Gamma_respiration_shaded = (R_d_shaded * (K_c_shaded * pressure * 1000) * (1 + O_i / K_o_shaded) + V_c_max_shaded * Gamma_shaded * (pressure / 1000)) / (V_c_max_shaded - R_d_shaded) / (pressure / 1000)
            d_1 = a_1 * f_w / (C_s_shaded - Gamma_respiration_shaded) / (1 + VPD_s_shaded / VPD_0)

            if V_c_shaded < V_j_shaded:
                d_2 = V_c_max_shaded
                d_3 = K_c_shaded * (1 + O_i / K_o_shaded)
            else:
                d_2 = 0.25 * J_shaded
                d_3 = 2 * Gamma_shaded

            b_0 = -(1 - d_1 * C_s_shaded) * (d_2 * Gamma_shaded + d_3 * R_d_shaded) - G_0_shaded * d_3 * C_s_shaded
            b_1 = (1 - d_1 * C_s_shaded) * (d_2 - R_d_shaded) + G_0_shaded * (d_3 - C_s_shaded) - d_1 * (d_2 * Gamma_shaded + d_3 * R_d_shaded)
            b_2 = G_0_shaded + d_1 * (d_2 - R_d_shaded)

            if b_1 ** 2 > 4 * b_0 * b_2 and b_2 != 0:
                C_i_shaded = (-b_1 + math.sqrt(b_1 ** 2 - 4 * b_0 * b_2)) / (2 * b_2)
            else:
                C_i_shaded = C_a_out * 0.7

            if math.isnan(C_i_shaded) or C_i_shaded + d_3 == 0:
                C_i_shaded = C_a_out * 0.7
            else:
                A_c_shaded = d_2 * (C_i_shaded - Gamma_shaded) / (C_i_shaded + d_3) - R_d_sunlit
        else:
            Gamma_respiration_shaded = Gamma_respiration_shaded_pre
    else:
        Gamma_respiration_sunlit, Gamma_respiration_shaded = Gamma_respiration_sunlit_pre, Gamma_respiration_shaded_pre
    
    Gamma_respiration_sunlit_pre = Gamma_respiration_sunlit
    Gamma_respiration_shaded_pre = Gamma_respiration_shaded
            
    return k_n, Gamma_respiration_sunlit, Gamma_respiration_shaded, C_i_sunlit, C_i_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, f_w, Gamma_respiration_sunlit_pre, Gamma_respiration_shaded_pre, W_retention

def Photosynthesis(vegetation_type: int, leaf_type: int, transmissivity: list, reflectance: list, scattering_coefficient: list, reflectance_b: float, reflectance_d: float, k_b: float, k_d: float, k_b_black: float, absorptance: float, V_c_max_org: list, Ha_V: list, J_m_25: list, Ha_J: list, Hd_J: list, Delta_S: list, Gamma_25: list, T_cold: list, R: float, f: float, a: float, R_s_b: list, R_s_d: list, R_s_total: float, pressure: float, S_max: float, T_a_C_mean: float, LAI: float, LAI_g: float, LAI_sunlit: float, LAI_shaded: float, Q_sunlit: list, C_i_sunlit: float, C_i_shaded: float, T_c_C_sunlit: float, T_c_K_sunlit: float, T_c_C_shaded: float, T_c_K_shaded: float, J_sunlit_pre: float, J_shaded_pre: float, V_j_sunlit_pre: float, V_j_shaded_pre: float, called_from: str):
    L = 0 # PAR
        
    K_c_sunlit = 260 * np.exp((T_c_C_sunlit - 25) * 59366 / (298 * R * T_c_K_sunlit)) / 1000 / pressure  # mol/mol
    K_c_shaded = 260 * np.exp((T_c_C_shaded - 25) * 59366 / (298 * R * T_c_K_shaded)) / 1000 / pressure  # mol/mol
    K_o_sunlit = 179 * np.exp((T_c_C_sunlit - 25) * 35948 / (298 * R * T_c_K_sunlit)) / pressure  # mol/mol
    K_o_shaded = 179 * np.exp((T_c_C_shaded - 25) * 35948 / (298 * R * T_c_K_shaded)) / pressure  # mol/mol
    Gamma_sunlit = Gamma_25[vegetation_type][leaf_type] * np.exp((T_c_C_sunlit - 25) * 29000 / (298 * R * T_c_K_sunlit)) / (pressure / 1000)  # CO2 compensation point without dark respiration (μmol/mol)
    Gamma_shaded = Gamma_25[vegetation_type][leaf_type] * np.exp((T_c_C_shaded - 25) * 29000 / (298 * R * T_c_K_shaded)) / (pressure / 1000)  # CO2 compensation point without dark respiration (μmol/mol)
    O_i = 209 / pressure  # mol/mol

    k_n, extinction_coeffcient_b_black_plus_n, extinction_coefficient_n, V_c_max_sunlit, V_c_max_shaded, V_c_sunlit, V_c_shaded = CalVc(vegetation_type, leaf_type,R, V_c_max_org, Ha_V, k_b_black, T_c_C_sunlit, T_c_K_sunlit, T_c_C_shaded, T_c_K_shaded, LAI, LAI_g, LAI_sunlit, LAI_shaded, C_i_sunlit, C_i_shaded, K_c_sunlit, K_c_shaded, K_o_sunlit, K_o_shaded, Gamma_sunlit, Gamma_shaded, O_i, called_from)

    if R_s_total > 0 and S_max > 0:
        J_sunlit, J_shaded, V_j_sunlit, V_j_shaded = CalVj(vegetation_type, leaf_type, transmissivity, reflectance, scattering_coefficient, reflectance_b, reflectance_d, k_b, k_d, k_b_black, absorptance, J_m_25, Ha_J, Hd_J, Delta_S, T_cold, R, f, a, R_s_b, R_s_d, T_a_C_mean, T_c_C_sunlit, T_c_K_sunlit, T_c_C_shaded, T_c_K_shaded, LAI, LAI_g, Q_sunlit, C_i_sunlit, C_i_shaded, Gamma_sunlit, Gamma_shaded, L, extinction_coeffcient_b_black_plus_n, extinction_coefficient_n, called_from)

        if V_c_sunlit < V_j_sunlit:
            V_n_sunlit = V_c_sunlit
        else:
            V_n_sunlit = V_j_sunlit

        if V_c_shaded < V_j_shaded:
            V_n_shaded = V_c_shaded
        else:
            V_n_shaded = V_j_shaded
    else:    
        V_n_sunlit = 0
        V_n_shaded = 0
        J_sunlit = J_sunlit_pre
        J_shaded = J_shaded_pre
        V_j_sunlit = V_j_sunlit_pre
        V_j_shaded = V_j_shaded_pre

    R_d_sunlit = 0.015 * V_c_max_sunlit
    R_d_shaded = 0.015 * V_c_max_shaded
    A_c_sunlit = V_n_sunlit - R_d_sunlit
    A_c_shaded = V_n_shaded - R_d_shaded

    J_sunlit_pre = J_sunlit
    J_shaded_pre = J_shaded
    V_j_sunlit_pre = V_j_sunlit
    V_j_shaded_pre = V_j_shaded

    return k_n, K_c_sunlit, K_c_shaded, K_o_sunlit, K_o_shaded, Gamma_sunlit, Gamma_shaded, O_i, V_c_max_sunlit, V_c_max_shaded, V_c_sunlit, V_c_shaded, J_sunlit, J_shaded, V_j_sunlit, V_j_shaded, V_n_sunlit, V_n_shaded, R_d_sunlit, R_d_shaded, A_c_sunlit, A_c_shaded, J_sunlit_pre, J_shaded_pre, V_j_sunlit_pre, V_j_shaded_pre

def CalVc(vegetation_type: int, leaf_type: int, R: float, V_c_max_org: list, Ha_V: list, k_b_black: float, T_c_C_sunlit: float, T_c_K_sunlit: float, T_c_C_shaded: float, T_c_K_shaded: float, LAI: float, LAI_g: float, LAI_sunlit: float, LAI_shaded: float, C_i_sunlit: float, C_i_shaded: float, K_c_sunlit: float, K_c_shaded: float, K_o_sunlit: float, K_o_shaded: float, Gamma_sunlit: float, Gamma_shaded: float, O_i: float, called_from: str):
    V_c_max_top_sunlit = V_c_max_org[vegetation_type][leaf_type]*np.exp((T_c_C_sunlit-25)*Ha_V[vegetation_type][leaf_type]/(298*R*T_c_K_sunlit))	# μmol/m^2/s
    V_c_max_top_shaded = V_c_max_org[vegetation_type][leaf_type]*np.exp((T_c_C_shaded-25)*Ha_V[vegetation_type][leaf_type]/(298*R*T_c_K_shaded))	# μmol/m^2/s
    k_n = 1.5*np.exp(0.00963*(V_c_max_top_sunlit*LAI_sunlit+V_c_max_top_shaded*LAI_shaded)/LAI-2.43)

    extinction_coeffcient_b_black_plus_n = ExtinctionCoeffcientFunction(LAI,k_b_black+k_n)
    V_c_max_sunlit = V_c_max_top_sunlit*extinction_coeffcient_b_black_plus_n

    extinction_coefficient_n = ExtinctionCoeffcientFunction(LAI,k_n)
    V_c_max_shaded = V_c_max_top_shaded*(extinction_coefficient_n - extinction_coeffcient_b_black_plus_n)
    
    if called_from == "NoConverged":
        V_c_sunlit = V_c_max_sunlit*(C_i_sunlit-Gamma_sunlit)/(C_i_sunlit+K_c_sunlit*(1+O_i/K_o_sunlit))
        V_c_shaded = V_c_max_shaded*(C_i_shaded-Gamma_shaded)/(C_i_shaded+K_c_shaded*(1+O_i/K_o_shaded))
    else:
        V_c_sunlit = V_c_max_sunlit*(C_i_sunlit-Gamma_sunlit)/(C_i_sunlit+K_c_sunlit*(1+O_i/K_o_sunlit))*LAI_g/LAI
        V_c_shaded = V_c_max_shaded*(C_i_shaded-Gamma_shaded)/(C_i_shaded+K_c_shaded*(1+O_i/K_o_shaded))*LAI_g/LAI
    
    return k_n, extinction_coeffcient_b_black_plus_n, extinction_coefficient_n, V_c_max_sunlit, V_c_max_shaded, V_c_sunlit, V_c_shaded

def CalVj(vegetation_type: int, leaf_type: int, transmissivity: list, reflectance: list, scattering_coefficient: list, reflectance_b: float, reflectance_d: float, k_b: float, k_d: float, k_b_black: float, absorptance: float, J_m_25: list, Ha_J: list, Hd_J: list, Delta_S: list, T_cold: list, R: float, f: float, a: float, R_s_b: list, R_s_d: list, T_a_C_mean: float, T_c_C_sunlit: float, T_c_K_sunlit: float, T_c_C_shaded: float, T_c_K_shaded: float, LAI: float, LAI_g: int, Q_sunlit: list, C_i_sunlit: float, C_i_shaded: float, Gamma_sunlit: float, Gamma_shaded: float, L: int, extinction_coeffcient_b_black_plus_n: float, extinction_coefficient_n: float, called_from: str):

    q_sunlit = Q_sunlit[0]/(1-transmissivity[L])/(1-transmissivity[L]-reflectance[L])*4.6 # 4.6 μmol/J

    J_max_sunlit = J_m_25[vegetation_type][leaf_type]*np.exp((T_c_C_sunlit-25)*Ha_J[vegetation_type][leaf_type]/(R*T_c_K_sunlit*298))*(1+np.exp((Delta_S[leaf_type]*298-Hd_J[vegetation_type][leaf_type])/(R*298)))/(1+np.exp((Delta_S[leaf_type]*T_c_K_sunlit-Hd_J[vegetation_type][leaf_type])/(R*T_c_K_sunlit)))
    J_max_shaded = J_m_25[vegetation_type][leaf_type]*np.exp((T_c_C_shaded-25)*Ha_J[vegetation_type][leaf_type]/(R*T_c_K_shaded*298))*(1+np.exp((Delta_S[leaf_type]*298-Hd_J[vegetation_type][leaf_type])/(R*298)))/(1+np.exp((Delta_S[leaf_type]*T_c_K_shaded-Hd_J[vegetation_type][leaf_type])/(R*T_c_K_shaded)))

    J_max_sunlit *= extinction_coeffcient_b_black_plus_n
    J_max_shaded *= extinction_coefficient_n - extinction_coeffcient_b_black_plus_n

    if called_from == None and leaf_type > 0 and T_a_C_mean < T_cold[vegetation_type][leaf_type]:
        J_max_sunlit /= 10
        J_max_shaded /= 10

    APAR_su = q_sunlit * absorptance * (1 - f) / 2
    b_su = -(APAR_su + J_max_sunlit)
    c_su = APAR_su * J_max_sunlit
    J_sunlit = (-b_su - math.sqrt(b_su ** 2 - 4 * a * c_su)) / (2 * a)
    J_sunlit = max(0, J_sunlit)
	
    if called_from == 'NoConverged':
        j_max_top = J_max_shaded
        q_shaded_top = Q_sunlit[0]/(1-transmissivity[L])/(1-transmissivity[L]-reflectance[L])*4.6 # 4.6 μmol/J
        APAR_sh = q_shaded_top * absorptance * (1 - f) / 2
        b_sh = -(APAR_sh + j_max_top)
        c_sh = APAR_sh * j_max_top
        J_shaded = (-b_sh - math.sqrt(b_sh ** 2 - 4 * a * c_sh)) / (2 * a)
        J_shaded = max(J_shaded, 0)

        V_j_sunlit = (C_i_sunlit - Gamma_sunlit) * J_sunlit / (4 * (C_i_sunlit + 2 * Gamma_sunlit))
        V_j_shaded = (C_i_shaded - Gamma_shaded) * J_shaded / (4 * (C_i_shaded + 2 * Gamma_shaded))

    else:
        extinction_coeffcient_d = ExtinctionCoeffcientFunction(LAI, k_d)
        extinction_coeffcient_b_black_plus_d = ExtinctionCoeffcientFunction(LAI, k_b_black + k_d)
        j_max_top = J_max_shaded / (extinction_coeffcient_d - extinction_coeffcient_b_black_plus_d)

        extinction_coeffcient_b = ExtinctionCoeffcientFunction(LAI, k_b)
        extinction_coeffcient_b_plus_b_black = ExtinctionCoeffcientFunction(LAI, k_b + k_b_black)
        extinction_coeffcient_d_plus_b_black = ExtinctionCoeffcientFunction(LAI, k_d + k_b_black)
        extinction_coeffcient_b_black = ExtinctionCoeffcientFunction(LAI, k_b_black)
        extinction_coeffcient_2b_black = ExtinctionCoeffcientFunction(LAI, 2 * k_b_black)
        c_2 = (R_s_b[L] * ((1 - reflectance_b) * k_b * (extinction_coeffcient_b - extinction_coeffcient_b_plus_b_black)) - (1 - scattering_coefficient[L]) * k_b_black * (extinction_coeffcient_b_black - extinction_coeffcient_2b_black)) / (R_s_d[L] * (1 - reflectance_d) * k_d * (extinction_coeffcient_d - extinction_coeffcient_d_plus_b_black))
        q_shaded_top = (1 + c_2) * R_s_d[L] * k_d * (1 - reflectance_d) * 4.2  # 4.2 μmol/J

        APAR_sh = q_shaded_top * absorptance * (1 - f) / 2
        b_sh = -(APAR_sh + j_max_top)
        c_sh = APAR_sh * j_max_top
        J_shaded = (extinction_coeffcient_d - extinction_coeffcient_d_plus_b_black) * (-b_sh - math.sqrt(b_sh ** 2 - 4 * a * c_sh)) / (2 * a)
        J_shaded = max(J_shaded, 0)

        V_j_sunlit = (C_i_sunlit - Gamma_sunlit) * J_sunlit / (4 * (C_i_sunlit + 2 * Gamma_sunlit)) * LAI_g / LAI
        if C_i_shaded + 2 * Gamma_shaded == 0:
            V_j_shaded = 0
        else:
            V_j_shaded = (C_i_shaded - Gamma_shaded) * J_shaded / (4 * (C_i_shaded + 2 * Gamma_shaded)) * LAI_g / LAI

    return J_sunlit, J_shaded, V_j_sunlit, V_j_shaded

def CalET_SW(emissivity_c: float, emissivity_s: float, extinction_coef: float, albedo_soil: list, SteBol_const: float, T_a_K: float, R_s_total: float, VPD_a: float, rho_a: float, C_p: float, c_p: float, gamma: float, Delta: float, S_max: float, L_down: float, COSTHETA: float, LAI: float, LAI_g: float, T_c_K_sunlit: float, T_c_K_shaded: float, Q_n_isothermal_sunlit: float, Q_n_isothermal_shaded: float, r_a: float, r_a_soil: float, g_r_sunlit: float, g_r_shaded: float, r_a_c: float, r_s_c: float, r_s: float, Q_n_sunlit_pre: float, Q_n_shaded_pre: float): # calculate evapotranspiration rate using Shuttleworth-Wallace approach
    lambda_mol_sunlit = 56780.3 - 42.84 * T_a_K # latent heat of vaporization (J/kg)
    lambda_mol_shaded = 56780.3 - 42.84 * T_a_K # latent heat of vaporization (J/kg)
    lambda_kg = (56780.3 - 42.84 * T_a_K) / 0.018 # latent heat of vaporization (J/kg)

    #Delta_soil = a0 * b0 * c0 / (T_c_C_soil + b0) ** 2 * np.exp(a0 * T_c_C_soil / (T_c_C_soil + b0)) # slope of the saturated vapor pressure (hPa/°C)
    Q_n_sunlit, Q_n_shaded, R_n_soil, G, R_n_sum, Q_n_sunlit_pre, Q_n_shaded_pre = CalRnSoil(emissivity_c, emissivity_s, extinction_coef, albedo_soil, SteBol_const, T_a_K, R_s_total, c_p, S_max, L_down, COSTHETA, LAI, T_c_K_sunlit, T_c_K_shaded, Q_n_isothermal_sunlit, Q_n_isothermal_shaded, g_r_sunlit, g_r_shaded, Q_n_sunlit_pre, Q_n_shaded_pre)
    rr_a = (Delta + gamma) * r_a
    rr_s = (Delta + gamma) * r_a_soil + gamma * r_s

    if LAI > 0:
        PM_c = (Delta * (R_n_sum - G) + (rho_a * C_p * VPD_a - Delta * r_a_c * (R_n_soil - G)) / (r_a + r_a_c)) / (Delta + gamma * (1 + r_s_c / (r_a + r_a_c))) * LAI_g / LAI
        rr_c = (Delta + gamma) * r_a_c + gamma * r_s_c
        Cc = 1 / (1 + rr_c * rr_a / rr_s / (rr_c + rr_a))
        Cs = 1 / (1 + rr_s * rr_a / rr_c / (rr_s + rr_a))
    else:
        PM_c = 0
        Cc = 0
        Cs = 1
    PM_s = (Delta * (R_n_sum - G) + (rho_a * C_p * VPD_a - Delta * r_a_soil * (R_n_sum - R_n_soil)) / (r_a + r_a_soil)) / (Delta + gamma * (1 + r_s / (r_a + r_a_soil)))

    if S_max <= 0 or R_s_total <= 0:
        PM_c = 0
        PM_s = 0
    lE_c = max(0, Cc * PM_c)
    lE_soil = max(0, Cs * PM_s)

    lE_sum = lE_c + lE_soil

    if lE_sum > abs(R_n_sum * 2):
        lE_sum = R_n_sum
        lE_c = lE_sum
        lE_soil = 0
    ET_c_sum = lE_c / lambda_kg
    ET_soil = lE_soil / lambda_kg
    ET_sum = lE_sum / lambda_kg

    H_sum = R_n_sum - lE_sum - G
    
    return Q_n_sunlit, Q_n_shaded, Q_n_sunlit_pre, Q_n_shaded_pre, R_n_soil, G, R_n_sum, lE_soil, lE_sum, ET_c_sum, ET_soil, ET_sum, H_sum

def CalRnSoil(emissivity_c: float, emissivity_s: float, extinction_coef: float, albedo_soil: list, SteBol_const: float, T_a_K: float, R_s_total: float, c_p: float, S_max: float, L_down: float, COSTHETA: float, LAI: float, T_c_K_sunlit: float, T_c_K_shaded: float, Q_n_isothermal_sunlit: float, Q_n_isothermal_shaded: float, g_r_sunlit: float, g_r_shaded: float, Q_n_sunlit_pre: float, Q_n_shaded_pre: float):
    if LAI > 0:
        Q_n_sunlit = Q_n_isothermal_sunlit - g_r_sunlit * c_p * (T_c_K_sunlit - T_a_K)
        Q_n_shaded = Q_n_isothermal_shaded - g_r_shaded * c_p * (T_c_K_shaded - T_a_K)
        
        if S_max > 0 and R_s_total > 0:
            k_r = math.sqrt(0.5) * 0.5 / COSTHETA  # k_r: canopy extinction coefficient of net radiation
            
            if k_r * LAI > 1:
                k_r_exp = np.exp(-k_r * LAI)
                k_r_exp = max(0.01, min(k_r_exp, 0.99))                    
                R_n_soil = (Q_n_sunlit + Q_n_shaded) / (1 / k_r_exp - 1)
                R_n_sum = Q_n_sunlit + Q_n_shaded + R_n_soil
                
                if abs(R_n_sum) > abs(R_s_total * (1 - albedo_soil[1]) + L_down - emissivity_c * SteBol_const * T_a_K ** 4) * 1.5:
                    R_n_sum = R_s_total * (1 - albedo_soil[1]) + L_down - emissivity_c * SteBol_const * T_a_K ** 4
                    R_n_soil = R_n_sum - (Q_n_sunlit + Q_n_shaded)
            elif LAI < 1:
                R_n_sum = (Q_n_sunlit + Q_n_shaded) / LAI
                R_n_soil = R_n_sum - (Q_n_sunlit + Q_n_shaded)
            else:
                R_n_sum = R_s_total * (1 - albedo_soil[1]) + L_down - emissivity_c * SteBol_const * T_a_K**4
                R_n_soil = R_n_sum - (Q_n_sunlit + Q_n_shaded)
            G = R_n_sum * 0.036
        else:
            R_n_soil = (L_down - emissivity_c * SteBol_const * T_a_K**4) * np.exp(-extinction_coef * LAI)
            G = R_n_soil * 0.036
            R_n_sum = L_down - emissivity_c * SteBol_const * T_a_K**4
    else:
        Q_n_sunlit = Q_n_sunlit_pre
        Q_n_shaded = Q_n_shaded_pre
        R_n_soil = R_s_total * (1 - albedo_soil[1]) + L_down - emissivity_s * SteBol_const * T_a_K ** 4
        G = R_n_soil * 0.036
        R_n_sum = R_n_soil

    Q_n_sunlit_pre = Q_n_sunlit
    Q_n_shaded_pre = Q_n_shaded

    return Q_n_sunlit, Q_n_shaded, R_n_soil, G, R_n_sum, Q_n_sunlit_pre, Q_n_shaded_pre
    
"""
def InitialTemp():
    T_c_K[N_layer + 1] = T_a_K[N_layer + 1]
    if solar_elevation > 0:
        if LAI >= 2:
            T_a_K[0] = T_a_K[N_layer + 1] - 2
        elif LAI >= 1:
            T_a_K[0] = T_a_K[N_layer + 1] - 1
        else:
            T_a_K[0] = T_a_K[N_layer + 1]
    else:
        if LAI >= 2:
            T_a_K[0] = T_a_K[N_layer + 1]
        elif LAI >= 1:
            T_a_K[0] = T_a_K[N_layer + 1]
        else:
            T_a_K[0] = T_a_K[N_layer + 1]
    
    T_c_K[0] =	T_a_K[0]
    L_up[0] = emissivity_s * SteBol_const * T_c_K[0] ** 4

    for i_l in range(N_layer, 0, -1):
        T_c_K[i_l] = T_a_K[N_layer + 1] + (T_c_K[0] - T_c_K[N_layer + 1]) / N_layer * (N_layer - i_l) + 1
        if solar_elevation > 0:
            if LAI >= 2:
                T_a_K[i_l] = T_a_K[N_layer + 1]
            elif LAI >= 1:
                T_a_K[i_l] = T_a_K[N_layer + 1]
            else:
                T_a_K[i_l] = T_a_K[N_layer + 1]
        else:
            if LAI >= 2:
                T_a_K[i_l] = T_a_K[N_layer + 1]
            elif LAI >= 1:
                T_a_K[i_l] = T_a_K[N_layer + 1]
            else:
                T_a_K[i_l] = T_a_K[N_layer + 1]

    for i_l in range(N_layer, 0, -1):
        T_c_C[i_l] = T_c_K[i_l] - 273.15
        T_a_C[i_l] = T_a_K[i_l] - 273.15
        L_down[i_l - 1] = L_down[i_l] * I_d[i_l] + emissivity_c * SteBol_const * T_c_K[i_l] ** 4 * (1 - I_d[i_l])
        L_up[N_layer - (i_l - 1)] = L_up[N_layer - 1 - (i_l - 1)] * I_d[N_layer - (i_l - 1)] + emissivity_c * SteBol_const * T_c_K[N_layer - (i_l - 1)] ** 4 * (1 - I_d[N_layer - (i_l - 1)])

    T_c_C[0] = T_c_K[0] - 273.15
    T_a_C[0] = T_a_K[0] - 273.15
"""

def WaterBalance(t_step: int, W: float, W_pre: float, W_capacity: float, W_wilting: float, rainfall: float, ET_sum: float):
    W += rainfall
    if W > W_capacity:
        W = W_capacity

    W += -ET_sum * 60 * t_step
    if W < 0:
        ET_sum = W_pre
        W = 0
        
    if W > W_capacity:
        W = W_capacity

    if W == W_capacity:
        W_retention = 1
    elif W >= W_wilting and W < W_capacity:
        W_retention = (W - W_wilting) / (W_capacity - W_wilting)
    else:
        W_retention = 0

    if W_retention < 0:
        beta_water = 0
    else:
        beta_water = W_retention

    G_water = 1 - (1 - beta_water) ** 2
    if G_water < 0.001:
        G_water = 0.001

    SW_retention = 1

    W_pre = W
    
    return W, W_retention, G_water, SW_retention, W_pre

"""
def NoConverged(vegetation_type: int, leaf_type: int, emissivity_c: float, reflectance: list, transmissivity: list, albedo_soil: list, k_b_black: float, k_d_black: float, k_n: float, absorptance: float, C_a_out: int, V_c_max_org: list, Ha_V: list, J_m_25: list, Ha_J: list, Hd_J: list, Delta_S: list, Gamma_25: list, a_1: int, g_0_sunlit: float, g_0_shaded: float, VPD_0: int, SteBol_const: float, R: float, f: float, a: float, W: float, W_wilting: float, W_capacity: float, LAI: float, LAI_g: float, T_a_C: float, T_a_K: float, pressure: float, R_s_total: float, VPD_a: float, rho_a: float, C_p: float, L_down: float, gamma: float, Delta: float, COSTHETA: float, LAI_sunlit: float, LAI_shaded: float, L: int, g_s_sunlit: float, g_s_shaded: float, C_i_sunlit: float, C_i_shaded: float, T_c_C_soil: float, Q_sunlit: list, Q_shaded: list, VPD_s_sunlit: float, VPD_s_shaded: float, f_w: float, r_a: float, r_a_soil: float, r_a_c: float, r_s: float, Gamma_respiration_sunlit: float, Gamma_respiration_shaded: float, Q_n_sunlit: float, Q_n_shaded: float):
    g_s_sunlit = 1
    g_s_shaded = 0.5
    C_s_sunlit = C_a_out * 0.8
    C_s_shaded = C_a_out * 0.8
    C_i_sunlit = C_a_out * 0.7
    C_i_shaded = C_a_out * 0.7
    diff_T_sunlit = 0
    diff_T_shaded = 0
    T_c_C_soil = 0
    T_c_K_soil = T_c_C_soil + 273.15
    Q_long_isothermal_sunlit, Q_long_isothermal_shaded = L_balance(emissivity_c , k_b_black , k_d_black , SteBol_const , LAI , R_s_total , T_a_K , L_down , LAI_sunlit , LAI_shaded)
    Q_n_isothermal_sunlit, Q_n_isothermal_shaded = CalRn(Q_sunlit , Q_shaded , Q_long_isothermal_sunlit , Q_long_isothermal_shaded)
    if (k_r*LAI > 1):
        k_r_exp = math.exp(-k_r*LAI)
        if k_r_exp < 0.01:
            k_r_exp = 0.01
        elif k_r_exp >= 1:
            k_r_exp = 0.99
        R_n_soil = (Q_n_sunlit + Q_n_shaded) / (1 / k_r_exp - 1)
        R_n_sum = Q_n_sunlit + Q_n_shaded + R_n_soil
        if abs(R_n_sum) > abs(R_s_total * (1 - albedo_soil[1]) + L_down - emissivity_c * SteBol_const * T_a_K ** 4) * 1.5:
            R_n_sum = R_s_total * (1 - albedo_soil[1]) + L_down - emissivity_c * SteBol_const * T_a_K ** 4
            R_n_soil = R_n_sum - (Q_n_sunlit + Q_n_shaded)
    elif LAI < 1:
        R_n_sum = (Q_n_sunlit + Q_n_shaded) / LAI
        R_n_soil = R_n_sum - (Q_n_sunlit + Q_n_shaded)
    else:
        R_n_sum = R_s_total * (1 - albedo_soil[1]) + L_down - emissivity_c * SteBol_const * T_a_K ** 4
        R_n_soil = R_n_sum - (Q_n_sunlit + Q_n_shaded)
    G = R_n_sum * 0.036

    W_retention = (W - W_wilting) / (W_capacity - W_wilting)
    beta_water = max(0, W_retention)
    G_water = max(0.001, 1 - (1-beta_water) ** 2)
    
    K_c_sunlit = 260 * np.exp((T_c_C_sunlit - 25) * 59366 / (298 * R * T_c_K_sunlit)) / 1000 / pressure  # mol/mol
    K_c_shaded = 260 * np.exp((T_c_C_shaded - 25) * 59366 / (298 * R * T_c_K_shaded)) / 1000 / pressure  # mol/mol
    K_o_sunlit = 179 * np.exp((T_c_C_sunlit - 25) * 35948 / (298 * R * T_c_K_sunlit)) / pressure  # mol/mol
    K_o_shaded = 179 * np.exp((T_c_C_shaded - 25) * 35948 / (298 * R * T_c_K_shaded)) / pressure  # mol/mol
    Gamma_sunlit = Gamma_25[vegetation_type][leaf_type] * np.exp((T_c_C_sunlit - 25) * 29000 / (298 * R * T_c_K_sunlit)) / (pressure / 1000)  # CO2 compensation point without dark respiration (μmol/mol)
    Gamma_shaded = Gamma_25[vegetation_type][leaf_type] * np.exp((T_c_C_shaded - 25) * 29000 / (298 * R * T_c_K_shaded)) / (pressure / 1000)  # CO2 compensation point without dark respiration (μmol/mol)
    O_i = 209 / pressure  # mol/mol
    extinction_coeffcient_b_black_plus_n, extinction_coefficient_n, V_c_max_sunlit, V_c_max_shaded, V_c_sunlit, V_c_shaded = CalVc(vegetation_type, leaf_type,R, V_c_max_org, Ha_V, k_b_black, T_c_C_sunlit, T_c_K_sunlit, T_c_C_shaded, T_c_K_shaded, LAI, LAI_g, LAI_sunlit, LAI_shaded, C_i_sunlit, C_i_shaded, K_c_sunlit, K_c_shaded, K_o_sunlit, K_o_shaded, Gamma_sunlit, Gamma_shaded, O_i, called_from = "NoConverged")
    J_sunlit, J_shaded, V_j_sunlit, V_j_shaded = CalVj(vegetation_type, leaf_type, transmissivity, reflectance, scattering_coefficient, reflectance_b, reflectance_d, k_b, k_d, k_b_black, absorptance, J_m_25, Ha_J, Hd_J, Delta_S, T_cold, R, f, a, R_s_b, R_s_d, T_a_C_mean, T_c_C_sunlit, T_c_K_sunlit, T_c_C_shaded, T_c_K_shaded, LAI, LAI_g, Q_sunlit, C_i_sunlit, C_i_shaded, Gamma_sunlit, Gamma_shaded, L, extinction_coeffcient_b_black_plus_n, extinction_coefficient_n, called_from = "NoConverged")
    V_n_sunlit = min(V_c_sunlit, V_j_sunlit)
    V_n_shaded = min(V_c_shaded, V_j_shaded)
    R_d_sunlit = 0.015 * V_c_max_sunlit
    R_d_shaded = 0.015 * V_c_max_shaded
    A_c_sunlit = V_n_sunlit - R_d_sunlit
    A_c_shaded = V_n_shaded - R_d_shaded
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
"""

def CalLAI(vegetation_type: int, leaf_type: int, k_n: float, R_g_parameter: list, allocation_parameter: list, allocation_parameter_index: list, r_00: list, s_00: list, sensitivity_allocation: list, R_m_base_stem: list, R_m_base_root: list, death_rate_leaf: float, push_down_rate: float, loss_rate_W_max: list, loss_rate_T_max: float, b_W: list, b_T: float, T_cold: list, loss_rate_stem: list, loss_rate_root: list, loss_rate_standby_leaf: list, loss_rate_standby_stem: list, loss_rate_standby_root: list, virtual_LAI_succeed_day: int, dormancy_terminate_day: int, leaf_onset: int, leaf_normal: int, leaf_dormant: int, C_increase_dy: int, phenophase: int, C_leaf: float, Cg_leaf: float, Cd_leaf: float, C_stem: float, C_root: float, C_all: float, phase2to3_dy: int, phase3to2_dy: int, phase3to2_dy2: int, phaseto4_dy: int, W: float, W_wilting: float, W_retention: float, virtual_LAI: int, virtual_LAI_day: int, dormancy_dys: int, pass200: int, dy_month: list, DOY_max: int, growing_dys: int, loss_rate: float, y2: int, month: int, day: int, DOY: int, virtual_DOY: list, virtual_T_a_C_mean: list, virtual_R_s_sum: list, virtual_pressure_mean: list, virtual_rainfall_sum: list, virtual_rh_mean: list, virtual_u_z_mean: list, virtual_W_mean: list, virtual_A_sum_daily: list, virtual_R_a_c: list, virtual_R_m_leaf_sum_daily: list, virtual_ET_daily: list, virtual_ET_c_daily: list, virtual_ET_eq_daily: list, virtual_LAI_list: list, virtual_LAI_g: list, virtual_C_leaf: list, virtual_Cg_leaf: list, virtual_Cd_leaf: list, virtual_C_stem: list, virtual_C_root: list, virtual_C_all: list, virtual_A_n: float, LAI: float, LAI_g: float, A_sum_daily: float, A_n_sum_daily: float, R_m_leaf_sum_daily: float, ET_daily: float, ET_c_daily: float, ET_eq_daily: float, R_s_sum: float, rainfall_sum: float, sun_duration_pos: float, T_a_C_mean: float, pressure_mean: float, rh_mean: float, u_z_mean: float, W_mean: float):
    continueLoop = False
    if leaf_dormant == 1:
        dormancy_dys += 1
    if dormancy_dys == dormancy_terminate_day:
        leaf_dormant = 0

    if virtual_LAI == 1:
        A_sum_daily *= (12 / 1_000_000_000)  # (μmol/m^2/dy–>kgC/m^2/dy)
        A_n_sum_daily *= (12 / 1_000_000_000)  # (μmol/m^2/dy–>kgC/m^2/dy)
        R_m_leaf_sum_daily *= (12 / 1_000_000_000)  # (μmol/m^2/dy–>kgC/m^2/dy)
        virtual_A_n = A_sum_daily - R_m_leaf_sum_daily
        print(LAI, "\t", A_sum_daily)
        A_sum_daily /= (12 / 1_000_000_000)  # (kgC/m^2/dy–>μmol/m^2/dy)
        A_n_sum_daily /= (12 / 1_000_000_000)  # (kgC/m^2/dy–>μmol/m^2/dy)
        R_m_leaf_sum_daily /= (12 / 1_000_000_000)  # (kgC/m^2/dy–>μmol/m^2/dy)

        if virtual_A_n > 0 and T_a_C_mean >= T_cold[vegetation_type][leaf_type] and ((leaf_type == 0 and sun_duration_pos >= 11) or vegetation_type == 1) and W > W_wilting:
            virtual_LAI_day += 1
            if virtual_LAI_day == virtual_LAI_succeed_day:
                virtual_LAI = 10
            print("~~~\n", virtual_A_n, "\n", "~~~\n")
        else:
            virtual_LAI_day += 1
            print("~\n", virtual_A_n, "\n", "~\n")
            virtual_LAI = -9

    if virtual_LAI == -9:
        y2 -= 24 * virtual_LAI_day
        day -= virtual_LAI_day
        if y2 < -1:
            y2 += 24 * DOY_max
        if day < 0:
            month -= 1
            if pass200 == 0 and month == 0:
                month = 1
                day = 0
            else:
                if month == 0:
                    month = 12
                day = dy_month[month - 1] + day
        virtual_LAI_day = 0
        continueLoop = True
        result = (continueLoop, dormancy_dys, growing_dys, y2, month, day, LAI, virtual_LAI, virtual_LAI_day, A_sum_daily, A_n_sum_daily, R_m_leaf_sum_daily, virtual_A_n)
        return result
    else:
        A_sum_daily, A_n_sum_daily, A_stem_daily, A_root_daily, R_m_leaf_sum_daily, R_g_leaf_daily, leaf_onset, leaf_normal, C_increase_dy, dormancy_dys, R_a_c, L_leaf, L_leaf_d, L_all, Cg_leaf, Cd_leaf, C_leaf, C_stem, C_root, C_all, phenophase, phase2to3_dy, phase3to2_dy, phase3to2_dy2, phaseto4_dy, leaf_dormant, dormancy_dys = AllocationModel(vegetation_type, leaf_type, k_n, R_g_parameter, allocation_parameter, allocation_parameter_index, r_00, s_00, sensitivity_allocation, R_m_base_stem, R_m_base_root, death_rate_leaf, push_down_rate, loss_rate_W_max, loss_rate_T_max, b_W, b_T, T_cold, loss_rate_stem, loss_rate_root, loss_rate_standby_leaf, loss_rate_standby_stem, loss_rate_standby_root, leaf_onset, leaf_normal, leaf_dormant, C_increase_dy, dormancy_dys, phenophase, C_leaf, Cg_leaf, Cd_leaf, C_stem, C_root, C_all, phase2to3_dy, phase3to2_dy, phase3to2_dy2, phaseto4_dy, W, W_wilting, W_retention, virtual_LAI, loss_rate, DOY, LAI, LAI_g, A_sum_daily, A_n_sum_daily, R_m_leaf_sum_daily, sun_duration_pos, T_a_C_mean, pressure_mean, rh_mean, u_z_mean, W_mean)
        if leaf_type == 0 and growing_dys >= virtual_LAI_succeed_day:
            growing_dys += 1
    
    if virtual_LAI==1 or virtual_LAI==10:
        virtual_DOY[virtual_LAI_day-1] = DOY
        virtual_T_a_C_mean[virtual_LAI_day-1] = T_a_C_mean
        virtual_R_s_sum[virtual_LAI_day-1] = R_s_sum
        virtual_pressure_mean[virtual_LAI_day-1] = pressure_mean
        virtual_rainfall_sum[virtual_LAI_day-1] = rainfall_sum
        virtual_rh_mean[virtual_LAI_day-1] = rh_mean
        virtual_u_z_mean[virtual_LAI_day-1] = u_z_mean
        virtual_W_mean[virtual_LAI_day-1] = W_mean
        virtual_A_sum_daily[virtual_LAI_day-1] = A_sum_daily
        virtual_R_a_c[virtual_LAI_day-1] = R_a_c
        virtual_R_m_leaf_sum_daily[virtual_LAI_day-1] = R_m_leaf_sum_daily
        virtual_ET_daily[virtual_LAI_day-1] = ET_daily
        virtual_ET_c_daily[virtual_LAI_day-1] = ET_c_daily
        virtual_ET_eq_daily[virtual_LAI_day-1] = ET_eq_daily
        virtual_LAI_list[virtual_LAI_day-1] = LAI
        virtual_LAI_g[virtual_LAI_day-1] = LAI_g
        virtual_C_leaf[virtual_LAI_day-1] = C_leaf
        virtual_Cg_leaf[virtual_LAI_day-1] = Cg_leaf
        virtual_Cd_leaf[virtual_LAI_day-1] = Cd_leaf
        virtual_C_stem[virtual_LAI_day-1] = C_stem
        virtual_C_root[virtual_LAI_day-1] = C_root
        virtual_C_all[virtual_LAI_day-1] = C_all 

    result = (continueLoop, dormancy_dys, growing_dys, y2, month, day, LAI, virtual_LAI, virtual_LAI_day, A_sum_daily, A_n_sum_daily, A_stem_daily, A_root_daily, R_m_leaf_sum_daily, R_g_leaf_daily, virtual_A_n, leaf_onset, leaf_normal, C_increase_dy, R_a_c, L_leaf, L_leaf_d, L_all, Cg_leaf, Cd_leaf, C_leaf, C_stem, C_root, C_all, phenophase, phase2to3_dy, phase3to2_dy, phase3to2_dy2, phaseto4_dy, leaf_dormant, virtual_DOY, virtual_T_a_C_mean, virtual_R_s_sum, virtual_pressure_mean, virtual_rainfall_sum, virtual_rh_mean, virtual_u_z_mean, virtual_W_mean, virtual_A_sum_daily, virtual_R_a_c, virtual_R_m_leaf_sum_daily, virtual_ET_daily, virtual_ET_c_daily, virtual_ET_eq_daily, virtual_LAI_list, virtual_LAI_g, virtual_C_leaf, virtual_Cg_leaf, virtual_Cd_leaf, virtual_C_stem, virtual_C_root, virtual_C_all)
    return result

@jit(nopython=True, cache=True)
def AllocationModel(vegetation_type: int, leaf_type: int, k_n: float, R_g_parameter: list, allocation_parameter: list, allocation_parameter_index: list, r_00: list, s_00: list, sensitivity_allocation: list, R_m_base_stem: list, R_m_base_root: list, death_rate_leaf: float, push_down_rate: float, loss_rate_W_max: list, loss_rate_T_max: float, b_W: list, b_T: float, T_cold: list, loss_rate_stem: list, loss_rate_root: list, loss_rate_standby_leaf: list, loss_rate_standby_stem: list, loss_rate_standby_root: list, leaf_onset: int, leaf_normal: int, leaf_dormant: int, C_increase_dy: int, dormancy_dys: int, phenophase: int, C_leaf: float, Cg_leaf: float, Cd_leaf: float, C_stem: float, C_root: float, C_all: float, phase2to3_dy: int, phase3to2_dy: int, phase3to2_dy2: int, phaseto4_dy: int, W: float, W_wilting: float, W_retention: float, virtual_LAI: int, loss_rate: float, DOY: int, LAI: float, LAI_g: float, A_sum_daily: float, A_n_sum_daily: float, R_m_leaf_sum_daily: float, sun_duration_pos: float, T_a_C_mean: float, pressure_mean: float, rh_mean: float, u_z_mean: float, W_mean: float):
    A_root_daily = 0
    A_stem_daily = 0
    A_sum_daily *= 12/1_000_000_000  # (μmol/m^2/dy–>kgC/m^2/dy)
    A_n_sum_daily *= 12/1_000_000_000  # (μmol/m^2/dy–>kgC/m^2/dy)
    R_m_leaf_sum_daily *= 12/1_000_000_000  # (μmol/m^2/dy–>kgC/m^2/dy)
    Q_10 = 3.22 - 0.046*T_a_C_mean  # Q_10: Q10 value for maintenance respiration
    f_20 = Q_10**((T_a_C_mean-20)/10)  # f_20: temperature dependent function
    f_15 = Q_10**((T_a_C_mean-15)/10)  # f_15: temperature dependent function
    lf_stem = max(0.05, min(np.exp(-0.2835*C_stem), 1))
    lf_root = max(0.05, min(np.exp(-0.2835*C_root), 1))
    R_m_stem_daily = R_m_base_stem[leaf_type]*lf_stem*C_stem*f_20
    R_m_root_daily = R_m_base_root[vegetation_type][leaf_type]*lf_root*C_root*f_20
    R_m_c = R_m_leaf_sum_daily + R_m_stem_daily + R_m_root_daily  # R_m: maintenance respiration (kgC/m^2/dy)
    light_availability = np.exp(-k_n*LAI)  # for trees and crops
    water_availability = W_retention
    Cg_max = ((C_stem + C_root) / allocation_parameter[vegetation_type][leaf_type]) ** (1 / allocation_parameter_index[vegetation_type][leaf_type])
    if vegetation_type == 1 and leaf_type == 0 and ((LAI < 1.5 and phenophase == 1) or virtual_LAI == 1):
        A_root_frac = 0.8
    elif (vegetation_type != 1 or leaf_type != 0) and ((Cg_leaf < Cg_max * 0.45 and phenophase == 1) or virtual_LAI == 1):
        A_root_frac = 0
        A_stem_frac = 0
    else:
        if phenophase == 1:
            leaf_onset = 0
            leaf_normal = 1
            phenophase = 2
        elif leaf_type > 0 and phenophase == 3 and Cg_leaf < Cg_max * 0.45:
            phenophase = 2
        elif phenophase < 3 and Cg_leaf > Cg_max:
            print("Cg_leaf reached the maximum value")
            phenophase = 3
            C_increase_dy = 0

        if vegetation_type == 0:
            A_root_frac = (r_00[vegetation_type][leaf_type] + sensitivity_allocation[vegetation_type][leaf_type] * (1 - water_availability)) / (1 + sensitivity_allocation[vegetation_type][leaf_type] * (2 - light_availability - water_availability))
            A_stem_frac = (s_00[leaf_type] + sensitivity_allocation[vegetation_type][leaf_type] * (1 - light_availability)) / (1 + sensitivity_allocation[vegetation_type][leaf_type] * (2 - light_availability - water_availability))
            if phenophase >= 3:
                A_root_frac = A_root_frac / (A_root_frac + A_stem_frac)
                A_stem_frac = 1 - A_root_frac
        elif vegetation_type == 1:
            A_root_frac = (r_00[vegetation_type][leaf_type] + sensitivity_allocation[vegetation_type][leaf_type] * (1 - water_availability)) / (1 + sensitivity_allocation[vegetation_type][leaf_type] * (1 + light_availability - water_availability))
            A_stem_frac = 0
            if phenophase >= 3:
                A_root_frac = 0.33
                A_stem_frac = 0

        if leaf_type == 0 and phenophase > 0 and (A_sum_daily - R_m_leaf_sum_daily - R_m_stem_daily - R_m_root_daily < 0 or (A_sum_daily - R_m_c > 0 and A_sum_daily - R_m_leaf_sum_daily - R_m_stem_daily - R_m_root_daily - R_g_parameter[vegetation_type] * (A_sum_daily - R_m_c) < 0) or T_a_C_mean < T_cold[vegetation_type][leaf_type] or (sun_duration_pos < 11 and vegetation_type == 0)):
            C_increase_dy = 0
    A_leaf_frac = 1 - (A_root_frac + A_stem_frac)

    # R_g: growth respiration (kgC/m^2/dy)
    if A_sum_daily - R_m_c > 0:
        R_g_c = R_g_parameter[vegetation_type] * (A_sum_daily - R_m_c)
        # R_g_[X]: growth respiration for [X]
        R_g_leaf_daily = A_leaf_frac * R_g_c
        R_g_stem_daily = A_stem_frac * R_g_c
        R_g_root_daily = A_root_frac * R_g_c
    else:
        R_g_c = 0
        R_g_leaf_daily = 0
        R_g_stem_daily = 0
        R_g_root_daily = 0

    R_a_c = R_m_c + R_g_c  # R_a_c: autotrophic respiration (kgC/m^2/dy)

    # A_[X]: allocation to [X] (kgC/m^dy)
    if phenophase == 1 and (vegetation_type == 0 or (vegetation_type == 1 and leaf_type == 1)):
        A_stem_daily = 0
        A_root_daily = 0
    elif phenophase >= 2 or (vegetation_type == 1 and leaf_type == 0 and phenophase == 1):
        if A_sum_daily - R_a_c >= 0:
            A_stem_daily = (A_sum_daily - R_a_c) * A_stem_frac + R_m_stem_daily + R_g_stem_daily
            A_root_daily = (A_sum_daily - R_a_c) * A_root_frac + R_m_root_daily + R_g_root_daily
        else:
            A_stem_daily = A_sum_daily * A_stem_frac
            A_root_daily = A_sum_daily * A_root_frac

    if vegetation_type == 0 and leaf_type == 2:
        T_retention = 1
    else:
        if T_a_C_mean > T_cold[vegetation_type][leaf_type]:
            T_retention = 1
        elif T_a_C_mean > T_cold[vegetation_type][leaf_type] - 5:
            T_retention = (T_a_C_mean - (T_cold[vegetation_type][leaf_type] - 5)) / 5
        else:
            T_retention = 0
    
    # loss_rate_[X]: leaf loss rate because of [X] stress (/dy)
    if phenophase > 0:
        loss_rate_W = loss_rate_W_max[vegetation_type][leaf_type] * (1 - W_retention) ** b_W[vegetation_type][leaf_type]
        loss_rate_T = loss_rate_T_max * (1 - T_retention) ** b_T
    else:
        loss_rate_W = 0
        loss_rate_T = 0

    S_leaf = death_rate_leaf * Cg_leaf
    # L_[X]: amount of C for [X] litter (kgC/m^2/dy)
    L_leaf_d = push_down_rate * Cd_leaf

    if phenophase >= 2 and LAI > 0 and virtual_LAI != 1:
        L_leaf = (loss_rate + loss_rate_W + loss_rate_T) * Cg_leaf
        if phenophase == 4 and leaf_type == 0:
            L_leaf += ((C_stem + C_root) / allocation_parameter[vegetation_type][leaf_type]) ** (1 / allocation_parameter_index[vegetation_type][leaf_type]) * (1 - 0.075) / 30
        if vegetation_type == 0 and u_z_mean >= 6.0:
            if leaf_type == 0:
                L_leaf += 0.001 * u_z_mean ** 2 * Cg_leaf
                L_leaf_d += 0.003 * u_z_mean ** 2 * Cd_leaf
            elif leaf_type == 1 or leaf_type == 2:
                L_leaf += 0.0003 * u_z_mean ** 2 * Cg_leaf
                L_leaf_d += 0.001 * u_z_mean ** 2 * Cd_leaf
            elif leaf_type == 3:
                L_leaf += 0.0003 * u_z_mean * Cg_leaf
                L_leaf_d += 0.001 * u_z_mean * Cd_leaf
    else:
        L_leaf = 0

    if phenophase == 1:
        if vegetation_type == 0:
            L_stem = loss_rate_stem[vegetation_type][leaf_type] * C_stem
            L_root = loss_rate_root[vegetation_type][leaf_type] * C_root
        elif vegetation_type == 1:
            L_stem = 0
            L_root = 2 * (loss_rate_root[vegetation_type][leaf_type] * C_root)
    else:
        if vegetation_type == 0:
            L_stem = loss_rate_stem[vegetation_type][leaf_type] * C_stem * 2
            L_root = loss_rate_root[vegetation_type][leaf_type] * C_root * 2
        elif vegetation_type == 1:
            L_stem = 0
            L_root = loss_rate_root[vegetation_type][leaf_type] * C_root * 2

    L_all = L_leaf + L_stem + L_root

    # C_[X]: amount of C for [X] (kgC/m^2)
    if phenophase >= 3:
        if A_sum_daily - A_stem_daily - A_root_daily - R_m_leaf_sum_daily - R_g_leaf_daily - L_leaf > 0:
            L_leaf += A_sum_daily - A_stem_daily - A_root_daily - R_m_leaf_sum_daily - R_g_leaf_daily - L_leaf

    Cg_leaf += A_sum_daily - A_stem_daily - A_root_daily - R_m_leaf_sum_daily - R_g_leaf_daily - L_leaf
    Cg_leaf = max(0, Cg_leaf)
    Cd_leaf += L_leaf - L_leaf_d
    Cd_leaf = max(0, Cd_leaf)
    C_leaf = Cg_leaf + Cd_leaf
    C_leaf = max(0, C_leaf)
    C_stem += A_stem_daily - R_m_stem_daily - R_g_stem_daily - L_stem
    C_stem = max(0, C_stem)
    C_root += A_root_daily - R_m_root_daily - R_g_root_daily - L_root
    C_root = max(0, C_root)
    C_all += A_sum_daily - R_a_c - L_all
    C_all = max(0, C_all)

    if leaf_type > 0 and phenophase >= 4 and A_sum_daily - R_a_c - L_all >0:
        C_increase_dy += 1
        if C_increase_dy == 5 and W >= W_wilting:
            phenophase = 2
    else:
        C_increase_dy = 0

    if phenophase == 2 and (A_sum_daily - A_stem_daily - A_root_daily - R_m_leaf_sum_daily - R_g_leaf_daily - L_leaf < 0 or T_a_C_mean < T_cold[vegetation_type][leaf_type] or (sun_duration_pos < 11 and vegetation_type == 0 and leaf_type == 0)):
        phase2to3_dy += 1
        if phase2to3_dy == 5:
            phenophase = 3
    else:
        phase2to3_dy = 0

    if leaf_type > 0 and phenophase == 3 and A_sum_daily - A_stem_daily - A_root_daily - R_m_leaf_sum_daily - R_g_leaf_daily - L_leaf >= 0 and T_a_C_mean > T_cold[vegetation_type][leaf_type] and W >= W_wilting:
        phase3to2_dy += 1
        if phase3to2_dy == 7:
            phenophase = 2
    else:
        phase3to2_dy = 0

    if leaf_type > 0 and phenophase == 3 and DOY >= 1 and DOY <= 180 and T_a_C_mean > T_cold[vegetation_type][leaf_type] and W >= W_wilting:
        phase3to2_dy2 += 1
        if phase3to2_dy2 == 7:
            phenophase = 2
    else:
        phase3to2_dy2 = 0

    if vegetation_type == 0 and leaf_type == 0 and (phenophase == 2 or phenophase == 3) and (T_a_C_mean < T_cold[vegetation_type][leaf_type] or sun_duration_pos < 11):
        phaseto4_dy += 1
        if phaseto4_dy == 5:
            phenophase = 4
    else:
        phaseto4_dy = 0

    if leaf_type == 0 and leaf_dormant == 0 and phenophase != 0 and ((vegetation_type == 0 and LAI_g < 0.3) or (vegetation_type == 1 and LAI_g < 0.2) or Cg_leaf < Cg_max * 0.075):
        if phenophase >= 2:
            Cg_leaf = 0
            Cd_leaf = 0
            C_leaf = 0
            C_all = C_stem + C_root
            leaf_dormant = 1
            dormancy_dys = 0
            leaf_normal = 0
            phenophase = 0
            leaf_fall = 0

    return A_sum_daily, A_n_sum_daily, A_stem_daily, A_root_daily, R_m_leaf_sum_daily, R_g_leaf_daily, leaf_onset, leaf_normal, C_increase_dy, dormancy_dys, R_a_c, L_leaf, L_leaf_d, L_all, Cg_leaf, Cd_leaf, C_leaf, C_stem, C_root, C_all, phenophase, phase2to3_dy, phase3to2_dy, phase3to2_dy2, phaseto4_dy, leaf_dormant, dormancy_dys

def DailyMeanClimate(t: int, t_start: int, t_step: int, T_a_C: float, rainfall: float, pressure: float, R_s_total: float, A_n_obs: float, rh: float, u_z: float, W: float, T_a_C_sum: float, R_s_sum: float, pressure_sum: float, rainfall_sum: float, rh_sum: float, u_z_sum: float, W_sum: float, A_n_obs_sum: float, T_a_C_mean: float, pressure_mean: float, rh_mean: float, u_z_mean: float, W_mean: float):
    T_a_C_sum += T_a_C
    R_s_sum += R_s_total * 60 * t_step / 1_000_000 #Jm^-2s^-1 -> MJm^-2dy^-1
    pressure_sum += pressure
    rainfall_sum += rainfall
    rh_sum += rh
    u_z_sum += u_z
    W_sum += W
    A_n_obs_sum += A_n_obs * 12 * 60 * t_step / 1_000_000 #μmolm^-2s^-1 -> gCm^-2dy^-1

    if t == (1440 - t_start):
        T_a_C_mean = T_a_C_sum / (1440/t_step)
        pressure_mean = pressure_sum / (1440/t_step)
        rh_mean = rh_sum / (1440/t_step)
        u_z_mean = u_z_sum / (1440/t_step)
        W_mean = W_sum / (1440/t_step)

    return T_a_C_sum, R_s_sum, pressure_sum, rainfall_sum, rh_sum, u_z_sum, W_sum, A_n_obs_sum, T_a_C_mean, pressure_mean, rh_mean, u_z_mean, W_mean
