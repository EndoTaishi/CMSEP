#!/usr/bin/perl

use Math::Trig;

# $site = 'CA-Oas';
$site = 'KWG';
# if ($site eq 'CA-Oas') {
if ($site eq 'KWG') {
	# $lat = 53.6289;	$lon = -106.198;	$t_offset_GMT = -6;	# Site local time offset from UTC (GMT) (hours)
	$lat = 35.8725;	$lon = 139.4869;	$t_offset_GMT = 9;	# Site local time offset from UTC (GMT) (hours)
	# $elv = 530;	# $elv: elevation (m)
	$elv = 26;	# $elv: elevation (m)
	$z = 30;	# $z: reference height (m)
	$PFT = 'DBF';
	$vegetation_type = 0;
	$leaf_type = 0;
	$t_step = 30;	# (min)
	$year_start = 2000;	$year_end = 2000;
}

$t_start = $t_step/2;
$lon_LST = 15*$t_offset_GMT;

## vegetatiom type ######################
# 0: forest
# 1: grass
#$vegetation_type = 0;
#########################################

## vegetatiom type ######################
# 0: deciduous
# 1: conifer-evergreen
#$leaf_type = 1;
#########################################

## input parameters ####################################
$emissivity_c = 0.97;	# $emissivity_c: emissivity of canopy
$emissivity_s = 0.97;	# $emissivity_s: emissivity of soil
$extinction_coef = 0.5;	# $extinction_coef: extinction coeffcient
$u_attenuation_coef = 2.5;	# $u_attenuation_coef: in-canopy extinction coefficient ~ (Shuttleworth & Wallace, 1985)
@d_leaf = (0.05,0.05,0.03,0.03);	# $d: characteristic dimension of leaf (m) ~ (Jones, 1992)
$albedo_c = 0.15;
@reflectance = (0.09,0.3);	# PAR, SR
@transmissivity = (0.06,0.2);	# PAR, SR
@albedo_soil = (0.06,0.26);	# PAR, SR
@scattering_coefficient = (0.20,0.85);	# scattering coefficient of the leaf (PAR, SR)
$reflectance_b = 0.09;	# effective canopy-soil reflectance for direct beam radiaiton
$reflectance_d = 0.09;	# effective canopy-soil reflectance for diffusive radiaiton
$k_b = 0.5;	# $k_b: extinction_coeffcient of a canopy for direct beam radiation
$k_d = 0.4;	# $k_d: extinction_coeffcient of a canopy for direct diffuse radiation
$k_b_black = 0.5;	# $k_b_black: extinction_coeffcient of a canopy of black leaves for direct beam radiation
$k_d_black = 0.4;	# $k_d_black: extinction_coeffcient of a canopy of black leaves for direct diffuse radiation
$k_n = 0.5;	# $k_n: extinction_coeffcient of a canopy for leaf nitrogen
$k_u = 0.5;	# $k_n: extinction_coeffcient of a canopy for wind speed
$absorptance = 1-$reflectance[0]-$transmissivity[0];	# PAR
@alpha = (5,2.5,2.5,2.5);
@beta = (2,4.5,4.5,4.5);
$C_a_out = 380;	# atmospheric CO2 (ppm = μmol/mol)
$O_a_out  = 209490;	# oxygen partial pressure in chloroplast (ppm) (205–209 μbar)
@g_s_max = ([7.5,8.3,8.3,8.3],[12,10]);	# maximum stomatal conductance (mm/s) ~ (IGBP)
@R_g_parameter = (0.25,0.25);	# growth respiration parameter ~ (Knorr, 2000)
@allocation_parameter = ([40,52,30,35],[0.8,3]);
@allocation_parameter_index = ([1.4,1.3,1.4,1.6],[1.2,1.2]);
$r_0 = 0.3;	# r_0: fractional carbon allocation to root for non-limiting conditions
$s_0 = 0.3;	# s_0: fractional carbon allocation to stem for non-limiting conditions
@r_00 = ([0.55,0.60,0.55,0.6],[0.8,0.8]);	# r_0: fractional carbon allocation to root for non-limiting conditions
@s_00 = (0.1,0.05,0.2,0.05);	# s_0: fractional carbon allocation to stem for non-limiting conditions
@sensitivity_allocation = ([0.8,0.8,0.8,0.5],[1,1]);	# $sensitivity_allocation: sensitivity of allocation to changes in W & L
$specified_R_rate = 0.218;	# $specified_R_rate: specified respiration rate based on N content (kgC/kgN/dy) ~ (Keyser et al., 2000)
$CN_ratio_stem = 50;	# $CN_ratio_stem: C:N ratio of stem ~ (Arora, 2003)
$CN_ratio_root = 50;	# $CN_ratio_root: C:N ratio of root ~ (Arora, 2003)
@V_c_max_org = ([80,96,50,94],[80,80]);
@Ha_V = ([62000,71000,71000,71000],[65000,70000]);	# activation_energy (J/mol)
@J_m_25 = ([80,84,60,82],[100,60]);	# μmol/m^2/s
@Ha_J = ([48000,59000,52000,52000],[48000,52000]);	# activation_energy (J/mol)
@Hd_J = ([219400,219400,219400,219400],[219400,219400]);	# deactivation_energy (J/mol)
@Delta_S = (635.2,650,650,650);	# entropy factor (J/mol/K)
@Gamma_25 = ([37,44.7,40,44.7],[45,42]);	# CO2 compensation point without dark respiration (μbar)
$quantum_efficiency = 0.385;	# quantum efficiency of RuBP regeneration (mol e mol-1 quanta)
$a_1 = 10;
$g_0_sunlit = 0.1;
$g_0_shaded = 0.1;
$VPD_0 = 15;	# Leuning, 1995
@R_m_base_stem = (0.00005,0.0008,0.00002,0.0008);
@R_m_base_root = ([0.0012,0.0008,0.0012,0.001],[0.0002,0.0015]);
$death_rate_leaf = 0.0045;	# $death_rate_leaf: death rate of leaf (/dy)
$push_down_rate = 0.05;	# $push_down_rate: rate of leaf pushed down (/dy)
@loss_rate_W_max = ([0.025,0.005,0.005,0.015],[0.025,0.015]);	# $loss_rate_W_max: maximum drought leaf loss rate (/dy)
$loss_rate_T_max = 3;	# $loss_rate_T_max: maximum colf leaf loss rate (/dy)
@b_W = ([1.1,3,3,2],[2,3]);	# $b_W: shape parameter for leaf loss (drought)
$b_T = 3;	# $b_T: shape parameter for leaf loss (cold)
@T_cold = ([7,0,5,-5],[5,8]);	# $T_cold: Temperature threshold for leaf loss because of cold stress (°C)
@loss_rate_stem = ([0.000047,0.000061,0.000061,0.000055],[0.00035,0.00035]);
@loss_rate_root = ([0.0001,0.000061,0.0003,0.000055],[0.00035,0.00035]);
@loss_rate_standby_leaf = ([0,0.00052,0.00052,0.00055],[0,0.0001]);
@loss_rate_standby_stem = (0.00012,0.00012,0.00022,0.000055);
@loss_rate_standby_root = ([0.00012,0.00012,0.00022,0.000055],[0.01,0.0006]);
$virtual_LAI_succeed_day = 7;	# $virtual_LAI_succeed_day: successive days neccessary for virtual leaf to survive at onset
$dormancy_terminate_day = 90;
# soil water characteristics
$porosity = 0.5;	# $porosity: porosity (i.e. saturated moisture content)
########################################################

## constants ###########################################
$pi = atan2(1,1)*4;
$R = 8.31447; # $R: gas constant (J/K/mol)
$rho = 1.20;	# $rho: air density (kg/m^3)
$M_d = 0.02897;	# $M_d: molecular mass of dry air (kg/mol)
$lapse_rate = 0.0065;	# $lapse_rate: lapse-rate of air temperature (K/m)
$refraction = 0.01;	# $refraction: horizontal refraction (rad)
$SteBol_const = 5.67*10**(-8);	# $SteBol_const: Stefan-Bolzman constatnt (W/m^2/K^4)
$C_pd = 1005;	# $C_pd: specific heat of dry air at constant pressure (J/kg/K)
$gamma = 0.68;	# $gamma: psychrometric constant (hPa=mbar)
$a0 = 17.67; $b0 = 243.5; $c0 = 6.112; # empirical parameters ~ (Tetens's equation, 1930)
$Karman_const = 0.41;	# $Karman_const: the von-Karman constant
$f = 0.15;	# to correct for spectral quality of the light (Evans, 1987)
$a = 0.7;	# empirical curvature factor (0.7 is a good average value, Evans, 1989)
$f_s = 0.7;	# Direct solar fraction of solar radiation incident on the canopy
########################################################

if ($vegetation_type == 0) {
	$N_layer = 20;
} elsif ($vegetation_type == 1) {
	$N_layer = 10;
}
$N_angle = 10;
$N2_angle = 40;
$diff_level = 0.01;

$latitude = $lat/180*$pi;

if (!-d './../out/'.$site) {
	mkdir './../out/'.$site;
}
if (!-d './../out/'.$site.'/yearly') {
	mkdir './../out/'.$site.'/yearly';
}

open (YEARLY,'>./../out/'.$site.'/yearly/yearly.csv') or die;
#open (FCAP,'<./../data/'.$site.'/soil/fieldcap.out') or die;
#$FCAP_line = <FCAP>;
#chomp $FCAP_line;
#@FCAP_data=split(/\t/,$FCAP_line);
#$W_capacity = $FCAP_data[2];
#close (FCAP);
#open (WILT,'<./../data/'.$site.'/soil/wiltpont.out') or die;
#$WILT_line = <WILT>;
#chomp $WILT_line;
#@WILT_data=split(/\t/,$WILT_line);
#$W_wilting = $WILT_data[2];
$W_capacity = 400; $W_wiling = 200;
$W_critical = $W_capacity*0.75;	# $W_critical: critical point (mm)
#close (WILT);


$leaf_onset = 0;
$leaf_fall = 0;
$leaf_dormant = 0;
$leaf_normal = 0;
$C_increase_dy = 0;
$C_decrease_dy = 0;
if ($leaf_type == 0) {
	$phenophase = 0;
} elsif ($leaf_type==1 or $leaf_type==2 or $leaf_type==3) {
	$phenophase = 2;
}
if ($leaf_type == 0) {
	$C_leaf=0.0; $Cg_leaf=0.0;	$Cd_leaf=0.0;	$C_stem=8.0;	$C_root=1.8; $C_all=$C_stem+$C_root;
} elsif ($leaf_type == 1) {
	$C_leaf=0.0; $Cg_leaf=0.0;	$Cd_leaf=0.0;	$C_stem=0.5;	$C_root=0.2; $C_all=$C_stem+$C_root;
} elsif ($leaf_type == 2) {
	$C_leaf=0.15; $Cg_leaf=0.13;	$Cd_leaf=0.02; $C_stem=6.5;	$C_root=1.8; $C_all=$C_leaf+$C_stem+$C_root;
} elsif ($leaf_type == 3) {
	$C_leaf=0.18; $Cg_leaf=0.16;	$Cd_leaf=0.02; $C_stem=(0.6-0.03*($lat-40))*10; $C_all=$C_leaf+$C_stem+$C_root;
}
if ($vegetation_type == 1) {
	$C_leaf=0.0; $Cg_leaf=0.0;	$Cd_leaf=0.0;	$C_stem=0;
	if ($leaf_type == 0) {
		$C_root = 0;
	} elsif ($leaf_type == 1) {
		$C_root = 0.55;
	}
	$C_all=$C_leaf+$C_root;
}
$sun_duration_sum = 0;
$sun_duration_count = -9;
$solar_elevation_pre = -30;


if ($vegetation_type == 0) {
	$W = $W_capacity*0.9;
} elsif ($vegetation_type == 1) {
	$W = $W_capacity*0.8;
}

$virtual_LAI = 0;
$virtual_LAI_day = 0;
$dormancy_dys = 0;
$below_minus5 = 0;
$above_minus5 = 0;
$below_0 = 0;
$above_0 = 0;

$itrn = 0;	$pass200 = 0;
$diff_count_yr = 0;
$A_sum_yearly_pre = 0;
$R_m_leaf_sum_yearly_pre = 0;
$R_a_c_yearly_pre = 0;
$NPP_yearly_pre = 0;
$ET_yearly_pre = 0;
$ET_c_yearly_pre = 0;
# iterate calclulations until steady conditions
LOOP: while (0<1) {

	for ($year=$year_start; $year<=$year_end; $year++) {
	 

		if ($year%4==0) {
			@month_days = (0,31,60,91,121,152,182,213,244,274,305,335);
			@dy_month = (31,29,31,30,31,30,31,31,30,31,30,31);
		} else {
			@month_days = (0,31,59,90,120,151,181,212,243,273,304,334);
			@dy_month = (31,28,31,30,31,30,31,31,30,31,30,31);
		}
		if ($year%4==0) {
			$DOY_max = 366;
		} else {
			$DOY_max = 365;
		}

		if ($leaf_type == 0) {
			$growing_dys = 200;	# $growing_dys: leaf life span (dy)
			$loss_rate = 1/$growing_dys;	# $loss_rate: normal turnover (/dy)
		} elsif ($leaf_type==1 or $leaf_type==2 or $leaf_type==3) {
			if ($vegetation_type == 0) {
				$loss_rate = 1/$DOY_max/1.75;
			} elsif ($vegetation_type == 1) {
				$loss_rate = 1/$DOY_max;
			}
		}

		if ($virtual_LAI != 1) {

			$diff_count_yr = 0;
			$A_sum_yearly = 0;
			$R_a_c_yearly = 0;
			$R_m_leaf_sum_yearly = 0;
			$NPP_yearly = 0;
			$ET_yearly = 0;
			$ET_c_yearly = 0;
			$ET_eq_yearly = 0;

			if ($vegetation_type == 0) {
				if ($leaf_type == 0) {
					$SLA = 27;	# $SLA: specidic leaf area (m^2/kgC)
				} elsif ($leaf_type == 1) {
					$SLA = 20;
				} elsif ($leaf_type==2 or $leaf_type==3) {
					$SLA = 12;	# $SLA: specidic leaf area (m^2/kgC)
				}
			} elsif ($vegetation_type == 1) {
				if ($leaf_type == 0) {
					$SLA = 30;	# $SLA: specidic leaf area (m^2/kgC)
				} elsif ($leaf_type==1) {
					$SLA = 25.0;	# $SLA: specidic leaf area (m^2/kgC)
				}
			}

		}

		&MakeDirectory;

		$y2 = -1;
		for ($month='01'; $month<='12'; $month++) {

			if ($month==1) {$month='01';}
			if ($month==2) {$month='02';}
			if ($month==3) {$month='03';}
			if ($month==4) {$month='04';}
			if ($month==5) {$month='05';}
			if ($month==6) {$month='06';}
			if ($month==7) {$month='07';}
			if ($month==8) {$month='08';}
			if ($month==9) {$month='09';}


			for ($day='01'; $day<=$dy_month[$month-1]; $day++) {
				if ($day==1) {$day='01';}
				if ($day==2) {$day='02';}
				if ($day==3) {$day='03';}
				if ($day==4) {$day='04';}
				if ($day==5) {$day='05';}
				if ($day==6) {$day='06';}
				if ($day==7) {$day='07';}
				if ($day==8) {$day='08';}
				if ($day==9) {$day='09';}

				$DOY = $month_days[$month-1]+$day;
				if ($DOY==200) {$pass200=1;}
				# initial conditions
				if ($virtual_LAI == -9) {
					$LAI = 0;
					$virtual_LAI = 0;
					$W = $W_bare;
					$A_sum_yearly = $A_sum_yearly_bare;
					$R_m_leaf_sum_yearly = $R_m_leaf_sum_yearly_bare;
					$ET_yearly = $ET_yearly_bare;
					$ET_c_yearly = $ET_c_yearly_bare;
					$ET_eq_yearly = $ET_eq_yearly_bare;
					$C_leaf = $C_leaf_bare;
					$Cg_leaf = $Cg_leaf_bare;
					$Cd_leaf = $Cd_leaf_bare;
					$C_stem = $C_stem_bare;
					$C_root = $C_root_bare;
					$C_all = $C_all_bare;
				} elsif ($leaf_dormant==0 && $virtual_LAI!=1 && (($vegetation_type==0 && $LAI_g<0.3) or ($vegetation_type==1 && $LAI_g<0.2) or $Cg_leaf<(($C_stem+$C_root)/$allocation_parameter[$vegetation_type][$leaf_type])**(1/$allocation_parameter_index[$vegetation_type][$leaf_type])*0.075)) {
					if ($leaf_type == 0) {
						$leaf_fall = 0;
						$W_bare = $W;
						$A_sum_yearly_bare = $A_sum_yearly;
						$R_m_leaf_sum_yearly_bare = $R_m_leaf_sum_yearly;
						$ET_yearly_bare = $ET_yearly;
						$ET_c_yearly_bare = $ET_c_yearly;
						$ET_eq_yearly_bare = $ET_eq_yearly;
						$C_leaf_bare = $C_leaf;
						$Cg_leaf_bare = $Cg_leaf;
						$Cd_leaf_bare = $Cd_leaf;
						$C_stem_bare = $C_stem;
						$C_root_bare = $C_root;
						$C_all_bare = $C_all;
						$C_leaf_0 = (($C_stem+$C_root)/$allocation_parameter[$vegetation_type][$leaf_type])**(1/$allocation_parameter_index[$vegetation_type][$leaf_type])*0.075;
							if ($vegetation_type == 0) {
								if ($C_leaf_0 < 0.3/$SLA) {$C_leaf_0 = 0.3/$SLA;	$Cg_leaf = $C_leaf_0;	$Cd_leaf = 0;}
							} elsif ($vegetation_type == 1) {
								if ($C_leaf_0 < 0.2/$SLA) {$C_leaf_0 = 0.2/$SLA;	$Cg_leaf = $C_leaf_0;	$Cd_leaf = 0;	$C_root=$C_leaf_0;}
							}
						$C_leaf = $C_leaf_0;	$Cg_leaf = $C_leaf_0;	$Cd_leaf = 0;	$C_all += $C_leaf+$C_root;
						$virtual_LAI = 1;
						$virtual_LAI_day = 0;
						$virtual_A_n = 0;
						if ($leaf_type==0 && $growing_dys>$virtual_LAI_succeed_day) {
							$growing_dys_pre = $growing_dys;
							$growing_dys = 0;
						}
						print ">", $growing_dys, "\t", $Cg_leaf, "\t", $Cd_leaf, "\t", $C_leaf, "\n";
					} elsif ($leaf_type==1 or $leaf_type==2 or $leaf_type==3) {
						if ($W>=$W_wilting) {
							$leaf_normal = 0;
							$phenophase = 1;
						}
					}
				}

				$LAI = $SLA*$C_leaf;
				$LAI_g = $SLA*$Cg_leaf;
				$LAI_d = $SLA*$Cd_leaf;
				$N_layer = int $LAI*4;
				if ($N_layer <= 4) {$N_layer = 4;}

				if ($vegetation_type == 0) {
					$a_ch = 6.5;	# for forest ~ (Lawrence & Slingo, 2004)
				} elsif ($vegetation_type == 1) {
					$a_ch = 0.5;	# for grass ~ (Lawrence & Slingo, 2004)
				}
	
				# $height_c: canopy height (m)
				if ($LAI == 0) {
						if ($vegetation_type == 0) {
							$height_c = 10.0*$C_stem**0.385;
						} elsif ($vegetation_type == 1) {
							$height_c = 0.1;
						}
				} else {
						if ($vegetation_type == 0) {
							$height_c = 10.0*$C_stem**0.385;
						} elsif ($vegetation_type == 1) {
							$height_c = ($a_ch*$LAI)**(2/3);
						}
						if ($height_c > $z) {$z = $height_c+1;}
				}
	
				$d_z_h = 1/$N_layer;
				$d_z = $d_z_h*$height_c;

				$A_sum_daily = 0;
				$A_n_sum_daily = 0;
				$R_m_leaf_sum_daily = 0;
				$ET_daily = 0;
				$ET_c_daily = 0;
				$ET_eq_daily = 0;
				$T_a_C_sum = 0;
				$R_s_sum = 0;
				$pressure_sum = 0;
				$rainfall_sum = 0;
				$rh_sum = 0;
				$u_z_sum = 0;
				$W_sum = 0;

				$t = -10;	# time of the day (min)
				open (IN, '<./../data/'.$site.'/forcing/'.$site.'_'.$year.'_'.$month.'_'.$day.'.csv') or die;
				open (CANOPY, '>./../out/'.$site.'/'.$year.'/canopy/canopy_'.$year.'_'.$month.'_'.$day.'.csv') or die;
				for ($t=$t_start; $t<=1440; $t+=$t_step) {
					$hour = $t/60;	# (hr)
					$hr = int $hour;
					$min = int $t%60;	# (min)

					#if ($t%30==10) {
					$line=<IN>;
					chomp $line;
					@flux_data = split(/,/, $line);
					#}
					&InputClimateData;

					### Main ##########
					&DailyMeanClimate;
					if ($LAI>0) {
						&LAdistribution;
						if ($S_max <= 0) {
							for ($i_l=$N_layer; $i_l>=0; $i_l--) {
								$S_b_down[$i_l] = 0;
								$S_d_down[$i_l] = 0;
								$S_d_up[$i_l] = 0;
							}
						}
						&CalI;
					}
					if ($S_max <= 0) {
						for ($i_l=$N_layer; $i_l>=0; $i_l--) {
							$S_b_down[$i_l] = 0;
							$S_d_down[$i_l] = 0;
							$S_d_up[$i_l] = 0;
						}
						$LAI_sunlit = 0;
						$LAI_shaded = $LAI;
					} else {
						if ($LAI>0) {
							$L = 1;	# near infra red radiation
							&S_balance;
							$LAI_sunlit = 0;
							if ($R_s_b_total > 0) {
								for ($i_l=$N_layer; $i_l>=1; $i_l--) {
									$LAI_sunlit += $d_LA_su[$i_l];
								}
							}
							$LAI_shaded = $LAI-$LAI_sunlit;
							$L = 0;	# PAR
							&S_balance;
						}
					}

					if ($LAI>0) {
						if ($R_s_total > 0 && $S_max > 0) {
						
							$g_s_sunlit = 1*$W_retention;
							$g_s_shaded = 0.5*$W_retention;
							$C_s_sunlit = $C_a_out*0.8;
							$C_s_shaded = $C_a_out*0.8;
							$C_i_sunlit = $C_a_out*0.7;
							$C_i_shaded = $C_a_out*0.7;
							$diff_T_sunlit = 0;
							$diff_T_shaded = 0;
							$T_c_C_soil = $T_a_C;
							$T_c_K_soil = $T_c_C_soil+273.15;
							$count_i = 0;
							$no_converged = 0;

							while () {

								$diff_T_sunlit_pre = $diff_T_sunlit;
								$diff_T_shaded_pre = $diff_T_shaded;
								$C_s_sunlit_pre = $C_s_sunlit;
								$C_s_shaded_pre = $C_s_shaded;
								$C_i_sunlit_pre = $C_i_sunlit;
								$C_i_shaded_pre = $C_i_shaded;

								&SimultaneousEquations;
								$count_i++;
								if ($S_max>0 && $R_s_total>0) {
									if (abs($diff_T_sunlit-$diff_T_sunlit_pre)<0.01 && abs($diff_T_shaded-$diff_T_shaded_pre)<0.01 && abs($C_s_sunlit-$C_s_sunlit_pre)<0.01 && abs($C_s_shaded-$C_s_shaded_pre)<0.01 && abs($C_i_sunlit-$C_i_sunlit_pre)<0.01 && abs($C_i_shaded-$C_i_shaded_pre)<0.01) {last;}
								} else {
									if (abs($diff_T_sunlit-$diff_T_sunlit_pre)<0.01 && abs($diff_T_shaded-$diff_T_shaded_pre)<0.01) {last;}
								}
								if ($count_i>100) {last;}
							}
							if ($no_converged == 1) {
								&NoConverged;
							}
						} else {
							$g_s_sunlit = 1;
							$g_s_shaded = 0.5;
							$C_s_sunlit = $C_a_out*0.8;
							$C_s_shaded = $C_a_out*0.8;
							$C_i_sunlit = $C_a_out*0.7;
							$C_i_shaded = $C_a_out*0.7;
							$diff_T_sunlit = 0;
							$diff_T_shaded = 0;
							$T_c_C_soil = $T_a_C;
							$T_c_K_soil = $T_c_C_soil+273.15;
							$count = 0;
							$no_converged = 0;
							while () {
								&SimultaneousEquationsNight;
								$diff_T_sunlit_pre = $diff_T_sunlit;
								$diff_T_shaded_pre = $diff_T_shaded;
								$C_s_sunlit_pre = $C_s_sunlit;
								$C_s_shaded_pre = $C_s_shaded;
								if (abs($diff_T_shaded-$diff_T_shaded_pre)<0.01) {last;}
								if ($count_i>10000) {last;}
							}
						}
					} else {
						$T_c_C_soil = $T_a_C;
						$T_c_K_soil = $T_c_C_soil+273.15;
						$Q_n_sunlit = 0;
						$Q_n_shaded = 0;
						&CalConductance;
						&CalET_SW;
					}
					$R_n_sum = $Q_n_sunlit+$Q_n_shaded+$R_n_soil;
					$A_sum = $V_n_sunlit+$V_n_shaded;
					$A_n_sum = $A_c_sunlit+$A_c_shaded;
					$R_m_leaf_sum = $R_d_sunlit+$R_d_shaded;

					$H_sum = $R_n_sum-$lE_sum-$G;
					&WaterBalance;
					
					#####################

					print CANOPY $hr, ":", $min, ",", $T_a_C, ",", $R_s_total, ",", $rainfall, ",", $rh, ",", $u_z, ",", $ET_soil, ",", $lE_soil, ",", $ET_sum, ",", $lE_sum, ",", $A_sum, ",", $R_m_leaf_sum, ",", $T_c_C[$N_layer], ",", $T_c_C[$N_layer/2], ",", $g_s_W_top, ",", $g_s_W_middle, ",", $H_sum, ",", $G, ",", $R_n_sum, ",", $ET_eq, ",", $g_s_W_sum, ",", $sun_duration, "\n";

					$A_n_sum_daily += $A_n_sum*60*$t_step;	# (μmol/m^2/dy)
					$R_m_leaf_sum_daily += $R_m_leaf_sum*60*$t_step;	# (μmol/m^2/dy)
					$A_sum_daily += $A_sum*60*$t_step;	# (μmol/m^2/dy)
					$ET_daily += $ET_sum*60*$t_step;
					$ET_c_daily += $ET_c_sum*60*$t_step;
					$ET_eq_daily += $ET_eq*60*$t_step;



				}

				&CalLAI;
				close (CANOPY);
				close (IN);
				print DEBUG $DOY, ",", $A_n_sum_daily, ",", $A_stem_daily, ",", $A_root_daily, ",", $R_m_leaf_sum_daily, ",", $R_g_leaf_daily, ",", $L_leaf, "\n";
				if ($T_a_C_mean <= -5) {
					$below_minus5++;
					$above_minus5 = 0;
				} else {
					$below_minus5 = 0;
					$above_minus5++;
					if ($T_a_C_mean <=0) {
						$below_0++;
						$above_0 = 0;
					} else {
						$below_0 = 0;
						$above_0++;
					}
				}
				if ($DOY == 1) {
					$T_a_C_mean_7[5] = $T_a_C_mean;
					$T_a_C_mean_7[4] = $T_a_C_mean;
					$T_a_C_mean_7[3] = $T_a_C_mean;
					$T_a_C_mean_7[2] = $T_a_C_mean;
					$T_a_C_mean_7[1] = $T_a_C_mean;
					$T_a_C_mean_7[0] = $T_a_C_mean;
				}
				$T_a_C_mean_7[6] = $T_a_C_mean_7[5];
				$T_a_C_mean_7[5] = $T_a_C_mean_7[4];
				$T_a_C_mean_7[4] = $T_a_C_mean_7[3];
				$T_a_C_mean_7[3] = $T_a_C_mean_7[2];
				$T_a_C_mean_7[2] = $T_a_C_mean_7[1];
				$T_a_C_mean_7[1] = $T_a_C_mean_7[0];
				$T_a_C_mean_7[0] = $T_a_C_mean;
				$T_a_C_mean_7_mean = ($T_a_C_mean_7[0]+$T_a_C_mean_7[1]+$T_a_C_mean_7[2]+$T_a_C_mean_7[3]+$T_a_C_mean_7[4]+$T_a_C_mean_7[5]+$T_a_C_mean_7[6])/7;
				$T_a_C_mean_5_mean = ($T_a_C_mean_7[0]+$T_a_C_mean_7[1]+$T_a_C_mean_7[2]+$T_a_C_mean_7[3]+$T_a_C_mean_7[4])/5;
				$T_a_C_mean_3_mean = ($T_a_C_mean_7[0]+$T_a_C_mean_7[1]+$T_a_C_mean_7[2])/3;
				
				if ($virtual_LAI==0) {
					print $year, "-", $itrn, "-", $DOY, "\t", $LAI, "[", $phenophase, "]";
					if ($leaf_onset==1) {print "*";}
					print "\t", $A_sum_daily*1000, "\t", $R_a_c*1000, "\t", $ET_daily, "\t", $ET_c_daily, "\t", $W_mean, "\n";
					print MEANCLIM $DOY, ",", $T_a_C_mean, ",", $R_s_sum, ",", $pressure_mean, ",", $rainfall_sum, ",", $rh_mean, ",", $u_z_mean, ",", $W_mean, "\n";
					print DAILY $DOY, ",", $A_sum_daily*1000, ",", $R_a_c*1000, ",", $R_m_leaf_sum_daily*1000, ",", $ET_daily, ",", $ET_c_daily, ",", $ET_eq_daily, "\n";
					print ALLOCATION $DOY, ",", $LAI, ",", $C_leaf, ",", $C_stem, ",", $C_root, ",", $C_all, ",", $LAI_g, ",", $phenophase, "\n";

					$A_sum_yearly += $A_sum_daily;
					$R_m_leaf_sum_yearly += $R_m_leaf_sum_daily;
					$R_a_c_yearly += $R_a_c;
					$ET_yearly += $ET_daily;
					$ET_c_yearly += $ET_c_daily;
					$ET_eq_yearly += $ET_eq_daily;
				} elsif ($virtual_LAI==10) {
					for ($virtual_LAI_day=1; $virtual_LAI_day<=$virtual_LAI_succeed_day; $virtual_LAI_day++) {
						print $virtual_DOY[$virtual_LAI_day-1], ",", $virtual_A_sum_daily[$virtual_LAI_day-1]*1000, ",", $virtual_R_a_c[$virtual_LAI_day-1]*1000, ",", $virtual_ET_daily[$virtual_LAI_day-1], ",", $virtual_ET_c_daily[$virtual_LAI_day-1], ",", $virtual_ET_eq_daily[$virtual_LAI_day-1], "\n";
						print MEANCLIM $virtual_DOY[$virtual_LAI_day-1], ",", $virtual_T_a_C_mean[$virtual_LAI_day-1], ",", $virtual_R_s_sum[$virtual_LAI_day-1], ",", $virtual_pressure_mean[$virtual_LAI_day-1], ",", $virtual_rainfall_sum[$virtual_LAI_day-1], ",", $virtual_rh_mean[$virtual_LAI_day-1], ",", $virtual_u_z_mean[$virtual_LAI_day-1], ",", $virtual_W_mean[$virtual_LAI_day-1], "\n";
						print DAILY $virtual_DOY[$virtual_LAI_day-1], ",", $virtual_A_sum_daily[$virtual_LAI_day-1]*1000, ",", $virtual_R_a_c[$virtual_LAI_day-1]*1000, ",", $virtual_R_m_leaf_sum_daily[$virtual_LAI_day-1]*1000, ",", $virtual_ET_daily[$virtual_LAI_day-1], ",", $virtual_ET_c_daily[$virtual_LAI_day-1], ",", $virtual_ET_eq_daily[$virtual_LAI_day-1], "\n";
						print ALLOCATION $virtual_DOY[$virtual_LAI_day-1], ",", $virtual_LAI[$virtual_LAI_day-1], ",", $virtual_C_leaf[$virtual_LAI_day-1], ",", $virtual_C_stem[$virtual_LAI_day-1], ",", $virtual_C_root[$virtual_LAI_day-1], ",", $virtual_C_all[$virtual_LAI_day-1], ",", $virtual_LAI_g[$virtual_LAI_day-1], ",", 1, "\n";

						$A_sum_yearly += $virtual_A_sum_daily[$virtual_LAI_day-1];
						$R_m_leaf_sum_yearly += $virtual_R_m_leaf_sum_daily[$virtual_LAI_day-1];
						$R_a_c_yearly += $virtual_R_a_c[$virtual_LAI_day-1];
						$ET_yearly += $virtual_ET_daily[$virtual_LAI_day-1];
						$ET_c_yearly += $virtual_ET_c_daily[$virtual_LAI_day-1];
						$ET_eq_yearly += $virtual_ET_eq_daily[$virtual_LAI_day-1];
						
						if ($virtual_DOY[$virtual_LAI_day-1]==$DOY_max) {
							print $DOY, ",", $A_sum_yearly, ",", $R_m_leaf_sum_yearly, ",", $R_a_c_yearly, ",", $ET_yearly, ",", $ET_c_yearly, ",", $W_mean, "\n";
							$NPP_yearly = ($A_sum_yearly-$R_a_c_yearly)*10;	# (tC/ha)
							$itrn++;
							if ($itrn>=3 && $DOY==$DOY_max) {last LOOP;}
							$A_sum_yearly = 0;
							$R_a_c_yearly = 0;
							$R_m_leaf_sum_yearly = 0;
							$NPP_yearly = print0;
							$ET_yearly = 0;
							$ET_c_yearly = 0;
							$ET_eq_yearly = 0;
						}
					}
					$virtual_LAI=0;
					$leaf_onset = 1;
					$phenophase = 1;
					$growing_dys = $virtual_LAI_succeed_day;
				}

				print ">", $growing_dys, "\t", $Cg_leaf, "(", $LAI_g, ")", "\t", $Cd_leaf, "(", $LAI_d, ")", "\t", $C_leaf, "(", $LAI, ")", "\n";
			}

		}

		print DAILY "#SUM", ",", $A_sum_yearly, ",", $R_m_leaf_sum_yearly, ",", $R_a_c_yearly, ",", $ET_yearly, ",", $ET_c_yearly, ",", $ET_eq_yearly, "\n";
		print YEARLY $year, ",", $A_sum_yearly*10, ",", $R_a_c_yearly*10, ",", $NPP_yearly, ",", $ET_yearly, ",", $ET_c_yearly, ",", $ET_eq_yearly, "\n";
		
		close (ALLOCATION);
		close (DAILY);
		close (MEANCLIM);
		close (DEBUG);

	}

	$diff_C_stem = abs($C_stem_pre-$C_stem);
	$diff_C_root = abs($C_root_pre-$C_root);
	if ($diff_C_stem > 0.01) {$diff_count_yr++;}
	if ($diff_C_root > 0.01) {$diff_count_yr++;}
	print $itrn, "\t", $diff_C_stem, "\t", $diff_C_root, "\n";
	if ($diff_count_yr==0 && $DOY==$DOY_max) {last;}

	$C_stem_pre = $C_stem;
	$C_root_pre = $C_root;

}
close (YEARLY);

########################################################

sub MakeDirectory {
	
	if (!-d './../out/'.$site.'/'.$year) {
		mkdir './../out/'.$site.'/'.$year;
	}
	if (!-d './../out/'.$site.'/'.$year.'/canopy') {
		mkdir './../out/'.$site.'/'.$year.'/canopy';
	}
	if (!-d './../out/'.$site.'/'.$year.'/daily') {
		mkdir './../out/'.$site.'/'.$year.'/daily';
	}
	if (!-d './../out/'.$site.'/'.$year.'/meanclim') {
		mkdir './../out/'.$site.'/'.$year.'/meanclim';
	}
	if (!-d './../out/'.$site.'/'.$year.'/debug') {
		mkdir './../out/'.$site.'/'.$year.'/debug';
	}

	if (!-f './../out/'.$site.'/'.$year.'/daily/allocation_'.$year.'.csv') {
		open (ALLOCATION,'>./../out/'.$site.'/'.$year.'/daily/allocation_'.$year.'.csv') or die;
	} else {
		open (ALLOCATION,'>>./../out/'.$site.'/'.$year.'/daily/allocation_'.$year.'.csv') or die;
	}
	if (!-f './../out/'.$site.'/'.$year.'/daily/daily_'.$year.'.csv') {
		open (DAILY,'>./../out/'.$site.'/'.$year.'/daily/daily_'.$year.'.csv') or die;
	} else {
		open (DAILY,'>>./../out/'.$site.'/'.$year.'/daily/daily_'.$year.'.csv') or die;
	}
	if (!-f './../out/'.$site.'/'.$year.'/meanclim/meanclim_'.$year.'.csv') {
		open (MEANCLIM, '>./../out/'.$site.'/'.$year.'/meanclim/meanclim_'.$year.'.csv') or die;
	} else {
		open (MEANCLIM, '>>./../out/'.$site.'/'.$year.'/meanclim/meanclim_'.$year.'.csv') or die;
	}
	if (!-f './../out/'.$site.'/'.$year.'/debug/debug_'.$year.'.csv') {
		open (DEBUG, '>./../out/'.$site.'/'.$year.'/debug/debug_'.$year.'.csv') or die;
	} else {
		open (DEBUG, '>>./../out/'.$site.'/'.$year.'/debug/debug_'.$year.'.csv') or die;
	}

}

sub InputClimateData {
	
	if ($t==$t_start) {
		$rainfall_pre = 0;
	} else {
		$rainfall_pre = $rainfall;
	}
	
	if ($t==$t_start) {
		$T_a_C_pre = 20;
	} else {
		$T_a_C_pre = $T_a_C;
	}
	$T_a_K = $T_a_C+273.15;	# air temperature (°K)
	if ($t==$t_start) {
		$pressure_pre = 1013.25*(1-0.0065/$T_a_K*$elv)**5.2552;
	} else {
		$pressure_pre = $pressure;
	}

	if ($t==$t_start) {
		$rh_pre = 60;
	} else {
		$rh_pre = $rh;	# (%)
	}

	if ($t==$t_start) {
		$u_z_pre = 2;
	} else {
		$u_z_pre = $u_z;	# (m/s)
	}

	if ($t==15) {
		$sun_duration_pre = 0;
	} else {
		$sun_duration_pre = $sun_duration;	# (min/10min)
	}

	if ($t==$t_start) {
		$R_s_total_pre = 0;
	} else {
		$R_s_total_pre = $R_s_total;	# (min/10min)
	}

	$solar_elevation_pre = $solar_elevation;
	if ($flux_data[2] == -9999) {
		$T_a_C = $T_a_C_pre;
	} else {
		$T_a_C = $flux_data[2]; # air temperature (℃)
	}
	$T_a_K = $T_a_C + 273.15; # air temperature (K)

	if ($flux_data[1] == -9999) {
		$rainfall = $rainfall_pre;
	} else {
		$rainfall = $flux_data[1]; # * 60 * $t_step; # (mm/30min)
	}

	if ($flux_data[3] == -9999) {
		$rh = $rh_pre;
	} else {
		$rh = $flux_data[3] > 0.0001 ? $flux_data[3] : 0.0001;
	}

	if ($flux_data[5] == -9999) {
		$u_z = $u_z_pre;
	} else { 
		$u_z = $flux_data[5] > 0.1 ? $flux_data[5] : 0.1;
	}

	# $pressure = $flux_data[9] / 100; # (hPa)
	$pressure = 1013.25 * (1 - 0.0065 / $T_a_K * $elv) ** 5.2552;

	if ($flux_data[6] == -9999) {
		$R_s_total = $R_s_total_pre;
	} else {
		$R_s_total = $flux_data[6]; # (W/m^2)
	}
	
	$T_a_K = $T_a_C+273.15;	# air temperature (°K)
	&CalPressVPD;
	&CalSunDuration;
	$I_max = $I*sin($solar_elevation);
	$R_s_b_total = $R_s_total-$R_s_d_total;
	if ($solar_elevation > 0) {
		if ($solar_elevation_pre < 0) {$sun_duration_sum = 0;	$sun_duration_count = 0;}
		$sun_duration_sum += $sun_duration;
		$sun_duration_count++;
	} else {
		if ($sun_duration_count==-9) {
			$sun_duration = 0.8;
		} else {
			$sun_duration = $sun_duration_sum/$sun_duration_count;
		}
	}
	if ($S_max > 0) {
		$sun_duration = $R_s_b_total/($S_max*(1-0.156));
	}
	if ($sun_duration < 0) {
		$sun_duration = 0;
	} elsif ($sun_duration > 1) {
		$sun_duration = 1;
	}
	$C = 0.826*$sun_duration**3-1.234*$sun_duration**2+1.135*$sun_duration+0.298;
	if ($C > 1) {$C = 1;}
	$L_down = $SteBol_const*$T_a_K**4*(1-(1-$R_l_f/($SteBol_const*$T_a_K**4))*$C);	# downward longwave radiation ~ (Kondo & Xu, 1997)

	$T_a_K = $T_a_C+273.15;	# air temperature (°K)

	$R_s_b[0] = $R_s_b_total*0.43;	# PAR content: 43%;	4.6 μmol/J * 43 %	(~ Larcher, 1995)
	$R_s_d[0] = $R_s_d_total*0.57;	# PAR content: 57%;	4.2 μmol/J * 57 %
	$R_s_b[1] = $R_s_b_total*0.57;	# 4.6 μmol/J * 43 %	(~ Larcher, 1995)
	$R_s_d[1] = $R_s_d_total*0.43;	# 4.2 μmol/J * 57 %

}

sub CalPressVPD {

	$e_s_Ta = $c0*exp($a0*$T_a_C/($T_a_C+$b0));	# saturation vapor pressure (hPa)
	$e_a_Ta = $e_s_Ta * $rh / 100;	# vapor pressure (hPa)
	$q = 0.622*$e_a_Ta/($pressure-0.379 * $e_a_Ta);	# specific humidity (kg/kg)
	# vapor pressure deficit (hPa)
	$VPD_a = $e_s_Ta-$e_a_Ta;
	# mole fraction of vapor
	$e_frac_s = $e_s_Ta/$pressure;
	$e_frac_a = $e_a_Ta/$pressure;
	$VPD_frac = $e_frac_s-$e_frac_a;

	# $rho_m: molar density (mol/m^3)
	$rho_m = $pressure*100/$R/($T_a_K);
	# $rho_a: air density (kg/m^3)
	$rho_a = $rho_m*$M_d*(1-0.378*$e_a_Ta/$pressure);
	# $M_a: molecular mass of air (kg/mol)
	$M_a = $rho_a/$rho_m;
	# $C_p: specific heat of air at constant pressure (J/kg/K)
	$C_p = $C_pd*(1+0.84*$q);
	# $c_p (J/mol/K)
	$c_p = $C_pd*(1+0.84*$q)*$M_a;
	# latent heat of vaporization (J/mol)
	$lambda_mol = 56780.3-42.84*$T_a_K;
	$lambda_kg = $lambda_mol/0.018;
	# gamma: psychrometric constant (hPa/K)
	$gamma = $C_p*$pressure/0.622/$lambda_kg;
	
	$Delta = $a0*$b0*$c0/($T_a_C+$b0)**2*exp($a0*$T_a_C/($T_a_C+$b0));	# slope of the saturated vapor pressure (hPa/°C)
	
}

sub CalSunDuration {	# to calculate possible sunshine duration

	&vapor;
	
	$eta = (2*$pi/$DOY_max)*$DOY;
	$SA2 = 4.871+$eta+0.033*sin($eta);
	$A22 = 0.398*sin($SA2);
	$declination = atan($A22/sqrt(1-$A22**2));	# $declination (rad)
	$DD2 = 1.00011+0.034221*cos($eta)+0.00128*sin($eta)+0.000719*cos(2*$eta)+0.000077*sin(2*$eta);
	$j1 = (0.066+0.34*sqrt($dust))*($albedo-0.15);
	$C1 = 0.21-0.2*$dust;
	$C2 = 0.15-0.2*$dust;
	if ($dust>0.3) {$C1=0.15; $C2=0.09;}
	$Io = 1365*$DD2;
	$F1 = 0.056+0.16*sqrt($dust);
	$F2 = 0.075+0.65*$dust;
	
	$hour_angle = asin(sqrt(sin($pi/4+($latitude-$declination+$refraction)/2)*sin($pi/4-($latitude-$declination-$refraction)/2)/(cos($latitude)*cos($declination))));
	# $sun_duration_pos: possible sunshine duration (hr)
	$sun_duration_pos = 4*$hour_angle/0.2618;

	$t_local = $t+4*($lon-$lon_LST);
	$H = $pi*($t_local/60/12-1);
	$COSTHETA = sin($latitude)*sin($declination)+cos($latitude)*cos($declination)*cos($H);
	$solar_elevation = asin($COSTHETA);	# (rad)

	$So = $Io*$COSTHETA;
	$M = $pressure/1013.25/$COSTHETA;
	if ($COSTHETA<=0) {
		$I=0;	$So=0;	$S=0;	$S_max= 0;	$R_s_total = 0;	$R_s_d_total = 0;
		$sun_duration = $sun_duration_pre;
	} else {
		$I1 = 0.014*($M+7+2*$LW10)*$LW10;
		$I2 = 0.02*($M+5.5+1.5*$LW10)*$LW10;
		$I = $Io*($C2+0.75*10**(-$M*$F2))*(1-$I2);
		$S_max = $So*($C1+0.7*10**(-$M*$F1))*(1-$I1)*(1+$j1);	# 全天日射量 (W/m^2)
		$I_max = $I*sin($solar_elevation);
		$sun_duration = $R_s_total/$S_max;
		$K_Tt = $R_s_total/$So;
		if ($K_Tt < 0) {$K_Tt = 0;}
		if ($K_Tt <= 0.15) {
			$R_s_d_total = $R_s_total;
		} elsif ($K_Tt <= 0.22) {
			$R_s_d_total = $R_s_total*(1-0.09*$K_Tt);
		} elsif ($K_Tt <= 0.80) {
			$R_s_d_total = $R_s_total*(0.9511-0.1604*$K_Tt+4.388*$K_Tt**2-16.638*$K_Tt**3+12.336*$K_Tt**4);
		} else {
			$R_s_d_total = $R_s_total*0.156;
		}
	}

	if ($R_s_d_total > $R_s_total) {
		$R_s_d_total = $R_s_total;
	} elsif ($R_s_total>0 && $R_s_d_total<=0) {
		$R_s_d_total = 0.01;
	}

	# $R_l_f: downward longwave radiation under clear sky (W/m^2) ~ (Zuo et al., 1991)
	$xi = 46.5*$e_a_Ta/$T_a_K;
	$emissivity_a = 1-(1+$xi)*exp(-sqrt(1.2+3*$xi));
	$R_l_f = $SteBol_const*$T_a_K**4*$emissivity_a;	# Prata (1996)
	
	if ($solar_elevation > 0) {
		$k_b_black = 0.5/sin($solar_elevation);
	}
	if ($k_b_black>1 or $solar_elevation<=0) {
		$k_b_black = 1;
	}
	$k_b = sqrt(0.5)*$k_b_black;
	$k_d_black = 0.8;
	$k_d = sqrt(0.5)*$k_d_black;

}

sub vapor {
	$vapor_pressure = $e_a_Ta;	# (hPa)
	$albedo = 0.15;
	$dust = 0.05;	# 混濁係数


	$dew_point = (27*log($vapor_pressure/6.11))/(17.27-log($vapor_pressure/6.11));	# (°C)
	$Bo = 0;
	if ($dew_point<-5) {
		$LNW = 0.0622*$dew_point+1.958-$Bo;
	} elsif ($dew_point<23) {
		$LNW = 0.0714*$dew_point+2.003-$Bo;
	} else {
		$LNW = 0.0345*$dew_point+2.851-$Bo;
	}
	$WSTOP = exp($LNW);
	$LW10 = 0.4343*log((1.234*$WSTOP-0.21)/10);

}

sub LAdistribution {	# describe leaf area distribution for the multi-layer canopy model

	$sum = 0;
	for ($i_l=$N_layer; $i_l>=1; $i_l--) {
		$x = $d_z_h*(2*$i_l-1)/2;
		$sum += $x**($alpha[$leaf_type]-1)*(1-$x)**($beta[$leaf_type]-1)*$d_z_h;
	}

	$f_sum = 0;	$LA_z_max = 0;
	for ($i_l=$N_layer; $i_l>=1; $i_l--) {
		$x = $d_z_h*(2*$i_l-1)/2;
		$f[$i_l] = $x**($alpha[$leaf_type]-1)*(1-$x)**($beta[$leaf_type]-1)/$sum*$LAI;	# (m^2/m^2)
		$z_c[$i_l] = $x*$height_c;	# (m)
		$LA_z[$i_l] = $f[$i_l]/$height_c;	# $LA_z: Leaf Area Density (m^2/m^3)
		$d_LA[$i_l] = $LA_z[$i_l]*$d_z;	# (m^2/m^2)
		$d_LA_g[$i_l] = $d_LA[$i_l]*$Cg_leaf/$C_leaf;	# (m^2/m^2)
		if ($LA_z[$i_l] > $LA_z_max) {
			$LA_z_max = $LA_z[$i_l];
		}
	}
	$z_c[0] = 0;
	$z_c[$N_layer+1] = $d_z_h*(2*($N_layer+1)-1)/2*$height_c;

}

sub ExtinctionCoeffcientFunction {

	$extinction_coeffcient = (1-exp(-$_[1]*$_[0]))/$_[1];

}

sub I_b {

	&G_layer;
	$I_b[$i_l] = 1-$G_layer*$d_LA[$i_l]/sin($x);	# the probability of no contact with direct radiation within a layer between z and z + dz

	if ($I_b[$i_l] < 0) {
		$I_b[$i_l] = 0;
	} elsif ($I_b[$i_l] > 1) {
		$I_b[$i_l] = 1;
	}

}

sub I_d {	# the probability of no contact with diffuse radiation within a layer between z and z + dz

	$I_d[$i_l] = 0;
	$d_angle = $pi/2/$N_angle;
	for ($i_a3=1; $i_a3<=$N_angle; $i_a3++) {
		$x = $d_angle*(2*$i_a3-1)/2;
		&I_b;
		$I_d[$i_l] += $I_b[$i_l]*sin($x)*cos($x)*$d_angle;
	}
	$I_d[$i_l] *= 2;
	if ($I_d[$i_l] < 0) {
		$I_d[$i_l] = 0;
	} elsif ($I_d[$i_l] > 1) {
		$I_d[$i_l] = 1;
	}

}

sub G_layer {

	$angle_ave = 50*$pi/180;	# $angle_ave: the mean of leaf inclination (degree)
	$angle_SD = 18*$pi/180;	# $angle_SD: the standard deviation of leaf inclination (degree)

	$G_layer = 0;
	$d_angle = $pi/2/$N_angle;
	$d2_angle = 2*$pi/$N2_angle;
	for ($i_a2=1; $i_a2<=$N2_angle; $i_a2++) {

		$leaf_orientation_angle = $d2_angle*(2*$i_a2-1)/2;

		$G_layer0 = 0;
		for ($i_a=1; $i_a<=$N_angle; $i_a++) {
			$leaf_inclination = $d_angle*(2*$i_a-1)/2;	# (rad)
			$g_alpha = exp(-(($leaf_inclination-$angle_ave)/$angle_SD)**2/2)/($angle_SD*sqrt(2*$pi));

			$G_leaf = abs(cos($leaf_inclination)*sin($x)+sin($leaf_inclination)*cos($leaf_orientation_angle)*cos($x));
			$G_layer0 += $g_alpha*$G_leaf*$d_angle;	# $G_layer: the ratio of the area of leaves, projected into a plane normal to the solar elevation, to the leaf area index within the layer
		}
		
		$G_layer += $G_layer0/(2*$pi)*$d2_angle;

	}

}

sub CalI {

	for ($i_l=$N_layer; $i_l>=1; $i_l--) {
		&I_d;
		$x = $solar_elevation;
		&I_b;
	}

}

sub S_b_balance {	# calculate direct radiation balance within each layer

	for ($i_l=$N_layer; $i_l>=1; $i_l--) {
		$S_b_down[$i_l-1] = $S_b_down[$i_l]*$I_b[$i_l];
	}

}

sub S_d_balance {	# calculate diffuse radiation balance within each layer

	$S_d_up[0] = ($S_b_down[0]+$S_d_down[0])*$albedo_soil[$L];
	for ($i_l=$N_layer; $i_l>=1; $i_l--) {
		
		$S_d_down[$i_l-1] = $S_d_down[$i_l]*($transmissivity[$L]*(1-$I_d[$i_l])+$I_d[$i_l])+$S_d_up[$i_l-1]*$reflectance[$L]*(1-$I_d[$i_l])+$S_b_down[$i_l]*$transmissivity[$L]*(1-$I_b[$i_l]);
		$S_d_up[$N_layer-($i_l-1)] = $S_d_up[$N_layer-1-($i_l-1)]*($transmissivity[$L]*(1-$I_d[$N_layer-($i_l-1)])+$I_d[$N_layer-($i_l-1)])+$S_d_down[$N_layer-($i_l-1)]*$reflectance[$L]*(1-$I_d[($N_layer-($i_l-1))])+$S_b_down[$N_layer-($i_l-1)]*$reflectance[$L]*(1-$I_b[$N_layer-($i_l-1)]);
	}

}

sub S_balance {	# calculate radiation balance within each layer

	$S_b_down[$N_layer] = $R_s_b[$L];
	$S_d_down[$N_layer] = $R_s_d[$L];
	&S_b_balance;

	if ($S_max>0 && $G_layer>0) {
		for ($i_l=$N_layer; $i_l>=1; $i_l--) {
			if ($S_b_down[$N_layer]<=0) {
				$d_LA_su[$i_l] = 0;
			} else {
				$d_LA_su[$i_l] = ($S_b_down[$i_l]-$S_b_down[$i_l-1])/$S_b_down[$N_layer]*sin($solar_elevation)/$G_layer;
				if ($d_LA_su[$i_l] > $d_LA[$i_l]) {$d_LA_su[$i_l] = $d_LA[$i_l];}
			}
			$d_LA_sh[$i_l] = $d_LA[$i_l]-$d_LA_su[$i_l];
		}
	}

	for ($i_l=$N_layer; $i_l>=1; $i_l--) {
		$S_d_down[$i_l-1] = $S_d_down[$N_layer]*$i_l/$N_layer;
		$S_d_up[$i_l-1] = $S_b_down[$N_layer]*$reflectance[$L]*$reflectance[$L]*$i_l/$N_layer;

	}
	$S_d_up[$N_layer] = $R_s_b[$L]*$reflectance[$L]*$reflectance[$L];

	$i = 0;
	while (0<1) {	# iterate until calculations are converged

		$diff_count = 0;

		for ($i_l=$N_layer; $i_l>=1; $i_l--) {
			$S_d_down_pre[$i_l-1] = $S_d_down[$i_l-1];
			$S_d_up_pre[$i_l-1] = $S_d_up[$i_l-1];
		}
		$S_d_up_pre[$N_layer] = $S_d_up[$N_layer];
		&S_d_balance;
		for ($i_l=$N_layer; $i_l>=1; $i_l--) {
			$diff_S_d_down[$i_l-1] = abs($S_d_down[$i_l-1]-$S_d_down_pre[$i_l-1]);
			$diff_S_d_up[$i_l-1] = abs($S_d_up[$i_l-1]-$S_d_up_pre[$i_l-1]);
		}
		$diff_S_d_up[$N_layer] = abs($S_d_up[$N_layer]-$S_d_up_pre[$N_layer]);

		for ($i_l=$N_layer; $i_l>=1; $i_l--) {
			if ($diff_S_d_down[$i_l-1]>$diff_level or $diff_S_d_up[$i_l-1]>$diff_level) {
				$diff_count++;
			}
			if ($diff_S_d_up[$N_layer]>$diff_level) {
				$diff_count++;
			}
		}

		if ($diff_count == 0) {last;}

		$i++;
	}

			$Q_shaded[$L] = 0;
			$Q_sunlit[$L] = 0;
			for ($i_l=$N_layer; $i_l>=1; $i_l--) {
				if ($S_b_down[$N_layer]<=0) {
					$d_LA_su[$i_l] = 0;
				} else {
					$d_LA_su[$i_l] = ($S_b_down[$i_l]-$S_b_down[$i_l-1])/$S_b_down[$N_layer]*sin($solar_elevation)/$G_layer;
					if ($d_LA_su[$i_l] > $d_LA[$i_l]) {$d_LA_su[$i_l] = $d_LA[$i_l];}
				}
				$d_LA_sh[$i_l] = $d_LA[$i_l]-$d_LA_su[$i_l];
				$Q_shaded[$L] += (1-$transmissivity[$L])*(1-$transmissivity[$L]-$reflectance[$L])*($S_d_down[$i_l]+$S_d_up[$i_l-1])*$d_LA_sh[$i_l];
				$Q_sunlit[$L] += (1-$transmissivity[$L])*(1-$transmissivity[$L]-$reflectance[$L])*($S_d_down[$i_l]+$S_d_up[$i_l-1]+$S_b_down[$N_layer]*$G_layer/sin($solar_elevation))*$d_LA_su[$i_l];
			}

}

sub L_balance {	# calculate longwave radiation absorbed by sunlit/sunshaed leaves

	&ExtinctionCoeffcientFunction($LAI,$k_d_black);
	$extinction_coeffcient_d_black = $extinction_coeffcient;
	&ExtinctionCoeffcientFunction($LAI,2*$k_d_black);
	$extinction_coeffcient_2d_black = $extinction_coeffcient;
	&ExtinctionCoeffcientFunction($LAI,$k_b_black+$k_d_black);
	$extinction_coeffcient_b_black_plus_d_black = $extinction_coeffcient;
	&ExtinctionCoeffcientFunction($LAI,$k_b_black-$k_d_black);
	$extinction_coeffcient_b_black_minus_d_black = $extinction_coeffcient;
	$Q_long_isothermal_sunlit = ($L_down-$emissivity_c*$SteBol_const*$T_a_K**4)*$LAI_sunlit/$LAI;
	$Q_long_isothermal_shaded = ($L_down-$emissivity_c*$SteBol_const*$T_a_K**4)*$LAI_shaded/$LAI;
	if ($LAI<1) {
		$Q_long_isothermal_sunlit *= $LAI;
		$Q_long_isothermal_shaded *= $LAI;
	}
	if (abs($Q_long_isothermal_sunlit) > $R_s_total) {$Q_long_isothermal_sunlit = 0;}
	if (abs($Q_long_isothermal_shaded) > $R_s_total) {$Q_long_isothermal_shaded = 0;}

}

sub CalRn {

	$Q_n_isothermal_sunlit = $Q_sunlit[0]+$Q_sunlit[1]+$Q_long_isothermal_sunlit;
	$Q_n_isothermal_shaded = $Q_shaded[0]+$Q_shaded[1]+$Q_long_isothermal_shaded;
	
}

sub ResConst {
	$g_w_sunlit = 1/(1/$g_a+1/$g_b_sunlit+1/$g_s_sunlit);
	$g_w_shaded = 1/(1/$g_a+1/$g_b_shaded+1/$g_s_shaded);
	$gamma_modified_sunlit = $gamma*($g_h_sunlit+$g_r_sunlit)/$g_w_sunlit;
	$gamma_modified_shaded = $gamma*($g_h_shaded+$g_r_shaded)/$g_w_shaded;

	$g_c_sunlit = 1/(1/$g_a+1/(0.11/0.147*$g_b_forced_sunlit+0.038/0.055*$g_b_free_sunlit)+1/0.64/$g_s_sunlit);
	$g_c_shaded = 1/(1/$g_a+1/(0.11/0.147*$g_b_forced_shaded+0.038/0.055*$g_b_free_shaded)+1/0.64/$g_s_shaded);
	$C_s_sunlit = $C_i_sunlit+($C_a_out-$C_i_sunlit)*$g_c_sunlit/0.64/$g_s_sunlit;
	$C_s_shaded = $C_i_shaded+($C_a_out-$C_i_shaded)*$g_c_shaded/0.64/$g_s_shaded;

	$VPD_s_sunlit = ($VPD_a+$Delta*$diff_T_sunlit)*$g_w_sunlit/$g_s_sunlit;
	$VPD_s_shaded = ($VPD_a+$Delta*$diff_T_shaded)*$g_w_shaded/$g_s_shaded;

}

sub CalDiffT {

	$diff_T_sunlit = $gamma_modified_sunlit/($Delta+$gamma_modified_sunlit)*$Q_n_isothermal_sunlit/$c_p/($g_h_sunlit+$g_r_sunlit)-$VPD_a/($Delta+$gamma_modified_sunlit);
	$diff_T_shaded = $gamma_modified_shaded/($Delta+$gamma_modified_shaded)*$Q_n_isothermal_shaded/$c_p/($g_h_shaded+$g_r_shaded)-$VPD_a/($Delta+$gamma_modified_shaded);

}

sub CalConductance {
	if ($vegetation_type==1 && $leaf_type==0 && $LAI==0) {
		$d = 0.01;
		$z_0_m = 0.01;
	} else {
		if ($vegetation_type == 0) {
			$z_0_m = 0.10*$height_c;	# for forest ~ (Verseghy et al., 1993)
			$z_0_h = $z_0_m/2.0;	# for forest
			$d = 0.70*$height_c;	# for forest ~ (Verseghy et al., 1993)
		} elsif ($vegetation_type == 1) {
			$z_0_m = 0.123*$height_c;	# for crop and grass ~ (Monteith, 1981)
			$z_0_h = $z_0_m/12.0;	# for grass
			$d = 0.67*$height_c;	# for crop and grass ~ (Monteith, 1981)
		}
		$u_h = $u_z*log(($height_c-$d)/$z_0_m)*log(($height_c-$d)/$z_0_h)/log(($z-$d)/$z_0_m)/log(($z-$d)/$z_0_h);
		$u = $u_h*exp(-$u_attenuation_coef*0.5);	# within-canopy wind speed ~ (Cionco, 1972)
		if ($u<0.0001) {$u = 0.0001;}
		$r_a_alpha = log(($z-$d)/$z_0_m)/($Karman_const**2*$u_z)*(log(($z-$d)/($height_c-$d))+$height_c/(2.5*($height_c-$d))*(exp(2.5*(1-($d+$z_0_m)/$height_c))-1));
		$r_s_a_alpha = log(($z-$d)/$z_0_m)/($Karman_const**2*$u_z)*$height_c/(2.5*($height_c-$d))*(exp(2.5)-exp(2.5*(1-($d+$z_0_m)/$height_c)));
	}
		$z_0_e = 0.01;
	
		$u_soil = $u_h*exp(-$u_attenuation_coef);
		if ($u_soil<0.0001) {$u_soil = 0.0001;}

		$r_s_a_0 = log($z/$z_0_e)*log(($d+$z_0_m)/$z_0_e)/($Karman_const**2*$u_z);
		if ($r_s_a_0 < 0) {
			$r_s_a_0 = 0;
		}
	$r_a_0 = (log($z/$z_0_e))**2/($Karman_const**2*$u_z)-$r_s_a_0;
	if ($r_a_0 <= 0) {$r_a_0 = 0.1}
	if ($LAI <= 4) {
		$r_a = $r_a_alpha*$LAI/4+$r_a_0*(1-$LAI/4);
		$r_a_soil = $r_s_a_alpha*$LAI/4+$r_s_a_0*(1-$LAI/4);
	} else {
		$r_a = $r_a_alpha;
		$r_a_soil = $r_s_a_alpha;
	}
		$g_a = 1/$r_a*$pressure*100/$R/$T_a_K;

	if ($LAI > 0) {
		&ExtinctionCoeffcientFunction($LAI,0.5*$k_u+$k_b_black);
		$extinction_coeffcient_u_plus_b_black = $extinction_coeffcient;
		&ExtinctionCoeffcientFunction($LAI,0.5*$k_u);
		$extinction_coeffcient_u = $extinction_coeffcient;
		$g_b_forced = 0.147*sqrt($u_h/$d_leaf[$leaf_type]);	# $g_b_forced: boundary conductance due to forced convection for H2O (mol/m2/s)
		$g_b_forced_sunlit = $g_b_forced*$extinction_coeffcient_u_plus_b_black;
		if ($extinction_coeffcient_u-$extinction_coeffcient_u_plus_b_black < 0.0001) {
			$g_b_forced_shaded = $g_b_forced*0.0001;
		} else {
			$g_b_forced_shaded = $g_b_forced*($extinction_coeffcient_u-$extinction_coeffcient_u_plus_b_black);
		}
		if ($g_b_free_sunlit eq 'nan') {$g_b_free_sunlit = 0;}
		if ($g_b_free_shaded eq 'nan') {$g_b_free_shaded = 0;}
		$g_b_sunlit = $g_b_free_sunlit+$g_b_forced_sunlit;
		$g_b_shaded = $g_b_free_shaded+$g_b_forced_shaded;
		$g_h_sunlit = 1/(1/$g_a+1/(0.135/0.147*$g_b_forced_sunlit+0.05/0.055*$g_b_free_sunlit));
		$g_h_shaded = 1/(1/$g_a+1/(0.135/0.147*$g_b_forced_shaded+0.05/0.055*$g_b_free_shaded));
		$r_a_c = 2/($g_b_sunlit+$g_b_shaded)*$pressure*100/$R/$T_a_K;	# (s/m)
		$r_s_c = 2/($g_s_sunlit+$g_s_shaded)*$pressure*100/$R/$T_a_K;	# (s/m)
	}

	if ($g_b_free_soil eq 'nan') {$g_b_free_soil = 0;}
	
	$h_beta = 1/4*(1-cos($W/$W_capacity*$pi))**2;
}

sub SimultaneousEquations {
	$T_c_C_sunlit = $T_a_C+$diff_T_sunlit;
	$T_c_K_sunlit = $T_c_C_sunlit+273.15;
	$T_c_C_shaded = $T_a_C+$diff_T_shaded;
	$T_c_K_shaded = $T_c_C_shaded+273.15;
	&L_balance;
	&CalRn;
	&CalConductance;
	&ResConst;
	&SimultaneousPhotosynthesis;
	&CalConductance;
	&CalET_SW;
	&CalDiffT;
	if (abs($diff_T_sunlit) > 10) {$diff_T_sunlit = 0;	$no_converged = 1;	last;}
	if (abs($diff_T_shaded) > 10) {$diff_T_shaded = 0;	$no_converged = 1;	last;}
	
}

sub SimultaneousEquationsNight {
	$T_c_C_sunlit = $T_a_C+$diff_T_sunlit;
	$T_c_K_sunlit = $T_c_C_sunlit+273.15;
	$T_c_C_shaded = $T_a_C+$diff_T_shaded;
	$T_c_K_shaded = $T_c_C_shaded+273.15;
			&L_balance;
			$Q_n_isothermal_sunlit = $Q_long_isothermal_sunlit;
			$Q_n_isothermal_shaded = $Q_long_isothermal_shaded;
			&CalConductance;
			&ResConst;
			&SimultaneousPhotosynthesis;
			&CalConductance;
			&CalET_SW;
			&CalDiffT;
	
}

sub SimultaneousPhotosynthesis {
	
	&Photosynthesis;
		$W_retention = ($W-$W_wilting)/($W_capacity-$W_wilting);
	if ($vegetation_type == 0) {
		$f_w = 10*($W-$W_wilting)/3/($W_capacity-$W_wilting);
	} elsif ($vegetation_type == 1) {
		$f_w = ($W-$W_wilting)/($W_capacity-$W_wilting);
	}
	if ($f_w < 0) {
		$f_w = 0;
	}
	$f_w = $h_beta;

	if ($R_s_total>0 && $S_max>0) {	# analytical solution of photosynthesis rate and incellular CO2 concenteration

		if ($R_s_b_total>0 && $A_c_sunlit >= 0) {

			$Gamma_respiration_sunlit = ($R_d_sunlit*($K_c_sunlit*$pressure*1000)*(1+$O_i/$K_o_sunlit)+$V_c_max_sunlit*$Gamma_sunlit*($pressure/1000))/($V_c_max_sunlit-$R_d_sunlit)/($pressure/1000);
			$d_1 = $a_1*$f_w/($C_s_sunlit-$Gamma_respiration_sunlit)/(1+$VPD_s_sunlit/$VPD_0);
			if ($V_c_sunlit < $V_j_sunlit) {
				$d_2 = $V_c_max_sunlit;
				$d_3 = $K_c_sunlit*(1+$O_i/$K_o_sunlit);
			} else {
				$d_2 = 0.25*$J_sunlit;
				$d_3 = 2*$Gamma_sunlit;
			}
			$b_0 = -(1-$d_1*$C_s_sunlit)*($d_2*$Gamma_sunlit+$d_3*$R_d_sunlit)-$G_0_sunlit*$d_3*$C_s_sunlit;
			$b_1 = (1-$d_1*$C_s_sunlit)*($d_2-$R_d_sunlit)+$G_0_sunlit*($d_3-$C_s_sunlit)-$d_1*($d_2*$Gamma_sunlit+$d_3*$R_d_sunlit);
			$b_2 = $G_0_sunlit+$d_1*($d_2-$R_d_sunlit);
			if ($b_1**2>4*$b_0*$b_2 && $b_2!=0) {
				$C_i_sunlit = (-$b_1+sqrt($b_1**2-4*$b_0*$b_2))/(2*$b_2);
			} else {
				$C_i_sunlit = $C_a_out*0.7;
			}
			if (($C_i_sunlit eq 'nan') or ($C_i_sunlit+$d_3==0)) {
				$C_i_sunlit = $C_a_out*0.7;
			} else {
				$A_c_sunlit = $d_2*($C_i_sunlit-$Gamma_sunlit)/($C_i_sunlit+$d_3)-$R_d_sunlit;
			}

		}
	
		if ($R_s_d_total>0 && $A_c_shaded >= 0) {
	
			$Gamma_respiration_shaded = ($R_d_shaded*($K_c_shaded*$pressure*1000)*(1+$O_i/$K_o_shaded)+$V_c_max_shaded*$Gamma_shaded*($pressure/1000))/($V_c_max_shaded-$R_d_shaded)/($pressure/1000);
			$d_1 = $a_1*$f_w/($C_s_shaded-$Gamma_respiration_shaded)/(1+$VPD_s_shaded/$VPD_0);
			if ($V_c_shaded < $V_j_shaded) {
				$d_2 = $V_c_max_shaded;
				$d_3 = $K_c_shaded*(1+$O_i/$K_o_shaded);
			} else {
				$d_2 = 0.25*$J_shaded;
				$d_3 = 2*$Gamma_shaded;
			}
			$b_0 = -(1-$d_1*$C_s_shaded)*($d_2*$Gamma_shaded+$d_3*$R_d_shaded)-$G_0_shaded*$d_3*$C_s_shaded;
			$b_1 = (1-$d_1*$C_s_shaded)*($d_2-$R_d_shaded)+$G_0_shaded*($d_3-$C_s_shaded)-$d_1*($d_2*$Gamma_shaded+$d_3*$R_d_shaded);
			$b_2 = $G_0_shaded+$d_1*($d_2-$R_d_shaded);
			if ($b_1**2>4*$b_0*$b_2 && $b_2!=0) {
				$C_i_shaded = (-$b_1+sqrt($b_1**2-4*$b_0*$b_2))/(2*$b_2);
			} else {
				$C_i_shaded = $C_a_out*0.7;
			}
			if (($C_i_shaded eq 'nan') or ($C_i_shaded+$d_3==0)) {
				$C_i_shaded = $C_a_out*0.7;
			} else  {
				$A_c_shaded = $d_2*($C_i_shaded-$Gamma_shaded)/($C_i_shaded+$d_3)-$R_d_sunlit;
			}

		}
	
	}
	
}

sub Photosynthesis {	# calculate photosynthesis rate using Farquhar approach

	$L = 0;	# PAR

	$K_c_sunlit = 260*exp(($T_c_C_sunlit-25)*59366/(298*$R*$T_c_K_sunlit))/1000/$pressure;	# mol/mol
	$K_c_shaded = 260*exp(($T_c_C_shaded-25)*59366/(298*$R*$T_c_K_shaded))/1000/$pressure;	# mol/mol
	$K_o_sunlit = 179*exp(($T_c_C_sunlit-25)*35948/(298*$R*$T_c_K_sunlit))/$pressure;	# mol/mol
	$K_o_shaded = 179*exp(($T_c_C_shaded-25)*35948/(298*$R*$T_c_K_shaded))/$pressure;	# mol/mol
	$Gamma_sunlit = $Gamma_25[$vegetation_type][$leaf_type]*exp(($T_c_C_sunlit-25)*29000/(298*$R*$T_c_K_sunlit))/($pressure/1000);	# CO2 compensation point without dark respiration (μmol/mol)
	$Gamma_shaded = $Gamma_25[$vegetation_type][$leaf_type]*exp(($T_c_C_shaded-25)*29000/(298*$R*$T_c_K_shaded))/($pressure/1000);	# CO2 compensation point without dark respiration (μmol/mol)
	$O_i = 209/$pressure;	# mol/mol
	
	&CalVc;

	if ($R_s_total>0 && $S_max>0) {

		&CalVj;


		if ($V_c_sunlit < $V_j_sunlit) {
			$V_n_sunlit = $V_c_sunlit;
		} else {
			$V_n_sunlit = $V_j_sunlit;
		}
		if ($V_c_shaded < $V_j_shaded) {
			$V_n_shaded = $V_c_shaded;
		} else {
			$V_n_shaded = $V_j_shaded;
		}
	} else {
		$V_n_sunlit = 0;
		$V_n_shaded = 0;
	}

	$R_d_sunlit = 0.015*$V_c_max_sunlit;
	$R_d_shaded = 0.015*$V_c_max_shaded;
	$A_c_sunlit = $V_n_sunlit-$R_d_sunlit;
	$A_c_shaded = $V_n_shaded-$R_d_shaded;
}

sub CalVc {

	$V_c_max_top_sunlit = $V_c_max_org[$vegetation_type][$leaf_type]*exp(($T_c_C_sunlit-25)*$Ha_V[$vegetation_type][$leaf_type]/(298*$R*$T_c_K_sunlit));	# μmol/m^2/s
	$V_c_max_top_shaded = $V_c_max_org[$vegetation_type][$leaf_type]*exp(($T_c_C_shaded-25)*$Ha_V[$vegetation_type][$leaf_type]/(298*$R*$T_c_K_shaded));	# μmol/m^2/s
	$k_n = 1.5*exp(0.00963*($V_c_max_top_sunlit*$LAI_sunlit+$V_c_max_top_shaded*$LAI_shaded)/$LAI-2.43);

	&ExtinctionCoeffcientFunction($LAI,$k_b_black+$k_n);
	$extinction_coeffcient_b_black_plus_n = $extinction_coeffcient;
	$V_c_max_sunlit = $V_c_max_top_sunlit*$extinction_coeffcient_b_black_plus_n;

	&ExtinctionCoeffcientFunction($LAI,$k_n);
	$extinction_coeffcient_n = $extinction_coeffcient;
	$V_c_max_shaded = $V_c_max_top_shaded*($extinction_coeffcient_n-$extinction_coeffcient_b_black_plus_n);

	$V_c_sunlit = $V_c_max_sunlit*($C_i_sunlit-$Gamma_sunlit)/($C_i_sunlit+$K_c_sunlit*(1+$O_i/$K_o_sunlit))*$LAI_g/$LAI;
	$V_c_shaded = $V_c_max_shaded*($C_i_shaded-$Gamma_shaded)/($C_i_shaded+$K_c_shaded*(1+$O_i/$K_o_shaded))*$LAI_g/$LAI;

}

sub CalVj {
	
	$q_sunlit = $Q_sunlit[0]/(1-$transmissivity[$L])/(1-$transmissivity[$L]-$reflectance[$L])*4.6;	# 4.6 μmol/J

	$J_max_sunlit = $J_m_25[$vegetation_type][$leaf_type]*exp(($T_c_C_sunlit-25)*$Ha_J[$vegetation_type][$leaf_type]/($R*$T_c_K_sunlit*298))*(1+exp(($Delta_S[$leaf_type]*298-$Hd_J[$vegetation_type][$leaf_type])/($R*298)))/(1+exp(($Delta_S[$leaf_type]*$T_c_K_sunlit-$Hd_J[$vegetation_type][$leaf_type])/($R*$T_c_K_sunlit)));
	$J_max_shaded = $J_m_25[$vegetation_type][$leaf_type]*exp(($T_c_C_shaded-25)*$Ha_J[$vegetation_type][$leaf_type]/($R*$T_c_K_shaded*298))*(1+exp(($Delta_S[$leaf_type]*298-$Hd_J[$vegetation_type][$leaf_type])/($R*298)))/(1+exp(($Delta_S[$leaf_type]*$T_c_K_shaded-$Hd_J[$vegetation_type][$leaf_type])/($R*$T_c_K_shaded)));

	$J_max_sunlit *= $extinction_coeffcient_b_black_plus_n;
	$J_max_shaded *= $extinction_coeffcient_n-$extinction_coeffcient_b_black_plus_n;

	if ($leaf_type>0 && $T_a_C_mean<$T_cold[$vegetation_type][$leaf_type]) {
		$J_max_sunlit /= 10;
		$J_max_shaded /= 10;
	}

	$APAR_su = $q_sunlit*$absorptance*(1-$f)/2;

	$b_su = -($APAR_su+$J_max_sunlit);
	$c_su = $APAR_su*$J_max_sunlit;
	$J_sunlit = (-$b_su-sqrt($b_su**2-4*$a*$c_su))/(2*$a);
	if ($J_sunlit < 0) {$J_sunlit = 0;}
	
	&ExtinctionCoeffcientFunction($LAI,$k_d);
	$extinction_coeffcient_d = $extinction_coeffcient;
	&ExtinctionCoeffcientFunction($LAI,$k_b_black+$k_d);
	$extinction_coeffcient_b_black_plus_d = $extinction_coeffcient;
	$j_max_top = $J_max_shaded/($extinction_coeffcient_d-$extinction_coeffcient_b_black_plus_d);

	&ExtinctionCoeffcientFunction($LAI,$k_b);
	$extinction_coeffcient_b = $extinction_coeffcient;
	&ExtinctionCoeffcientFunction($LAI,$k_b+$k_b_black);
	$extinction_coeffcient_b_plus_b_black = $extinction_coeffcient;
	&ExtinctionCoeffcientFunction($LAI,$k_d+$k_b_black);
	$extinction_coeffcient_d_plus_b_black = $extinction_coeffcient;
	&ExtinctionCoeffcientFunction($LAI,$k_b_black);
	$extinction_coeffcient_b_black = $extinction_coeffcient;
	&ExtinctionCoeffcientFunction($LAI,2*$k_b_black);
	$extinction_coeffcient_2b_black = $extinction_coeffcient;
	$c_2 = ($R_s_b[$L]*((1-$reflectance_b)*$k_b*($extinction_coeffcient_b-$extinction_coeffcient_b_plus_b_black))-(1-$scattering_coefficient[$L])*$k_b_black*($extinction_coeffcient_b_black-$extinction_coeffcient_2b_black))/($R_s_d[$L]*(1-$reflectance_d)*$k_d*($extinction_coeffcient_d-$extinction_coeffcient_d_plus_b_black));
	$q_shaded_top = (1+$c_2)*$R_s_d[$L]*$k_d*(1-$reflectance_d)*4.2;	# 4.2 μmol/J

	$APAR_sh = $q_shaded_top*$absorptance*(1-$f)/2;
	$b_sh = -($APAR_sh+$j_max_top);
	$c_sh = $APAR_sh*$j_max_top;
	$J_shaded = ($extinction_coeffcient_d-$extinction_coeffcient_d_plus_b_black)*(-$b_sh-sqrt($b_sh**2-4*$a*$c_sh))/(2*$a);
	if ($J_shaded < 0) {$J_shaded = 0;}
	
	$V_j_sunlit = ($C_i_sunlit-$Gamma_sunlit)*$J_sunlit/(4*($C_i_sunlit+2*$Gamma_sunlit))*$LAI_g/$LAI;
	
	if ($C_i_shaded+2*$Gamma_shaded==0) {
		$V_j_shaded = 0;
	} else {
	$V_j_shaded = ($C_i_shaded-$Gamma_shaded)*$J_shaded/(4*($C_i_shaded+2*$Gamma_shaded))*$LAI_g/$LAI;
	}

}

sub CalET_SW {	# calculate evapotranspiration rate using Shuttleworth-Wallace approach
	
			$lambda_mol_sunlit = 56780.3-42.84*$T_a_K;	# latent heat of vaporization (J/kg)
			$lambda_mol_shaded = 56780.3-42.84*$T_a_K;	# latent heat of vaporization (J/kg)
			$lambda_kg = (56780.3-42.84*$T_a_K)/0.018;	# latent heat of vaporization (J/kg)
	$Delta_soil = $a0*$b0*$c0/($T_c_C_soil+$b0)**2*exp($a0*$T_c_C_soil/($T_c_C_soil+$b0));	# slope of the saturated vapor pressure (hPa/°C)
	if ($LAI > 0) {
		$Q_n_sunlit = $Q_n_isothermal_sunlit-$g_r_sunlit*$c_p*($T_c_K_sunlit-$T_a_K);
		$Q_n_shaded = $Q_n_isothermal_shaded-$g_r_shaded*$c_p*($T_c_K_shaded-$T_a_K);
		if ($S_max>0 && $R_s_total>0) {
			$k_r = sqrt(0.5)*0.5/$COSTHETA;	# $k_r: canopy extinction coefficient of net raditation
			if ($k_r*$LAI > 1) {
				$k_r_exp = exp(-$k_r*$LAI);
				if ($k_r_exp < 0.01) {
					$k_r_exp = 0.01;
				} elsif ($k_r_exp >= 1) {
					$k_r_exp = 0.99;
				}
				$R_n_soil = ($Q_n_sunlit+$Q_n_shaded)/(1/$k_r_exp-1);
				$R_n_sum = $Q_n_sunlit+$Q_n_shaded+$R_n_soil;
				if (abs($R_n_sum) > abs($R_s_total*(1-$albedo_soil[1])+$L_down-$emissivity_c*$SteBol_const*$T_a_K**4)*1.5) {
					$R_n_sum = $R_s_total*(1-$albedo_soil[1])+$L_down-$emissivity_c*$SteBol_const*$T_a_K**4;
					$R_n_soil = $R_n_sum-($Q_n_sunlit+$Q_n_shaded);
				}
			} elsif ($LAI < 1) {
				$R_n_sum = ($Q_n_sunlit+$Q_n_shaded)/$LAI;
				$R_n_soil = $R_n_sum-($Q_n_sunlit+$Q_n_shaded);
			} else {
				$R_n_sum = $R_s_total*(1-$albedo_soil[1])+$L_down-$emissivity_c*$SteBol_const*$T_a_K**4;
				$R_n_soil = $R_n_sum-($Q_n_sunlit+$Q_n_shaded);
			}
			if ($LAI <= 4) {
				$G = $R_n_sum*0.036;
			} else {
				$G = $R_n_sum*0.036;
			}
		} else {
			$R_n_soil = ($L_down-$emissivity_c*$SteBol_const*$T_a_K**4)*exp(-$extinction_coef*$LAI);
			$G = $R_n_sum*0.036;
			$R_n_sum = $L_down-$emissivity_c*$SteBol_const*$T_a_K**4;
		}
	} else {
		$R_n_soil = $R_s_total*(1-$albedo_soil[1])+$L_down-$emissivity_s*$SteBol_const*$T_a_K**4;
		$G = $R_n_soil*0.036;
		$R_n_sum = $R_n_soil;
	}

	$rr_a = ($Delta+$gamma)*$r_a;
	$rr_s = ($Delta+$gamma)*$r_a_soil+$gamma*$r_s;

	if ($LAI > 0) {
		$PM_c = ($Delta*($R_n_sum-$G)+($rho_a*$C_p*$VPD_a-$Delta*$r_a_c*($R_n_soil-$G))/($r_a+$r_a_c))/($Delta+$gamma*(1+$r_s_c/($r_a+$r_a_c)))*$LAI_g/$LAI;
		$rr_c = ($Delta+$gamma)*$r_a_c+$gamma*$r_s_c;
		$Cc = 1/(1+$rr_c*$rr_a/$rr_s/($rr_c+$rr_a));
		$Cs = 1/(1+$rr_s*$rr_a/$rr_c/($rr_s+$rr_a));
	} else {
		$PM_c = 0;
		$Cc = 0;
		$Cs = 1;
	}
	$PM_s = ($Delta*($R_n_sum-$G)+($rho_a*$C_p*$VPD_a-$Delta*$r_a_soil*($R_n_sum-$R_n_soil))/($r_a+$r_a_soil))/($Delta+$gamma*(1+$r_s/($r_a+$r_a_soil)));

	if ($S_max<=0 or $R_s_total<=0) {
		$PM_c = 0;
		$PM_s = 0;
	}
	$lE_c = $Cc*$PM_c;
	if ($lE_c < 0) {
		$lE_c = 0;
	}
	$lE_soil = $Cs*$PM_s;
	if ($lE_soil< 0) {
		$lE_soil = 0;
	}
	$lE_sum = $lE_c+$lE_soil;

	if ($lE_sum > abs($R_n_sum*2)) {
		$lE_sum = $R_n_sum;
		$lE_c = $lE_sum;
		$lE_soil = 0;
	}
	$ET_c_sum = $lE_c/$lambda_kg;
	$ET_soil = $lE_soil/$lambda_kg;
	$ET_sum = $lE_sum/$lambda_kg;
	
	$H_sum = $R_n_sum-$lE_sum-$G;

}

sub InitialTemp {

	$T_c_K[$N_layer+1] = $T_a_K[$N_layer+1];
		if ($solar_elevation > 0) {
			if ($LAI >= 2) {
				$T_a_K[0] = $T_a_K[$N_layer+1]-2;
			} elsif ($LAI >= 1) {
				$T_a_K[0] = $T_a_K[$N_layer+1]-1;	
			} else {
				$T_a_K[0] = $T_a_K[$N_layer+1];	
			}
		} else {
			if ($LAI >= 2) {
				$T_a_K[0] = $T_a_K[$N_layer+1];
			} elsif ($LAI >= 1) {
				$T_a_K[0] = $T_a_K[$N_layer+1];	
			} else {
				$T_a_K[0] = $T_a_K[$N_layer+1];	
			}
		}
	$T_c_K[0] =	$T_a_K[0];
	$L_up[0] = $emissivity_s*$SteBol_const*$T_c_K[0]**4;

	for ($i_l=$N_layer; $i_l>=1; $i_l--) {
		$T_c_K[$i_l] = $T_a_K[$N_layer+1]+($T_c_K[0]-$T_c_K[$N_layer+1])/$N_layer*($N_layer-$i_l)+1;
		if ($solar_elevation > 0) {
			if ($LAI >= 2) {
				$T_a_K[$i_l] = $T_a_K[$N_layer+1];
			} elsif ($LAI >= 1) {
				$T_a_K[$i_l] = $T_a_K[$N_layer+1];
			} else {
				$T_a_K[$i_l] = $T_a_K[$N_layer+1];
			}
		} else {
			if ($LAI >= 2) {
				$T_a_K[$i_l] = $T_a_K[$N_layer+1];
			} elsif ($LAI >= 1) {
				$T_a_K[$i_l] = $T_a_K[$N_layer+1];
			} else {
				$T_a_K[$i_l] = $T_a_K[$N_layer+1];
			}
		}
	}
	for ($i_l=$N_layer; $i_l>=1; $i_l--) {
		$T_c_C[$i_l] = $T_c_K[$i_l]-273.15;
		$T_a_C[$i_l] = $T_a_K[$i_l]-273.15;
		$L_down[$i_l-1] = $L_down[$i_l]*$I_d[$i_l]+$emissivity_c*$SteBol_const*$T_c_K[$i_l]**4*(1-$I_d[$i_l]);
		$L_up[$N_layer-($i_l-1)] = $L_up[$N_layer-1-($i_l-1)]*$I_d[$N_layer-($i_l-1)]+$emissivity_c*$SteBol_const*$T_c_K[$N_layer-($i_l-1)]**4*(1-$I_d[$N_layer-($i_l-1)]);
	}

	$T_c_C[0] = $T_c_K[0]-273.15;
	$T_a_C[0] = $T_a_K[0]-273.15;

}

sub WaterBalance {

	$W += $rainfall;
	if ($W > $W_capacity) {
		$W = $W_capacity;
	}

	$W += -$ET_sum*60*$t_step;
	if ($W < 0) {
		$ET_sum = $W_pre;
		$W = 0;
	}
	if ($W > $W_capacity) {
		$W = $W_capacity;
	}
		$W_retention = ($W-$W_wilting)/($W_capacity-$W_wilting);
	if ($W_retention < 0) {
		$beta_water = 0;
	} else {
		$beta_water = $W_retention;
	}
	$G_water = 1-(1-$beta_water)**2;
	if ($G_water < 0.001) {
		$G_water = 0.001;
	}
	$SW_retention = 1;
	
}

sub NoConverged {
	
			$g_s_sunlit = 1;
			$g_s_shaded = 0.5;
			$C_s_sunlit = $C_a_out*0.8;
			$C_s_shaded = $C_a_out*0.8;
			$C_i_sunlit = $C_a_out*0.7;
			$C_i_shaded = $C_a_out*0.7;
			$diff_T_sunlit = 0;
			$diff_T_shaded = 0;
			$T_c_C_soil = $T_a_C;
			$T_c_K_soil = $T_c_C_soil+273.15;
			&L_balance;
			&CalRn;
			if ($k_r*$LAI > 1) {
				$k_r_exp = exp(-$k_r*$LAI);
				if ($k_r_exp < 0.01) {
					$k_r_exp = 0.01;
				} elsif ($k_r_exp >= 1) {
					$k_r_exp = 0.99;
				}
				$R_n_soil = ($Q_n_sunlit+$Q_n_shaded)/(1/$k_r_exp-1);
				$R_n_sum = $Q_n_sunlit+$Q_n_shaded+$R_n_soil;
				if (abs($R_n_sum) > abs($R_s_total*(1-$albedo_soil[1])+$L_down-$emissivity_c*$SteBol_const*$T_a_K**4)*1.5) {
					$R_n_sum = $R_s_total*(1-$albedo_soil[1])+$L_down-$emissivity_c*$SteBol_const*$T_a_K**4;
					$R_n_soil = $R_n_sum-($Q_n_sunlit+$Q_n_shaded);
				}
			} elsif ($LAI < 1) {
				$R_n_sum = ($Q_n_sunlit+$Q_n_shaded)/$LAI;
				$R_n_soil = $R_n_sum-($Q_n_sunlit+$Q_n_shaded);
			} else {
				$R_n_sum = $R_s_total*(1-$albedo_soil[1])+$L_down-$emissivity_c*$SteBol_const*$T_a_K**4;
				$R_n_soil = $R_n_sum-($Q_n_sunlit+$Q_n_shaded);
			}
			$G = $R_n_sum*0.036;
		$W_retention = ($W-$W_wilting)/($W_capacity-$W_wilting);
	if ($W_retention < 0) {
		$beta_water = 0;
	} else {
		$beta_water = $W_retention;
	}
	$G_water = 1-(1-$beta_water)**2;
	if ($G_water < 0.001) {
		$G_water = 0.001;
	}
			$C_i_sunlit = $C_a_out*0.7;
			$C_i_shaded = $C_a_out*0.7;
			$K_c_sunlit = 260*exp(($T_c_C_sunlit-25)*59366/(298*$R*$T_c_K_sunlit))/1000/$pressure;	# mol/mol
			$K_c_shaded = 260*exp(($T_c_C_shaded-25)*59366/(298*$R*$T_c_K_shaded))/1000/$pressure;	# mol/mol
			$K_o_sunlit = 179*exp(($T_c_C_sunlit-25)*35948/(298*$R*$T_c_K_sunlit))/$pressure;	# mol/mol
			$K_o_shaded = 179*exp(($T_c_C_shaded-25)*35948/(298*$R*$T_c_K_shaded))/$pressure;	# mol/mol
			$Gamma_sunlit = $Gamma_25[$vegetation_type][$leaf_type]*exp(($T_c_C_sunlit-25)*29000/(298*$R*$T_c_K_sunlit))/($pressure/1000);	# CO2 compensation point without dark respiration (μmol/mol)
			$Gamma_shaded = $Gamma_25[$vegetation_type][$leaf_type]*exp(($T_c_C_shaded-25)*29000/(298*$R*$T_c_K_shaded))/($pressure/1000);	# CO2 compensation point without dark respiration (μmol/mol)
			$O_i = 209/$pressure;	# mol/mol
			$V_c_max_top_sunlit = $V_c_max_org[$vegetation_type][$leaf_type]*exp(($T_c_C_sunlit-25)*$Ha_V[$vegetation_type][$leaf_type]/(298*$R*$T_c_K_sunlit));	# μmol/m^2/s
			$V_c_max_top_shaded = $V_c_max_org[$vegetation_type][$leaf_type]*exp(($T_c_C_shaded-25)*$Ha_V[$vegetation_type][$leaf_type]/(298*$R*$T_c_K_shaded));	# μmol/m^2/s
			$V_c_max_sunlit = $V_c_max_top_sunlit*$extinction_coeffcient_b_black_plus_n;
			$V_c_max_shaded = $V_c_max_top_shaded*($extinction_coeffcient_n-$extinction_coeffcient_b_black_plus_n);
			$V_c_sunlit = $V_c_max_sunlit*($C_i_sunlit-$Gamma_sunlit)/($C_i_sunlit+$K_c_sunlit*(1+$O_i/$K_o_sunlit));
			$V_c_shaded = $V_c_max_shaded*($C_i_shaded-$Gamma_shaded)/($C_i_shaded+$K_c_shaded*(1+$O_i/$K_o_shaded));
			$J_max_sunlit = $J_m_25[$vegetation_type][$leaf_type]*exp(($T_c_C_sunlit-25)*$Ha_J[$vegetation_type][$leaf_type]/($R*$T_c_K_sunlit*298))*(1+exp(($Delta_S[$leaf_type]*298-$Hd_J[$vegetation_type][$leaf_type])/($R*298)))/(1+exp(($Delta_S[$leaf_type]*$T_c_K_sunlit-$Hd_J[$vegetation_type][$leaf_type])/($R*$T_c_K_sunlit)));
			$J_max_shaded = $J_m_25[$vegetation_type][$leaf_type]*exp(($T_c_C_shaded-25)*$Ha_J[$vegetation_type][$leaf_type]/($R*$T_c_K_shaded*298))*(1+exp(($Delta_S[$leaf_type]*298-$Hd_J[$vegetation_type][$leaf_type])/($R*298)))/(1+exp(($Delta_S[$leaf_type]*$T_c_K_shaded-$Hd_J[$vegetation_type][$leaf_type])/($R*$T_c_K_shaded)));
			$J_max_sunlit *= $extinction_coeffcient_b_black_plus_n;
			$J_max_shaded *= $extinction_coeffcient_n-$extinction_coeffcient_b_black_plus_n;
			$APAR_su = $q_sunlit*$absorptance*(1-$f)/2;
			$b_su = -($APAR_su+$J_max_sunlit);
			$c_su = $APAR_su*$J_max_sunlit;
			$J_sunlit = (-$b_su-sqrt($b_su**2-4*$a*$c_su))/(2*$a);
			if ($J_sunlit < 0) {$J_sunlit = 0;}
			$q_shaded = $Q_sunlit[0]/(1-$transmissivity[$L])/(1-$transmissivity[$L]-$reflectance[$L])*4.6;
			$APAR_sh = $q_shaded*$absorptance*(1-$f)/2;
			$b_sh = -($APAR_sh+$J_max_shaded);
			$c_sh = $APAR_sh*$J_max_shaded;
			$J_shaded = (-$b_sh-sqrt($b_sh**2-4*$a*$c_sh))/(2*$a);
			if ($J_shaded < 0) {$J_shaded = 0;}
			$V_j_sunlit = ($C_i_sunlit-$Gamma_sunlit)*$J_sunlit/(4*($C_i_sunlit+2*$Gamma_sunlit));
			$V_j_shaded = ($C_i_shaded-$Gamma_shaded)*$J_shaded/(4*($C_i_shaded+2*$Gamma_shaded));
			if ($V_c_sunlit < $V_j_sunlit) {
				$V_n_sunlit = $V_c_sunlit;
			} else {
				$V_n_sunlit = $V_j_sunlit;
			}
			if ($V_c_shaded < $V_j_shaded) {
				$V_n_shaded = $V_c_shaded;
			} else {
				$V_n_shaded = $V_j_shaded;
			}
			$R_d_sunlit = 0.015*$V_c_max_sunlit;
			$R_d_shaded = 0.015*$V_c_max_shaded;
			$A_c_sunlit = $V_n_sunlit-$R_d_sunlit;
			$A_c_shaded = $V_n_shaded-$R_d_shaded;
			$g_s_sunlit = $G_0_sunlit+$a_1*$f_w*$A_c_sunlit/($C_s_sunlit-$Gamma_respiration_sunlit)/(1+$VPD_s_sunlit/$VPD_0);	# mol m^2/s
			$g_s_shaded = $G_0_shaded+$a_1*$f_w*$A_c_shaded/($C_s_shaded-$Gamma_respiration_shaded)/(1+$VPD_s_shaded/$VPD_0);	# mol m^2/s
			$r_s_c = 2/($g_s_sunlit+$g_s_shaded)*$pressure*100/$R/$T_a_K;	# (s/m)
			$lE_c = ($Delta*($R_n_sum-$R_n_soil)+($rho_a*$C_p*$VPD_a)/($r_a+$r_a_c))/($Delta+$gamma*(1+$r_s_c/($r_a+$r_a_c)))*$LAI_g/$LAI;
			$lE_soil = ($Delta*($R_n_soil-$G)+($rho_a*$C_p*$VPD_a)/($r_a+$r_a_soil))/($Delta+$gamma*(1+$r_s/($r_a+$r_a_soil)))*$h_beta;
			if (abs($lE_sum) > abs($R_n_sum*1.5)) {
				$lE_sum = $R_n_sum;
				$lE_c = $lE_sum;
				$lE_soil = 0;
			}
			$ET_c_sum = $lE_c/$lambda_kg;
			$ET_soil = $lE_soil/$lambda_kg;
			$ET_sum = $lE_sum/$lambda_kg;
			$H_sum = $R_n_sum-$lE_sum-$G;
	
}

sub CalLAI {

	if ($leaf_dormant == 1) {
		$dormancy_dys++;
	}
	if ($dormancy_dys == $dormancy_terminate_day) {
		$leaf_dormant = 0;
	}

	if ($virtual_LAI==1) {
		$A_sum_daily *= 12/1000000000;	# (μmol/m^2/dy–>kgC/m^2/dy)
		$A_n_sum_daily *= 12/1000000000;	# (μmol/m^2/dy–>kgC/m^2/dy)
		$R_m_leaf_sum_daily *= 12/1000000000;	# (μmol/m^2/dy–>kgC/m^2/dy)
		$virtual_A_n = $A_sum_daily-$R_m_leaf_sum_daily;
		print $LAI, "\t", $A_sum_daily, "\n";
		$A_sum_daily /= 12/1000000000;	# (kgC/m^2/dy–>μmol/m^2/dy)
		$A_n_sum_daily /= 12/1000000000;	# (kgC/m^2/dy–>μmol/m^2/dy)
		$R_m_leaf_sum_daily /= 12/1000000000;	# (kgC/m^2/dy–>μmol/m^2/dy)


				if ($virtual_A_n > 0 && $T_a_C_mean>=$T_cold[$vegetation_type][$leaf_type] && (($leaf_type==0 && $sun_duration_pos>=11) or $vegetation_type==1) && $W>$W_wilting) {
					$virtual_LAI_day++;
					if ($virtual_LAI_day==$virtual_LAI_succeed_day) {
						$virtual_LAI = 10;
					}
					print "~~~\n", $virtual_A_n, "\n", "~~~\n";
				} else {
					$virtual_LAI_day++;
					print "~\n", $virtual_A_n, "\n", "~\n";
					$virtual_LAI = -9;
				}
	}
			if ($virtual_LAI == -9) {
				$y2 -= 24*$virtual_LAI_day;
				$day -= $virtual_LAI_day;
				if ($y2 < -1) {
					$y2 += 24*$DOY_max;
				}
				if ($day < 0) {
					$month--;
					if ($pass200==0 && $month==0) {
						$month = '01';
						$day = 0;
					} else {
						if ($month==0) {$month=12;}
						if ($month==1) {$month='01';}
						if ($month==2) {$month='02';}
						if ($month==3) {$month='03';}
						if ($month==4) {$month='04';}
						if ($month==5) {$month='05';}
						if ($month==6) {$month='06';}
						if ($month==7) {$month='07';}
						if ($month==8) {$month='08';}
						if ($month==9) {$month='09';}
						$day = $dy_month[$month-1]+$day;
					}
				}
				$virtual_LAI_day = 0;
				next;
			} else {
				&AllocationModel;
				if ($leaf_type==0 && $growing_dys>=$virtual_LAI_succeed_day) {
					$growing_dys++;
				}
			}
		if ($virtual_LAI==1 or $virtual_LAI==10) {
					$virtual_DOY[$virtual_LAI_day-1] = $DOY;
					$virtual_T_a_C_mean[$virtual_LAI_day-1] = $T_a_C_mean;
					$virtual_R_s_sum[$virtual_LAI_day-1] = $R_s_sum;
					$virtual_pressure_mean[$virtual_LAI_day-1] = $pressure_mean;
					$virtual_rainfall_sum[$virtual_LAI_day-1] = $rainfall_sum;
					$virtual_rh_mean[$virtual_LAI_day-1] = $rh_mean;
					$virtual_u_z_mean[$virtual_LAI_day-1] = $u_z_mean;
					$virtual_W_mean[$virtual_LAI_day-1] = $W_mean;
					$virtual_A_sum_daily[$virtual_LAI_day-1] = $A_sum_daily;
					$virtual_R_a_c[$virtual_LAI_day-1] = $R_a_c;
					$virtual_R_m_leaf_sum_daily[$virtual_LAI_day-1] = $R_m_leaf_sum_daily;
					$virtual_ET_daily[$virtual_LAI_day-1] = $ET_daily;
					$virtual_ET_c_daily[$virtual_LAI_day-1] = $ET_c_daily;
					$virtual_ET_eq_daily[$virtual_LAI_day-1] = $ET_eq_daily;
					$virtual_LAI[$virtual_LAI_day-1] = $LAI;
					$virtual_LAI_g[$virtual_LAI_day-1] = $LAI_g;
					$virtual_C_leaf[$virtual_LAI_day-1] = $C_leaf;
					$virtual_Cg_leaf[$virtual_LAI_day-1] = $Cg_leaf;
					$virtual_Cd_leaf[$virtual_LAI_day-1] = $Cd_leaf;
					$virtual_C_stem[$virtual_LAI_day-1] = $C_stem;
					$virtual_C_root[$virtual_LAI_day-1] = $C_root;
					$virtual_C_all[$virtual_LAI_day-1] = $C_all;
		}

}

sub AllocationModel {
	
		$A_root_daily = 0;
		$A_stem_daily = 0;
		$A_sum_daily *= 12/1000000000;	# (μmol/m^2/dy–>kgC/m^2/dy)
		$A_n_sum_daily *= 12/1000000000;	# (μmol/m^2/dy–>kgC/m^2/dy)
		$R_m_leaf_sum_daily *= 12/1000000000;	# (μmol/m^2/dy–>kgC/m^2/dy)
		$Q_10 = 3.22-0.046*$T_a_C_mean;	# $Q_10: Q10 value for maintenance respiration
		$f_20 = $Q_10**(($T_a_C_mean-20)/10);	# $f_20: temperature dependent function
		$f_15 = $Q_10**(($T_a_C_mean-15)/10);	# $f_15: temperature dependent function
		$lf_stem = exp(-0.2835*$C_stem);
		if ($lf_stem > 1) {$lf_stem = 1;}
		if ($lf_stem < 0.05) {$lf_stem = 0.05;}
		$lf_root = exp(-0.2835*$C_root);
		if ($lf_root > 1) {$lf_root = 1;}
		if ($lf_root < 0.05) {$lf_root = 0.05;}
		$R_m_stem_daily = $R_m_base_stem[$leaf_type]*$lf_stem*$C_stem*$f_20;
		$R_m_root_daily = $R_m_base_root[$vegetation_type][$leaf_type]*$lf_root*$C_root*$f_20;
		$R_m_c = $R_m_leaf_sum_daily+$R_m_stem_daily+$R_m_root_daily;	# $R_m: maintenance respiration (kgC/m^2/dy)
		$light_availability = exp(-$k_n*$LAI);	# for trees and crops
		$water_availability = $W_retention;
		if ($vegetation_type==1 && $leaf_type==0 && (($LAI<1.5 && $phenophase==1) or $virtual_LAI==1)) {
				$A_root_frac = 0.8;
		} elsif (($vegetation_type!=1 or $leaf_type!=0) && (($Cg_leaf<(($C_stem+$C_root)/$allocation_parameter[$vegetation_type][$leaf_type])**(1/$allocation_parameter_index[$vegetation_type][$leaf_type])*0.45 && $phenophase==1) or $virtual_LAI==1)) {
				$A_root_frac = 0;
				$A_stem_frac = 0;
		} else {
			if ($phenophase == 1) {
				$leaf_onset = 0;
				$leaf_normal = 1;
				$phenophase = 2;
			}
			if ($leaf_type>0 && $phenophase==3 && $Cg_leaf<(($C_stem+$C_root)/$allocation_parameter[$vegetation_type][$leaf_type])**(1/$allocation_parameter_index[$vegetation_type][$leaf_type])*0.45) {
				$phenophase = 2;
			}
			if ($vegetation_type == 0) {
				$A_root_frac = ($r_00[$vegetation_type][$leaf_type]+$sensitivity_allocation[$vegetation_type][$leaf_type]*(1-$water_availability))/(1+$sensitivity_allocation[$vegetation_type][$leaf_type]*(2-$light_availability-$water_availability));
				$A_stem_frac = ($s_00[$leaf_type]+$sensitivity_allocation[$vegetation_type][$leaf_type]*(1-$light_availability))/(1+$sensitivity_allocation[$vegetation_type][$leaf_type]*(2-$light_availability-$water_availability));
			} elsif ($vegetation_type == 1) {
				$A_root_frac = ($r_00[$vegetation_type][$leaf_type]+$sensitivity_allocation[$vegetation_type][$leaf_type]*(1-$water_availability))/(1+$sensitivity_allocation[$vegetation_type][$leaf_type]*(1+$light_availability-$water_availability));
				$A_stem_frac = 0;
			}
			if ($Cg_leaf>(($C_stem+$C_root)/$allocation_parameter[$vegetation_type][$leaf_type])**(1/$allocation_parameter_index[$vegetation_type][$leaf_type])) {
				if ($phenophase < 3) {
						$phenophase = 3;
					$C_increase_dy = 0;
				}
			}
			if ($leaf_type==0 && $phenophase>0 && ($A_sum_daily-$R_m_leaf_sum_daily-$R_m_stem_daily-$R_m_root_daily<0 or ($A_sum_daily-$R_m_c>0 && $A_sum_daily-$R_m_leaf_sum_daily-$R_m_stem_daily-$R_m_root_daily-$R_g_parameter[$vegetation_type]*($A_sum_daily-$R_m_c)<0) or $T_a_C_mean<$T_cold[$vegetation_type][$leaf_type] or ($sun_duration_pos<11 && $vegetation_type==0))) {
					$C_increase_dy = 0;
			}
			if ($phenophase >= 3) {
				if ($vegetation_type == 0) {
					$A_root_frac = $A_root_frac/($A_root_frac+$A_stem_frac);
					$A_stem_frac = 1-$A_root_frac;
				} elsif ($vegetation_type == 1) {
					$A_root_frac = 0.33;
					$A_stem_frac = 0;
				}
			}
		}
		$A_leaf_frac = 1-($A_root_frac+$A_stem_frac);
		# R_g: growth respiration (kgC/m^2/dy)
		if ($A_sum_daily-$R_m_c > 0) {
			$R_g_c = $R_g_parameter[$vegetation_type]*($A_sum_daily-$R_m_c);
			# $R_g_[X]: growth respirtation for [X]
			$R_g_leaf_daily = $A_leaf_frac*$R_g_c;
			$R_g_stem_daily = $A_stem_frac*$R_g_c;
			$R_g_root_daily = $A_root_frac*$R_g_c;
		} else {
			$R_g_c = 0;
			$R_g_leaf_daily = 0;
			$R_g_stem_daily = 0;
			$R_g_root_daily = 0;
		}

	$R_a_c = $R_m_c+$R_g_c;	# $R_a_c: autotrophic respiration (kgC/m^2/dy)
		
	# $A_[X]: allocation to [X] (kgC/m^dy)
	if ($phenophase==1 && ($vegetation_type==0 or ($vegetation_type==1 && $leaf_type==1))) {
		$A_stem_daily = 0;
		$A_root_daily = 0;
	} elsif ($phenophase>=2 or ($vegetation_type==1 && $leaf_type==0 && $phenophase==1)) {
		if ($A_sum_daily-$R_a_c > 0) {
			$A_stem_daily = ($A_sum_daily-$R_a_c)*$A_stem_frac+$R_m_stem_daily+$R_g_stem_daily;
			$A_root_daily = ($A_sum_daily-$R_a_c)*$A_root_frac+$R_m_root_daily+$R_g_root_daliy;
		} else {
			$A_stem_daily = $A_sum_daily*$A_stem_frac;
			$A_root_daily = $A_sum_daily*$A_root_frac;
		}
	}

	if ($vegetation_type==0 && $leaf_type==2) {
		$T_retention = 1;
	} else {
		if ($T_a_C_mean > $T_cold[$vegetation_type][$leaf_type]) {
			$T_retention = 1;
		} elsif ($T_a_C_mean > $T_cold[$vegetation_type][$leaf_type]-5) {
			$T_retention = ($T_a_C_mean-($T_cold[$vegetation_type][$leaf_type]-5))/5;
		} else {
			$T_retention = 0;
		}
	}

	# $loss_rate_[X]: leaf loss rate because of [X] stress (/dy)
	if ($phenophase > 0) {
			$loss_rate_W = $loss_rate_W_max[$vegetation_type][$leaf_type]*(1-$W_retention)**$b_W[$vegetation_type][$leaf_type];
		$loss_rate_T = $loss_rate_T_max*(1-$T_retention)**$b_T;
	} else {
		$loss_rate_W = 0;
		$loss_rate_T = 0;
	}
	$S_leaf = $death_rate_leaf*$Cg_leaf;
	# $L_[X]: amount of C for [X] litter (kgC/m^2/dy)
	$L_leaf_d = $push_down_rate*$Cd_leaf;
	if ($phenophase>=2 && $LAI>0 && $virtual_LAI!=1) {
			$L_leaf = ($loss_rate+$loss_rate_W+$loss_rate_T)*$Cg_leaf;
		if ($phenophase==4 && $leat_type==0) {
			$L_leaf += (($C_stem+$C_root)/$allocation_parameter[$vegetation_type][$leaf_type])**(1/$allocation_parameter_index[$vegetation_type][$leaf_type])*(1-0.075)/30;
		}
		if ($vegetation_type==0 && $u_z_mean>=5) {
			if ($leaf_type ==0) {
				$L_leaf += 0.001*$u_z_mean**2*$Cg_leaf;
				$L_leaf_d += 0.003*$u_z_mean**2*$Cd_leaf;
			} elsif ($leaf_type==1 or $leaf_type==2) {
				$L_leaf += 0.0003*$u_z_mean**2*$Cg_leaf;
				$L_leaf_d += 0.001*$u_z_mean**2*$Cd_leaf;
			} elsif ($leaf_type == 3) {
				$L_leaf += 0.0003*$u_z_mean*$Cg_leaf;
				$L_leaf_d += 0.001*$u_z_mean*$Cd_leaf;
			}
		}
	} else {
		$L_leaf = 0;
	}

	if ($vegetation_type == 0) {
		$L_stem = $loss_rate_stem[$vegetation_type][$leaf_type]*$C_stem;
		if ($phenophase>=2 or $phenophase==0) {
			$L_stem += $loss_rate_standby_stem[$vegetation_type]*$C_stem;
		}
	} elsif ($vegetation_type == 1) {
		$L_stem = 0;
	}
		$L_root = $loss_rate_root[$vegetation_type][$leaf_type]*$C_root;
		if ($phenophase>=2 or $phenophase==0) {
			$L_root += $loss_rate_standby_root[$vegetation_type][$leaf_type]*$C_root;
		}
	$L_all = $L_leaf+$L_stem+$L_root;
	
	# $C_[X]: amount of C for [X] (kgC/m^2)
	if ($phenophase>=3) {
		if ($A_sum_daily-$A_stem_daily-$A_root_daily-$R_m_leaf_sum_daily-$R_g_leaf_daily-$L_leaf > 0) {
			$L_leaf += $A_sum_daily-$A_stem_daily-$A_root_daily-$R_m_leaf_sum_daily-$R_g_leaf_daily-$L_leaf;
		}
	}
	$Cg_leaf += $A_sum_daily-$A_stem_daily-$A_root_daily-$R_m_leaf_sum_daily-$R_g_leaf_daily-$L_leaf;
	if ($Cg_leaf < 0) {$Cg_leaf = 0;}
	$Cd_leaf += $L_leaf-$L_leaf_d;
	if ($Cd_leaf < 0) {$Cd_leaf = 0;}
	$C_leaf = $Cg_leaf+$Cd_leaf;
	if ($C_leaf < 0) {$C_leaf = 0;}
	$C_stem += $A_stem_daily-$R_m_stem_daily-$R_g_stem_daily-$L_stem;
	if ($C_stem < 0) {$C_stem = 0;}
	$C_root += $A_root_daily-$R_m_root_daily-$R_g_root_daily-$L_root;
	if ($C_root < 0) {$C_root = 0;}
	$C_all += $A_sum_daily-$R_a_c-$L;
	if ($C_all < 0) {$C_all = 0;}
	
	if ($leaf_type>0 && $phenophase>=4) {
		if ($A_sum_daily-$R_a_c-$L > 0) {
			$C_increase_dy++;
		} else {
			$C_increase_dy = 0;
		}
		if ($C_increase_dy == 5 && $W>=$W_wilting) {
			$phenophase = 2;
		}
	}
	
	if ($phenophase==2 && ($A_sum_daily-$A_stem_daily-$A_root_daily-$R_m_leaf_sum_daily-$R_g_leaf_daily-$L_leaf<0 or $T_a_C_mean<$T_cold[$vegetation_type][$leaf_type] or ($sun_duration_pos<11 && $vegetation_type==0 && $leaf_type==0))) {
		$phase2to3_dy++;
		if ($phase2to3_dy == 7) {
			$phenophase = 3;
		}
	} else {
		$phase2to3_dy = 0;
	}
	if ($leaf_type>0 && $phenophase==3 && $A_sum_daily-$A_stem_daily-$A_root_daily-$R_m_leaf_sum_daily-$R_g_leaf_daily-$L_leaf>=0 && $T_a_C_mean>$T_cold[$vegetation_type][$leaf_type] && $W>=$W_wilting) {
		$phase3to2_dy++;
		if ($phase3to2_dy == 7) {
			$phenophase = 2;
		}
	} else {
		$phase3to2_dy = 0;
	}
	if ($leaf_type>0 && $phenophase==3 && $DOY>=1 && $DOY<=180 && $T_a_C_mean>$T_cold[$vegetation_type][$leaf_type] && $W>=$W_wilting) {
		$phase3to2_dy2++;
		if ($phase3to2_dy2 == 7) {
			$phenophase = 2;
		}
	} else {
		$phase3to2_dy2 = 0;
	}
	if ($vegetation_type==0 && $leaf_type==0 && ($phenophase==2 or $phenophase==3) && ($T_a_C_mean<$T_cold[$vegetation_type][$leaf_type] or $sun_duration_pos<11)) {
		$phaseto4_dy++;
		if ($phaseto4_dy == 5) {
			$phenophase = 4;
		}
	} else {
		$phaseto4_dy = 0;
	}

	if ($leaf_type==0 && $leaf_dormant==0 && $phenophase!=0 && (($vegetation_type==0 && $LAI_g<0.3) or ($vegetation_type==1 && $LAI_g<0.2) or $Cg_leaf<(($C_stem+$C_root)/$allocation_parameter[$vegetation_type][$leaf_type])**(1/$allocation_parameter_index[$vegetation_type][$leaf_type])*0.075)) {
			if ($phenophase>=2) {
				$Cg_leaf = 0;	$Cd_leaf = 0;	$C_leaf = 0;	$C_all = $C_stem+$C_root;
				$leaf_dormant = 1;
				$dormancy_dys = 0;
				$leaf_normal = 0;
				$phenophase = 0;
				$leaf_fall = 0;
			}
	}

}

sub DailyMeanClimate {

	$T_a_C_sum += $T_a_C;
	$R_s_sum += $R_s_total*60*$t_step/1000000;
	$pressure_sum += $pressure;
	$rainfall_sum += $rainfall;
	$rh_sum += $rh;
	$u_z_sum += $u_z;
	$W_sum += $W;

	if ($t == 1440-$t_start) {
		$T_a_C_mean = $T_a_C_sum/(1440/$t_step);
		$pressure_mean  = $pressure_sum/(1440/$t_step);
		$rh_mean = $rh_sum/(1440/$t_step);
		$u_z_mean = $u_z_sum/(1440/$t_step);
		$W_mean = $W_sum/(1440/$t_step);
	}
	
}