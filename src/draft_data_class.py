#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 13:47:51 2018

@author: mathiasberger
"""

import pandas as pd
from numpy import isnan, asarray, arange, split, repeat, tile, multiply
from datetime import datetime
import itertools
import os

class Data(object):

    def __init__(self, path, scenario="RES_GW"):

        if path[-1]=="/":
            path=path[:-1]
        self.path = path
        self.solar_path = path + "/" + "generation/solar" + "/"
        self.w_on_path = path + "/" + "generation/wind/onshore" + "/"
        self.w_off_path = path + "/" + "generation/wind/offshore" + "/"
        self.load_path = path + "/" + "load" + "/"
        self.time_path = path + "/" + "time" + "/"
        self.capacities_path = path + "/" + "capacities" + "/"
        self.efficiencies_path = path + "/" + "efficiencies" + "/"
        self.costs_path = path + "/" + "costs" + "/"
        self.other_path = path + "/" + "other" + "/"
        self.elprice_path = path + "/" + "el_price" + "/"
        self.gas_path = path + "/" + "gas" + "/"
        self.timeseries_path = path + "/" + "timeseries" + "/"

        self.scenario = scenario

        self.time_df = fetch_file(self.time_path)
        self.capacities_df = fetch_file(self.capacities_path)
        self.efficiencies_df = fetch_file(self.efficiencies_path)
        self.costs_df = fetch_file(self.costs_path)
        self.other_df = fetch_file(self.other_path)

    @property
    def gamma_S(self):
        dict_tmp = fetch_elia_generation_data(self.solar_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output = self.n_y)
        return dict(zip(self.time, _to_float(list(asarray(list(dict_tmp.values()))))))


    @property
    def gamma_W_on(self):
        dict_tmp = fetch_elia_generation_data(self.w_on_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output = self.n_y)
        return dict(zip(self.time, _to_float(list(asarray(list(dict_tmp.values()))))))

    
    @property
    def gamma_W_off(self):
        dict_tmp = fetch_elia_generation_data(self.w_off_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output = self.n_y)
        return dict(zip(self.time, _to_float(list(asarray(list(dict_tmp.values()))))))
    
    
    @property
    def gamma_L(self):
        dict_tmp = fetch_elia_load_data(self.load_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output = self.n_y)
        return dict(zip(self.time, _to_float(list(asarray(list(dict_tmp.values()))))))

    @property
    def n_y(self):
        return float(self.time_df.loc["value", "n_y"])
    
    @property
    def n_input(self):
        return float(self.time_df.loc["value", "n_input"])

    @property
    def T(self):
        return self.time_df.loc["value", "T"]

    @property
    def t_max(self):
        return self.time_df.loc["value", "t_max"]

    @property
    def delta_t(self):
        return float(self.time_df.loc["value", "delta_t"])

    @property
    def time(self):
        t_vec = [t for t in range(self.T)]
        return t_vec

    @property
    def kappa_L(self):
        return float(self.capacities_df.loc["load",self.scenario])
    
    @property
    def growth_rate(self):
        return float(self.other_df.loc['growth', 'load'])

    @property
    def pi_L(self):
        dict_tmp = build_peakload_timeseries(self.n_y, self.kappa_L, self.growth_rate)
        return dict(zip(self.time, _to_float(list(multiply(asarray(list(self.gamma_L.values())), asarray(list(dict_tmp.values())))))))

    @property
    def kappa_W_on_0(self):
        return float(self.capacities_df.loc["won_0",self.scenario])

    @property
    def kappa_W_on_max(self):
        return float(self.capacities_df.loc["won_max",self.scenario])

    @property
    def kappa_W_off_0(self):
        return float(self.capacities_df.loc["woff_0",self.scenario])

    @property
    def kappa_W_off_max(self):
        return float(self.capacities_df.loc["woff_max",self.scenario])

    @property
    def kappa_S_0(self):
        return float(self.capacities_df.loc["pv_0",self.scenario])

    @property
    def kappa_S_max(self):
        return float(self.capacities_df.loc["pv_max",self.scenario])

    @property
    def kappa_NGNet(self):
        return float(self.capacities_df.loc["NGNet",self.scenario])

    @property
    def kappa_NG_0(self):
        y_list, p_list = fetch_timeseries_data(self.timeseries_path+"gas_capacity/")
        return build_capacity_timeseries(y_list, p_list, self.n_y)

    @property
    def kappa_NG_max(self):
        return float(self.capacities_df.loc["NG_max",self.scenario])

    @property
    def kappa_NK(self):
        y_list, p_list = fetch_timeseries_data(self.timeseries_path+"nuclear_capacity/")
        return build_capacity_timeseries(y_list, p_list, self.n_y)

    @property
    def pi_NG(self):
        dict_tmp = fetch_fluxys_demand_data(self.gas_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output=self.n_y)
        return dict(zip(self.time, _to_float(list(asarray(list(dict_tmp.values()))))))
    
#    @property
#    def pi_W_on(self):
#        return dict(zip(self.time, _to_float(list(asarray(list(self.gamma_W_on.values()))*self.kappa_W_on_max))))
#    
#    @property
#    def pi_W_off(self):
#        return dict(zip(self.time, _to_float(list(asarray(list(self.gamma_W_off.values()))*self.kappa_W_off_max))))
#    
#    @property
#    def pi_S(self):
#        return dict(zip(self.time, _to_float(list(asarray(list(self.gamma_S.values()))*self.kappa_S_max))))
#
#    @property
#    def pi_C(self):
#        prod_W_on, prod_W_off, prod_S, cons_L = self.pi_W_on, self.pi_W_off, self.pi_S, self.pi_L
#        surplus, C_list = [prod_W_on[t]+prod_W_off[t]+prod_S[t]-cons_L[t] for t in range(len(prod_W_on))], []
#        for v in surplus:
#            if v >= 0:
#                C_list.append(v)
#            else:
#                C_list.append(0)
#        return dict(zip(self.time, C_list))

    @property
    def kappa_PtG(self):
        return float(self.capacities_df.loc["PtG",self.scenario])

    @property
    def xi_H2(self):
        return float(self.capacities_df.loc["H2_s",self.scenario])

    @property
    def kappa_H2(self):
        return float(self.capacities_df.loc["H2",self.scenario])

    @property
    def kappa_H2tCH4(self):
        return float(self.capacities_df.loc["H2tCH4",self.scenario])

    @property
    def xi_CH4(self):
        return float(self.capacities_df.loc["CH4_s",self.scenario])

    @property
    def kappa_CH4tNG(self):
        return float(self.capacities_df.loc["CH4tNG",self.scenario])

    @property
    def kappa_PH(self):
        return float(self.capacities_df.loc["PH_p",self.scenario])

    @property
    def kappa_PtPH(self):
        return float(self.capacities_df.loc["PH_p",self.scenario])

    @property
    def kappa_PHtP(self):
        return float(self.capacities_df.loc["PH_p",self.scenario])

    @property
    def xi_PH(self):
        return float(self.capacities_df.loc["PH_e", self.scenario])
    
    @property
    def xi_B(self):
        return float(self.capacities_df.loc["batt_e", self.scenario])
    
    @property
    def kappa_B(self):
        return float(self.capacities_df.loc["batt_p", self.scenario])

    @property
    def kappa_trs(self):
        return float(self.capacities_df.loc["trs",self.scenario])

    @property
    def kappa_disp(self):
        return float(self.capacities_df.loc["disp",self.scenario])

    @property
    def psi_CO2(self):
        return float(self.capacities_df.loc["co2_budget",self.scenario])

    @property
    def eta_NK(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "NK"])

    @property
    def eta_H2tP(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "H2"])

    @property
    def eta_NGtP(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "ocgt"])
    
    @property
    def eta_B(self):
        return float(self.efficiencies_df.loc["EFF_SELF", "batt_s"])
    
    @property
    def eta_PtB(self):
        return float(self.efficiencies_df.loc["EFF_IN", "batt_s"])
    
    @property
    def eta_BtP(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "batt_s"])

    @property
    def eta_PtG(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "PtG"])

    @property
    def eta_H2tCH4(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "H2tCH4"])

    @property
    def eta_PHtP(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "PH"])

    @property
    def eta_PtPH(self):
        return float(self.efficiencies_df.loc["EFF_IN", "PH"])

    @property
    def eta_disp(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "disp"])

    @property
    def zeta_W_on(self):
        if self.costs_df.loc["ANNUITY", "won"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "won"], self.other_df.loc["lifetime", "won"], self.other_df.loc["wacc", "won"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "won"])
        
    @property
    def zeta_W_off(self):
        if self.costs_df.loc["ANNUITY", "woff"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "woff"], self.other_df.loc["lifetime", "woff"], self.other_df.loc["wacc", "woff"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "woff"])

    @property
    def zeta_S(self):
        if self.costs_df.loc["ANNUITY", "pv"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "pv"], self.other_df.loc["lifetime", "pv"], self.other_df.loc["wacc", "pv"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "pv"])

    @property
    def zeta_PtG(self):
        if self.costs_df.loc["ANNUITY", "PtG"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "PtG"], self.other_df.loc["lifetime", "PtG"], self.other_df.loc["wacc", "PtG"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "PtG"])

    @property
    def zeta_H2_s(self):
        if self.costs_df.loc["ANNUITY", "H2_s"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "H2_s"], self.other_df.loc["lifetime", "H2_s"], self.other_df.loc["wacc", "H2_s"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "H2_s"])
        
    @property
    def zeta_B(self):
        if self.costs_df.loc["ANNUITY", "batt"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "batt"], self.other_df.loc["lifetime", "batt"], self.other_df.loc["wacc", "batt"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "batt"]) 

    @property
    def zeta_H2(self):
        if self.costs_df.loc["ANNUITY", "H2"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "H2"], self.other_df.loc["lifetime", "H2"], self.other_df.loc["wacc", "H2"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "H2"])

    @property
    def zeta_H2tCH4(self):
        if self.costs_df.loc["ANNUITY", "H2tCH4"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "H2tCH4"], self.other_df.loc["lifetime", "H2tCH4"], self.other_df.loc["wacc", "H2tCH4"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "H2tCH4"])

    @property
    def zeta_CH4(self):
        if self.costs_df.loc["ANNUITY", "CH4_s"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "CH4_s"], self.other_df.loc["lifetime", "CH4_s"], self.other_df.loc["wacc", "CH4_s"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "CH4_s"])

    @property
    def zeta_NG(self):
        if self.costs_df.loc["ANNUITY", "ccgt"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "ccgt"], self.other_df.loc["lifetime", "ccgt"], self.other_df.loc["wacc", "ccgt"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "ccgt"])
        
    @property
    def varsigma_ENS(self):
        return float(self.costs_df.loc["OPEX", "ENS"])

    @property
    def varsigma_C(self):
        return float(self.costs_df.loc['OPEX', 'C'])

    @property
    def theta_CO2(self):
        return float(self.costs_df.loc["CO2", "ccgt"])

    @property
    def theta_W_on_f(self):
        return float(self.costs_df.loc["FOM", "won"])

    @property
    def theta_W_on_v(self):
        return float(self.costs_df.loc["VOM", "won"])

    @property
    def theta_W_off_f(self):
        return float(self.costs_df.loc["FOM", "woff"])

    @property
    def theta_W_off_v(self):
        return float(self.costs_df.loc["VOM", "woff"])

    @property
    def theta_S_f(self):
        return float(self.costs_df.loc["FOM", "pv"])

    @property
    def theta_S_v(self):
        return float(self.costs_df.loc["VOM", "pv"])

    @property
    def theta_PtG_f(self):
        return float(self.costs_df.loc["FOM", "PtG"])

    @property
    def theta_H2_s_f(self):
        return float(self.costs_df.loc["FOM","H2_s"])
    
    @property
    def theta_B_f(self):
        return float(self.costs_df.loc["FOM","batt"])
    
    @property
    def theta_B_v(self):
        return float(self.costs_df.loc["VOM","batt"])

    @property
    def theta_H2_f(self):
        return float(self.costs_df.loc["FOM", "H2"])

    @property
    def theta_H2_v(self):
        return float(self.costs_df.loc["VOM", "H2"])

    @property
    def theta_H2tCH4_f(self):
        return float(self.costs_df.loc["FOM", "H2tCH4"])

    @property
    def theta_CH4_f(self):
        return float(self.costs_df.loc["FOM","CH4_s"])

    @property
    def theta_NG_f(self):
        return float(self.costs_df.loc["FOM", "ccgt"])

    @property
    def theta_NG_v(self):
        return float(self.costs_df.loc["VOM","ccgt"])

    @property
    def theta_NG_fuel(self):
        return float(self.costs_df.loc["FUEL", "ccgt"])

    @property
    def theta_NK_v(self):
        return float(self.costs_df.loc["VOM", "NK"])

    @property
    def theta_NK_fuel(self):
        return float(self.costs_df.loc["FUEL", "NK"])

    @property
    def theta_disp_v(self):
        return float(self.costs_df.loc["VOM", "disp"])

    @property
    def theta_disp_fuel(self):
        return float(self.costs_df.loc["FUEL", "disp"])

    @property
    def theta_PH_v(self):
        return float(self.costs_df.loc["VOM","PH"])

    @property
    def theta_el(self):
        dict_tmp = fetch_elix_index_data(self.elprice_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output=self.n_y)
        return dict(zip(self.time, _to_float(list(asarray(list(dict_tmp.values()))))))

    @property
    def delta_disp_inc(self):
        return float(self.other_df.loc["delta_p","disp"])

    @property
    def delta_disp_dec(self):
        return float(self.other_df.loc["delta_m", "disp"])

    @property
    def delta_NK_inc(self):
        return float(self.other_df.loc["delta_p", "NK"])

    @property
    def delta_NK_dec(self):
        return float(self.other_df.loc["delta_m", "NK"])

    @property
    def sigma_H2(self):
        return float(self.other_df.loc["sigma", "H2_s"])

    @property
    def sigma_PH(self):
        return float(self.other_df.loc["sigma", "PH"])

    @property
    def sigma_CH4(self):
        return float(self.other_df.loc["sigma", "CH4_s"])
    
    @property
    def sigma_B(self):
        return float(self.other_df.loc["sigma", "batt"])
    
    @property
    def chi_B(self):
        return float(self.other_df.loc["chi", "batt"])
    
    @property
    def rho_B(self):
        return float(self.other_df.loc["rho", "batt"])

    @property
    def mu_disp(self):
        return float(self.other_df.loc["mu", "disp"])

    @property
    def mu_NK(self):
        return float(self.other_df.loc["mu", "NK"])

    @property
    def mu_trs(self):
        return float(self.other_df.loc["mu","trs"])

    @property
    def nu_NG_CO2(self):
        return float(self.other_df.loc["nu", "ccgt"])

    @property
    def nu_CH4_CO2(self):
        return float(self.other_df.loc["nu", "CH4"])

    @property
    def nu_disp_CO2(self):
        return float(self.other_df.loc["nu", "disp"])
    
    @property
    def nu_trs_CO2(self):
        return float(self.other_df.loc["nu", "trs"])

def fetch_file(path_datafiles):
    dirs = os.listdir(path_datafiles)
    for dirr in dirs:
        if len(dirr)>=5:
            if ".xls" not in dirr[-5:]:
                dirs.remove(dirr)
            elif dirr[0] == "~" or dirr[0] == ".":
                dirs.remove(dirr)
        else:
            dirs.remove(dirr)
    path_files = [path_datafiles+"/"+dirs[index] for index in range(len(dirs))]
    df_list = [pd.read_excel(path_files[index],sheet_name=0) for index in range(len(path_files))]
    if len(df_list)==1:
        return df_list[0]

def fetch_elia_load_data(path_datafiles, data_years_input, year_no, data_years_output):
    
    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)

    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.xls')]
    
    if len(year_no) != data_years_input:
        raise InputError('Length of input sequence should match number of input years!')
        
    n_files_yearly = 12
    file_sequence = []
    for i in arange(0, len(year_no), 1):
        files_temp = files[(year_no[i]-1)*n_files_yearly:year_no[i]*n_files_yearly]
        file_sequence.extend(files_temp)
    
    df = pd.DataFrame()
    for f in file_sequence:
        data = pd.read_excel(path_datafiles+str(f), index_col=0, usecols=[0, 1])
        data.columns = ['load']
        data.index = pd.to_datetime(data.index, dayfirst=True)
        data = data.resample('H').mean()
        data.fillna(method='ffill', inplace=True)
        df = pd.concat([df, data], axis=0)

    for y in df.index.year.unique():
        df.loc[df.index.year == y, 'peak'] = df['load'][df.index.year == y].max()

    df['load_n'] = df['load']/df['peak']
    df = df[~((df.index.month == 2) & (df.index.day == 29))]
    df.index = arange(0, df.shape[0], 1)

    df.drop(['load', 'peak'], inplace=True, axis=1)
    df.clip(0.0, 1.0, inplace=True)

    if df.isnull().values.any():
        raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

    x = split(df['load_n'], data_years_input, axis=0)
    x_r = shuffle_ts(x, data_years_input, data_years_output)

    return x_r.to_dict()

def fetch_elia_generation_data(path_datafiles, data_years_input, year_no, data_years_output):
    
    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)

    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.xls')]
    
    if len(year_no) != data_years_input:
        raise InputError('Length of input sequence should match number of input years!')
        
    n_files_yearly = 12
    file_sequence = []
    for i in arange(0, len(year_no), 1):
        files_temp = files[(year_no[i]-1)*n_files_yearly:year_no[i]*n_files_yearly]
        file_sequence.extend(files_temp)
    
    df = pd.DataFrame()

    if files[0].startswith('Solar'):

        for f in file_sequence:
            data = pd.read_excel(path_datafiles+str(f), index_col=0, skiprows=[0,1,2], usecols=[0, 1, 5, 6])
            data.columns = ['forecast', 'prod', 'peak']
            data['prod'].fillna(data['forecast'], inplace=True)
            data.index = pd.to_datetime(data.index, dayfirst=True)
            data = data.resample('H').mean()
            data.fillna(method='ffill', inplace=True)
            df = pd.concat([df, data], axis=0)

    else:

        for f in file_sequence:
            data = pd.read_excel(path_datafiles+str(f), index_col=0, skiprows=[0,1,2], usecols=[0, 3, 4, 6])
            data.columns = ['forecast', 'prod', 'peak']
            data['prod'].fillna(data['forecast'], inplace=True)
            data.index = pd.to_datetime(data.index, dayfirst=True)
            data = data.resample('H').mean()
            data.fillna(method='ffill', inplace=True)
            df = pd.concat([df, data], axis=0)

    df['prod_n'] = df['prod']/df['peak']
    df = df[~((df.index.month == 2) & (df.index.day == 29))]
    df.index = arange(0, df.shape[0], 1)

    df.drop(['prod', 'peak', 'forecast'], inplace=True, axis=1)
    df.clip(0.0, 1.0, inplace=True)
    
    df[df < 0.01] = 0.0

    if df.isnull().values.any():
        raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

    x = split(df['prod_n'], data_years_input, axis=0)
    x_r = shuffle_ts(x, data_years_input, data_years_output)

    return x_r.to_dict()


def fetch_fluxys_demand_data(path_datafiles, data_years_input, year_no, data_years_output):
    
    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)

    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.csv')]
    
    if len(year_no) != data_years_input:
        raise InputError('Length of input sequence should match number of input years!')
    
    n_files_yearly = 4
    file_sequence = []
    for i in arange(0, len(year_no), 1):
        files_temp = files[(year_no[i]-1)*n_files_yearly:year_no[i]*n_files_yearly]
        file_sequence.extend(files_temp)

    df = pd.DataFrame()
    
    for f in file_sequence:
        data = pd.read_csv(path_datafiles+'/'+str(f), sep=';', usecols=['subGrid', 'gasDay', 'gasHour', 'clientType', 'physicalFlow'])
        if f == file_sequence[0]:
            d_start_str = str(data['gasDay'].iloc[0])+' '+str(int(data['gasHour'].iloc[0])+5)
            d_start = pd.to_datetime(d_start_str, format='%d/%m/%Y %H')
        else:
            d_start = df.index[-1] + pd.Timedelta(hours=1)
            
        for item in data['clientType'].unique():
            data_client = data[data['clientType'] == item]
            for net in data_client['subGrid'].unique():
                data_network = data_client[data_client['subGrid'] == net]
                data_network.index = pd.date_range(d_start, freq='H', periods=data_network.shape[0])
                df = pd.concat([df, data_network], axis=0, ignore_index=False)
    
    df = df[~((df.index.month == 2) & (df.index.day == 29))]            
    df = df[df['clientType'] != 'End User Domestic Exit Point PP']
    
    df_sum = df.groupby(df.index)['physicalFlow'].sum().reset_index(drop=True)
#    df_sum = df_sum * (-1) * 1e-3 
    df_sum = df_sum * (-1) * 1e-6    
    
    # ADJUST TIME SERIES FOR THE GAS DAY THAT STARTS AT 6 AM
    data_shift = df_sum[:6]
    df_shifted = df_sum.shift(-6)
    df_shifted[-6:] = data_shift

    df_shifted = df_shifted.astype(float)

    if df_shifted.isnull().values.any():
        raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

    x = split(df_shifted, data_years_input, axis=0)
    x_r = shuffle_ts(x, data_years_input, data_years_output)

    return x_r.to_dict()


def build_capacity_timeseries(y_list, p_list, data_years_output, start_year = 2019):
    
    data_years_output = int(data_years_output)
    
    end_year = start_year + data_years_output - 1
    x = pd.Series(index = arange(start_year, end_year + 1))

    if end_year != y_list[-1]:
        raise InputError('End year MUST match the last year in the list')

    for i in arange(0, len(p_list)):
        if i == 0:
            x.loc[start_year : y_list[0]] = p_list[0]
        else:
            x.loc[y_list[i-1] + 1 : y_list[i]] = p_list[i]

    ts_cap = pd.Series()

    for idx in x.index:
            x_r = pd.Series(repeat(x[idx], 8760))
            ts_cap = pd.concat([ts_cap, x_r], axis=0, ignore_index=True)

    if ts_cap.isnull().values.any():
        raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

    return ts_cap.to_dict()

def build_peakload_timeseries(data_years_output, base_peak, growth_rate):
    
        data_years_output = int(data_years_output)
        
        x = pd.Series(index=arange(1, data_years_output+1, 1))
        for i in range(1, data_years_output+1):
            x[i] = round(base_peak * (1+growth_rate)**i, 1)

        ts_cap = pd.Series()

        for idx in x.index:
            x_r = pd.Series(repeat(x[idx], 8760))
            ts_cap = pd.concat([ts_cap, x_r], axis=0, ignore_index=True)

        if ts_cap.isnull().values.any():
            raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

        return ts_cap.to_dict()

def fetch_elix_index_data(path_datafiles, data_years_input, year_no, data_years_output):
    # DATA ACQUISITION HERE MIGHT BE AN ISSUE...
    
    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)
    
    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.csv')]

    df = pd.DataFrame()
    df_sequence = pd.DataFrame()
    
    for f in files:
        ts = pd.read_csv(path_datafiles+'/'+str(f), sep=';', index_col=0)
        df = pd.concat([df, ts], axis=0, ignore_index=True)
    
    df.clip(lower=0.0, inplace=True)
#    df.drop(df.tail(24).index, inplace=True)
    
    for i in arange(0, len(year_no), 1):
        df_temp = df.loc[(year_no[i]-1)*8760:(year_no[i]*8760-1)]
        df_sequence = pd.concat([df_sequence, df_temp], axis=0, ignore_index=False)
    
    df = df_sequence

    df[df < 1.0] = 0.0
    df = df * 1e-3
        
    if df.isnull().values.any():
        raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

    x = split(df['price'], data_years_input, axis=0)
    x_r = shuffle_ts(x, data_years_input, data_years_output)

    return x_r.to_dict()


def capex_annuity(capex, lifetime, wacc, n_y):
    
    #x = (capex/(1+wacc))*(wacc/(1-(1+wacc)**(-lifetime)))
    x = n_y*capex/lifetime
    return x


def fetch_timeseries_data(path_datafiles, output_dict=False, output_list=True):

    files = [f for f in os.listdir(path_datafiles) if f.endswith('.xlsx')]
    for f in files:
        df = pd.read_excel(path_datafiles+str(f))
    if output_dict:
        return dict(zip(list(df["time"]), list(df["value"])))
    if output_list:
        return list(df["time"]), list(df["value"])


def shuffle_ts(data, input_years, output_years):

    if input_years == 1:
        seq = repeat(0, output_years)

    elif input_years == 2:
        seq = tile([0, 1], (output_years // input_years + 1))

    elif input_years == 3:
        seq = [0, 1, 2, 1, 0, 2, 0, 2, 1, 2, 1, 0, 2, 0, 1, 1, 2, 0, 1, 0, 2, 0, 2, 1, 2, 1, 0, 2, 0, 1, 1, 2, 0]

    elif input_years == 4:
        seq = [0, 1, 2, 3, 1, 3, 0, 2, 3, 0, 1, 2, 3, 2, 1, 0, 1, 0, 3, 2, 0, 3, 2, 1, 3, 2, 0, 1, 2, 1, 3, 0]

    elif input_years == 5:
        seq = [0, 1, 2, 3, 4, 1, 2, 4, 3, 0, 3, 0, 1, 4, 2, 4, 3, 1, 2, 0, 1, 4, 0, 2, 3, 2, 3, 4, 0, 1, 4, 1, 3, 0, 2]

    elif input_years == 6:
        seq = [0, 1, 2, 3, 4, 5, 1, 2, 4, 5, 0, 3, 2, 0, 3, 4, 5, 2, 4, 3, 5, 2, 0, 1, 3, 4, 0, 5, 1, 2, 5, 0, 1, 3, 2, 4]

    else:
        seq = tile(list(range(0, input_years)), (output_years // input_years + 1))

    if output_years > len(seq):
        raise InputError('Time horizon ("output_years" argument) larger than sequence length. Update accordingly.')

    ts = pd.Series()
    for i in arange(0, output_years):
        ts = pd.concat([ts, data[seq[i]]], axis=0, ignore_index=True)

    return ts

def year_selection(path_datafiles):
    
    data = pd.read_excel(path_datafiles + '/time.xlsx')
    
    str_list = data.loc['value', 'year_no']
    if isinstance(str_list, float):
        l = [int(str_list)]
    else:
        l = [int(i) for i in str_list.split(',')]

    return l 


def _sweep_list(input_list, value):
    cnt = 0
    for ind in input_list:
        if ind==value:
            cnt+=1
        else:
            break
    return cnt

def interpolate_nan(df):
    df_bool = list(isnan(df["Consumption"]))
    nan_indices = []
    for index in range(len(df)):
        if df_bool[index]:
            nan_indices.append(index)
    nan_indices_delta = [nan_indices[1+index]-nan_indices[index] for index in range(len(nan_indices)-1)]
    index=0
    safety = []
    while index < len(nan_indices):
        if index<len(nan_indices_delta):
            if nan_indices_delta[index]==1:
                cnt = _sweep_list(nan_indices_delta[index:], 1)
                interp_index = nan_indices[index] + cnt + 1
                interp_indices = [nan_indices[index] + ind for ind in range(cnt + 1)]
                delta_cons = (df["Consumption"][nan_indices[index]] - df["Consumption"][interp_index])/(cnt + 1)
                for ind in interp_indices:
                    # for some reason the next line does not work
                    df.loc[ind, "Consumption"] = df["Consumption"][nan_indices[index]] - delta_cons
                    index+=(cnt+1)
            else:
                df.loc[nan_indices[index], "Consumption"] = (df["Consumption"][nan_indices[index]-1]+df["Consumption"][nan_indices[index]+1])/2
                index+=1
        else:
            df.loc[nan_indices[index], "Consumption"] = df["Consumption"][nan_indices[index]-1]
            index+=1
        safety.append(index)
        if len(safety)>len(df):
            break
    return df

def _to_float(input_list):
    output_list = [float(input_list[index]) for index in range(len(input_list))]
    return output_list

class Error(Exception):
    pass

class InputError(Error):
    def __init__(self, message):
        self.message = message

#dirs = os.listdir(path)
#path_files = [path+"/"+dirs[index] for index in range(len(dirs))]
#for dirr in dirs:
#    if dirr[0]==".":
#        dirs.remove(dirr)
#    elif "." in dirr[-5:-1]:
#        dirs.remove(dirr)
#
#paths, branch_temp = [], []
#for dirr in dirs:
#    root_dir = path + dirr
#    dirs_temp = os.listdir(root_dir)
#    if all([".xls" in dirs_temp[index] for index in range(len(dirs_temp))]):
#        paths.append(root_dir)
#    else:
#        for dirrt in dirs_temp:
#            current_dir, cnt = root_dir +"/"+ dirrt, 0
#            dirs_t = os.listdir(current_dir)
#            for dirst in dirs_t:
#                if dirst[0]==".":
#                    dirs_t.remove(dirst)
#            print(dirs_t)
#            while not all([".xls" in dirs_t[index] for index in range(len(dirs_t))]) and cnt < 50:
#                current_dir = current_dir
#                branch_temp.append(current_dir)
#                dirs_t = os.listdir(current_dir)
#                #current_dir = current_dir +
#                cnt+=1
