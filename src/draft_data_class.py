#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 13:47:51 2018

@author: mathiasberger
"""

import pandas as pd
from numpy import isnan, asarray, arange, split, repeat, tile, multiply, where
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
        self.gas_HIDM_path = path + "/" + "gas_HIDM" + "/"
        self.timeseries_path = path + "/" + "timeseries" + "/"

        self.scenario = scenario

        self.time_df = fetch_file(self.time_path)
        self.capacities_df = fetch_file(self.capacities_path)
        self.efficiencies_df = fetch_file(self.efficiencies_path)
        self.costs_df = fetch_file(self.costs_path)
        self.other_df = fetch_file(self.other_path)

    @property
    def gamma_S(self):
        return fetch_elia_generation_data(self.solar_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output = self.n_y)

    @property
    def gamma_W_on(self):
        return fetch_elia_generation_data(self.w_on_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output = self.n_y)
    
    @property
    def gamma_W_off(self):
        return fetch_elia_generation_data(self.w_off_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output = self.n_y)
    
    @property
    def gamma_L(self):
        return fetch_elia_load_data(self.load_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output = self.n_y)

    @property
    def n_y(self):
        return float(self.time_df.loc["value", "n_y"])
    
    @property
    def n_input(self):
        return float(self.time_df.loc["value", "n_input"])

    @property
    def T(self):
        return int(self.time_df.loc["value", "T"])

    @property
    def t_max(self):
        return self.time_df.loc["value", "t_max"]

    @property
    def delta_t(self):
        return float(self.time_df.loc["value", "delta_t"])
    
    @property
    def n_delta_t(self):
        return int(24/self.delta_t)

    @property
    def n_weeks(self):
        return int(self.t_max/168)

    @property
    def time(self):
        t_vec = [t for t in range(self.T)]
        return t_vec
    
    @property
    def time_d(self):
        t_vec = self.time
        t_d_vec = set([24*int(t/24) for t in t_vec])
        #t_d_vec.sort()
        return t_d_vec

    @property
    def kappa_L(self):
        return float(self.capacities_df.loc["load",self.scenario])
    
    @property
    def growth_rate_E(self):
        return float(self.other_df.loc['growth', 'load'])

    @property
    def lambda_E(self):
        return dict(self.gamma_L.multiply(build_peakload_timeseries(self.n_y, self.kappa_L, self.growth_rate_E), axis=0))

    @property
    def COP(self):
        return float(self.other_df.loc['COP', 'HP'])

    @property
    def lambda_E_HT(self):
        return fetch_electric_heating_data(self.timeseries_path, self.n_y, self.COP)
    
    @property
    def lambda_E_TR(self):
        return fetch_electric_transport_data(self.timeseries_path+'ev_profile_electric/', self.n_input, year_selection(self.time_path), self.n_y, self.t_max, self.n_delta_t)

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
    def kappa_NG_I(self):
        return float(self.capacities_df.loc["NG_I",self.scenario])

    @property
    def kappa_OCGT_0(self):
        # y_list, p_list = fetch_timeseries_data(self.timeseries_path+"ocgt_capacity/")
        # return build_timeseries(y_list, p_list, self.n_y)

        return float(self.capacities_df.loc["OCGT", self.scenario])
    
    @property
    def kappa_CCGT_0(self):
        # y_list, p_list = fetch_timeseries_data(self.timeseries_path+"ccgt_capacity/")
        # return build_timeseries(y_list, p_list, self.n_y)

        return float(self.capacities_df.loc["CCGT", self.scenario])

    @property
    def kappa_NG_max(self):
        return float(self.capacities_df.loc["NG_max",self.scenario])

    # @property
    # def kappa_NK(self):
    #     y_list, p_list = fetch_timeseries_data(self.timeseries_path+"nuclear_capacity/")
    #     return build_timeseries(y_list, p_list, self.n_y)

    @property
    def kappa_NK(self):
        return float(self.capacities_df.loc["NK",self.scenario])

    @property
    def kappa_E(self):
        return sum(self.lambda_E.values()) + sum(self.lambda_E_TR.values()) + sum(self.lambda_E_HT.values())

    @property
    def lambda_NG_HT(self):
        return fetch_fluxys_demand_data(self.gas_HIDM_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output=self.n_y,
                                              client_type='residential')

    @property
    def lambda_NG_ID(self):
        return fetch_fluxys_demand_data(self.gas_HIDM_path, data_years_input = self.n_input,
                                              year_no = year_selection(self.time_path), data_years_output=self.n_y,
                                              client_type='industry')
        
    @property
    def lambda_NGtH2(self):
        return dict.fromkeys(arange(0, self.t_max), self.capacities_df.loc["NGtH2_0", self.scenario])
    
    @property
    def lambda_NG_TR(self):
        return fetch_gas_markets_data(self.timeseries_path+"NG_transport/", data_years_input = self.n_input,
                                      year_no = year_selection(self.time_path), data_years_output = self.n_y, 
                                      vol_multiple = self.capacities_df.loc["NGtransport", self.scenario],
                                      price_multiple = self.other_df.loc['price', 'NGtrs'])[0]
    
#    @property
#    def lambda_NG_TR(self):
#        return fetch_gas_transit_data(self.timeseries_path+"NG_transit/", data_years_input = self.n_input,
#                                      year_no = year_selection(self.time_path), data_years_output = self.n_y)

    @property
    def lambda_H2_TR(self):
        return fetch_gas_markets_data(self.timeseries_path+"H2_transport/", data_years_input = self.n_input,
                                      year_no = year_selection(self.time_path), data_years_output = self.n_y, 
                                      vol_multiple = self.capacities_df.loc["H2transport", self.scenario],
                                      price_multiple = self.other_df.loc['price', 'H2trs'])[0]
    
    @property
    def lambda_H2_ID(self):
        return fetch_gas_markets_data(self.timeseries_path+"H2_industry/", data_years_input = self.n_input,
                                      year_no = year_selection(self.time_path), data_years_output = self.n_y, 
                                      vol_multiple = self.capacities_df.loc["H2industry", self.scenario], 
                                      price_multiple = self.other_df.loc['price', 'H2ind'])[0]
    
    @property
    def kappa_EL(self):
        return float(self.capacities_df.loc["EL",self.scenario])

    @property
    def kappa_FC(self):
        return float(self.capacities_df.loc["FC",self.scenario])
    
    @property
    def xi_H2S(self):
        return float(self.capacities_df.loc["H2S",self.scenario])

    @property
    def kappa_H2S(self):
        return float(self.capacities_df.loc["P_H2S",self.scenario])

    @property
    def kappa_MT(self):
        return float(self.capacities_df.loc["MT",self.scenario])

    @property
    def xi_CH4(self):
        return float(self.capacities_df.loc["CH4S",self.scenario])

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
    def kappa_NGS(self):
        return float(self.capacities_df.loc["P_NGS", self.scenario])
    
    @property
    def xi_NGS(self):
        return float(self.capacities_df.loc["NGS", self.scenario])
    
    @property
    def xi_B(self):
        return float(self.capacities_df.loc["batt_e", self.scenario])
    
    @property
    def kappa_B(self):
        return float(self.capacities_df.loc["batt_p", self.scenario])

    @property
    def kappa_IE(self):
        return float(self.capacities_df.loc["trs",self.scenario])

    @property
    def kappa_H2_I_size(self):
        return float(self.capacities_df.loc["H2_I_size",self.scenario])

    @property
    def kappa_H2_I_freq(self):
        return float(self.capacities_df.loc["H2_I_freq",self.scenario])

    @property
    def kappa_H2_I(self):
        return build_hydrogen_supply_profile(self.n_input, self.n_y,
                                             tanker_size=self.kappa_H2_I_size, tankers_per_week=self.kappa_H2_I_freq, hours_to_discharge=24, start_year=2014)

    @property
    def kappa_CO2_E(self):
        return float(self.capacities_df.loc["CO2_trs",self.scenario])
    
    @property
    def kappa_SMR(self):
        return float(self.capacities_df.loc["SMR",self.scenario])

    @property
    def kappa_CHP(self):
        return float(self.capacities_df.loc["CHP",self.scenario])
    
    @property
    def kappa_BM(self):
        return float(self.capacities_df.loc["biomass",self.scenario])
    
    @property
    def kappa_WS(self):
        return float(self.capacities_df.loc["waste",self.scenario])

    @property
    def psi_CO2(self):
        return float(self.capacities_df.loc["co2_budget",self.scenario] * self.n_input)

    @property
    def xi_CO2S(self):
        return float(self.capacities_df.loc["CO2S", self.scenario])

    @property
    def kappa_CO2S(self):
        return float(self.capacities_df.loc["Q_CO2S",self.scenario])
    
    @property
    def kappa_L_E_ENS(self):
        return float(self.capacities_df.loc["ENS_E", self.scenario])
    
    @property
    def kappa_L_NG_ENS(self):
        return float(self.capacities_df.loc["ENS_NG", self.scenario])
    
    @property
    def kappa_L_H2_ENS(self):
        return float(self.capacities_df.loc["ENS_H2", self.scenario])
    
    @property
    def kappa_CO2_OCGT_CCS(self):
        return float(self.capacities_df.loc["OCGT_CCS", self.scenario])
    
    @property
    def kappa_CO2_CCGT_CCS(self):
        return float(self.capacities_df.loc["CCGT_CCS", self.scenario])
    
    @property
    def kappa_CO2_CHP_CCS(self):
        return float(self.capacities_df.loc["CHP_CCS", self.scenario])
    
    @property
    def kappa_CO2_BM_CCS(self):
        return float(self.capacities_df.loc["BM_CCS", self.scenario])
    
    @property
    def kappa_CO2_WS_CCS(self):
        return float(self.capacities_df.loc["WS_CCS", self.scenario])
    
    @property
    def kappa_CO2_SMR_CCS(self):
        return float(self.capacities_df.loc["SMR_CCS", self.scenario])

    @property
    def kappa_CO2_ACC(self):
        return float(self.capacities_df.loc["kappa_ACC",self.scenario])

    @property
    def eta_NK(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "NK"])

    @property
    def eta_FC(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "FC"])

    @property
    def eta_OCGT(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "OCGT"])
    
    @property
    def eta_CCGT(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "CCGT"])
    
#    @property
#    def eta_NGtPH(self):
#        return float(self.efficiencies_df.loc["EFF_TOT", "ccgt"])
    
    @property
    def eta_NG_CCS(self):
        return float(self.efficiencies_df.loc["EFF_CCS", "OCGT"])
    
    @property
    def eta_B(self):
        return float(self.efficiencies_df.loc["EFF_SELF", "B"])
    
    @property
    def eta_PtB(self):
        return float(self.efficiencies_df.loc["EFF_IN", "B"])
    
    @property
    def eta_BtP(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "B"])

    @property
    def eta_EL(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "EL"])

    @property
    def eta_MT(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "MT"])

    @property
    def eta_PHtP(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "PH"])

    @property
    def eta_PtPH(self):
        return float(self.efficiencies_df.loc["EFF_IN", "PH"])
    
    @property
    def eta_NGS(self):
        return float(self.efficiencies_df.loc["EFF_IN", "NGS"])

    @property
    def eta_CHP(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "CHP"])
    
    @property
    def eta_CHP_CCS(self):
        return float(self.efficiencies_df.loc["EFF_CCS", "CHP"])
    
    @property
    def eta_BM(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "BM"])
    
    @property
    def eta_BM_CCS(self):
        return float(self.efficiencies_df.loc["EFF_CCS", "BM"])
    
    @property
    def eta_WS(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "WS"])
    
    @property
    def eta_WS_CCS(self):
        return float(self.efficiencies_df.loc["EFF_CCS", "WS"])
    
    @property
    def eta_SMR(self):
        return float(self.efficiencies_df.loc["EFF_OUT", "SMR"])
    
    @property
    def eta_SMR_CCS(self):
        return float(self.efficiencies_df.loc["EFF_CCS", "SMR"])

    @property
    def zeta_W_on(self):
        if self.costs_df.loc["ANNUITY", "WON"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "WON"], self.other_df.loc["lifetime", "WON"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "WON"])
        
    @property
    def theta_W_on_f(self):
        return float(self.costs_df.loc["FOM", "WON"])

    @property
    def theta_W_on_v(self):
        return float(self.costs_df.loc["VOM", "WON"])
        
    @property
    def zeta_W_off(self):
        if self.costs_df.loc["ANNUITY", "WOFF"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "WOFF"], self.other_df.loc["lifetime", "WOFF"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "WOFF"])
        
    @property
    def theta_W_off_f(self):
        return float(self.costs_df.loc["FOM", "WOFF"])
    
    @property
    def theta_W_off_v(self):
        return float(self.costs_df.loc["VOM", "WOFF"])

    @property
    def zeta_S(self):
        if self.costs_df.loc["ANNUITY", "PV"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "PV"], self.other_df.loc["lifetime", "PV"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "PV"])
        
    @property
    def theta_S_f(self):
        return float(self.costs_df.loc["FOM", "PV"])

    @property
    def theta_S_v(self):
        return float(self.costs_df.loc["VOM", "PV"])

    @property
    def zeta_EL(self):
        if self.costs_df.loc["ANNUITY", "EL"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "EL"], self.other_df.loc["lifetime", "EL"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "EL"])
        
    @property
    def theta_EL_f(self):
        return float(self.costs_df.loc["FOM", "EL"])

    @property
    def zeta_H2S(self):
        if self.costs_df.loc["ANNUITY", "H2S"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "H2S"], self.other_df.loc["lifetime", "H2S"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "H2S"])
        
    @property
    def theta_H2S_f(self):
        return float(self.costs_df.loc["FOM","H2S"])
    
    @property
    def zeta_CO2S(self):
        if self.costs_df.loc["ANNUITY", "CO2S"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "CO2S"], self.other_df.loc["lifetime", "CO2S"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "CO2S"])
        
    @property
    def theta_CO2S_f(self):
        return float(self.costs_df.loc["FOM","CO2S"])
    
    @property
    def theta_CO2S_v(self):
        return float(self.costs_df.loc["VOM","CO2S"])
        
    @property
    def zeta_B_E(self):
        if self.costs_df.loc["ANNUITY", "B_E"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "B_E"], self.other_df.loc["lifetime", "B_E"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "B_E"])
        
    @property
    def zeta_B_P(self):
        if self.costs_df.loc["ANNUITY", "B_P"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "B_P"], self.other_df.loc["lifetime", "B_P"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "B_P"])
    
    @property
    def theta_B_E_f(self):
        return float(self.costs_df.loc["FOM","B_E"])
    
    @property
    def theta_B_P_f(self):
        return float(self.costs_df.loc["FOM","B_P"])

    @property
    def zeta_FC(self):
        if self.costs_df.loc["ANNUITY", "FC"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "FC"], self.other_df.loc["lifetime", "FC"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "FC"])
        
    @property
    def theta_FC_f(self):
        return float(self.costs_df.loc["FOM", "FC"])

    @property
    def theta_FC_v(self):
        return float(self.costs_df.loc["VOM", "FC"])

    @property
    def zeta_MT(self):
        if self.costs_df.loc["ANNUITY", "MT"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "MT"], self.other_df.loc["lifetime", "MT"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "MT"])
        
    @property
    def theta_MT_f(self):
        return float(self.costs_df.loc["FOM", "MT"])

    @property
    def zeta_OCGT(self):
        if self.costs_df.loc["ANNUITY", "OCGT"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "OCGT"], self.other_df.loc["lifetime", "OCGT"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "OCGT"])
        
    @property
    def theta_OCGT_f(self):
        return float(self.costs_df.loc["FOM", "OCGT"])

    @property
    def theta_OCGT_v(self):
        return float(self.costs_df.loc["VOM","OCGT"])
        
    @property
    def zeta_CCGT(self):
        if self.costs_df.loc["ANNUITY", "CCGT"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "CCGT"], self.other_df.loc["lifetime", "CCGT"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "CCGT"])
        
    @property
    def theta_CCGT_f(self):
        return float(self.costs_df.loc["FOM", "CCGT"])

    @property
    def theta_CCGT_v(self):
        return float(self.costs_df.loc["VOM","CCGT"])
        
    @property
    def zeta_NG_CCS(self):
        if self.costs_df.loc["ANNUITY", "CCGT"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX_CCS", "CCGT"], self.other_df.loc["lifetime", "CCS"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX_CCS", "CCGT"])
        
    @property
    def theta_NG_CCS_f(self):
        return float(self.costs_df.loc["FOM_CCS", "CCGT"])

    @property
    def theta_NG_CCS_v(self):
        return float(self.costs_df.loc["VOM_CCS","CCGT"])
        
    @property
    def zeta_CHP_CCS(self):
        if self.costs_df.loc["ANNUITY", "CHP"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX_CCS", "CHP"], self.other_df.loc["lifetime", "CCS"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX_CCS", "CHP"])
        
    @property
    def theta_CHP_CCS_f(self):
        return float(self.costs_df.loc["FOM_CCS", "CHP"])
    
    @property
    def theta_CHP_CCS_v(self):
        return float(self.costs_df.loc["VOM_CCS", "CHP"])
        
    @property
    def zeta_BM_CCS(self):
        if self.costs_df.loc["ANNUITY", "BM"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX_CCS", "BM"], self.other_df.loc["lifetime", "CCS"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX_CCS", "BM"])
        
    @property
    def theta_BM_CCS_f(self):
        return float(self.costs_df.loc["FOM_CCS", "BM"])
    
    @property
    def theta_BM_CCS_v(self):
        return float(self.costs_df.loc["VOM_CCS", "BM"])
        
    @property
    def zeta_WS_CCS(self):
        if self.costs_df.loc["ANNUITY", "WS"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX_CCS", "WS"], self.other_df.loc["lifetime", "CCS"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX_CCS", "WS"])
        
    @property
    def theta_WS_CCS_f(self):
        return float(self.costs_df.loc["FOM_CCS", "WS"])
    
    @property
    def theta_WS_CCS_v(self):
        return float(self.costs_df.loc["VOM_CCS", "WS"])
        
    @property
    def zeta_SMR(self):
        if self.costs_df.loc["ANNUITY", "SMR"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "SMR"], self.other_df.loc["lifetime", "SMR"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "SMR"])
        
    @property
    def theta_SMR_f(self):
        return float(self.costs_df.loc["FOM", "SMR"])

    @property
    def theta_SMR_v(self):
        return float(self.costs_df.loc["VOM", "SMR"])
        
    @property
    def zeta_SMR_CCS(self):
        if self.costs_df.loc["ANNUITY", "SMR"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX_CCS", "SMR"], self.other_df.loc["lifetime", "CCS"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX_CCS", "SMR"])
        
    @property
    def theta_SMR_CCS_f(self):
        return float(self.costs_df.loc["FOM_CCS", "SMR"])
    
    @property
    def theta_SMR_CCS_v(self):
        return float(self.costs_df.loc["VOM_CCS", "SMR"])

    @property
    def zeta_ACC(self):
        if self.costs_df.loc["ANNUITY", "ACC"] == 1:
            compute_annuity = False
        else:
            compute_annuity = True
        if compute_annuity:
            return float(capex_annuity(self.costs_df.loc["CAPEX", "ACC"], self.other_df.loc["lifetime", "ACC_E"], self.n_y))
        else:
            return float(self.costs_df.loc["CAPEX", "ACC"])
    
    @property
    def theta_ACC_f(self):
        return float(self.costs_df.loc["FOM", "ACC"])
        
    @property
    def varsigma_ENS_E(self):
        return float(self.costs_df.loc["OPEX", "ENS_E"])
    
    @property
    def varsigma_ENS_NG(self):
        return float(self.costs_df.loc["OPEX", "ENS_NG"])
    
    @property
    def varsigma_ENS_H2(self):
        return float(self.costs_df.loc["OPEX", "ENS_H2"])

    @property
    def theta_H2_IE(self):
        return dict.fromkeys(arange(0, self.t_max), self.costs_df.loc["OPEX", "H2_I"])
    
    @property
    def theta_PH_f(self):
        return float(self.costs_df.loc["FOM","PH"])
    
    @property
    def theta_PH_v(self):
        return float(self.costs_df.loc["VOM","PH"])
    
    @property
    def theta_NGS_f(self):
        return float(self.costs_df.loc["FOM", "NGS"])
    
    @property
    def theta_NK_f(self):
        return float(self.costs_df.loc["FOM", "NK"])
    
    @property
    def theta_NK_v(self):
        return float(self.costs_df.loc["VOM", "NK"])
    
    @property
    def theta_CHP_f(self):
        return float(self.costs_df.loc["FOM", "CHP"])

    @property
    def theta_CHP_v(self):
        return float(self.costs_df.loc["VOM", "CHP"])
    
    @property
    def theta_BM_f(self):
        return float(self.costs_df.loc["FOM", "BM"])
    
    @property
    def theta_BM_v(self):
        return float(self.costs_df.loc["VOM", "BM"])
    
    @property
    def theta_WS_f(self):
        return float(self.costs_df.loc["FOM", "WS"])
    
    @property
    def theta_WS_v(self):
        return float(self.costs_df.loc["VOM", "WS"])
    
    @property
    def theta_NG_I(self):
        return fetch_gas_price_data(self.timeseries_path+"gas_price", data_years_input = 1, 
                                    year_no = [1], data_years_output=self.n_y, units = 'MEuroGWh')
    
    @property
    def theta_NK_fuel(self):
        return float(self.costs_df.loc["FUEL", "NK"])
    
    @property
    def theta_BM_fuel(self):
        return float(self.costs_df.loc["FUEL", "BM"])

    @property
    def theta_WS_fuel(self):
        return float(self.costs_df.loc["FUEL", "WS"])

    @property
    def theta_CO2(self):
        return float(self.costs_df.loc["CO2", "CCGT"])

    @property
    def theta_IE(self):
        return fetch_elix_index_data(self.elprice_path, data_years_input = self.n_input, 
                                              year_no = year_selection(self.time_path), data_years_output=self.n_y)
    @property
    def theta_CO2_E(self):
        return float(self.costs_df.loc["OPEX", "CO2_E"])
    
#    @property
#    def varsigma_H2_TRM(self):
#        return fetch_gas_markets_data(self.timeseries_path+"H2_transport/", data_years_input = self.n_input,
#                                      year_no = year_selection(self.time_path), data_years_output = self.n_y, 
#                                      vol_multiple = self.capacities_df.loc["H2transport", self.scenario], 
#                                      price_multiple = self.other_df.loc['price', 'H2trs'])[1]
#    
#    @property
#    def varsigma_H2_IDM(self):
#        return fetch_gas_markets_data(self.timeseries_path+"H2_industry/", data_years_input = self.n_input,
#                                      year_no = year_selection(self.time_path), data_years_output = self.n_y, 
#                                      vol_multiple = self.capacities_df.loc["H2industry", self.scenario], 
#                                      price_multiple = self.other_df.loc['price', 'H2ind'])[1]
#    
#    @property
#    def varsigma_NG_TRM(self):
#        return fetch_gas_markets_data(self.timeseries_path+"NG_transport/", data_years_input = self.n_input,
#                                      year_no = year_selection(self.time_path), data_years_output = self.n_y, 
#                                      vol_multiple = self.capacities_df.loc["NGtransport", self.scenario], 
#                                      price_multiple = self.other_df.loc['price', 'NGtrs'])[1]

    @property
    def delta_CHP_inc(self):
        return float(self.other_df.loc["delta_p","CHP"])

    @property
    def delta_CHP_dec(self):
        return float(self.other_df.loc["delta_m", "CHP"])
    
    @property
    def delta_BM_inc(self):
        return float(self.other_df.loc["delta_p","BM"])

    @property
    def delta_BM_dec(self):
        return float(self.other_df.loc["delta_m", "BM"])
    
    @property
    def delta_WS_inc(self):
        return float(self.other_df.loc["delta_p","WS"])

    @property
    def delta_WS_dec(self):
        return float(self.other_df.loc["delta_m", "WS"])

    @property
    def delta_NK_inc(self):
        return float(self.other_df.loc["delta_p", "NK"])

    @property
    def delta_NK_dec(self):
        return float(self.other_df.loc["delta_m", "NK"])

    @property
    def delta_OCGT_inc(self):
        return float(self.other_df.loc["delta_p", "OCGT"])

    @property
    def delta_OCGT_dec(self):
        return float(self.other_df.loc["delta_m", "OCGT"])

    @property
    def delta_CCGT_inc(self):
        return float(self.other_df.loc["delta_p", "CCGT"])

    @property
    def delta_CCGT_dec(self):
        return float(self.other_df.loc["delta_m", "CCGT"])

    @property
    def delta_EL_inc(self):
        return float(self.other_df.loc["delta_p", "EL"])

    @property
    def delta_EL_dec(self):
        return float(self.other_df.loc["delta_m", "EL"])

    @property
    def delta_MT_inc(self):
        return float(self.other_df.loc["delta_p", "MT"])

    @property
    def delta_MT_dec(self):
        return float(self.other_df.loc["delta_m", "MT"])

    @property
    def delta_FC_inc(self):
        return float(self.other_df.loc["delta_p", "FC"])

    @property
    def delta_FC_dec(self):
        return float(self.other_df.loc["delta_m", "FC"])

    @property
    def sigma_H2S(self):
        return float(self.other_df.loc["sigma", "H2S"])

    @property
    def sigma_PH(self):
        return float(self.other_df.loc["sigma", "PH"])
    
    @property
    def sigma_NGS(self):
        return float(self.other_df.loc["sigma", "NGS"])
    
    @property
    def sigma_B(self):
        return float(self.other_df.loc["sigma", "B_E"])
    
    @property
    def chi_H2S(self):
       return float(self.other_df.loc["chi", "H2S"])
    
    @property
    def chi_CO2S(self):
       return float(self.other_df.loc["chi", "CO2S"])
    
    @property
    def rho_B(self):
        return float(self.other_df.loc["rho", "B_P"])
    
    @property
    def rho_NGS(self):
        return float(self.other_df.loc["rho", "NGS"])

    @property
    def rho_H2S(self):
        return float(self.other_df.loc["rho", "H2S"])

    @property
    def rho_CO2S(self):
        return float(self.other_df.loc["rho", "CO2S"])

    @property
    def mu_CHP(self):
        return float(self.other_df.loc["mu", "CHP"])
    
    @property
    def mu_BM(self):
        return float(self.other_df.loc["mu", "BM"])
    
    @property
    def mu_WS(self):
        return float(self.other_df.loc["mu", "WS"])

    @property
    def mu_NK(self):
        return float(self.other_df.loc["mu", "NK"])
    
    @property
    def mu_NG_I(self):
        return float(self.other_df.loc["mu", "I_NG"])

    @property
    def mu_E_I(self):
        return float(self.other_df.loc["mu","I_E"])
    
    @property
    def mu_H2_I(self):
        return float(self.other_df.loc["mu","I_H2"])

    @property
    def mu_EL(self):
        return float(self.other_df.loc["mu","EL"])

    @property
    def mu_MT(self):
        return float(self.other_df.loc["mu","MT"])

    @property
    def nu_NG_CO2(self):
        return float(self.other_df.loc["nu", "CCGT"])
    
    @property
    def nu_BM_CO2(self):
        return float(self.other_df.loc["nu", "BM"])
    
    @property
    def nu_WS_CO2(self):
        return float(self.other_df.loc["nu", "WS"])
    
    @property
    def nu_E_I_CO2(self):
        return float(self.other_df.loc["nu", "I_E"])
    
    @property
    def phi_NG_CCS(self):
        return float(self.other_df.loc["phi", "CCGT"])
    
    @property
    def phi_CHP_CCS(self):
        return float(self.other_df.loc["phi", "CHP"])
    
    @property
    def phi_BM_CCS(self):
        return float(self.other_df.loc["phi", "BM"])
    
    @property
    def phi_WS_CCS(self):
        return float(self.other_df.loc["phi", "WS"])
    
    @property
    def phi_SMR_CCS(self):
        return float(self.other_df.loc["phi", "SMR_CCS"])
    
    @property
    def phi_SMR(self):
        return float(self.other_df.loc["phi", "SMR"])

    @property
    def phi_NG_ACC(self):
        return float(self.other_df.loc["phi", "ACC_NG"])

    @property
    def phi_E_ACC(self):
        return float(self.other_df.loc["phi", "ACC_E"])
    
    @property
    def MM_H2(self):
        return float(self.other_df.loc["MM", "H2"])
    
    @property
    def MM_H2O(self):
        return float(self.other_df.loc["MM", "H2O"])
    
    @property
    def MM_O2(self):
        return float(self.other_df.loc["MM", "O2"])
    
    @property
    def MM_CO2(self):
        return float(self.other_df.loc["MM", "CO2"])
    
    @property
    def MM_CH4(self):
        return float(self.other_df.loc["MM", "CH4"])
    
    @property
    def LHV_CH4(self):
        return float(self.other_df.loc["LHV", "CH4"])
    
    @property
    def LHV_H2(self):
        return float(self.other_df.loc["LHV", "H2"])
    
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
    df_list = [pd.read_excel(path_files[index],sheet_name=0, index_col=0) for index in range(len(path_files))]
    if len(df_list)==1:
        return df_list[0]

def fetch_elia_load_data(path_datafiles, data_years_input, year_no, data_years_output, n_files_yearly = 12):
    
    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)

    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.xls')]
    
    if len(year_no) != data_years_input:
        raise InputError('Length of input sequence should match number of input years!')
        
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

#    return x_r.to_dict() # needed to keep output as pandas dataframe for lambda_E property to work
    return x_r


def fetch_elia_generation_data(path_datafiles, data_years_input, year_no, data_years_output, n_files_yearly=12):
    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)

    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.xls')]

    if len(year_no) != data_years_input:
        raise InputError('Length of input sequence should match number of input years!')

    file_sequence = []
    for i in arange(0, len(year_no), 1):
        files_temp = files[(year_no[i] - 1) * n_files_yearly:year_no[i] * n_files_yearly]
        file_sequence.extend(files_temp)

    df = pd.DataFrame()

    if files[0].startswith('Solar'):

        for f in file_sequence:
            data = pd.read_excel(path_datafiles + str(f), index_col=0, skiprows=[0, 1, 2], usecols=[0, 1, 5, 6])
            data.columns = ['forecast', 'prod', 'peak']
            data['prod'].fillna(data['forecast'], inplace=True)
            data.index = pd.to_datetime(data.index, dayfirst=True)
            data = data.resample('H').mean()
            data.fillna(method='ffill', inplace=True)
            df = pd.concat([df, data], axis=0)

    else:

        for f in file_sequence:
            if (('2019' in f) and ('2018' in f)) or (('2018' in f) and not ('2017' in f)):
                data = pd.read_excel(path_datafiles + str(f), index_col=0, skiprows=[0, 1, 2], usecols=[0, 3, 4, 5])
            else:
                data = pd.read_excel(path_datafiles + str(f), index_col=0, skiprows=[0, 1, 2], usecols=[0, 3, 4, 6])
            data.columns = ['forecast', 'prod', 'peak']
            data['prod'].fillna(data['forecast'], inplace=True)
            data.index = pd.to_datetime(data.index, dayfirst=True)
            data = data.resample('H').mean()
            data.fillna(method='ffill', inplace=True)
            df = pd.concat([df, data], axis=0)

    df['prod_n'] = df['prod'] / df['peak']
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


def fetch_fluxys_demand_data(path_datafiles, data_years_input, year_no, data_years_output, client_type, n_files_yearly = 4):

    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)

    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.csv')]

    if len(year_no) != data_years_input:
        raise InputError('Length of input sequence should match number of input years!')

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

    if client_type == 'residential':
        df = df[df['clientType'] == 'Distribution Domestic Exit Point']
    elif client_type == 'industry':
        df = df[df['clientType'] == 'End User Domestic Exit Point IC']

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



#def fetch_fluxys_demand_data(path_datafiles, data_years_input, year_no, data_years_output, client_type, n_files_yearly = 4):
#
#    data_years_output = int(data_years_output)
#    data_years_input = int(data_years_input)
#
#    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.csv')]
#
#    if len(year_no) != data_years_input:
#        raise InputError('Length of input sequence should match number of input years!')
#
#    file_sequence = []
#    for i in arange(0, len(year_no), 1):
#        files_temp = files[(year_no[i]-1)*n_files_yearly:year_no[i]*n_files_yearly]
#        file_sequence.extend(files_temp)
#
#    df = pd.DataFrame()
#
#    for f in file_sequence:
#        data = pd.read_csv(path_datafiles+'/'+str(f), sep=';', usecols=['subGrid', 'gasDay', 'gasHour', 'clientType', 'physicalFlow'])
#        if f == file_sequence[0]:
#            d_start_str = str(data['gasDay'].iloc[0])+' '+str(int(data['gasHour'].iloc[0])+5)
#            d_start = pd.to_datetime(d_start_str, format='%d/%m/%Y %H')
#        else:
#            d_start = df.index[-1] + pd.Timedelta(hours=1)
#
#        for item in data['clientType'].unique():
#            data_client = data[data['clientType'] == item]
#            for net in data_client['subGrid'].unique():
#                data_network = data_client[data_client['subGrid'] == net]
#                data_network.index = pd.date_range(d_start, freq='H', periods=data_network.shape[0])
#                df = pd.concat([df, data_network], axis=0, ignore_index=False)
#
#    df = df[~((df.index.month == 2) & (df.index.day == 29))]
#    df = df[df['clientType'] != 'End User Domestic Exit Point PP']
#
#    if client_type == 'residential':
#        df = df[df['clientType'] == 'Distribution Domestic Exit Point']
#    elif client_type == 'industry':
#        df = df[df['clientType'] == 'End User Domestic Exit Point IC']
#
#    df_sum = df.groupby(df.index)['physicalFlow'].sum().reset_index(drop=True)
##    df_sum = df_sum * (-1) * 1e-3
#    df_sum = df_sum * (-1) * 1e-6
#
#    # ADJUST TIME SERIES FOR THE GAS DAY THAT STARTS AT 6 AM
#    data_shift = df_sum[:6]
#    df_shifted = df_sum.shift(-6)
#    df_shifted[-6:] = data_shift
#
#    df_shifted = df_shifted.astype(float)
#
#    if df_shifted.isnull().values.any():
#        raise InputError('You have some NaN in your timeseries and Gurobi will complain.')
#
#    x = split(df_shifted, data_years_input, axis=0)
#    x_r = shuffle_ts(x, data_years_input, data_years_output)
#
#    return x_r.to_dict()


def fetch_electric_heating_data(path_datafiles, data_years_output, COP):

    data_years_output = int(data_years_output)

    data = pd.read_excel(path_datafiles + 'heating_profile_electric/profile.xlsx', index_col=0)
    data_el = data['value']/COP
    data_el = pd.concat([data_el]*data_years_output, axis=0, ignore_index=True)

    return data_el.to_dict()

def fetch_electric_transport_data(path_datafiles, data_years_input, year_no, data_years_output, t_max, delta, n_files_yearly = 1):

    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)

    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.xlsx')]

    if len(year_no) != data_years_input:
        raise InputError('Length of input sequence should match number of input years!')

    file_sequence = []
    for i in arange(0, len(year_no), 1):
        files_temp = files[(year_no[i]-1)*n_files_yearly:year_no[i]*n_files_yearly]
        file_sequence.extend(files_temp)

    df = pd.DataFrame()

    for f in file_sequence:
        data = pd.read_excel(path_datafiles+'/'+str(f))
        df = pd.concat([df, data], axis=0, ignore_index=True)

    df = df['value']
    idx = arange(0, t_max, delta)

    x = split(df, data_years_input, axis=0)
    x_r = shuffle_ts(x, data_years_input, data_years_output)

    return dict(zip(idx, x_r))
    

def build_timeseries(y_list, p_list, data_years_output, start_year = 2019):
    
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

def build_peakload_timeseries(data_years_output, base_peak, growth_rate_E):
    
        data_years_output = int(data_years_output)
        
        x = pd.Series(index=arange(1, data_years_output+1, 1))
        for i in range(1, data_years_output+1):
            x[i] = round(base_peak * (1+growth_rate_E)**i, 1)

        ts_cap = pd.Series()

        for idx in x.index:
            x_r = pd.Series(repeat(x[idx], 8760))
            ts_cap = pd.concat([ts_cap, x_r], axis=0, ignore_index=True)

        if ts_cap.isnull().values.any():
            raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

#        return ts_cap.to_dict()
        return ts_cap

def fetch_elix_index_data(path_datafiles, data_years_input, year_no, data_years_output, units = 'MEuroGWh'):
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
    
    if units == 'MEuroMWh':
        df = df * 1e-6
    elif units == 'MEuroGWh':
        df = df * 1e-3
        
    if df.isnull().values.any():
        raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

    x = split(df['price'], data_years_input, axis=0)
    x_r = shuffle_ts(x, data_years_input, data_years_output)

    return x_r.to_dict()

def fetch_gas_price_data(path_datafiles, data_years_input, year_no, data_years_output, units='MEuroGWh'):
    # DATA ACQUISITION HERE MIGHT BE AN ISSUE...

    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)

    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.xlsx')]

    df = pd.DataFrame()
    df_sequence = pd.DataFrame()

    for f in files:
        ts = pd.read_excel(path_datafiles + '/' + str(f), index_col=0)
        df = pd.concat([df, ts], axis=0, ignore_index=True)

    # orders data as specified by order of years in year_no list?
    for i in arange(0, len(year_no), 1):
        df_temp = df.loc[(year_no[i] - 1) * 8760:(year_no[i] * 8760 - 1)]
        df_sequence = pd.concat([df_sequence, df_temp], axis=0, ignore_index=False)

    df = df_sequence
    df = df.append(df.tail(1))

    if units == 'MEuroMWh':
        df = df * 1e-6
    elif units == 'MEuroGWh':
        df = df * 1e-3

    if df.isnull().values.any():
        raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

    start_date = datetime(2017, 10, 1)
    end_date = datetime(2018, 10, 1)
    date_index = pd.date_range(start_date, end_date, freq='D')

    df.index = date_index
    df = df.resample('H').pad()
    df = df[:-1]

    df['day'] = df.index.day
    df['month'] = df.index.month
    df = df.sort_values(['month', 'day'], ascending=[True, True])

    # shuffles randomly?
    x = split(df['price'], data_years_input, axis=0)
    x_r = shuffle_ts(x, data_years_input, data_years_output)

    return x_r.to_dict()


def fetch_gas_transit_data(path_datafiles, data_years_input, year_no, data_years_output):
    
    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)
    
    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.csv')]

    df = pd.DataFrame()
    df_sequence = pd.DataFrame()
    
    for f in files:
        ts = pd.read_csv(path_datafiles+'/'+str(f), sep=';', index_col=0)
        df = pd.concat([df, ts], axis=0, ignore_index=True)
       
    for i in arange(0, len(year_no), 1):
        df_temp = df.loc[(year_no[i]-1)*8760:(year_no[i]*8760-1)]
        df_sequence = pd.concat([df_sequence, df_temp], axis=0, ignore_index=False)
    
    df = df_sequence

    if df.isnull().values.any():
        raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

    x = split(df['volume'], data_years_input, axis=0)
    x_r = shuffle_ts(x, data_years_input, data_years_output)

    return x_r.to_dict()



def fetch_gas_markets_data(path_datafiles, data_years_input, year_no, data_years_output, vol_multiple, price_multiple):

    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)
    
    files = [f for f in sorted(os.listdir(path_datafiles)) if f.endswith('.csv')]

    df = pd.DataFrame()
    df_sequence = pd.DataFrame()
    
    for f in files:
        ts = pd.read_csv(path_datafiles+'/'+str(f), sep=';', index_col=0)
        df = pd.concat([df, ts], axis=0, ignore_index=True)
    
    for i in arange(0, len(year_no), 1):
        df_temp = df.loc[(year_no[i]-1)*8760:(year_no[i]*8760-1)]
        df_sequence = pd.concat([df_sequence, df_temp], axis=0, ignore_index=False)
    
    df = df_sequence
    df['volume'] = df['volume']*vol_multiple
    df['price'] = df['price']*price_multiple
      
    if df.isnull().values.any():
        raise InputError('You have some NaN in your timeseries and Gurobi will complain.')

    x_volume = split(df['volume'], data_years_input, axis=0)
    x_r_volume = shuffle_ts(x_volume, data_years_input, data_years_output)
    
    x_price = split(df['price'], data_years_input, axis=0)
    x_r_price = shuffle_ts(x_price, data_years_input, data_years_output)

    return x_r_volume.to_dict(), x_r_price.to_dict()


def build_hydrogen_supply_profile(data_years_input, data_years_output, tanker_size, tankers_per_week, hours_to_discharge, start_year):

    data_years_output = int(data_years_output)
    data_years_input = int(data_years_input)
    date_range = pd.date_range(datetime(start_year, 1, 1, 0, 0, 0), freq='H', periods=data_years_input*8760)
    h_output = tanker_size / hours_to_discharge

    data = pd.DataFrame(index=date_range)

    if tankers_per_week == 1:
        data['value'] = where(data.index.dayofweek == 0, h_output, 0.0)
    elif tankers_per_week == 2:
        data['value'] = where(data.index.dayofweek == 0, h_output, where(data.index.dayofweek == 3, h_output, 0.))
    elif tankers_per_week == 3:
        data['value'] = where(data.index.dayofweek == 0, h_output, where(data.index.dayofweek == 3, h_output, where(data.index.dayofweek == 6, h_output, 0.)))
    elif tankers_per_week == 4:
        data['value'] = where(data.index.dayofweek == 0, h_output, where(data.index.dayofweek == 2, h_output, where(data.index.dayofweek == 4, h_output, where(data.index.dayofweek == 6, h_output, 0.))))
    else:
        raise InputError('Maximum four charges per week')

    data.index = arange(data_years_input*8760)

    x = split(data['value'], data_years_input, axis=0)
    x_r = shuffle_ts(x, data_years_input, data_years_output)

    return x_r.to_dict()



def capex_annuity(capex, lifetime, no_years):
    
    x = capex*(no_years/lifetime)
    return x


def fetch_timeseries_data(path_datafiles, output_dict=False):

    files = [f for f in os.listdir(path_datafiles) if f.endswith('.xlsx')]
    for f in files:
        df = pd.read_excel(path_datafiles+str(f))
    if output_dict:
        return dict(zip(list(df["time"]), list(df["value"])))
    else:
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
    
    data = pd.read_excel(path_datafiles + '/time.xlsx', index_col=0)
    
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
