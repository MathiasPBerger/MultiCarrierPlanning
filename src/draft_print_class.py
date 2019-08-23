# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 15:49:17 2018

@author: david
"""

import pandas as pd
import numpy as np
from os import makedirs
from os.path import join
from time import strftime
from pyomo.environ import Param, Var, value
from draft_postprocess_class import PostProcess

class Print(object):
    
    def __init__(self, model, data):
        
        self.model = model
        self.pp = PostProcess(data, model)
	
    def print_output(self, folder):

        df_ts = pd.DataFrame(index = np.arange(0, self.model.t_max()-1))
        df_cap = pd.DataFrame()
        df_pp = pd.DataFrame()

        for v in self.model.component_objects(Var, active=True):
            varobject = getattr(self.model, str(v))
    
            if str(v) in ['P_W_on', 'P_W_off', 'P_S',
                          'P_NG_OCGT', 'P_NG_CCGT', 'P_OCGT', 'P_CCGT', 'P_E_OCGT', 'P_E_CCGT',
                          'P_E_NK', 'P_CHP', 'P_E_CHP', 'P_NG_CHP', 'P_BM', 'P_E_BM', 'P_WS', 'P_E_WS',
                          'P_E_EL', 'P_H2_EL', 'P_H2_FC', 'P_E_FC', 'P_H2_MT', 'P_CH4_MT',
                          'E_H2', 'P_StH2', 'P_H2tS', 'P_H2S',
                          'E_NGS', 'P_NGS', 'P_NGStNG', 'P_NGtNGS',
                          'P_B', 'P_PtB', 'P_BtP', 'E_B',
                          'P_PtPH', 'P_PHtP', 'P_PH', 'E_PH',
                          'P_OCGT_CCS', 'P_CCGT_CCS', 'P_CHP_CCS', 'P_BM_CCS', 'P_WS_CCS', 'P_SMR_CCS',
                          'P_E_ACC', 'P_NG_ACC',
                          'Q_CO2_OCGT_CCS', 'Q_CO2_OCGT', 'Q_CO2_CCGT_CCS', 'Q_CO2_CCGT', 'Q_CO2_CHP_CCS', 'Q_CO2_CHP',
                          'Q_CO2_BM_CCS', 'Q_CO2_BM', 'Q_CO2_WS_CCS', 'Q_CO2_WS', 'Q_CO2_MT', 'Q_CO2_StCH4', 'Q_CO2_SP',
                          'Q_CO2_SMR', 'Q_CO2_SMR_CCS', 'Q_CO2_ACC', 'Q_CO2_E', 'Q_CO2S', 'Q_CO2S_StG', 'Q_CO2S_GtS', 'M_CO2',
                          'P_E_SMR', 'P_H2_SMR', 'P_NG_SMR',
                          'P_IE', 'P_E_I', 'P_E_E', 'P_NG_I', 'P_H2_I',
                          'L_E_ENS', 'L_E_TR', 'L_NG_ENS', 'L_H2_ENS',
                          'm_O2', 'm_H2O']:
                for index in varobject:
                    df_ts.loc[index, str(v)] = varobject[index].value

        for p in self.model.component_objects(Param, active=True):
            paramobject = getattr(self.model, str(p))
    
            if str(p) in ['lambda_E', 'lambda_E_HT', 'lambda_E_TR', 'lambda_NG_HT', 'lambda_NG_ID', 'lambda_NGtH2', 'lambda_NG_TR', 'lambda_H2_TR', 'lambda_H2_ID', 'theta_IE', 'kappa_H2_I']:
                for index in paramobject:
                    df_ts.loc[index, str(p)] = paramobject[index]   
                
                
        for v in self.model.component_objects(Var, active=True):
            varobject = getattr(self.model, str(v))
    
            if str(v) in ['K_W_on', 'K_W_off', 'K_S',
                          'K_E_EL', 'S_H2', 'K_E_FC', 'K_CH4_MT',
                          'K_E_OCGT', 'K_E_CCGT', 'S_B', 'K_B', 'K_CO2_OCGT_CCS', 'K_CO2_CCGT_CCS', 'K_CO2_CHP_CCS',
                          'K_CO2_WS_CCS', 'K_CO2_BM_CCS', 'S_CO2', 'K_H2_SMR', 'K_CO2_SMR_CCS', 'K_CO2_ACC']:
                for index in varobject:
                    df_cap.loc[str(v), 'value'] = varobject[index].value

        for p in self.model.component_objects(Param, active=True):
            
            if str(p) not in ['gamma_S', 'gamma_W_on', 'gamma_W_off',
                              'lambda_E', 'lambda_E_HT', 'lambda_E_TR', 'lambda_NG_HT', 'lambda_NG_ID', 'lambda_NGtH2', 'kappa_H2_I',
                              'lambda_NG_TR', 'lambda_H2_TR', 'lambda_H2_ID', 'theta_IE', 'theta_NG_I', 'theta_H2_IE']:
                df_cap.loc[str(p), 'value'] = value(p)
                
        pp_attr = [str(strr) for strr in list(PostProcess.__dict__.keys()) if "_" not in strr[:1]]
        for p in pp_attr:
            df_pp.loc[p, 'value'] = getattr(self.pp, p)

        df_ts['dual_E_balance'] = np.asarray([self.model.dual[self.model.power_balance[t]] for t in self.model.T])
        df_ts['dual_H2_balance'] = np.asarray([self.model.dual[self.model.H2_energy_flows_balance[t]] for t in self.model.T])
        df_ts['dual_NG_balance'] = np.asarray([self.model.dual[self.model.NG_network_balance[t]] for t in self.model.T])

        df_cap.loc['dual_CO2', 'value'] = self.model.dual[self.model.yearly_CO2_budget]

        df_ts.round(4)
        df_cap.round(4)
        df_pp.round(4)
        
        makedirs(folder+"/files")
        dir_files = folder+"/files"

        df_ts.to_csv(join(dir_files, 'output_timeseries.csv'), sep=';')
        df_cap.to_csv(join(dir_files, 'output_capacities.csv'), sep=';')
        df_pp.to_csv(join(dir_files, 'output_postprocessing.csv'), sep=';')
        
        return None

