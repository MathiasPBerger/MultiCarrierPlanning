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
    
            if str(v) in ['P_W_on', 'P_W_off', 'P_S', 'P_C', 
                  'Delta_P', 'P_PtG', 'P_PtH2', 'E_H2', 'P_H2_out', 'P_H2tP', 
                  'P_H2tCH4', 'P_H2', 'P_CH4', 'E_CH4', 'P_CH4tNG', 'P_NG_PP', 'P_NGtP', 'P_NG', 
                  'P_NK', 'P_disp', 'P_PtPH', 'P_PHtP', 'P_PH', 'E_PH', 'P_trs', "P_B", "E_B", "P_PtB", "P_BtP"]:
                for index in varobject:
                    df_ts.loc[index, str(v)] = varobject[index].value

        for p in self.model.component_objects(Param, active=True):
            paramobject = getattr(self.model, str(p))
    
            if str(p) in ['pi_L', 'pi_NG', 'kappa_NK', 'kappa_NG_0', 'theta_el']:
                for index in paramobject:
                    df_ts.loc[index, str(p)] = paramobject[index]   
                
                
        for v in self.model.component_objects(Var, active=True):
            varobject = getattr(self.model, str(v))
    
            if str(v) in ['K_W_on', 'K_W_off', 'K_S', 'K_PtG', 'S_H2', 'K_H2', 'K_H2tCH4', 
                          'S_CH4', 'K_NG', 'S_B', 'K_B']:
                for index in varobject:
                    df_cap.loc[str(v), 'value'] = varobject[index].value

        for p in self.model.component_objects(Param, active=True):
            
            if str(p) not in ['gamma_S', 'gamma_W_on', 'gamma_W_off', 'gamma_L', 'pi_L', 'pi_NG',
                              'kappa_NK', 'kappa_NG_0', 'theta_el']:
                df_cap.loc[str(p), 'value'] = value(p)
                
        pp_attr = [str(strr) for strr in list(PostProcess.__dict__.keys()) if "_" not in strr[:1]]
        for p in pp_attr:
            df_pp.loc[p, 'value'] = getattr(self.pp, p)
        
        df_ts.round(1)
        df_cap.round(1)
        
        makedirs(folder+"/files")
        dir_files = folder+"/files"
        
        df_ts.to_csv(join(dir_files, 'output_timeseries.csv'), sep=';')
        df_cap.to_csv(join(dir_files, 'output_capacities.csv'), sep=';')
        df_pp.to_csv(join(dir_files, 'output_postprocessing.csv'), sep=';')
        
        return None
         
