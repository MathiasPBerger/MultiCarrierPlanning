#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 13:47:51 2018

@author: mathiasberger
"""

from pyomo.environ import ConcreteModel, Param, RangeSet, Var, Constraint, Objective, NonNegativeReals, Reals, minimize, Suffix
from pyomo.opt import ProblemFormat
from os.path import join

def model_builder(data, folder, write_lp=False):

    model = ConcreteModel()

    ### Parameters Definition

    # Time Parameters
    
    model.t_max = Param(initialize=data.t_max)
    model.delta_t = Param(initialize=data.delta_t)
    model.n_delta_t = Param(initialize=data.n_delta_t)
    model.T_d = RangeSet(0, model.t_max-1, model.n_delta_t)
    model.T = RangeSet(0, model.t_max-1)
    model.n_y = Param(initialize=data.n_y)

    # Load factor parameters
    
    model.gamma_S = Param(model.T, initialize=data.gamma_S)
    model.gamma_W_on = Param(model.T, initialize=data.gamma_W_on)
    model.gamma_W_off = Param(model.T, initialize=data.gamma_W_off)

    # Capacity parameters
    
    model.kappa_W_on_max = Param(initialize=data.kappa_W_on_max)
    model.kappa_W_on_0 = Param(initialize=data.kappa_W_on_0)
    
    model.kappa_W_off_max = Param(initialize=data.kappa_W_off_max)
    model.kappa_W_off_0 = Param(initialize=data.kappa_W_off_0)
    
    model.kappa_S_max = Param(initialize=data.kappa_S_max)
    model.kappa_S_0 = Param(initialize=data.kappa_S_0)

    model.kappa_OCGT_0 = Param(initialize=data.kappa_OCGT_0)
    model.kappa_CCGT_0 = Param(initialize=data.kappa_CCGT_0)
    model.kappa_NK = Param(initialize=data.kappa_NK)
    model.kappa_CHP = Param(initialize=data.kappa_CHP)
    model.kappa_BM = Param(initialize=data.kappa_BM)
    model.kappa_WS = Param(initialize=data.kappa_WS)
    
    model.kappa_EL = Param(initialize=data.kappa_EL)
    model.kappa_FC = Param(initialize=data.kappa_FC)
    model.kappa_MT = Param(initialize=data.kappa_MT)
    model.kappa_SMR = Param(initialize=data.kappa_SMR)
    
    model.kappa_NGS = Param(initialize=data.kappa_NGS)
    model.xi_NGS = Param(initialize=data.xi_NGS)
    
    model.xi_H2S = Param(initialize=data.xi_H2S)
    
    model.xi_PH = Param(initialize=data.xi_PH)
    model.kappa_PtPH = Param(initialize=data.kappa_PtPH)
    model.kappa_PHtP = Param(initialize=data.kappa_PHtP)

    model.xi_B = Param(initialize=data.xi_B)
    model.kappa_B = Param(initialize=data.kappa_B)
    
    model.xi_CO2S = Param(initialize=data.xi_CO2S)
    
    model.kappa_NGNet = Param(initialize=data.kappa_NGNet)
    model.kappa_IE = Param(initialize=data.kappa_IE)
    model.kappa_H2_I = Param(model.T, initialize=data.kappa_H2_I)
    model.kappa_CO2_E = Param(initialize=data.kappa_CO2_E)

    model.kappa_CO2_OCGT_CCS = Param(initialize=data.kappa_CO2_OCGT_CCS)
    model.kappa_CO2_CCGT_CCS = Param(initialize=data.kappa_CO2_CCGT_CCS)
    model.kappa_CO2_CHP_CCS = Param(initialize=data.kappa_CO2_CHP_CCS)
    model.kappa_CO2_BM_CCS = Param(initialize=data.kappa_CO2_BM_CCS)
    model.kappa_CO2_WS_CCS = Param(initialize=data.kappa_CO2_WS_CCS)
    model.kappa_CO2_SMR_CCS = Param(initialize=data.kappa_CO2_SMR_CCS)
    model.kappa_CO2_ACC = Param(initialize=data.kappa_CO2_ACC)
    
    model.kappa_L_E_ENS = Param(initialize=data.kappa_L_E_ENS)
    model.kappa_L_NG_ENS = Param(initialize=data.kappa_L_NG_ENS)
    model.kappa_L_H2_ENS = Param(initialize=data.kappa_L_H2_ENS)
    
    model.psi_CO2 = Param(initialize=data.psi_CO2)
    model.psi_E = Param(initialize=data.mu_E_I * data.kappa_E)
    model.psi_NG = Param(initialize=data.mu_NG_I * data.kappa_NG_I * data.n_y)
    model.psi_H2 = Param(initialize=data.mu_H2_I * data.kappa_H2_I_size * data.kappa_H2_I_freq * data.n_weeks)

    # Demand & Curtailment Parameters
     
    model.lambda_E = Param(model.T, initialize=data.lambda_E) # tertiary, industry, residential and railway electricity demand
    model.lambda_E_HT = Param(model.T, initialize=data.lambda_E_HT) # electricity heating demand
    model.lambda_E_TR = Param(model.T_d, initialize=data.lambda_E_TR) # electricity demand for EV charging
    model.lambda_NG_HT = Param(model.T, initialize=data.lambda_NG_HT) # heating demand for natural gas
    model.lambda_NG_ID = Param(model.T, initialize=data.lambda_NG_ID) # industry demand for natural
    model.lambda_NG_TR = Param(model.T, initialize=data.lambda_NG_TR) # natural gas transportation demand
    model.lambda_NGtH2 = Param(model.T, initialize=data.lambda_NGtH2) # assumed historical natural gas demand to feed SMRs
    model.lambda_H2_TR = Param(model.T, initialize=data.lambda_H2_TR) # hydrogen transportation demand
    model.lambda_H2_ID = Param(model.T, initialize=data.lambda_H2_ID) # hydrogen industry demand

    # Efficiency parameters
    
    model.eta_OCGT = Param(initialize=data.eta_OCGT)
    model.eta_CCGT = Param(initialize=data.eta_CCGT)
    
    model.eta_NK = Param(initialize=data.eta_NK)
    model.eta_CHP = Param(initialize=data.eta_CHP)
    model.eta_BM = Param(initialize=data.eta_BM)
    model.eta_WS = Param(initialize=data.eta_WS)
    
    model.eta_FC = Param(initialize=data.eta_FC)
    model.eta_EL = Param(initialize=data.eta_EL)
    model.eta_MT = Param(initialize=data.eta_MT)
    model.eta_SMR = Param(initialize=data.eta_SMR)
    
    model.eta_BtP = Param(initialize=data.eta_BtP)
    model.eta_PtB = Param(initialize=data.eta_PtB)
    model.eta_B = Param(initialize=data.eta_B)
    
    model.eta_PHtP = Param(initialize=data.eta_PHtP)
    model.eta_PtPH = Param(initialize=data.eta_PtPH)

    model.eta_NGS = Param(initialize=data.eta_NGS)
    
    model.eta_NG_CCS = Param(initialize=data.eta_NG_CCS)
    model.eta_CHP_CCS = Param(initialize=data.eta_CHP_CCS)
    model.eta_BM_CCS = Param(initialize=data.eta_BM_CCS)
    model.eta_WS_CCS = Param(initialize=data.eta_WS_CCS)
    model.eta_SMR_CCS = Param(initialize=data.eta_SMR_CCS)

    # Cost parameters
    
    model.zeta_W_on = Param(initialize=data.zeta_W_on)
    model.theta_W_on_f = Param(initialize=data.theta_W_on_f)
    model.theta_W_on_v = Param(initialize=data.theta_W_on_v)
    
    model.zeta_W_off = Param(initialize=data.zeta_W_off)
    model.theta_W_off_f = Param(initialize=data.theta_W_off_f)
    model.theta_W_off_v = Param(initialize=data.theta_W_off_v)
    
    model.zeta_S = Param(initialize=data.zeta_S)
    model.theta_S_f = Param(initialize=data.theta_S_f)
    model.theta_S_v = Param(initialize=data.theta_S_v)
    
    model.zeta_EL = Param(initialize=data.zeta_EL)
    model.theta_EL_f = Param(initialize=data.theta_EL_f)
    
    model.zeta_H2S = Param(initialize=data.zeta_H2S)
    model.theta_H2S_f = Param(initialize=data.theta_H2S_f)
    
    model.zeta_FC = Param(initialize=data.zeta_FC)
    model.theta_FC_f = Param(initialize=data.theta_FC_f)
    model.theta_FC_v = Param(initialize=data.theta_FC_v)
    
    model.zeta_MT = Param(initialize=data.zeta_MT)
    model.theta_MT_f = Param(initialize=data.theta_MT_f)
    
    model.zeta_OCGT = Param(initialize=data.zeta_OCGT)
    model.theta_OCGT_f = Param(initialize=data.theta_OCGT_f)
    model.theta_OCGT_v = Param(initialize=data.theta_OCGT_v)
    
    model.zeta_CCGT = Param(initialize=data.zeta_CCGT)
    model.theta_CCGT_f = Param(initialize=data.theta_CCGT_f)
    model.theta_CCGT_v = Param(initialize=data.theta_CCGT_v)
    
    model.zeta_SMR = Param(initialize=data.zeta_SMR)
    model.theta_SMR_f = Param(initialize=data.theta_SMR_f)
    model.theta_SMR_v = Param(initialize=data.theta_SMR_v)
    
    model.zeta_CO2S = Param(initialize=data.zeta_CO2S)
    model.theta_CO2S_f = Param(initialize=data.theta_CO2S_f)
    model.theta_CO2S_v = Param(initialize=data.theta_CO2S_v)
    
    model.zeta_NG_CCS = Param(initialize=data.zeta_NG_CCS)
    model.theta_NG_CCS_f = Param(initialize=data.theta_NG_CCS_f)
    model.theta_NG_CCS_v = Param(initialize=data.theta_NG_CCS_v)
    
    model.zeta_CHP_CCS = Param(initialize=data.zeta_CHP_CCS)
    model.theta_CHP_CCS_f = Param(initialize=data.theta_CHP_CCS_f)
    model.theta_CHP_CCS_v = Param(initialize=data.theta_CHP_CCS_v)
    
    model.zeta_BM_CCS = Param(initialize=data.zeta_BM_CCS)
    model.theta_BM_CCS_f = Param(initialize=data.theta_BM_CCS_f)
    model.theta_BM_CCS_v = Param(initialize=data.theta_BM_CCS_v)
    
    model.zeta_WS_CCS = Param(initialize=data.zeta_WS_CCS)
    model.theta_WS_CCS_f = Param(initialize=data.theta_WS_CCS_f)
    model.theta_WS_CCS_v = Param(initialize=data.theta_WS_CCS_v)
    
    model.zeta_B_E = Param(initialize=data.zeta_B_E)
    model.zeta_B_P = Param(initialize=data.zeta_B_P)
    model.theta_B_E_f = Param(initialize=data.theta_B_E_f)
    model.theta_B_P_f = Param(initialize=data.theta_B_P_f)
    
    model.zeta_SMR_CCS = Param(initialize=data.zeta_SMR_CCS)
    model.theta_SMR_CCS_f = Param(initialize=data.theta_SMR_CCS_f)
    model.theta_SMR_CCS_v = Param(initialize=data.theta_SMR_CCS_v)

    model.zeta_ACC = Param(initialize=data.zeta_ACC)
    model.theta_ACC_f = Param(initialize=data.theta_ACC_f)
    
    model.theta_NK_v = Param(initialize=data.theta_NK_v)
    model.theta_CHP_v = Param(initialize=data.theta_CHP_v)
    model.theta_BM_v = Param(initialize=data.theta_BM_v)
    model.theta_WS_v = Param(initialize=data.theta_WS_v)
    model.theta_PH_v = Param(initialize=data.theta_PH_v)
    model.theta_NGS_f = Param(initialize=data.theta_NGS_f)
    
    model.theta_NK_fuel = Param(initialize=data.theta_NK_fuel)
    model.theta_BM_fuel = Param(initialize=data.theta_BM_fuel)
    model.theta_WS_fuel = Param(initialize=data.theta_WS_fuel)
    
    model.varsigma_ENS_E = Param(initialize=data.varsigma_ENS_E)
    model.varsigma_ENS_NG = Param(initialize=data.varsigma_ENS_NG)
    model.varsigma_ENS_H2 = Param(initialize=data.varsigma_ENS_H2)
    model.theta_CO2 = Param(initialize=data.theta_CO2)
    model.theta_IE = Param(model.T, initialize=data.theta_IE)
    model.theta_NG_I = Param(model.T, initialize=data.theta_NG_I)
    model.theta_H2_IE = Param(model.T, initialize=data.theta_H2_IE)
    model.theta_CO2_E = Param(initialize=data.theta_CO2_E)

    # Other parameters
    
    model.delta_CHP_inc = Param(initialize=data.delta_CHP_inc)
    model.delta_CHP_dec = Param(initialize=data.delta_CHP_dec)
    model.delta_BM_inc = Param(initialize=data.delta_BM_inc)
    model.delta_BM_dec = Param(initialize=data.delta_BM_dec)
    model.delta_WS_inc = Param(initialize=data.delta_WS_inc)
    model.delta_WS_dec = Param(initialize=data.delta_WS_dec)
    model.delta_NK_inc = Param(initialize=data.delta_NK_inc)
    model.delta_NK_dec = Param(initialize=data.delta_NK_dec)
    model.delta_EL_inc = Param(initialize=data.delta_EL_inc)
    model.delta_EL_dec = Param(initialize=data.delta_EL_dec)
    model.delta_MT_inc = Param(initialize=data.delta_MT_inc)
    model.delta_MT_dec = Param(initialize=data.delta_MT_dec)
    model.delta_OCGT_inc = Param(initialize=data.delta_OCGT_inc)
    model.delta_OCGT_dec = Param(initialize=data.delta_OCGT_dec)
    model.delta_CCGT_inc = Param(initialize=data.delta_CCGT_inc)
    model.delta_CCGT_dec = Param(initialize=data.delta_CCGT_dec)
    model.delta_FC_inc = Param(initialize=data.delta_FC_inc)
    model.delta_FC_dec = Param(initialize=data.delta_FC_dec)
    
    model.sigma_H2S = Param(initialize=data.sigma_H2S)
    model.sigma_PH = Param(initialize=data.sigma_PH)
    model.sigma_B = Param(initialize=data.sigma_B)
    model.sigma_NGS = Param(initialize=data.sigma_NGS)
    
    model.rho_B = Param(initialize=data.rho_B)
    model.rho_NGS = Param(initialize=data.rho_NGS)
    model.rho_H2S = Param(initialize=data.rho_H2S)
    model.rho_CO2S = Param(initialize=data.rho_CO2S)
    
    model.chi_H2S = Param(initialize=data.chi_H2S)
    model.chi_CO2S = Param(initialize=data.chi_CO2S)
    
    model.mu_CHP = Param(initialize=data.mu_CHP)
    model.mu_BM = Param(initialize=data.mu_BM)
    model.mu_WS = Param(initialize=data.mu_WS)
    model.mu_NK = Param(initialize=data.mu_NK)
    model.mu_EL = Param(initialize=data.mu_EL)
    model.mu_MT = Param(initialize=data.mu_MT)
    model.mu_E_I = Param(initialize=data.mu_E_I)
    model.mu_NG_I = Param(initialize=data.mu_NG_I)
    model.mu_H2_I = Param(initialize=data.mu_H2_I)
    
    model.nu_NG_CO2 = Param(initialize=data.nu_NG_CO2)
    model.nu_BM_CO2 = Param(initialize=data.nu_BM_CO2)
    model.nu_WS_CO2 = Param(initialize=data.nu_WS_CO2)
    model.nu_E_I_CO2 = Param(initialize=data.nu_E_I_CO2)
    
    model.phi_NG_CCS = Param(initialize=data.phi_NG_CCS)
    model.phi_CHP_CCS = Param(initialize=data.phi_CHP_CCS)
    model.phi_BM_CCS = Param(initialize=data.phi_BM_CCS)
    model.phi_WS_CCS = Param(initialize=data.phi_WS_CCS)
    model.phi_SMR_CCS = Param(initialize=data.phi_SMR_CCS)
    model.phi_SMR = Param(initialize=data.phi_SMR)
    model.phi_E_ACC = Param(initialize=data.phi_E_ACC)
    model.phi_NG_ACC = Param(initialize=data.phi_NG_ACC)
    
    model.MM_H2 = Param(initialize=data.MM_H2)
    model.MM_H2O = Param(initialize=data.MM_H2O)
    model.MM_O2 = Param(initialize=data.MM_O2)
    model.MM_CO2 = Param(initialize=data.MM_CO2)
    model.MM_CH4 = Param(initialize=data.MM_CH4)
    
    model.LHV_CH4 = Param(initialize=data.LHV_CH4)
    model.LHV_H2 = Param(initialize=data.LHV_H2)

    ### Variables Definition

    # Operational Variables
    
    model.P_W_on = Var(model.T, within=NonNegativeReals)
    model.P_W_off = Var(model.T, within=NonNegativeReals)
    model.P_S = Var(model.T, within=NonNegativeReals)
    model.P_NG_OCGT = Var(model.T, within=NonNegativeReals)
    model.P_NG_CCGT = Var(model.T, within=NonNegativeReals)
    model.P_OCGT = Var(model.T, within=NonNegativeReals)
    model.P_CCGT = Var(model.T, within=NonNegativeReals)
    model.P_E_OCGT = Var(model.T, within=NonNegativeReals)
    model.P_E_CCGT = Var(model.T, within=NonNegativeReals)
    model.P_E_NK = Var(model.T, within=NonNegativeReals)
    model.P_CHP = Var(model.T, within=NonNegativeReals)
    model.P_NG_CHP = Var(model.T, within=NonNegativeReals)
    model.P_BM = Var(model.T, within=NonNegativeReals)
    model.P_WS = Var(model.T, within=NonNegativeReals)
    model.P_E_CHP = Var(model.T, within=NonNegativeReals)
    model.P_E_BM = Var(model.T, within=NonNegativeReals)
    model.P_E_WS = Var(model.T, within=NonNegativeReals)
    model.P_E_SMR = Var(model.T, within=NonNegativeReals)
    model.P_H2_SMR = Var(model.T, within=NonNegativeReals)
    model.P_NG_SMR = Var(model.T, within=NonNegativeReals)
    
    model.P_E_EL = Var(model.T, within=NonNegativeReals)
    model.P_H2_EL = Var(model.T, within=NonNegativeReals)
    model.P_H2_FC = Var(model.T, within=NonNegativeReals)
    model.P_E_FC = Var(model.T, within=NonNegativeReals)
    model.P_H2_MT = Var(model.T, within=NonNegativeReals)
    model.P_CH4_MT = Var(model.T, within=NonNegativeReals)
    
    model.E_H2 = Var(model.T, within=NonNegativeReals)
    model.P_StH2 = Var(model.T, within=NonNegativeReals)
    model.P_H2tS = Var(model.T, within=NonNegativeReals)
    model.P_H2S = Var(model.T, within=Reals)
    
    model.E_NGS = Var(model.T, within=NonNegativeReals)
    model.P_NGS = Var(model.T, within=Reals)
    model.P_NGStNG = Var(model.T, within=NonNegativeReals)
    model.P_NGtNGS = Var(model.T, within=NonNegativeReals)
    
    model.P_B = Var(model.T, within=Reals)
    model.P_PtB = Var(model.T, within=NonNegativeReals)
    model.P_BtP = Var(model.T, within=NonNegativeReals)
    model.E_B = Var(model.T, within=NonNegativeReals)
    
    model.P_PtPH = Var(model.T, within=NonNegativeReals)
    model.P_PHtP = Var(model.T, within=NonNegativeReals)
    model.P_PH = Var(model.T, within=Reals)
    model.E_PH = Var(model.T, within=NonNegativeReals)
    
    model.P_OCGT_CCS = Var(model.T, within=NonNegativeReals)
    model.P_CCGT_CCS = Var(model.T, within=NonNegativeReals)
    model.P_CHP_CCS = Var(model.T, within=NonNegativeReals)
    model.P_BM_CCS = Var(model.T, within=NonNegativeReals)
    model.P_WS_CCS = Var(model.T, within=NonNegativeReals)
    model.P_SMR_CCS = Var(model.T, within=NonNegativeReals)
    
    model.P_E_ACC = Var(model.T, within=NonNegativeReals)
    model.P_NG_ACC = Var(model.T, within=NonNegativeReals)
    
    model.Q_CO2_OCGT_CCS = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_OCGT = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_CCGT_CCS = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_CCGT = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_CHP_CCS = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_CHP = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_BM_CCS = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_BM = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_WS_CCS = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_WS = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_MT = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_SMR = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_SMR_CCS = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_ACC = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_ACC_A = Var(model.T, within=NonNegativeReals)
    
    model.P_E_IE = Var(model.T, within=Reals)
    model.P_E_I = Var(model.T, within=NonNegativeReals)
    model.P_E_E = Var(model.T, within=NonNegativeReals)
    model.P_NG_I = Var(model.T, within=NonNegativeReals)
    model.P_H2_I = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_E = Var(model.T, within=NonNegativeReals)
    
    model.M_CO2 = Var(model.T, within=NonNegativeReals)
    model.Q_CO2S = Var(model.T, within=Reals)
    model.Q_CO2_StG = Var(model.T, within=NonNegativeReals)
    model.Q_CO2_GtS = Var(model.T, within=NonNegativeReals)

    model.L_E_TR = Var(model.T, within=NonNegativeReals)
    model.L_E_ENS = Var(model.T, within=NonNegativeReals)
    model.L_NG_ENS = Var(model.T, within=NonNegativeReals)
    model.L_H2_ENS = Var(model.T, within=NonNegativeReals)
    
    model.m_O2 = Var(model.T, within=NonNegativeReals)
    model.m_H2O = Var(model.T, within=NonNegativeReals)
    
    # Sizing Variables
    
    model.K_W_on = Var(within=NonNegativeReals)
    model.K_W_off = Var(within=NonNegativeReals)
    model.K_S = Var(within=NonNegativeReals)
    model.K_E_EL = Var(within=NonNegativeReals)
    model.S_H2 = Var(within=NonNegativeReals)
    model.K_H2S = Var(within=NonNegativeReals)
    model.K_E_FC = Var(within=NonNegativeReals)
    model.K_CH4_MT = Var(within=NonNegativeReals)
    model.K_E_OCGT = Var(within=NonNegativeReals)
    model.K_E_CCGT = Var(within=NonNegativeReals)
    model.K_CO2_OCGT_CCS = Var(within=NonNegativeReals)
    model.K_CO2_CCGT_CCS = Var(within=NonNegativeReals)
    model.K_CO2_CHP_CCS = Var(within=NonNegativeReals)
    model.K_CO2_BM_CCS = Var(within=NonNegativeReals)
    model.K_CO2_WS_CCS = Var(within=NonNegativeReals)
    model.S_B = Var(within=NonNegativeReals)
    model.K_B = Var(within=NonNegativeReals)
    model.S_CO2 = Var(within=NonNegativeReals)
    model.K_CO2S  = Var(within=NonNegativeReals)
    model.K_H2_SMR = Var(within=NonNegativeReals)
    model.K_CO2_SMR_CCS = Var(within=NonNegativeReals)
    model.K_CO2_ACC = Var(within=NonNegativeReals)


    ### Constraints Definition

    # Net Power Balance

    def power_balance_rule(model, t):
        return model.P_S[t] + model.P_W_on[t] + model.P_W_off[t] + model.P_PH[t] + model.P_E_FC[t] + model.P_E_OCGT[t] + model.P_E_CCGT[t] + model.P_E_CHP[t] +\
             model.P_E_BM[t] + model.P_E_WS[t] + model.P_E_NK[t] + model.P_B[t] + model.P_E_IE[t] + model.L_E_ENS[t] == model.lambda_E[t] + model.lambda_E_HT[t] + model.L_E_TR[t] + model.P_E_EL[t] + model.P_E_SMR[t] + model.P_SMR_CCS[t] + model.P_E_ACC[t]

    def electricity_load_shedding_upper_bound_rule(model, t):
        return model.L_E_ENS[t] <= model.kappa_L_E_ENS
    
    model.power_balance = Constraint(model.T, rule=power_balance_rule)
    model.electricity_load_shedding_upper_bound = Constraint(model.T, rule=electricity_load_shedding_upper_bound_rule)
    
    
    # EV Charging
    
    def ev_charging_rule(model, d):
        return sum(model.L_E_TR[d+t] for t in range(model.n_delta_t.value)) == model.lambda_E_TR[d]
    
    model.ev_charging = Constraint(model.T_d, rule=ev_charging_rule)


    # Variable Renewable Production Technologies Production & Sizing

    def solar_PV_power_output_definition_rule(model, t):
        return model.P_S[t] <= model.gamma_S[t] * (model.kappa_S_0 + model.K_S)

    def solar_PV_sizing_upper_bound_rule(model):
        return model.K_S <= model.kappa_S_max - model.kappa_S_0

    def wind_onshore_power_output_definition_rule(model, t):
        return model.P_W_on[t] <= model.gamma_W_on[t] * (model.kappa_W_on_0 + model.K_W_on)

    def wind_onshore_sizing_upper_bound_rule(model):
        return model.K_W_on <= model.kappa_W_on_max - model.kappa_W_on_0

    def wind_offshore_power_output_definition_rule(model, t):
        return model.P_W_off[t] <= model.gamma_W_off[t] * (model.kappa_W_off_0 + model.K_W_off)

    def wind_offshore_sizing_upper_bound_rule(model):
        return model.K_W_off <= model.kappa_W_off_max - model.kappa_W_off_0

    model.solar_PV_power_output_definition = Constraint(model.T, rule=solar_PV_power_output_definition_rule)
    model.solar_PV_sizing_upper_bound = Constraint(rule=solar_PV_sizing_upper_bound_rule)
    model.wind_onshore_power_output_definition = Constraint(model.T, rule=wind_onshore_power_output_definition_rule)
    model.wind_onshore_sizing_upper_bound = Constraint(rule=wind_onshore_sizing_upper_bound_rule)
    model.wind_offshore_power_output_definition = Constraint(model.T, rule=wind_offshore_power_output_definition_rule)
    model.wind_offshore_sizing_upper_bound = Constraint(rule=wind_offshore_sizing_upper_bound_rule)
 
    
    # Nuclear Power Plants

    def nuclear_power_lower_bound_rule(model, t):
        return model.mu_NK * model.kappa_NK <= model.P_E_NK[t]

    def nuclear_power_upper_bound_rule(model, t):
        return model.P_E_NK[t] <= model.kappa_NK

    def nuclear_units_inc_ramp_rate_rule(model, t):
        if t>0:
            return model.P_E_NK[t] <= model.P_E_NK[t - 1] + model.delta_NK_inc * model.kappa_NK
        else:
            return Constraint.Skip

    def nuclear_units_dec_ramp_rate_rule(model, t):
        if t>0:
            return model.P_E_NK[t] >= model.P_E_NK[t - 1] - model.delta_NK_dec * model.kappa_NK
        else:
            return Constraint.Skip
        
    model.nuclear_power_lower_bound = Constraint(model.T, rule=nuclear_power_lower_bound_rule)
    model.nuclear_power_upper_bound = Constraint(model.T, rule=nuclear_power_upper_bound_rule)
    model.nuclear_units_inc_ramp_rate = Constraint(model.T, rule=nuclear_units_inc_ramp_rate_rule)
    model.nuclear_units_dec_ramp_rate = Constraint(model.T, rule=nuclear_units_dec_ramp_rate_rule)


    # Other Dispatchable Power Plant (CHP)
    
    def CHP_power_gas_consumption_rule(model, t):
        return model.P_CHP[t] == model.eta_CHP * model.P_NG_CHP[t]

    def CHP_power_lower_bound_rule(model, t):
        return model.mu_CHP * model.kappa_CHP <= model.P_CHP[t]

    def CHP_power_upper_bound_rule(model, t):
        return model.P_CHP[t] <= model.kappa_CHP

    def CHP_units_inc_ramp_rate_rule(model, t):
        if t>0:
            return model.P_CHP[t] <= model.P_CHP[t-1] + model.delta_CHP_inc * model.kappa_CHP
        else:
            return Constraint.Skip

    def CHP_units_dec_ramp_rate_rule(model, t):
        if t>0:
            return model.P_CHP[t] >= model.P_CHP[t-1] - model.delta_CHP_dec * model.kappa_CHP
        else:
            return Constraint.Skip
        
    def CHP_plant_output_power_balance_rule(model, t):
        return model.P_CHP[t] == model.P_CHP_CCS[t] + model.P_E_CHP[t]

    model.CHP_power_gas_consumption = Constraint(model.T, rule=CHP_power_gas_consumption_rule)
    model.CHP_power_lower_bound = Constraint(model.T, rule=CHP_power_lower_bound_rule)
    model.CHP_power_upper_bound = Constraint(model.T, rule=CHP_power_upper_bound_rule)
    model.CHP_units_inc_ramp_rate = Constraint(model.T, rule=CHP_units_inc_ramp_rate_rule)
    model.CHP_units_dec_ramp_rate = Constraint(model.T, rule=CHP_units_dec_ramp_rate_rule)
    model.CHP_plant_output_power_balance = Constraint(model.T, rule=CHP_plant_output_power_balance_rule)
    
    
    # Other Dispatchable Power Plant (Biomass)

    def biomass_power_lower_bound_rule(model, t):
        return model.mu_BM * model.kappa_BM <= model.P_BM[t]

    def biomass_power_upper_bound_rule(model, t):
        return model.P_BM[t] <= model.kappa_BM

    def biomass_units_inc_ramp_rate_rule(model, t):
        if t>0:
            return model.P_BM[t] <= model.P_BM[t-1] + model.delta_BM_inc * model.kappa_BM
        else:
            return Constraint.Skip

    def biomass_units_dec_ramp_rate_rule(model, t):
        if t>0:
            return model.P_BM[t] >= model.P_BM[t-1] - model.delta_BM_dec * model.kappa_BM
        else:
            return Constraint.Skip
        
    def biomass_plant_output_power_balance_rule(model, t):
        return model.P_BM[t] == model.P_BM_CCS[t] + model.P_E_BM[t]

    model.biomass_power_lower_bound = Constraint(model.T, rule=biomass_power_lower_bound_rule)
    model.biomass_power_upper_bound = Constraint(model.T, rule=biomass_power_upper_bound_rule)
    model.biomass_units_inc_ramp_rate = Constraint(model.T, rule=biomass_units_inc_ramp_rate_rule)
    model.biomass_units_dec_ramp_rate = Constraint(model.T, rule=biomass_units_dec_ramp_rate_rule)
    model.biomass_plant_output_power_balance = Constraint(model.T, rule=biomass_plant_output_power_balance_rule)
    
    
    # Other Dispatchable Power Plant (Waste)

    def waste_power_lower_bound_rule(model, t):
        return model.mu_WS * model.kappa_WS <= model.P_WS[t]

    def waste_power_upper_bound_rule(model, t):
        return model.P_WS[t] <= model.kappa_WS

    def waste_units_inc_ramp_rate_rule(model, t):
        if t>0:
            return model.P_WS[t] <= model.P_WS[t-1] + model.delta_WS_inc * model.kappa_WS
        else:
            return Constraint.Skip

    def waste_units_dec_ramp_rate_rule(model, t):
        if t>0:
            return model.P_WS[t] >= model.P_WS[t-1] - model.delta_WS_dec * model.kappa_WS
        else:
            return Constraint.Skip
        
    def waste_plant_output_power_balance_rule(model, t):
        return model.P_WS[t] == model.P_WS_CCS[t] + model.P_E_WS[t]

    model.waste_power_lower_bound = Constraint(model.T, rule=waste_power_lower_bound_rule)
    model.waste_power_upper_bound = Constraint(model.T, rule=waste_power_upper_bound_rule)
    model.waste_units_inc_ramp_rate = Constraint(model.T, rule=waste_units_inc_ramp_rate_rule)
    model.waste_units_dec_ramp_rate = Constraint(model.T, rule=waste_units_dec_ramp_rate_rule)
    model.waste_plant_output_power_balance = Constraint(model.T, rule=waste_plant_output_power_balance_rule)
    
    
    # Electrical Interconnection
    
    def transmission_definition_rule(model, t):
        return model.P_E_IE[t] == model.P_E_I[t] - model.P_E_E[t]
    
    def imports_upper_bound_rule(model, t):
        return model.P_E_I[t] <= model.kappa_IE
    
    def exports_upper_bound_rule(model, t):
        return model.P_E_E[t] <= model.kappa_IE

    def transmission_budget_bound_rule(model):
        return sum(model.P_E_I[:]) <= model.psi_E
    
    model.transmission_definition = Constraint(model.T, rule=transmission_definition_rule)
    model.imports_upper_bound = Constraint(model.T, rule=imports_upper_bound_rule)
    model.exports_upper_bound = Constraint(model.T, rule=exports_upper_bound_rule)
    model.transmission_budget_bound = Constraint(rule=transmission_budget_bound_rule)


    # Pumped-Hydro Plant

    def net_pumped_hydro_power_definition_rule(model, t):
        return model.P_PH[t] == -model.P_PtPH[t]+ model.P_PHtP[t]

    def pumped_hydro_SOC_definition_rule(model, t):
        if t==0:
            return model.E_PH[t] == model.E_PH[model.t_max-1]
        else:
            return model.E_PH[t] == model.E_PH[t-1] + model.eta_PtPH * model.P_PtPH[t] * model.delta_t - model.P_PHtP[t] * model.delta_t / model.eta_PHtP

    def pumped_hydro_SOC_lower_bound_rule(model, t):
        return model.sigma_PH * model.xi_PH <= model.E_PH[t]

    def pumped_hydro_SOC_upper_bound_rule(model, t):
        return model.E_PH[t] <= model.xi_PH
        
    def pumped_hydro_input_power_upper_bound_rule(model, t):
       return model.P_PtPH[t] <= model.kappa_PtPH
       
    def pumped_hydro_output_power_upper_bound_rule(model, t):
        return model.P_PHtP[t] <= model.kappa_PHtP

    model.net_pumped_hydro_power_definition = Constraint(model.T, rule=net_pumped_hydro_power_definition_rule)
    model.pumped_hydro_SOC_definition = Constraint(model.T, rule=pumped_hydro_SOC_definition_rule)
    model.pumped_hydro_SOC_lower_bound = Constraint(model.T, rule=pumped_hydro_SOC_lower_bound_rule)
    model.pumped_hydro_SOC_upper_bound = Constraint(model.T, rule=pumped_hydro_SOC_upper_bound_rule)
    model.pumped_hydro_input_power_upper_bound = Constraint(model.T, rule=pumped_hydro_input_power_upper_bound_rule)
    model.pumped_hydro_output_power_upper_bound = Constraint(model.T, rule=pumped_hydro_output_power_upper_bound_rule)
    
    
    # Battery Plant

    def net_battery_power_definition_rule(model, t):
        return model.P_B[t] == -model.P_PtB[t] + model.P_BtP[t]

    def battery_SOC_definition_rule(model, t):
        if t == 0:
            return model.E_B[t] == model.E_B[model.t_max-1]
        else:
            return model.E_B[t] == model.eta_B * model.E_B[t-1] + model.eta_PtB * model.P_PtB[t] * model.delta_t - model.P_BtP[t] * model.delta_t / model.eta_BtP

    def battery_SOC_lower_bound_rule(model, t):
        return model.sigma_B * model.S_B <= model.E_B[t]

    def battery_SOC_upper_bound_rule(model, t):
        return model.E_B[t] <= model.S_B

    def battery_input_power_upper_bound_rule(model, t):
        return model.eta_PtB * model.P_PtB[t] <= model.rho_B * model.K_B
    
    def battery_output_power_upper_bound_rule(model, t):
        return model.P_BtP[t] <= model.K_B
    
    def battery_energy_capacity_upper_bound_rule(model):
        return model.S_B <= model.xi_B

    model.net_battery_power_definition = Constraint(model.T, rule=net_battery_power_definition_rule)
    model.battery_SOC_definition = Constraint(model.T, rule=battery_SOC_definition_rule)
    model.battery_SOC_lower_bound = Constraint(model.T, rule=battery_SOC_lower_bound_rule)
    model.battery_SOC_upper_bound = Constraint(model.T, rule=battery_SOC_upper_bound_rule)
    model.battery_input_power_upper_bound = Constraint(model.T, rule=battery_input_power_upper_bound_rule)
    model.battery_output_power_upper_bound = Constraint(model.T, rule=battery_output_power_upper_bound_rule)
    model.battery_energy_capacity_upper_bound = Constraint(rule=battery_energy_capacity_upper_bound_rule)


    # H2 Energy Flows Balance
    
    def H2_energy_flows_balance_rule(model, t):
        return model.P_H2_EL[t] + model.P_H2_SMR[t] + model.P_H2S[t] + model.L_H2_ENS[t] + model.P_H2_I[t] == model.P_H2_FC[t] + model.P_H2_MT[t] + model.lambda_H2_TR[t] + model.lambda_H2_ID[t]

    def H2_load_shedding_upper_bound_rule(model, t):
        return model.L_H2_ENS[t] <= model.kappa_L_H2_ENS
    
    model.H2_energy_flows_balance = Constraint(model.T, rule=H2_energy_flows_balance_rule)
    model.H2_load_shedding_upper_bound = Constraint(model.T, rule=H2_load_shedding_upper_bound_rule)


    # H2 Imports

    def H2_yearly_imports_budget_rule(model):
        return sum(model.P_H2_I[:]) <= model.psi_H2
                   
    def H2_imports_upper_bound_rule(model, t):
        return model.P_H2_I[t] <= model.kappa_H2_I[t]

    model.H2_yearly_imports_budget = Constraint(rule=H2_yearly_imports_budget_rule)
    model.H2_imports_upper_bound = Constraint(model.T, rule=H2_imports_upper_bound_rule)


    # Electrolysis Plant

    def electrolysis_power_definition_rule(model, t):
        return model.P_H2_EL[t] == model.eta_EL * model.P_E_EL[t]

    def electrolysis_power_lower_bound_rule(model, t):
        return model.mu_EL * model.K_E_EL <= model.P_E_EL[t]

    def electrolysis_power_upper_bound_rule(model, t):
        return model.P_E_EL[t] <= model.K_E_EL

    def electrolysis_plant_sizing_upper_bound_rule(model):
        return model.K_E_EL <= model.kappa_EL
    
    def electrolysis_water_consumption_rule(model, t):
        return model.m_H2O[t] == (model.MM_H2O * model.P_H2_EL[t]) / (model.MM_H2 * model.LHV_H2)
    
    def electrolysis_oxygen_production_rule(model, t):
        return model.m_O2[t] == (model.MM_O2 * model.P_H2_EL[t]) / (model.MM_H2 * model.LHV_H2)

    def electrolysis_plant_inc_ramp_rate_rule(model, t):
        if t>0:
            return model.P_H2_EL[t] <= model.P_H2_EL[t-1] + model.delta_EL_inc * model.K_E_EL
        else:
            return Constraint.Skip

    def electrolysis_plant_dec_ramp_rate_rule(model, t):
        if t>0:
            return model.P_H2_EL[t] >= model.P_H2_EL[t-1] - model.delta_EL_inc * model.K_E_EL
        else:
            return Constraint.Skip

    model.electrolysis_power_definition = Constraint(model.T, rule=electrolysis_power_definition_rule)
    model.electrolysis_power_lower_bound = Constraint(model.T, rule=electrolysis_power_lower_bound_rule)
    model.electrolysis_power_upper_bound = Constraint(model.T, rule=electrolysis_power_upper_bound_rule)
    model.electrolysis_plant_sizing_upper_bound = Constraint(rule=electrolysis_plant_sizing_upper_bound_rule)
    model.electrolysis_water_consumption = Constraint(model.T, rule=electrolysis_water_consumption_rule)
    model.electrolysis_oxygen_production = Constraint(model.T, rule=electrolysis_oxygen_production_rule)
    model.electrolysis_plant_inc_ramp_rate = Constraint(model.T, rule=electrolysis_plant_inc_ramp_rate_rule)
    model.electrolysis_plant_dec_ramp_rate = Constraint(model.T, rule=electrolysis_plant_dec_ramp_rate_rule)

    
    # Steam Methane Reformation Plant
    
    def SMR_plant_hydrogen_production_rule(model, t):
        return model.P_H2_SMR[t] == model.eta_SMR * model.P_NG_SMR[t]
    
    def SMR_plant_electricity_consumption_rule(model, t):
        return model.P_E_SMR[t] == model.phi_SMR * model.P_H2_SMR[t]
    
    def SMR_plant_sizing_rule(model, t):
        return model.P_H2_SMR[t] <= model.K_H2_SMR
    
    def SMR_plant_max_size_rule(model):
        return model.K_H2_SMR <= model.kappa_SMR
    
    model.SMR_plant_hydrogen_production = Constraint(model.T, rule=SMR_plant_hydrogen_production_rule)
    model.SMR_plant_electricity_consumption = Constraint(model.T, rule=SMR_plant_electricity_consumption_rule)
    model.SMR_plant_sizing = Constraint(model.T, rule=SMR_plant_sizing_rule)
    model.SMR_plant_max_size = Constraint(rule=SMR_plant_max_size_rule)


    # H2 Storage System
    
    def H2_storage_power_definition_rule(model, t):
        return model.P_H2S[t] == model.P_StH2[t] - model.P_H2tS[t]

    def H2_storage_output_power_initialisation_rule(model):
        return model.P_StH2[0] * model.delta_t <= model.E_H2[0] - model.sigma_H2S * model.S_H2

    def H2_storage_SOC_definition_rule(model, t):
        if t==0:
            return model.E_H2[t] == model.E_H2[model.t_max-1]
        else:
            return model.E_H2[t] == model.E_H2[t-1] + model.P_H2tS[t] * model.delta_t - model.P_StH2[t] * model.delta_t

    def H2_storage_SOC_lower_bound_rule(model, t):
        return model.sigma_H2S * model.S_H2 <= model.E_H2[t]

    def H2_storage_SOC_upper_bound_rule(model, t):
        return model.E_H2[t] <= model.S_H2

    def H2_storage_sizing_upper_bound_rule(model):
        return model.S_H2 <= model.xi_H2S
    
    def H2_storage_power_sizing_rule(model):
        return model.K_H2S == model.chi_H2S * model.S_H2

    def H2_storage_output_power_upper_bound_rule(model, t):
        return model.P_StH2[t] <= model.K_H2S

    def H2_storage_input_power_upper_bound_rule(model, t):
        return model.P_H2tS[t] <= model.rho_H2S * model.K_H2S
    
    model.H2_storage_power_definition = Constraint(model.T, rule=H2_storage_power_definition_rule)
    model.H2_storage_output_power_initialisation = Constraint(rule=H2_storage_output_power_initialisation_rule)
    model.H2_storage_SOC_definition = Constraint(model.T, rule=H2_storage_SOC_definition_rule)
    model.H2_storage_SOC_lower_bound = Constraint(model.T, rule=H2_storage_SOC_lower_bound_rule)
    model.H2_storage_SOC_upper_bound = Constraint(model.T, rule=H2_storage_SOC_upper_bound_rule)
    model.H2_storage_sizing_upper_bound = Constraint(rule=H2_storage_sizing_upper_bound_rule)
    model.H2_storage_power_sizing = Constraint(rule=H2_storage_power_sizing_rule)
    model.H2_storage_output_power_upper_bound = Constraint(model.T, rule=H2_storage_output_power_upper_bound_rule)
    model.H2_storage_input_power_upper_bound = Constraint(model.T, rule=H2_storage_input_power_upper_bound_rule)


     # Fuel Cells: H2 Repowering Plant

    def FC_plant_power_definition_rule(model, t):
        return model.P_E_FC[t] == model.eta_FC * model.P_H2_FC[t]

    def FC_plant_power_upper_bound_rule(model, t):
        return model.P_E_FC[t] <= model.K_E_FC

    def FC_plant_sizing_upper_bound_rule(model):
        return model.K_E_FC <= model.kappa_FC

    def FC_plant_inc_ramp_rate_rule(model, t):
        if t>0:
            return model.P_E_FC[t] <= model.P_E_FC[t-1] + model.delta_FC_inc * model.K_E_FC
        else:
            return Constraint.Skip

    def FC_plant_dec_ramp_rate_rule(model, t):
        if t>0:
            return model.P_E_FC[t] >= model.P_E_FC[t-1] - model.delta_FC_dec * model.K_E_FC
        else:
            return Constraint.Skip

    model.FC_plant_power_definition = Constraint(model.T, rule=FC_plant_power_definition_rule)
    model.FC_plant_power_upper_bound = Constraint(model.T, rule=FC_plant_power_upper_bound_rule)
    model.FC_plant_sizing_upper_bound = Constraint(rule=FC_plant_sizing_upper_bound_rule)
    model.FC_plant_inc_ramp_rate = Constraint(model.T, rule=FC_plant_inc_ramp_rate_rule)
    model.FC_plant_dec_ramp_rate = Constraint(model.T, rule=FC_plant_dec_ramp_rate_rule)


    # Methanation Plant

    def methanation_output_power_definition_rule(model, t):
        return model.P_CH4_MT[t] == model.eta_MT * model.P_H2_MT[t]

    def methanation_output_power_lower_bound_rule(model, t):
        return model.mu_MT * model.K_CH4_MT <= model.P_CH4_MT[t]

    def methanation_output_power_upper_bound_rule(model, t):
        return model.P_CH4_MT[t] <= model.K_CH4_MT

    def methanation_plant_sizing_upper_bound_rule(model):
        return model.K_CH4_MT <= model.kappa_MT
    
    def methanation_plant_CO2_consumption_rule(model, t):
        return model.Q_CO2_MT[t] == (model.MM_CO2 * model.P_CH4_MT[t]) / (model.MM_CH4 * model.LHV_CH4)
    
    def methanation_plant_inc_ramp_rate_rule(model, t):
        if t>0:
            return model.P_CH4_MT[t] <= model.P_CH4_MT[t-1] + model.delta_MT_inc * model.K_CH4_MT
        else:
            return Constraint.Skip

    def methanation_plant_dec_ramp_rate_rule(model, t):
        if t>0:
            return model.P_CH4_MT[t] >= model.P_CH4_MT[t-1] - model.delta_MT_inc * model.K_CH4_MT
        else:
            return Constraint.Skip

    model.methanation_output_power_definition = Constraint(model.T, rule=methanation_output_power_definition_rule)
    model.methanation_output_power_lower_bound = Constraint(model.T, rule=methanation_output_power_lower_bound_rule)
    model.methanation_output_power_upper_bound = Constraint(model.T, rule=methanation_output_power_upper_bound_rule)
    model.methanation_plant_sizing_upper_bound = Constraint(rule=methanation_plant_sizing_upper_bound_rule)
    model.methanation_plant_CO2_consumption = Constraint(model.T, rule=methanation_plant_CO2_consumption_rule)
    model.methanation_plant_inc_ramp_rate = Constraint(model.T, rule=methanation_plant_inc_ramp_rate_rule)
    model.methanation_plant_dec_ramp_rate = Constraint(model.T, rule=methanation_plant_dec_ramp_rate_rule)
    

    # Natural Gas Network Balance

    def NG_network_balance_rule(model, t):
        return model.P_NG_I[t] + model.P_CH4_MT[t] + model.P_NGS[t] + model.L_NG_ENS[t] == model.lambda_NG_HT[t] + model.lambda_NG_ID[t] - model.lambda_NGtH2[t] + model.lambda_NG_TR[t] + model.P_NG_OCGT[t] + model.P_NG_CCGT[t] + model.P_NG_CHP[t] + model.P_NG_SMR[t] + model.P_NG_ACC[t]

    def NG_load_shedding_upper_bound_rule(model, t):
        return model.L_NG_ENS[t] <= model.kappa_L_NG_ENS
    
    def NG_network_linepack_capacity_rule(model, t):
        return model.lambda_NG_HT[t] + model.lambda_NG_ID[t] - model.lambda_NGtH2[t] + model.lambda_NG_TR[t] + model.P_NG_OCGT[t] + model.P_NG_CCGT[t] + model.P_NG_CHP[t] + model.P_NG_SMR[t]  + model.P_NG_ACC[t] <= model.kappa_NGNet
    
    def NG_yearly_imports_budget_rule(model):
        return sum(model.P_NG_I[:]) <= model.psi_NG

    model.NG_network_balance = Constraint(model.T, rule=NG_network_balance_rule)
    model.NG_load_shedding_upper_bound = Constraint(model.T, rule=NG_load_shedding_upper_bound_rule)
    model.NG_network_linepack_capacity = Constraint(model.T, rule=NG_network_linepack_capacity_rule)
    model.NG_yearly_imports_budget = Constraint(rule=NG_yearly_imports_budget_rule)
    
    # NG Storage System

    def NG_storage_SOC_definition_rule(model, t):
        if t==0:
            return model.E_NGS[t] == model.E_NGS[model.t_max-1]
        else:
            return model.E_NGS[t] == model.E_NGS[t-1] + model.eta_NGS * model.P_NGtNGS[t] * model.delta_t - model.P_NGStNG[t] * model.delta_t

    def NG_storage_SOC_lower_bound_rule(model, t):
        return model.sigma_NGS * model.xi_NGS <= model.E_NGS[t]

    def NG_storage_SOC_upper_bound_rule(model, t):
        return model.E_NGS[t] <= model.xi_NGS

    def NG_storage_output_power_initialisation_rule(model):
        return model.P_NGStNG[0] <= model.E_NGS[0] - model.sigma_NGS * model.xi_NGS
    
    def NG_storage_power_definition_rule(model, t):
        return model.P_NGS[t] == - model.P_NGtNGS[t] + model.P_NGStNG[t]

    def NG_storage_output_power_upper_bound_rule(model, t):
        return model.P_NGStNG[t] <= model.kappa_NGS
    
    def NG_storage_input_power_upper_bound_rule(model, t):
        return model.P_NGtNGS[t] <= model.rho_NGS * model.kappa_NGS

    model.NG_storage_SOC_definition = Constraint(model.T, rule=NG_storage_SOC_definition_rule)
    model.NG_storage_SOC_lower_bound = Constraint(model.T, rule=NG_storage_SOC_lower_bound_rule)
    model.NG_storage_SOC_upper_bound = Constraint(model.T, rule=NG_storage_SOC_upper_bound_rule)
    model.NG_storage_output_power_initialisation = Constraint(rule=NG_storage_output_power_initialisation_rule)
    model.NG_storage_power_definition = Constraint(model.T, rule=NG_storage_power_definition_rule)
    model.NG_storage_output_power_upper_bound = Constraint(model.T, rule=NG_storage_output_power_upper_bound_rule)
    model.NG_storage_input_power_upper_bound = Constraint(model.T, rule=NG_storage_input_power_upper_bound_rule)

    # Gas-Fired Power Plants: OCGT

    def ocgt_plant_output_power_definition_rule(model, t):
        return model.P_OCGT[t] == model.eta_OCGT * model.P_NG_OCGT[t]

    def ocgt_plant_upper_power_bound_rule(model, t):
        return model.P_OCGT[t] <= model.kappa_OCGT_0 + model.K_E_OCGT

    
    def ocgt_plant_output_power_balance_rule(model, t):
        return model.P_OCGT[t] == model.P_OCGT_CCS[t] + model.P_E_OCGT[t]

    def ocgt_plant_inc_ramp_rate_rule(model, t):
        if t>0:
            return model.P_OCGT[t] <= model.P_OCGT[t-1] + model.delta_OCGT_inc * model.K_E_OCGT
        else:
            return Constraint.Skip

    def ocgt_plant_dec_ramp_rate_rule(model, t):
        if t>0:
            return model.P_OCGT[t] >= model.P_OCGT[t-1] - model.delta_OCGT_inc * model.K_E_OCGT
        else:
            return Constraint.Skip

    model.ocgt_plant_output_power_definition = Constraint(model.T, rule=ocgt_plant_output_power_definition_rule)
    model.ocgt_plant_upper_power_bound = Constraint(model.T, rule=ocgt_plant_upper_power_bound_rule)
    model.ocgt_plant_output_power_balance = Constraint(model.T, rule=ocgt_plant_output_power_balance_rule)
    model.ocgt_plant_inc_ramp_rate = Constraint(model.T, rule=ocgt_plant_inc_ramp_rate_rule)
    model.ocgt_plant_dec_ramp_rate = Constraint(model.T, rule=ocgt_plant_dec_ramp_rate_rule)
    
    
    # Gas-Fired Power Plants: CCGT

    def ccgt_plant_output_power_definition_rule(model, t):
        return model.P_CCGT[t] == model.eta_CCGT * model.P_NG_CCGT[t]

    def ccgt_plant_upper_power_bound_rule(model, t):
        return model.P_CCGT[t] <= model.kappa_CCGT_0 + model.K_E_CCGT

    
    def ccgt_plant_output_power_balance_rule(model, t):
        return model.P_CCGT[t] == model.P_CCGT_CCS[t] + model.P_E_CCGT[t]

    def ccgt_plant_inc_ramp_rate_rule(model, t):
        if t>0:
            return model.P_CCGT[t] <= model.P_CCGT[t-1] + model.delta_CCGT_inc * model.K_E_CCGT
        else:
            return Constraint.Skip

    def ccgt_plant_dec_ramp_rate_rule(model, t):
        if t>0:
            return model.P_CCGT[t] >= model.P_CCGT[t-1] - model.delta_CCGT_inc * model.K_E_CCGT
        else:
            return Constraint.Skip

    model.ccgt_plant_output_power_definition = Constraint(model.T, rule=ccgt_plant_output_power_definition_rule)
    model.ccgt_plant_upper_power_bound = Constraint(model.T, rule=ccgt_plant_upper_power_bound_rule)
    model.ccgt_plant_output_power_balance = Constraint(model.T, rule=ccgt_plant_output_power_balance_rule)
    model.ccgt_plant_inc_ramp_rate = Constraint(model.T, rule=ccgt_plant_inc_ramp_rate_rule)
    model.ccgt_plant_dec_ramp_rate = Constraint(model.T, rule=ccgt_plant_dec_ramp_rate_rule)


    # Gas-Fired Power Plants CCS: OCGT
    
    def ocgt_plant_carbon_balance_rule(model, t):
        return model.nu_NG_CO2 * model.P_NG_OCGT[t] == model.Q_CO2_OCGT_CCS[t] + model.Q_CO2_OCGT[t]
    
    def ocgt_plant_CCS_efficiency_rule(model, t):
        return model.Q_CO2_OCGT_CCS[t] <= model.eta_NG_CCS * model.nu_NG_CO2 * model.P_NG_OCGT[t]
    
    def ocgt_plant_CCS_consumption_rule(model, t):
        return model.phi_NG_CCS * model.Q_CO2_OCGT_CCS[t] == model.P_OCGT_CCS[t]
    
    def ocgt_plant_CCS_sizing_rule(model, t):
        return model.Q_CO2_OCGT_CCS[t] <= model.K_CO2_OCGT_CCS
    
    def ocgt_plant_CCS_sizing_upper_bound_rule(model):
         return model.K_CO2_OCGT_CCS <= model.kappa_CO2_OCGT_CCS
    
    model.ocgt_plant_carbon_balance = Constraint(model.T, rule=ocgt_plant_carbon_balance_rule)
    model.ocgt_plant_CCS_efficiency = Constraint(model.T, rule=ocgt_plant_CCS_efficiency_rule)
    model.ocgt_plant_CCS_consumption = Constraint(model.T, rule=ocgt_plant_CCS_consumption_rule)
    model.ocgt_plant_CCS_sizing = Constraint(model.T, rule=ocgt_plant_CCS_sizing_rule)
    model.ocgt_plant_CCS_sizing_upper_bound = Constraint(rule=ocgt_plant_CCS_sizing_upper_bound_rule)
    
    
    # Gas-Fired Power Plants CCS: CCGT
    
    def ccgt_plant_carbon_balance_rule(model, t):
        return model.nu_NG_CO2 * model.P_NG_CCGT[t] == model.Q_CO2_CCGT_CCS[t] + model.Q_CO2_CCGT[t]
    
    def ccgt_plant_CCS_efficiency_rule(model, t):
        return model.Q_CO2_CCGT_CCS[t] <= model.eta_NG_CCS * model.nu_NG_CO2 * model.P_NG_CCGT[t]
    
    def ccgt_plant_CCS_consumption_rule(model, t):
        return model.phi_NG_CCS * model.Q_CO2_CCGT_CCS[t] == model.P_CCGT_CCS[t]
    
    def ccgt_plant_CCS_sizing_rule(model, t):
        return model.Q_CO2_CCGT_CCS[t] <= model.K_CO2_CCGT_CCS
    
    def ccgt_plant_CCS_sizing_upper_bound_rule(model):
         return model.K_CO2_CCGT_CCS <= model.kappa_CO2_CCGT_CCS
    
    model.ccgt_plant_carbon_balance = Constraint(model.T, rule=ccgt_plant_carbon_balance_rule)
    model.ccgt_plant_CCS_efficiency = Constraint(model.T, rule=ccgt_plant_CCS_efficiency_rule)
    model.ccgt_plant_CCS_consumption = Constraint(model.T, rule=ccgt_plant_CCS_consumption_rule)
    model.ccgt_plant_CCS_sizing = Constraint(model.T, rule=ccgt_plant_CCS_sizing_rule)
    model.ccgt_plant_CCS_sizing_upper_bound = Constraint(rule=ccgt_plant_CCS_sizing_upper_bound_rule)
    
    
    # CHP Plants CCS
    
    def CHP_plant_carbon_balance_rule(model, t):
        return model.nu_NG_CO2 * model.P_NG_CHP[t] == model.Q_CO2_CHP_CCS[t] + model.Q_CO2_CHP[t]
    
    def CHP_plant_CCS_efficiency_rule(model, t):
        return model.Q_CO2_CHP_CCS[t] <= model.eta_CHP_CCS * model.nu_NG_CO2 * model.P_NG_CHP[t]
    
    def CHP_plant_CCS_consumption_rule(model, t):
        return model.phi_CHP_CCS * model.Q_CO2_CHP_CCS[t] == model.P_CHP_CCS[t]
    
    def CHP_plant_CCS_sizing_rule(model, t):
        return model.Q_CO2_CHP_CCS[t] <= model.K_CO2_CHP_CCS
    
    def CHP_plant_CCS_sizing_upper_bound_rule(model):
         return model.K_CO2_CHP_CCS <= model.kappa_CO2_CHP_CCS
    
    model.CHP_plant_carbon_balance = Constraint(model.T, rule=CHP_plant_carbon_balance_rule)
    model.CHP_plant_CCS_efficiency = Constraint(model.T, rule=CHP_plant_CCS_efficiency_rule)
    model.CHP_plant_CCS_consumption = Constraint(model.T, rule=CHP_plant_CCS_consumption_rule)
    model.CHP_plant_CCS_sizing = Constraint(model.T, rule=CHP_plant_CCS_sizing_rule)
    model.CHP_plant_CCS_sizing_upper_bound = Constraint(rule=CHP_plant_CCS_sizing_upper_bound_rule)
    
    
    # Biomass Plants CCS
    
    def biomass_plant_carbon_balance_rule(model, t):
        return model.nu_BM_CO2 * model.P_BM[t] / model.eta_BM == model.Q_CO2_BM_CCS[t] + model.Q_CO2_BM[t]
    
    def biomass_plant_CCS_efficiency_rule(model, t):
        return model.Q_CO2_BM_CCS[t] <= model.eta_BM_CCS * model.nu_BM_CO2 * model.P_BM[t] / model.eta_BM
    
    def biomass_plant_CCS_consumption_rule(model, t):
        return model.phi_BM_CCS * model.Q_CO2_BM_CCS[t] == model.P_BM_CCS[t]
    
    def biomass_plant_CCS_sizing_rule(model, t):
        return model.Q_CO2_BM_CCS[t] <= model.K_CO2_BM_CCS
    
    def biomass_plant_CCS_sizing_upper_bound_rule(model):
         return model.K_CO2_BM_CCS <= model.kappa_CO2_BM_CCS
    
    model.biomass_plant_carbon_balance = Constraint(model.T, rule=biomass_plant_carbon_balance_rule)
    model.biomass_plant_CCS_efficiency = Constraint(model.T, rule=biomass_plant_CCS_efficiency_rule)
    model.biomass_plant_CCS_consumption = Constraint(model.T, rule=biomass_plant_CCS_consumption_rule)
    model.biomass_plant_CCS_sizing = Constraint(model.T, rule=biomass_plant_CCS_sizing_rule)
    model.biomass_plant_CCS_sizing_upper_bound = Constraint(rule=biomass_plant_CCS_sizing_upper_bound_rule)
    
    
    # Waste Plants CCS
    
    def waste_plant_carbon_balance_rule(model, t):
        return model.nu_WS_CO2 * model.P_WS[t] / model.eta_WS == model.Q_CO2_WS_CCS[t] + model.Q_CO2_WS[t]
    
    def waste_plant_CCS_efficiency_rule(model, t):
        return model.Q_CO2_WS_CCS[t] <= model.eta_WS_CCS * model.nu_WS_CO2 * model.P_WS[t] / model.eta_WS
    
    def waste_plant_CCS_consumption_rule(model, t):
        return model.phi_WS_CCS * model.Q_CO2_WS_CCS[t] == model.P_WS_CCS[t]
    
    def waste_plant_CCS_sizing_rule(model, t):
        return model.Q_CO2_WS_CCS[t] <= model.K_CO2_WS_CCS
    
    def waste_plant_CCS_sizing_upper_bound_rule(model):
         return model.K_CO2_WS_CCS <= model.kappa_CO2_WS_CCS
    
    model.waste_plant_carbon_balance = Constraint(model.T, rule=waste_plant_carbon_balance_rule)
    model.waste_plant_CCS_efficiency = Constraint(model.T, rule=waste_plant_CCS_efficiency_rule)
    model.waste_plant_CCS_consumption = Constraint(model.T, rule=waste_plant_CCS_consumption_rule)
    model.waste_plant_CCS_sizing = Constraint(model.T, rule=waste_plant_CCS_sizing_rule)
    model.waste_plant_CCS_sizing_upper_bound = Constraint(rule=waste_plant_CCS_sizing_upper_bound_rule)
    
    
    # SMR Plants CCS
    
    def SMR_plant_carbon_balance_rule(model, t):
        return model.nu_NG_CO2 * model.P_NG_SMR[t] == model.Q_CO2_SMR_CCS[t] + model.Q_CO2_SMR[t]
    
    def SMR_plant_CCS_efficiency_rule(model, t):
        return model.Q_CO2_SMR_CCS[t] <= model.eta_SMR_CCS * model.nu_NG_CO2 * model.P_NG_SMR[t]
    
    def SMR_plant_CCS_consumption_rule(model, t):
        return model.phi_SMR_CCS * model.Q_CO2_SMR_CCS[t] == model.P_SMR_CCS[t]
    
    def SMR_plant_CCS_sizing_rule(model, t):
        return model.Q_CO2_SMR_CCS[t] <= model.K_CO2_SMR_CCS
    
    def SMR_plant_CCS_sizing_upper_bound_rule(model):
         return model.K_CO2_SMR_CCS <= model.kappa_CO2_SMR_CCS
    
    model.SMR_plant_carbon_balance = Constraint(model.T, rule=SMR_plant_carbon_balance_rule)
    model.SMR_plant_CCS_efficiency = Constraint(model.T, rule=SMR_plant_CCS_efficiency_rule)
    model.SMR_plant_CCS_consumption = Constraint(model.T, rule=SMR_plant_CCS_consumption_rule)
    model.SMR_plant_CCS_sizing = Constraint(model.T, rule=SMR_plant_CCS_sizing_rule)
    model.SMR_plant_CCS_sizing_upper_bound = Constraint(rule=SMR_plant_CCS_sizing_upper_bound_rule)

    # Atmospheric Carbon Capture

    def ACC_E_consumption_rule(model, t):
        return model.P_E_ACC[t] == model.phi_E_ACC * model.Q_CO2_ACC_A[t]

    def ACC_NG_consumption_rule(model, t):
        return model.P_NG_ACC[t] == model.phi_NG_ACC * model.Q_CO2_ACC_A[t]
    
    def ACC_CO2_production_rule(model, t):
        return model.Q_CO2_ACC[t] == model.Q_CO2_ACC_A[t] + model.nu_NG_CO2*model.P_NG_ACC[t]

    def ACC_sizing_rule(model, t):
        return model.Q_CO2_ACC_A[t] <= model.K_CO2_ACC

    def ACC_sizing_upper_bound_rule(model):
        return model.K_CO2_ACC <= model.kappa_CO2_ACC

    model.ACC_E_consumption = Constraint(model.T, rule=ACC_E_consumption_rule)
    model.ACC_NG_consumption = Constraint(model.T, rule=ACC_NG_consumption_rule)
    model.ACC_CO2_production = Constraint(model.T, rule=ACC_CO2_production_rule)
    model.ACC_sizing = Constraint(model.T, rule=ACC_sizing_rule)
    model.ACC_sizing_upper_bound = Constraint(rule=ACC_sizing_upper_bound_rule)


    # CO2 Storage System

    def CO2_storage_SOC_definition_rule(model, t):
        if t == 0:
            return model.M_CO2[t] == model.M_CO2[model.t_max-1]
        else:
            return model.M_CO2[t] == model.M_CO2[t - 1] + model.Q_CO2_GtS[t] * model.delta_t - model.Q_CO2_StG[t] * model.delta_t

    def CO2_storage_net_flow_rule(model, t):
        return model.Q_CO2S[t] == model.Q_CO2_StG[t] - model.Q_CO2_GtS[t]

    def CO2_storage_sizing_rule(model, t):
        return model.M_CO2[t] <= model.S_CO2
    
    def CO2_storage_sizing_upper_bound_rule(model):
        return model.S_CO2 <= model.xi_CO2S
    
    def CO2_storage_flow_sizing_rule(model):
        return model.K_CO2S == model.chi_CO2S * model.S_CO2

    def CO2_storage_output_flow_upper_bound_rule(model, t):
        return model.Q_CO2_StG[t] <= model.K_CO2S

    def CO2_storage_input_flow_upper_bound_rule(model, t):
        return model.Q_CO2_GtS[t] <= model.rho_CO2S * model.K_CO2S
    
    model.CO2_storage_SOC_definition = Constraint(model.T, rule=CO2_storage_SOC_definition_rule)
    model.CO2_storage_net_flow = Constraint(model.T, rule=CO2_storage_net_flow_rule)
    model.CO2_storage_sizing = Constraint(model.T, rule=CO2_storage_sizing_rule)
    model.CO2_storage_flow_sizing = Constraint(rule=CO2_storage_flow_sizing_rule)
    model.CO2_storage_sizing_upper_bound = Constraint(model.T, rule=CO2_storage_sizing_upper_bound_rule)
    model.CO2_storage_output_flow_upper_bound = Constraint(model.T, rule=CO2_storage_output_flow_upper_bound_rule)
    model.CO2_storage_input_flow_upper_bound = Constraint(model.T, rule=CO2_storage_input_flow_upper_bound_rule)


    # CO2 Exports

    def CO2_exports_upper_bound_rule(model, t):
        return model.Q_CO2_E[t] <= model.kappa_CO2_E
    
    model.CO2_exports_upper_bound = Constraint(model.T, rule=CO2_exports_upper_bound_rule)


    # CO2 Balance

    def CO2_balance_rule(model, t):
        return model.Q_CO2_ACC[t] - model.Q_CO2_E[t] + model.Q_CO2S[t] + model.Q_CO2_OCGT_CCS[t] + model.Q_CO2_CCGT_CCS[t] + model.Q_CO2_CHP_CCS[t] + \
                        model.Q_CO2_BM_CCS[t] + model.Q_CO2_WS_CCS[t] + model.Q_CO2_SMR_CCS[t] == model.Q_CO2_MT[t]
                        
    model.CO2_balance = Constraint(model.T, rule=CO2_balance_rule)


    # CO2 Budget

    def yearly_CO2_budget_rule(model):
        return model.nu_NG_CO2 * (sum(model.lambda_NG_HT[:]) + sum(model.lambda_NG_ID[:]) - sum(model.lambda_NGtH2[:]) - sum(model.L_NG_ENS[:]) + sum(model.lambda_NG_TR[:])) * model.delta_t + \
               sum(model.Q_CO2_OCGT[:]) * model.delta_t + sum(model.Q_CO2_CCGT[:]) * model.delta_t + \
               sum(model.Q_CO2_CHP[:]) * model.delta_t + sum(model.Q_CO2_SMR[:]) * model.delta_t + \
               sum(model.Q_CO2_BM[:]) * model.delta_t + sum(model.Q_CO2_WS[:]) * model.delta_t + \
               model.nu_E_I_CO2 * sum(model.P_E_I[:]) * model.delta_t - sum(model.Q_CO2_ACC_A[:]) * model.delta_t  <= model.psi_CO2

    model.yearly_CO2_budget = Constraint(rule=yearly_CO2_budget_rule)


    ### Objective Function

    def cost_rule(model):
        return (model.zeta_W_on + model.theta_W_on_f * model.n_y) * model.K_W_on +\
            model.theta_W_on_v * sum(model.P_W_on[:]) * model.delta_t +\
            (model.zeta_W_off + model.theta_W_off_f * model.n_y) * model.K_W_off +\
            model.theta_W_off_v * sum(model.P_W_off[:]) * model.delta_t +\
            (model.zeta_S + model.theta_S_f * model.n_y) * model.K_S +\
            model.theta_S_v * sum(model.P_S[:]) * model.delta_t +\
            model.zeta_B_E * model.S_B + model.zeta_B_P * model.K_B +\
            model.theta_B_E_f * model.n_y * model.S_B + model.theta_B_P_f * model.n_y * model.K_B +\
            (model.zeta_EL + model.theta_EL_f * model.n_y) * model.K_E_EL +\
            (model.zeta_H2S + model.theta_H2S_f * model.n_y) * model.S_H2 +\
            (model.zeta_FC + model.theta_FC_f * model.n_y) * model.K_E_FC + model.theta_FC_v * sum(model.P_E_FC[:]) * model.delta_t +\
            (model.zeta_MT + model.theta_MT_f * model.n_y) * model.K_CH4_MT +\
            (model.zeta_OCGT + model.theta_OCGT_f * model.n_y) * model.K_E_OCGT +\
            (model.zeta_CCGT + model.theta_CCGT_f * model.n_y) * model.K_E_CCGT +\
            (model.zeta_CO2S + model.theta_CO2S_f * model.n_y) * model.S_CO2 +\
            (model.zeta_NG_CCS + model.theta_NG_CCS_f * model.n_y) * model.K_CO2_OCGT_CCS +\
            (model.zeta_NG_CCS + model.theta_NG_CCS_f * model.n_y) * model.K_CO2_CCGT_CCS +\
            (model.zeta_CHP_CCS + model.theta_NG_CCS_f * model.n_y) * model.K_CO2_CHP_CCS +\
            (model.zeta_BM_CCS + model.theta_BM_CCS_f * model.n_y) * model.K_CO2_BM_CCS +\
            (model.zeta_WS_CCS + model.theta_WS_CCS_f * model.n_y) * model.K_CO2_WS_CCS +\
            (model.zeta_ACC + model.theta_ACC_f * model.n_y) * model.K_CO2_ACC +\
            model.theta_NG_CCS_v * sum(model.Q_CO2_OCGT_CCS[:]) * model.delta_t +\
            model.theta_NG_CCS_v * sum(model.Q_CO2_CCGT_CCS[:]) * model.delta_t +\
            model.theta_CHP_CCS_v * sum(model.Q_CO2_CHP_CCS[:]) * model.delta_t +\
            model.theta_BM_CCS_v * sum(model.Q_CO2_BM_CCS[:]) * model.delta_t +\
            model.theta_WS_CCS_v * sum(model.Q_CO2_WS_CCS[:]) * model.delta_t +\
            (model.zeta_SMR + model.theta_SMR_f * model.n_y) * model.K_H2_SMR +\
            model.theta_SMR_v * sum(model.P_H2_SMR[:]) * model.delta_t +\
            (model.zeta_SMR_CCS + model.theta_SMR_CCS_f * model.n_y) * model.K_CO2_SMR_CCS +\
            model.theta_OCGT_v * sum(model.P_OCGT[:]) * model.delta_t +\
            model.theta_CCGT_v * sum(model.P_CCGT[:]) * model.delta_t +\
            sum(model.theta_NG_I[t] * model.P_NG_I[t] * model.delta_t for t in model.T) +\
            sum(model.theta_H2_IE[t] * model.P_H2_I[t] for t in model.T) * model.delta_t +\
            model.theta_CO2 * sum(model.Q_CO2_OCGT[:]) * model.delta_t +\
            model.theta_CO2 * sum(model.Q_CO2_CCGT[:]) * model.delta_t +\
            model.theta_NK_v * sum(model.P_E_NK[:]) * model.delta_t +\
            model.theta_NK_fuel * sum(model.P_E_NK[:]) * model.delta_t / model.eta_NK +\
            model.theta_CHP_v * sum(model.P_CHP[:]) * model.delta_t +\
            model.theta_BM_v * sum(model.P_BM[:]) * model.delta_t +\
            model.theta_WS_v * sum(model.P_WS[:]) * model.delta_t +\
            model.theta_BM_fuel * sum(model.P_BM[:]) * model.delta_t / model.eta_BM +\
            model.theta_WS_fuel * sum(model.P_WS[:]) * model.delta_t / model.eta_WS +\
            model.theta_CO2 * sum(model.Q_CO2_CHP[:]) * model.delta_t +\
            model.theta_CO2 * sum(model.Q_CO2_BM[:]) * model.delta_t +\
            model.theta_CO2 * sum(model.Q_CO2_WS[:]) * model.delta_t +\
            model.theta_PH_v * sum(model.P_PtPH[:]) * model.delta_t +\
            model.varsigma_ENS_E * sum(model.L_E_ENS[:]) * model.delta_t +\
            sum(model.theta_IE[t] * model.P_E_IE[t] * model.delta_t for t in model.T) +\
            model.varsigma_ENS_NG * sum(model.L_NG_ENS[:]) * model.delta_t +\
            model.varsigma_ENS_H2 * sum(model.L_H2_ENS[:]) * model.delta_t +\
            model.theta_CO2_E * sum(model.Q_CO2_E[:]) * model.delta_t

    model.cost = Objective(rule=cost_rule, sense=minimize)
    
    model.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
    model.rc = Suffix(direction=Suffix.IMPORT_EXPORT)
    model.slack = Suffix(direction=Suffix.IMPORT_EXPORT)

    print('RES limits:', model.kappa_W_off_max.value, model.kappa_W_on_max.value, model.kappa_S_max.value)
    print('CO2 budget:', model.psi_CO2.value)
    print('Nuclear capacity:', model.kappa_NK.value)

    dir_run	= folder
	
    if write_lp:
        model.write(filename=join(dir_run, "model.lp"),
            format=ProblemFormat.cpxlp,
            io_options={"symbolic_solver_labels": True})
        model.write(filename=join(dir_run, 'model.mps'),
            format=ProblemFormat.mps)

    return model