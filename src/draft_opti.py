#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 13:47:51 2018

@author: mathiasberger
"""

from pyomo.environ import ConcreteModel, Param, RangeSet, Var, Constraint, Objective, NonNegativeReals, Reals, Binary, minimize
from pyomo.opt import ProblemFormat
from os.path import join

def model_builder(data, folder, write_lp=False):

    model = ConcreteModel()


    ### Parameters Definition


    # Time Parameters
    model.t_max = Param(initialize=data.t_max)
    model.delta_t = Param(initialize=data.delta_t)
    model.T = RangeSet(0, model.t_max-1)
    model.n_y = Param(initialize=data.n_y)

    # Load factor parameters
    model.gamma_S = Param(model.T, initialize=data.gamma_S)
    model.gamma_W_on = Param(model.T, initialize=data.gamma_W_on)
    model.gamma_W_off = Param(model.T, initialize=data.gamma_W_off)
    model.gamma_L = Param(model.T, initialize=data.gamma_L)

    # Capacity parameters
    model.kappa_L = Param(initialize=data.kappa_L)
    model.kappa_W_on_max = Param(initialize=data.kappa_W_on_max)
    model.kappa_W_on_0 = Param(initialize=data.kappa_W_on_0)
    model.kappa_W_off_max = Param(initialize=data.kappa_W_off_max)
    model.kappa_W_off_0 = Param(initialize=data.kappa_W_off_0)
    model.kappa_S_max = Param(initialize=data.kappa_S_max)
    model.kappa_S_0 = Param(initialize=data.kappa_S_0)
    model.kappa_NG_0 = Param(model.T, initialize=data.kappa_NG_0)
    model.kappa_NG_max = Param(initialize=data.kappa_NG_max)
    model.kappa_NK = Param(model.T, initialize=data.kappa_NK)
    model.kappa_PtG = Param(initialize=data.kappa_PtG)
    model.xi_H2 = Param(initialize=data.xi_H2)
    model.kappa_H2 = Param(initialize=data.kappa_H2)
    model.kappa_H2tCH4 = Param(initialize=data.kappa_H2tCH4)
    model.xi_CH4 = Param(initialize=data.xi_CH4)
    model.kappa_CH4tNG = Param(initialize=data.kappa_CH4tNG)
    model.xi_PH = Param(initialize=data.xi_PH)
    model.kappa_PtPH = Param(initialize=data.kappa_PtPH)
    model.kappa_PHtP = Param(initialize=data.kappa_PHtP)
    model.xi_B = Param(initialize=data.xi_B)
    model.kappa_B = Param(initialize=data.kappa_B)
    model.kappa_NGNet = Param(initialize=data.kappa_NGNet)
    model.kappa_disp = Param(initialize=data.kappa_disp)
    model.kappa_trs = Param(initialize=data.kappa_trs)
    model.psi_CO2 = Param(initialize=data.psi_CO2)

    # Demand & Curtailment Parameters
    model.pi_L = Param(model.T, initialize=data.pi_L)
    model.pi_NG = Param(model.T, initialize=data.pi_NG)

    # Efficiency parameters
    model.eta_NK = Param(initialize=data.eta_NK)
    model.eta_H2tP = Param(initialize=data.eta_H2tP)
    model.eta_NGtP = Param(initialize=data.eta_NGtP)
    model.eta_PtG = Param(initialize=data.eta_PtG)
    model.eta_H2tCH4 = Param(initialize=data.eta_H2tCH4)
    model.eta_PHtP = Param(initialize=data.eta_PHtP)
    model.eta_PtPH = Param(initialize=data.eta_PtPH)
    model.eta_BtP = Param(initialize=data.eta_BtP)
    model.eta_PtB = Param(initialize=data.eta_PtB)
    model.eta_B = Param(initialize=data.eta_B)
    model.eta_disp = Param(initialize=data.eta_disp)

    # Cost parameters
    model.zeta_W_on = Param(initialize=data.zeta_W_on)
    model.zeta_W_off = Param(initialize=data.zeta_W_off)
    model.zeta_S = Param(initialize=data.zeta_S)
    model.zeta_PtG = Param(initialize=data.zeta_PtG)
    model.zeta_H2_s = Param(initialize=data.zeta_H2_s)
    model.zeta_H2 = Param(initialize=data.zeta_H2)
    model.zeta_H2tCH4 = Param(initialize=data.zeta_H2tCH4)
    model.zeta_CH4 = Param(initialize=data.zeta_CH4)
    model.zeta_NG = Param(initialize=data.zeta_NG)
    model.zeta_B = Param(initialize=data.zeta_B)
    model.varsigma_ENS = Param(initialize=data.varsigma_ENS)
    model.varsigma_C = Param(initialize=data.varsigma_C)
    model.theta_CO2 = Param(initialize=data.theta_CO2)
    model.theta_NG_f = Param(initialize=data.theta_NG_f)
    model.theta_NG_v = Param(initialize=data.theta_NG_v)
    model.theta_NG_fuel = Param(initialize=data.theta_NG_fuel)
    model.theta_NK_v = Param(initialize=data.theta_NK_v)
    model.theta_NK_fuel = Param(initialize=data.theta_NK_fuel)
    model.theta_disp_v = Param(initialize=data.theta_disp_v)
    model.theta_disp_fuel = Param(initialize=data.theta_disp_fuel)
    model.theta_W_on_f = Param(initialize=data.theta_W_on_f)
    model.theta_W_on_v = Param(initialize=data.theta_W_on_v)
    model.theta_W_off_f = Param(initialize=data.theta_W_off_f)
    model.theta_W_off_v = Param(initialize=data.theta_W_off_v)
    model.theta_S_f = Param(initialize=data.theta_S_f)
    model.theta_S_v = Param(initialize=data.theta_S_v)
    model.theta_PtG_f = Param(initialize=data.theta_PtG_f)
    model.theta_H2_s_f = Param(initialize=data.theta_H2_s_f)
    model.theta_B_f = Param(initialize=data.theta_B_f)
    model.theta_H2_f = Param(initialize=data.theta_H2_f)
    model.theta_H2_v = Param(initialize=data.theta_H2_v)
    model.theta_H2tCH4_f = Param(initialize=data.theta_H2tCH4_f)
    model.theta_CH4_f = Param(initialize=data.theta_CH4_f)
    model.theta_PH_v = Param(initialize=data.theta_PH_v)
    model.theta_el = Param(model.T, initialize=data.theta_el)

    # Other parameters
    model.delta_disp_inc = Param(initialize=data.delta_disp_inc)
    model.delta_disp_dec = Param(initialize=data.delta_disp_dec)
    model.delta_NK_inc = Param(initialize=data.delta_NK_inc)
    model.delta_NK_dec = Param(initialize=data.delta_NK_dec)
    model.sigma_H2 = Param(initialize=data.sigma_H2)
    model.sigma_PH = Param(initialize=data.sigma_PH)
    model.sigma_CH4 = Param(initialize=data.sigma_CH4)
    model.sigma_B = Param(initialize=data.sigma_B)
    model.rho_B = Param(initialize=data.rho_B)
    model.chi_B = Param(initialize=data.chi_B)
    model.mu_disp = Param(initialize=data.mu_disp)
    model.mu_trs = Param(initialize=data.mu_trs)
    model.mu_NK = Param(initialize=data.mu_NK)
    model.nu_NG_CO2 = Param(initialize=data.nu_NG_CO2)
    model.nu_CH4_CO2 = Param(initialize=data.nu_CH4_CO2)
    model.nu_disp_CO2 = Param(initialize=data.nu_disp_CO2)
    model.nu_trs_CO2 = Param(initialize=data.nu_trs_CO2)


    ### Variables Definition


    # Operational Variables
    model.P_W_on = Var(model.T, within=NonNegativeReals)
    model.P_W_off = Var(model.T, within=NonNegativeReals)
    model.P_S = Var(model.T, within=NonNegativeReals)
    model.P_C = Var(model.T, within=NonNegativeReals)
#    model.B_C = Var(model.T, within=Binary)
#    model.P_C = Var(model.T, within=Reals)
#    model.P_C_p = Var(model.T, within=NonNegativeReals)
#    model.P_C_n = Var(model.T, within=NonNegativeReals)
    model.Delta_P = Var(model.T, within=NonNegativeReals)
    model.P_PtG = Var(model.T, within=NonNegativeReals)
    model.P_PtH2 = Var(model.T, within=NonNegativeReals)
    model.E_H2 = Var(model.T, within=NonNegativeReals)
    model.P_H2_out = Var(model.T, within=NonNegativeReals)
    model.P_H2tP = Var(model.T, within=NonNegativeReals)
    model.P_H2tCH4 = Var(model.T, within=NonNegativeReals)
    model.P_H2 = Var(model.T, within=NonNegativeReals)
    model.P_CH4 = Var(model.T, within=NonNegativeReals)
    model.E_CH4 = Var(model.T, within=NonNegativeReals)
    model.P_CH4tNG = Var(model.T, within=NonNegativeReals)
    model.P_NG_PP = Var(model.T, within=NonNegativeReals)
    model.P_NGtP = Var(model.T, within=NonNegativeReals)
    model.P_NG = Var(model.T, within=NonNegativeReals)
    model.P_NK = Var(model.T, within=NonNegativeReals)
    model.P_disp = Var(model.T, within=NonNegativeReals)
    model.P_B = Var(model.T, within=Reals)
    model.P_PtB = Var(model.T, within=NonNegativeReals)
    model.P_BtP = Var(model.T, within=NonNegativeReals)
    model.E_B = Var(model.T, within=NonNegativeReals)
    model.B_B = Var(model.T, within=Binary)
    model.P_PtPH = Var(model.T, within=NonNegativeReals)
    model.P_PHtP = Var(model.T, within=NonNegativeReals)
    model.P_PH = Var(model.T, within=Reals)
    model.B_PH = Var(model.T, within=Binary)
    model.E_PH = Var(model.T, within=NonNegativeReals)
    model.P_trs = Var(model.T, within=Reals)
    model.P_trs_im = Var(model.T, within=NonNegativeReals)
    model.P_trs_ex = Var(model.T, within=NonNegativeReals)
    
    model.V_CO2 = Var(model.T, within = NonNegativeReals)
    
    # Sizing Variables
    model.K_W_on = Var(within=NonNegativeReals)
    model.K_W_off = Var(within=NonNegativeReals)
    model.K_S = Var(within=NonNegativeReals)
    model.K_PtG = Var(within=NonNegativeReals)
    model.S_H2 = Var(within=NonNegativeReals)
    model.K_H2 = Var(within=NonNegativeReals)
    model.K_H2tCH4 = Var(within=NonNegativeReals)
    model.S_CH4 = Var(within=NonNegativeReals)
    model.K_NG = Var(within=NonNegativeReals)
    model.S_B = Var(within=NonNegativeReals)
    model.K_B= Var(within=NonNegativeReals)


    ### Constraints Definition


    # Net Power Balance

    def power_balance_rule(model, t):
        return model.P_S[t] + model.P_W_on[t] + model.P_W_off[t] + model.P_PH[t] + model.P_H2[t] + model.P_NG[t] + model.P_disp[t] +\
             model.P_NK[t] + model.P_B[t] - model.P_C[t] - model.P_PtG[t] + model.Delta_P[t] + model.P_trs[t] == model.pi_L[t]

    model.power_balance = Constraint(model.T, rule=power_balance_rule)

    def unserved_demand_upper_bound_rule(model, t):
        return model.Delta_P[t] <= model.pi_L[t]

    model.unserved_demand_upper_bound = Constraint(model.T, rule=unserved_demand_upper_bound_rule)

    # Variable Renewable Production Technologies Production & Sizing

    def solar_PV_power_output_definition_rule(model, t):
        return model.P_S[t] == model.gamma_S[t] * (model.kappa_S_0 + model.K_S)

    def solar_PV_sizing_upper_bound_rule(model):
        return model.K_S <= model.kappa_S_max - model.kappa_S_0

    def wind_onshore_power_output_definition_rule(model, t):
        return model.P_W_on[t] == model.gamma_W_on[t] * (model.kappa_W_on_0 + model.K_W_on)

    def wind_onshore_sizing_upper_bound_rule(model):
        return model.K_W_on <= model.kappa_W_on_max - model.kappa_W_on_0

    def wind_offshore_power_output_definition_rule(model, t):
        return model.P_W_off[t] == model.gamma_W_off[t] * (model.kappa_W_off_0 + model.K_W_off)

    def wind_offshore_sizing_upper_bound_rule(model):
        return model.K_W_off <= model.kappa_W_off_max - model.kappa_W_off_0

    model.solar_PV_power_output_definition = Constraint(model.T, rule=solar_PV_power_output_definition_rule)
    model.solar_PV_sizing_upper_bound = Constraint(rule=solar_PV_sizing_upper_bound_rule)
    model.wind_onshore_power_output_definition = Constraint(model.T, rule=wind_onshore_power_output_definition_rule)
    model.wind_onshore_sizing_upper_bound = Constraint(rule=wind_onshore_sizing_upper_bound_rule)
    model.wind_offshore_power_output_definition = Constraint(model.T, rule=wind_offshore_power_output_definition_rule)
    model.wind_offshore_sizing_upper_bound = Constraint(rule=wind_offshore_sizing_upper_bound_rule)

    # Variable Renewable Technologies and Curtailment

    def total_curtailment_upper_bound_rule(model, t):
        return model.P_C[t] <= model.P_W_on[t] + model.P_W_off[t] + model.P_S[t]
	
#    def total_curtailment_upper_bound_rule(model, t):
#        return model.P_C[t] <= model.P_W_on[t] + model.P_W_off[t] + model.P_S[t] - model.pi_L[t]
#
#    def total_curtailment_definition_rule(model, t):
#        return model.P_C[t] == model.P_C_p[t] - model.P_C_n[t]
#    
#    def positive_curtailment_upper_bound_rule(model, t):
#        return model.P_C_p[t] <= (model.kappa_S_max + model.kappa_W_on_max + model.kappa_W_off_max) * model.B_C[t]
#    
#    def negative_curtailment_upper_bound_rule(model, t):
#        return model.P_C_n[t] <= (model.kappa_S_max + model.kappa_W_on_max + model.kappa_W_off_max) * (1 - model.B_C[t])

    model.total_curtailment_upper_bound = Constraint(model.T, rule=total_curtailment_upper_bound_rule)
#    model.total_curtailment_definition = Constraint(model.T, rule=total_curtailment_definition_rule)
#    model.positive_curtailment_upper_bound = Constraint(model.T, rule=positive_curtailment_upper_bound_rule)
#    model.negative_curtailment_upper_bound = Constraint(model.T, rule=negative_curtailment_upper_bound_rule)

    # Electrical Interconnection

    def transmission_upper_bound_rule(model, t):
        return model.P_trs[t] <= model.kappa_trs
    
    def transmission_lower_bound_rule(model, t):
        return -model.kappa_trs <= model.P_trs[t]

    def transmission_budget_bound_rule(model):
        return sum(model.P_trs_im[:]) <= model.mu_trs * sum(model.pi_L[:])
    
    def transmission_definition_rule(model, t):
        return model.P_trs[t] == model.P_trs_im[t] - model.P_trs_ex[t]
    
    model.transmission_upper_bound = Constraint(model.T, rule=transmission_upper_bound_rule)
    model.transmission_lower_bound = Constraint(model.T, rule=transmission_lower_bound_rule)
    model.transmission_budget_bound = Constraint(rule=transmission_budget_bound_rule)
    model.transmission_definition = Constraint(model.T, rule=transmission_definition_rule)

    # Pumped-Hydro Plant

    def net_pumped_hydro_power_definition_rule(model, t):
        return model.P_PH[t] == -model.P_PtPH[t] + model.eta_PHtP * model.P_PHtP[t]

    def pumped_hydro_SOC_definition_rule(model, t):
        if t==0:
            return model.E_PH[t] == model.xi_PH
        else:
            return model.E_PH[t] == model.E_PH[t-1] + model.eta_PtPH * model.P_PtPH[t] * model.delta_t - model.P_PHtP[t] * model.delta_t

    def pumped_hydro_SOC_lower_bound_rule(model, t):
        return model.sigma_PH * model.xi_PH <= model.E_PH[t]

    def pumped_hydro_SOC_upper_bound_rule(model, t):
        return model.E_PH[t] <= model.xi_PH

    def pumped_hydro_input_power_upper_bound_rule(model, t):
        return model.P_PtPH[t] <= model.kappa_PtPH * model.B_PH[t]

    def pumped_hydro_output_power_upper_bound_rule(model, t):
        return model.P_PHtP[t] <= model.kappa_PHtP * (1 - model.B_PH[t])

    model.net_pumped_hydro_power_definition = Constraint(model.T, rule=net_pumped_hydro_power_definition_rule)
    model.pumped_hydro_SOC_definition = Constraint(model.T, rule=pumped_hydro_SOC_definition_rule)
    model.pumped_hydro_SOC_lower_bound = Constraint(model.T, rule=pumped_hydro_SOC_lower_bound_rule)
    model.pumped_hydro_SOC_upper_bound = Constraint(model.T, rule=pumped_hydro_SOC_upper_bound_rule)
    model.pumped_hydro_input_power_upper_bound = Constraint(model.T, rule=pumped_hydro_input_power_upper_bound_rule)
    model.pumped_hydro_output_power_upper_bound = Constraint(model.T, rule=pumped_hydro_output_power_upper_bound_rule)
    
    # Battery Plant

    def net_battery_power_definition_rule(model, t):
        return model.P_B[t] == -model.P_PtB[t] + model.eta_BtP * model.P_BtP[t]

    def battery_SOC_definition_rule(model, t):
        if t == 0:
            return model.E_B[t] == 0
        else:
            return model.E_B[t] == model.eta_B * model.E_B[t-1] + model.eta_PtB * model.P_PtB[t] * model.delta_t - model.P_BtP[t] * model.delta_t

    def battery_SOC_lower_bound_rule(model, t):
        return model.sigma_B * model.S_B <= model.E_B[t]

    def battery_SOC_upper_bound_rule(model, t):
        return model.E_B[t] <= model.S_B

    def battery_input_power_upper_bound_rule(model, t):
        return model.eta_PtB * model.P_PtB[t] <= model.rho_B * model.K_B
    
    def battery_output_power_upper_bound_rule(model, t):
        return model.P_BtP[t] <= model.K_B
    
    def battery_input_power_activation_rule(model, t):
        return model.P_PtB[t] <= model.kappa_B * model.B_B[t]
    
    def battery_output_power_activation_rule(model, t):
        return model.P_BtP[t] <= model.kappa_B * (1 - model.B_B[t])
    
    def battery_power_bound_definition_rule(model):
        return model.K_B == model.S_B / model.chi_B
    
    def battery_energy_capacity_upper_bound_rule(model):
        return model.S_B <= model.xi_B

    model.net_battery_power_definition = Constraint(model.T, rule=net_battery_power_definition_rule)
    model.battery_SOC_definition = Constraint(model.T, rule=battery_SOC_definition_rule)
    model.battery_SOC_lower_bound = Constraint(model.T, rule=battery_SOC_lower_bound_rule)
    model.battery_SOC_upper_bound = Constraint(model.T, rule=battery_SOC_upper_bound_rule)
    model.battery_input_power_upper_bound = Constraint(model.T, rule=battery_input_power_upper_bound_rule)
    model.battery_output_power_upper_bound = Constraint(model.T, rule=battery_output_power_upper_bound_rule)
    model.battery_input_power_activation = Constraint(model.T, rule=battery_output_power_activation_rule)
    model.battery_output_power_activation = Constraint(model.T, rule=battery_output_power_activation_rule)
    model.battery_power_bound_definition = Constraint(rule=battery_power_bound_definition_rule)
    model.battery_energy_capacity_upper_bound = Constraint(rule=battery_energy_capacity_upper_bound_rule)

    # Electrolysis Plant

    def electrolysis_power_definition_rule(model, t):
        return model.P_PtH2[t] == model.eta_PtG * model.P_PtG[t]

    def electrolysis_power_upper_bound_rule(model, t):
        return model.P_PtG[t] <= model.K_PtG

    def electrolysis_plant_sizing_upper_bound_rule(model):
        return model.K_PtG <= model.kappa_PtG

    model.electrolysis_power_definition = Constraint(model.T, rule=electrolysis_power_definition_rule)
    model.electrolysis_power_upper_bound = Constraint(model.T, rule=electrolysis_power_upper_bound_rule)
    model.electrolysis_plant_sizing_upper_bound = Constraint(rule=electrolysis_plant_sizing_upper_bound_rule)

    # H2 Storage System

    def H2_storage_output_power_definition_rule(model, t):
        return model.P_H2_out[t] == model.P_H2tP[t] + model.P_H2tCH4[t]

    def H2_storage_SOC_definition_rule(model, t):
        if t==0:
            return model.E_H2[t] == model.sigma_H2 * model.S_H2
        else:
            return model.E_H2[t] == model.E_H2[t-1] + model.P_PtH2[t] * model.delta_t - model.P_H2_out[t] * model.delta_t

    def H2_storage_SOC_lower_bound_rule(model, t):
        return model.sigma_H2 * model.S_H2 <= model.E_H2[t]

    def H2_storage_SOC_upper_bound_rule(model, t):
        return model.E_H2[t] <= model.S_H2

    def H2_storage_sizing_upper_bound_rule(model):
        return model.S_H2 <= model.xi_H2

    model.H2_storage_output_power_definition = Constraint(model.T, rule=H2_storage_output_power_definition_rule)
    model.H2_storage_SOC_definition = Constraint(model.T, rule=H2_storage_SOC_definition_rule)
    model.H2_storage_SOC_lower_bound = Constraint(model.T, rule=H2_storage_SOC_lower_bound_rule)
    model.H2_storage_SOC_upper_bound = Constraint(model.T, rule=H2_storage_SOC_upper_bound_rule)
    model.H2_storage_sizing_upper_bound = Constraint(rule=H2_storage_sizing_upper_bound_rule)

     # H2 Repowering Plant

    def repowered_H2_power_definition_rule(model, t):
        return model.P_H2[t] == model.eta_H2tP * model.P_H2tP[t]
    
    def repowered_H2_power_initialisation_rule(model):
        return model.P_H2[0] == 0

    def repowered_H2_power_upper_bound_rule(model, t):
        return model.P_H2[t] <= model.K_H2

    def H2_repowering_plant_sizing_upper_bound_rule(model):
        return model.K_H2 <= model.kappa_H2

    model.repowered_H2_power_definition = Constraint(model.T, rule=repowered_H2_power_definition_rule)
    model.repowered_H2_power_initialisation = Constraint(rule=repowered_H2_power_initialisation_rule)
    model.repowered_H2_power_upper_bound = Constraint(model.T, rule=repowered_H2_power_upper_bound_rule)
    model.H2_repowering_plant_sizing_upper_bound = Constraint(rule=H2_repowering_plant_sizing_upper_bound_rule)

    # Methanation Plant

    def methanation_output_power_definition_rule(model, t):
        return model.P_CH4[t] == model.eta_H2tCH4 * model.P_H2tCH4[t]

    def methanation_output_power_upper_bound_rule(model, t):
        return model.P_H2tCH4[t] <= model.K_H2tCH4

    def methanation_plant_sizing_upper_bound_rule(model):
        return model.K_H2tCH4 <= model.kappa_H2tCH4

    model.methanation_output_power_definition = Constraint(model.T, rule=methanation_output_power_definition_rule)
    model.methanation_output_power_upper_bound = Constraint(model.T, rule=methanation_output_power_upper_bound_rule)
    model.methanation_plant_sizing_upper_bound = Constraint(rule=methanation_plant_sizing_upper_bound_rule)

    # CH4 Storage System

    def CH4_storage_SOC_definition_rule(model, t):
        if t==0:
            return model.E_CH4[t] == model.sigma_CH4 * model.S_CH4
        else:
            return model.E_CH4[t] == model.E_CH4[t-1] + model.P_CH4[t] * model.delta_t - model.P_CH4tNG[t] * model.delta_t

    def CH4_storage_SOC_lower_bound_rule(model, t):
        return model.sigma_CH4 * model.S_CH4 <= model.E_CH4[t]

    def CH4_storage_SOC_upper_bound_rule(model, t):
        return model.E_CH4[t] <= model.S_CH4

    def CH4_storage_sizing_upper_bound_rule(model):
        return model.S_CH4 <= model.xi_CH4

    def CH4_storage_output_power_initialisation_rule(model):
        return model.P_CH4tNG[0] == 0

    def CH4_storage_output_power_upper_bound_rule(model, t):
        return model.P_CH4tNG[t] <= model.kappa_CH4tNG

    model.CH4_storage_SOC_definition = Constraint(model.T, rule=CH4_storage_SOC_definition_rule)
    model.CH4_storage_SOC_lower_bound = Constraint(model.T, rule=CH4_storage_SOC_lower_bound_rule)
    model.CH4_storage_SOC_upper_bound = Constraint(model.T, rule=CH4_storage_SOC_upper_bound_rule)
    model.CH4_storage_sizing_upper_bound = Constraint(rule=CH4_storage_sizing_upper_bound_rule)
    model.CH4_storage_output_power_initialisation = Constraint(rule=CH4_storage_output_power_initialisation_rule)
    model.CH4_storage_output_power_upper_bound = Constraint(model.T, rule=CH4_storage_output_power_upper_bound_rule)

    # Natural Gas Network (and Storage) System

    def NG_network_balance_rule(model, t):
        return model.pi_NG[t] + model.P_NGtP[t] <= model.kappa_NGNet

    def NG_network_output_power_rule(model, t):
        return model.P_NGtP[t] == model.P_NG_PP[t] + model.P_CH4tNG[t]

    model.NG_network_balance = Constraint(model.T, rule=NG_network_balance_rule)
    model.NG_network_output_power = Constraint(model.T, rule=NG_network_output_power_rule)

    # Gas-Fired Power Plants

    def gas_plant_output_power_definition_rule(model, t):
        return model.P_NG[t] == model.eta_NGtP * model.P_NGtP[t]

    def gas_plant_upper_power_bound_rule(model, t):
        return model.P_NG[t] <= model.kappa_NG_0[t] + model.K_NG

    def gas_plant_sizing_upper_bound_rule(model):
        return model.K_NG <= model.kappa_NG_max

    model.gas_plant_output_power_definition = Constraint(model.T, rule=gas_plant_output_power_definition_rule)
    model.gas_plant_upper_power_bound = Constraint(model.T, rule=gas_plant_upper_power_bound_rule)
    model.gas_plant_sizing_upper_bound = Constraint(rule=gas_plant_sizing_upper_bound_rule)

    # Nuclear Power Plants

    def nuclear_power_lower_bound_rule(model, t):
        return model.mu_NK * model.kappa_NK[t] <= model.P_NK[t]

    def nuclear_power_upper_bound_rule(model, t):
        return model.P_NK[t] <= model.kappa_NK[t]

    def nuclear_units_inc_ramp_rate_rule(model, t):
        if t==0:
            return 0 <= model.P_NK[t]
        else:
            return model.P_NK[t] <= model.P_NK[t-1] + model.delta_NK_inc * model.kappa_NK[t]

    def nuclear_units_dec_ramp_rate_rule(model, t):
        if t==0:
            return 0 <= model.P_NK[t]
        else:
            return model.P_NK[t] >= model.P_NK[t-1] - model.delta_NK_dec * model.kappa_NK[t]

    model.nuclear_power_lower_bound = Constraint(model.T, rule=nuclear_power_lower_bound_rule)
    model.nuclear_power_upper_bound = Constraint(model.T, rule=nuclear_power_upper_bound_rule)
    model.nuclear_units_inc_ramp_rate = Constraint(model.T, rule=nuclear_units_inc_ramp_rate_rule)
    model.nuclear_units_dec_ramp_rate = Constraint(model.T, rule=nuclear_units_dec_ramp_rate_rule)

    # Other Dispatchable Power Plant (Biomass, Waste, CHP)

    def dispatchable_power_lower_bound_rule(model, t):
        return model.mu_disp * model.kappa_disp <= model.P_disp[t]

    def dispatchable_power_upper_bound_rule(model, t):
        return model.P_disp[t] <= model.kappa_disp

    def dispatchable_units_inc_ramp_rate_rule(model, t):
        if t==0:
            return 0 <= model.P_disp[t]
        else:
            return model.P_disp[t] <= model.P_disp[t-1] + model.delta_disp_inc * model.kappa_disp

    def dispatchable_units_dec_ramp_rate_rule(model, t):
        if t==0:
            return 0 <= model.P_disp[t]
        else:
            return model.P_disp[t] >= model.P_disp[t-1] - model.delta_disp_dec * model.kappa_disp

    model.dispatchable_power_lower_bound = Constraint(model.T, rule=dispatchable_power_lower_bound_rule)
    model.dispatchable_power_upper_bound = Constraint(model.T, rule=dispatchable_power_upper_bound_rule)
    model.dispatchable_units_inc_ramp_rate = Constraint(model.T, rule=dispatchable_units_inc_ramp_rate_rule)
    model.dispatchable_units_dec_ramp_rate = Constraint(model.T, rule=dispatchable_units_dec_ramp_rate_rule)

    # CO2 Budget
            
    def yearly_CO2_consumption_rule(model):
        return model.nu_NG_CO2 * sum(model.P_NG_PP[:]) +\
            model.nu_CH4_CO2 * sum(model.P_CH4tNG[:]) +\
            model.nu_trs_CO2 * sum(model.P_trs[:]) +\
            model.nu_disp_CO2 * sum(model.P_disp[:]) <= model.psi_CO2

    model.yearly_CO2_consumption = Constraint(rule=yearly_CO2_consumption_rule)

    ### Objective Function


    def cost_rule(model):
        return (model.zeta_W_on + model.theta_W_on_f * model.n_y) * model.K_W_on +\
            (model.zeta_W_off + model.theta_W_off_f * model.n_y) * model.K_W_off +\
            (model.zeta_S + model.theta_S_f * model.n_y) * model.K_S +\
            model.theta_W_on_v * sum(model.P_W_on[:]) * model.delta_t +\
            model.theta_W_off_v * sum(model.P_W_off[:]) * model.delta_t +\
            model.theta_S_v * sum(model.P_S[:]) * model.delta_t +\
            (model.zeta_PtG + model.theta_PtG_f * model.n_y) * model.K_PtG +\
            (model.zeta_H2_s + model.theta_H2_s_f * model.n_y) * model.S_H2 +\
            (model.zeta_H2 + model.theta_H2_f * model.n_y) * model.K_H2 + model.theta_H2_v * sum(model.P_H2[:]) * model.delta_t +\
            (model.zeta_H2tCH4 + model.theta_H2tCH4_f * model.n_y) * model.K_H2tCH4 +\
            (model.zeta_CH4 + model.theta_CH4_f * model.n_y) * model.S_CH4 +\
            (model.zeta_NG + model.theta_NG_f * model.n_y) * model.K_NG +\
            model.theta_NG_v * sum(model.P_NG[:]) * model.delta_t +\
            model.theta_NG_fuel * sum(model.P_NG_PP[:]) * model.delta_t +\
            model.theta_CO2 * model.nu_NG_CO2 * sum(model.P_NG_PP[:]) * model.delta_t +\
            model.theta_CO2 * model.nu_CH4_CO2 * sum(model.P_CH4tNG[:]) * model.delta_t +\
            model.theta_NK_v * sum(model.P_NK[:]) * model.delta_t +\
            model.theta_NK_fuel * sum(model.P_NK[:]) * model.delta_t / model.eta_NK +\
            model.theta_disp_v * sum(model.P_disp[:]) * model.delta_t +\
            model.theta_disp_fuel * sum(model.P_disp[:]) * model.delta_t / model.eta_disp +\
            model.theta_CO2 * model.nu_disp_CO2 * sum(model.P_disp[:]) * model.delta_t +\
            model.theta_PH_v * sum(model.P_PtPH[:]) * model.delta_t +\
            model.varsigma_ENS * sum(model.Delta_P[:]) * model.delta_t + model.varsigma_C * sum(model.P_C[:]) * model.delta_t +\
            sum(model.theta_el[t] * model.P_trs_im[t] * model.delta_t for t in model.T) -\
            sum(model.theta_el[t] * model.P_trs_ex[t] for t in model.T) +\
            model.zeta_B * model.S_B + model.theta_B_f * model.n_y * model.K_B
            

    model.cost = Objective(rule=cost_rule, sense=minimize)

    dir_run	= folder
	
    if write_lp:
        model.write(filename=join(dir_run, "model.lp"),
            format=ProblemFormat.cpxlp,
            io_options={"symbolic_solver_labels": True})
        model.write(filename=join(dir_run, 'model.mps'),
            format=ProblemFormat.mps)

    return model


def _index(keys, values):
    """
    Builds a dictionary from the lists of keys and values.
    Params:
        keys: the list of keys
        values: the list of values
    Returns:
        A dictionary with corresponding keys and values.
    """
    assert len(keys) == len(values)
    return dict(zip(keys, values))
