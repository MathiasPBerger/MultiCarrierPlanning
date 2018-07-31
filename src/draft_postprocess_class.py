#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:09:50 2018

@author: mathiasberger
"""

class PostProcess(object):
    def __init__(self, data, model):
        self.data = data
        self.model = model
        
    @property
    def CO2_emissions(self):
        model, data = self.model, self.data
        return model.nu_NG_CO2.value * model.eta_NGtP.value * sum([model.P_NG_PP[i].value for i in data.time]) +\
            model.nu_CH4_CO2.value * model.eta_NGtP.value * sum([model.P_CH4tNG[i].value for i in data.time]) +\
            model.nu_trs_CO2.value * sum([model.P_trs[i].value for i in data.time]) +\
            model.nu_disp_CO2.value * sum([model.P_disp[i].value for i in data.time])
    @property
    def CO2_emissions_disp(self):
        model, data = self.model, self.data
        return model.nu_disp_CO2.value * sum([model.P_disp[i].value for i in data.time])
    
    @property
    def CO2_emissions_NG(self):
        model, data = self.model, self.data
        return model.nu_NG_CO2.value * model.eta_NGtP.value * sum([model.P_NG_PP[i].value for i in data.time])
    
    @property
    def CO2_emissions_CH4(self):
        model, data = self.model, self.data
        return model.nu_CH4_CO2.value * model.eta_NGtP.value * sum([model.P_CH4tNG[i].value for i in data.time])
    
    @property
    def CO2_emissions_trs(self):
        model, data = self.model, self.data
        return model.nu_trs_CO2.value * sum([model.P_trs[i].value for i in data.time])
            
    @property
    def capacity_factor_S(self):
        model, time = self.model, self.data.time
        if (model.K_S.value + model.kappa_S_0.value) != 0:
            return sum([model.P_S[t].value for t in time]) / (len(time) * (model.K_S.value + model.kappa_S_0.value))
        else:
            print("No solar capacity installed, no capacity factor can be computed.")
            return 0
            
    @property
    def capacity_factor_W_off(self):
        model, time = self.model, self.data.time
        if (model.K_W_off.value + model.kappa_W_off_0.value):
            return sum([model.P_W_off[t].value for t in time]) / (len(time) * (model.K_W_off.value + model.kappa_W_off_0.value))
        else:
            print("No offshore wind capacity built, no capacity factor can be computed.")
            return 0
    
    @property
    def capacity_factor_W_on(self):
        model, time = self.model, self.data.time
        if (model.K_W_on.value + model.kappa_W_on_0.value):
            return sum([model.P_W_on[t].value for t in time]) / (len(time) * (model.K_W_on.value + model.kappa_W_on_0.value))
        else:
            print("No onshore wind capacity built, no capacity factor can be computed.")
            return 0
    
    @property
    def capacity_factor_PtG(self):
        model, time = self.model, self.data.time
        if model.K_PtG.value != 0:
            return sum([model.P_PtG[t].value for t in time]) / (len(time) * model.K_PtG.value)
        else:
            print("%s GW of power-to-gas (electrolysis) capacity built, no capacity factor can be computed." % model.K_PtG.value)
            return 0
    
    @property
    def capacity_factor_H2(self):
        model, time = self.model, self.data.time
        if model.K_H2.value != 0:
            return sum([model.P_H2[t].value for t in time]) / (len(time) * model.K_H2.value)
        else: 
            print("%s GW of H2-repowering capacity built, no capacity factor can be computed." % model.K_H2.value)
            return 0
            
    @property
    def capacity_factor_H2tCH4(self):
        model, time = self.model, self.data.time
        if model.K_H2tCH4.value != 0:
            return sum([model.P_H2tCH4[t].value for t in time]) / (len(time) * model.K_H2tCH4.value)
        else:
            print("%s GW of methanation capacity built, no capacity factor can be computed." % model.K_H2tCH4.value)
            return 0
    
    @property
    def capacity_factor_CH4tNG(self):
        model, time = self.model, self.data.time
        if max([model.P_CH4tNG[t].value for t in time]) != 0:
            return sum([model.P_CH4tNG[t].value for t in time]) / (len(time) * max([model.P_CH4tNG[t].value for t in time]))
        else:
            print("%s GW of methane to natural gas network pipe capacity built, no capacity factor can be computed." % max([model.P_CH4tNG[t].value for t in time]))
            return 0
        
    @property
    def capacity_factor_H2s(self):
        model, time = self.model, self.data.time
        if model.S_H2.value != 0:
            return sum([model.E_H2[t].value for t in time]) / (len(time) * model.S_H2.value)
        else:
            print("%s GWh of hydrogen storage capacity built, no capacity factor can be computed." % model.S_H2.value)
            return 0
    
    @property
    def capacity_factor_CH4s(self):
        model, time = self.model, self.data.time
        if model.S_CH4.value != 0:
            return sum([model.E_CH4[t].value for t in time]) / (len(time) * model.S_CH4.value)
        else:
            print("%s GWh of methane storage capacity built, no capacity factor can be computed." % model.S_CH4.value)
            return 0
        
    @property
    def capacity_factor_NG(self):
        model, time = self.model, self.data.time
        if model.K_NG.value != 0:
            return sum([model.P_NG[t].value for t in time]) / (len(time) * model.K_NG.value)
        else:
            print("%s GW of natural gas-fired power plants capacity built, no capacity factor can be computed." % model.K_NG.value)
            return 0
        
    @property
    def capacity_factor_trs(self):
        model, time = self.model, self.data.time
        return sum([model.P_trs[t].value for t in time]) / (len(time) * model.kappa_trs.value)
        
    @property
    def n_cycles_H2s(self):
        model, time = self.model, self.data.time
        if model.S_H2.value != 0:
            return sum([model.P_PtH2[t].value for t in time]) / model.S_H2.value
        else:
            print("%s GWh of hydrogen storage capacity built, no number of cycles can be computed." % model.S_H2.value)
            return 0
        
    @property
    def n_cycles_CH4s(self):
        model, time = self.model, self.data.time
        if model.S_CH4.value != 0:
            return sum([model.P_CH4[t].value for t in time]) / model.S_CH4.value
        else:
            print("%s GWh of methane storage capacity built, no number of cycles can be computed." % model.S_CH4.value)
            return 0
        
    @property
    def energy_production_total(self):
        model, time = self.model, self.data.time
        return sum([model.P_S[t].value + model.P_W_on[t].value + model.P_W_off[t].value + model.P_H2[t].value + model.P_disp[t].value + model.P_trs[t].value + model.P_NG[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_production_S(self):
        model, time = self.model, self.data.time
        return sum([model.P_S[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_production_W_on(self):
        model, time = self.model, self.data.time
        return sum([model.P_W_on[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_production_W_off(self):
        model, time = self.model, self.data.time
        return sum([model.P_W_off[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_production_PtG(self):
        model, time = self.model, self.data.time
        return sum([model.P_PtG[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_production_H2(self):
        model, time = self.model, self.data.time
        return sum([model.P_H2[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_production_H2tCH4(self):
        model, time = self.model, self.data.time
        return sum([model.P_H2tCH4[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_production_disp(self):
        model, time = self.model, self.data.time
        return sum([model.P_disp[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_production_trs(self):
        model, time = self.model, self.data.time
        return sum([model.P_trs[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_production_NG(self):
        model, time = self.model, self.data.time
        return sum([model.P_NG[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_not_served(self):
        model, time = self.model, self.data.time
        return sum([model.Delta_P[t].value for t in time]) * model.delta_t.value
    
    @property
    def energy_curtailed(self):
        model, time = self.model, self.data.time
        return sum([model.P_C[t].value for t in time]) * model.delta_t.value
    
    @property
    def CAPEX_W_on(self):
        model = self.model
        return model.zeta_W_on.value * model.K_W_on.value
    
    @property
    def FOM_W_on(self):
        model = self.model
        return model.theta_W_on_f.value * model.n_y.value * model.K_W_on.value
    
    @property
    def OPEX_W_on(self):
        model, time = self.model, self.data.time
        return model.theta_W_on_v.value * sum([model.P_W_on[t].value for t in time]) * model.delta_t.value
    
    @property
    def total_cost_W_on(self):
        return self.CAPEX_W_on + self.FOM_W_on + self.OPEX_W_on
    
    @property
    def CAPEX_W_off(self):
        model = self.model
        return model.zeta_W_off.value * model.K_W_off.value
    
    @property
    def FOM_W_off(self):
        model = self.model
        return model.theta_W_off_f.value * model.n_y.value * model.K_W_off.value
    
    @property
    def OPEX_W_off(self):
        model, time = self.model, self.data.time
        return model.theta_W_off_v.value * sum([model.P_W_off[t].value for t in time]) * model.delta_t.value
    
    @property
    def total_cost_W_off(self):
        return self.CAPEX_W_off + self.FOM_W_off + self.OPEX_W_off
    
    @property
    def CAPEX_S(self):
        model = self.model
        return model.zeta_S.value * model.K_S.value
    
    @property
    def FOM_S(self):
        model = self.model
        return model.theta_S_f.value * model.n_y.value * model.K_S.value
    
    @property
    def OPEX_S(self):
        model, time = self.model, self.data.time
        return model.theta_S_v.value * sum([model.P_S[t].value for t in time]) * model.delta_t.value
    
    @property
    def total_cost_S(self):
        return self.CAPEX_S + self.FOM_S + self.OPEX_S
    
    @property
    def CAPEX_PtG(self):
        model = self.model
        return model.zeta_PtG.value * model.K_PtG.value
    
    @property
    def FOM_PtG(self):
        model = self.model
        return model.theta_PtG_f.value * model.K_PtG.value
    
    #@property
    #def OPEX_PtG(self):
    #    model, time = self.model, self.data.time
    #    return model.theta_PtG_v.value * sum([model.P_PtG[t].value for t in time]) * model.delta_t.value
    
    @property
    def total_cost_PtG(self):
        return self.CAPEX_PtG + self.FOM_PtG #+ self.OPEX_S
    
    @property
    def CAPEX_H2(self):
        model = self.model
        return model.zeta_H2.value * model.K_H2.value
    
    @property
    def FOM_H2(self):
        model = self.model
        return model.theta_H2_f.value * model.K_H2.value
    
    @property
    def OPEX_H2(self):
        model, time = self.model, self.data.time
        return model.theta_H2_v.value * sum(model.P_H2[t].value for t in time) * model.delta_t.value
    
    @property
    def total_cost_H2(self):
        return self.CAPEX_H2 + self.FOM_H2 + self.OPEX_H2
    
    @property
    def CAPEX_H2s(self):
        model = self.model
        return model.zeta_H2_s.value * model.S_H2.value
    
    @property
    def FOM_H2s(self):
        model = self.model
        return model.theta_H2_s_f.value * model.S_H2.value
    
    @property
    def total_cost_H2s(self):
        return self.CAPEX_H2s + self.FOM_H2s
    
    @property
    def CAPEX_H2tCH4(self):
        model = self.model
        return model.zeta_H2tCH4.value * model.K_H2tCH4.value
    
    @property
    def FOM_H2tCH4(self):
        model = self.model
        return model.theta_H2tCH4_f.value * model.K_H2tCH4.value
    
    @property
    def total_cost_H2tCH4(self):
        return self.CAPEX_H2tCH4 + self.FOM_H2tCH4
    
    @property
    def CAPEX_CH4s(self):
        model = self.model
        return model.zeta_CH4.value * model.S_CH4.value
    
    @property
    def FOM_CH4s(self):
        model = self.model
        return model.theta_CH4_f.value * model.S_CH4.value
    
    @property
    def total_cost_CH4s(self):
        return self.CAPEX_CH4s + self.FOM_CH4s
    
    @property
    def CAPEX_NG(self):
        model = self.model
        return  model.zeta_NG.value * model.K_NG.value
    
    @property
    def FOM_NG(self):
        model = self.model
        return model.theta_PtG_f.value * model.K_NG.value
    
    @property
    def OPEX_NG(self):
        model, time = self.model, self.data.time
        return model.theta_NG_v.value * sum([model.P_NG[t].value for t in time]) * model.delta_t.value +\
            model.theta_NG_fuel.value * sum([model.P_NG_PP[t].value for t in time]) * model.delta_t.value +\
            model.theta_CO2.value * model.nu_NG_CO2.value * model.eta_NGtP.value * sum([model.P_NG_PP[t].value for t in time]) * model.delta_t.value +\
            model.theta_CO2.value * model.nu_CH4_CO2.value * model.eta_NGtP.value * sum([model.P_CH4tNG[t].value for t in time]) * model.delta_t.value
            
    @property
    def CO2_costs_NG(self):
        model, time = self.model, self.data.time
        return model.theta_CO2.value * model.nu_NG_CO2.value * model.eta_NGtP.value * sum([model.P_NG_PP[t].value for t in time]) * model.delta_t.value
    
    @property
    def total_cost_NG(self):
        return self.CAPEX_NG + self.FOM_NG + self.OPEX_NG
    
    @property
    def OPEX_disp(self):
        model, time = self.model, self.data.time
        return (model.theta_disp_v.value + model.theta_disp_fuel.value / model.eta_disp +\
            model.theta_CO2.value * model.nu_disp_CO2.value) * sum([model.P_disp[t].value for t in time]) * model.delta_t.value
                
    @property
    def CO2_costs_disp(self):
        model, time = self.model, self.data.time
        return model.theta_CO2.value * model.nu_disp_CO2.value * sum([model.P_disp[t].value for t in time]) * model.delta_t.value
    
    @property
    def total_cost_disp(self):
        return self.OPEX_disp + self.CO2_costs_disp
    
    @property
    def OPEX_PH(self):
        model, time = self.model, self.data.time
        return model.theta_PH_v.value * sum([model.P_PtPH[t].value for t in time]) * model.delta_t.value
    
    @property
    def OPEX_ENS(self):
        model, time = self.model, self.data.time
        return model.varsigma_ENS.value * sum([model.Delta_P[t].value for t in time]) * model.delta_t.value
    
    @property
    def OPEX_C(self):
        model, time = self.model, self.data.time
        return model.varsigma_C.value * sum([model.P_C[t].value for t in time]) * model.delta_t.value
    
    @property
    def total_cost_trs(self):
        model, theta_el, time = self.model, self.data.theta_el, self.data.time
        return sum([theta_el[t] * model.P_trs[t].value for t in time]) * model.delta_t.value
    
    @property
    def total_cost(self):
        return self.total_cost_S + self.total_cost_W_on + self.total_cost_W_off + self.total_cost_PtG + self.total_cost_H2s + self.total_cost_H2tCH4 + self.total_cost_CH4s + self.total_cost_disp + self.total_cost_NG + self.OPEX_PH + self.OPEX_ENS + self.OPEX_C + self.total_cost_trs
    
    @property
    def budget_trs(self):
        model, time = self.model, self.data.time
        if sum([model.P_trs[t].value for t in time]) != 0:
            return sum([model.P_trs[t].value for t in time]) / sum([model.pi_L[t] for t in time])
        else:
            print("%s energy transmitted, no interconnection energy budget can be computed." % sum([model.P_trs[t].value for t in time]))
            return 0
        
    @property
    def total_cost_GWh(self):
        return self.total_cost / self.energy_production_total
    
    @property
    def cost_GWh_S(self):
        return self.total_cost_S / self.energy_production_S
    
    @property
    def cost_GWh_W_on(self):
        return self.total_cost_W_on / self.energy_production_W_on
    
    @property
    def cost_GWh_W_off(self):
        return self.total_cost_W_off / self.energy_production_W_off
    
    @property
    def cost_GWh_PtG(self):
        if self.energy_production_PtG !=0:
            return self.total_cost_PtG / self.energy_production_PtG
        else:
            print("No electrolyser built.")
            return 0
    
#    @property
#    def cost_GWh_H2s(self):
#        return self.total_cost_H2s / self.energy_production_H2s
    
    @property
    def cost_GWh_H2tCH4(self):
        if self.energy_production_H2tCH4 !=0:
            return self.total_cost_H2tCH4 / self.energy_production_H2tCH4
        else:
            print("No methanation built.")
            return 0
    
#    @property
#    def cost_GWH_CH4s(self):
#        return self.total_cost_CH4s / self.energy_production_CH4s
    
    @property
    def cost_GWh_H2(self):
        if self.energy_production_H2 !=0:
            return self.total_cost_H2 / self.energy_production_H2
        else:
            print("No hydrogen repowering plant built.")
            return 0
    
    @property
    def cost_GWh_disp(self):
        if self.energy_production_disp !=0:
            return self.total_cost_disp / self.energy_production_disp
        else:
            print("No energy produced by dispatchable units.")
            return 0
    
    @property
    def cost_GWh_NG(self):
        if self.energy_production_NG !=0:
            return self.total_cost_NG / self.energy_production_NG
        else:
            print("No natural gas-fired power plants built.")
            return 0
    
    @property
    def cost_GWh_trs(self):
        if self.energy_production_trs !=0:
            return self.total_cost_trs / self.energy_production_trs
        else:
            print("No electricity imported.")
            return 0
        