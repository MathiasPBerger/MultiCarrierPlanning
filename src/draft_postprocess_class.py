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
        model, time = self.model, self.data.time
        return model.nu_NG_CO2.value * sum(model.lambda_NG_HT[t] + model.lambda_NG_ID[t] + model.lambda_NG_TR[t] - model.lambda_NGtH2[t] - model.L_NG_ENS[t].value for t in time) * model.delta_t.value + \
               sum(model.Q_CO2_OCGT[t].value for t in time) * model.delta_t.value + sum(model.Q_CO2_CCGT[t].value for t in time) * model.delta_t.value + \
               sum(model.Q_CO2_CHP[t].value for t in time) * model.delta_t.value + sum(model.Q_CO2_SMR[t].value for t in time) * model.delta_t.value + \
               sum(model.Q_CO2_BM[t].value for t in time) * model.delta_t.value + sum(model.Q_CO2_WS[t].value for t in time) * model.delta_t.value + \
               model.nu_E_I_CO2.value * sum(model.P_E_I[t].value for t in time) * model.delta_t.value - sum(model.Q_CO2_ACC[t].value for t in time) * model.delta_t.value

    @property
    def CO2_emissions_BM(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_BM[t].value for t in data.time) * model.delta_t.value

    @property
    def CO2_emissions_BM_CCS(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_BM_CCS[t].value for t in data.time) * model.delta_t.value
    
    @property
    def CO2_emissions_CHP(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_CHP[t].value for t in data.time) * model.delta_t.value
    
    @property
    def CO2_emissions_CHP_CCS(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_CHP_CCS[t].value for t in data.time) * model.delta_t.value
    
    @property
    def CO2_emissions_WS(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_WS[t].value for t in data.time) * model.delta_t.value
    
    @property
    def CO2_emissions_WS_CCS(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_WS_CCS[t].value for t in data.time) * model.delta_t.value
    
    @property
    def CO2_emissions_OCGT(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_OCGT[t].value for t in data.time) * model.delta_t.value
    
    @property
    def CO2_emissions_OCGT_CCS(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_OCGT_CCS[t].value for t in data.time) * model.delta_t.value
    
    @property
    def CO2_emissions_CCGT(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_CCGT[t].value for t in data.time) * model.delta_t.value
    
    @property
    def CO2_emissions_CCGT_CCS(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_CCGT_CCS[t].value for t in data.time) * model.delta_t.value

    @property
    def CO2_emissions_SMR(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_SMR[t].value for t in data.time) * model.delta_t.value

    @property
    def CO2_emissions_SMR_CCS(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_SMR_CCS[t].value for t in data.time) * model.delta_t.value

    @property
    def CO2_emissions_NG_HT(self):
        model, data = self.model, self.data
        return model.nu_NG_CO2.value * sum([model.lambda_NG_HT[t] for t in data.time]) * model.delta_t.value

    @property
    def CO2_emissions_NG_ID(self):
        model, data = self.model, self.data
        return model.nu_NG_CO2.value * sum([model.lambda_NG_ID[t] for t in data.time]) * model.delta_t.value

    @property
    def CO2_emissions_NG_TR(self):
        model, data = self.model, self.data
        return model.nu_NG_CO2.value * sum([model.lambda_NG_TR[t] for t in data.time]) * model.delta_t.value
    
    @property
    def CO2_consumption_MT(self):
        model, data = self.model, self.data
        return sum([model.Q_CO2_MT[t].value for t in data.time]) * model.delta_t.value

    @property
    def CO2_absorbtion_ACC(self):
        model, data = self.model, self.data
        return sum(model.Q_CO2_ACC[t].value for t in data.time) * model.delta_t.value

    @property
    def CO2_capture_CCS(self):
        model, data = self.model, self.data
        return (sum(model.Q_CO2_OCGT_CCS[t].value + model.Q_CO2_CCGT_CCS[t].value + model.Q_CO2_WS_CCS[t].value + \
               model.Q_CO2_BM_CCS[t].value + model.Q_CO2_CHP_CCS[t].value + model.Q_CO2_SMR_CCS[t].value for t in data.time)) * model.delta_t.value

    @property
    def capacity_factor_S(self):
        model, time = self.model, self.data.time
        if (model.K_S.value + model.kappa_S_0.value) != 0:
            return sum([model.P_S[t].value for t in time]) / (len(time) * (model.K_S.value + model.kappa_S_0.value))
        else:
            return 0

    @property
    def capacity_factor_W_off(self):
        model, time = self.model, self.data.time
        if (model.K_W_off.value + model.kappa_W_off_0.value):
            return sum([model.P_W_off[t].value for t in time]) / (len(time) * (model.K_W_off.value + model.kappa_W_off_0.value))
        else:
            return 0
    
    @property
    def capacity_factor_W_on(self):
        model, time = self.model, self.data.time
        if (model.K_W_on.value + model.kappa_W_on_0.value):
            return sum([model.P_W_on[t].value for t in time]) / (len(time) * (model.K_W_on.value + model.kappa_W_on_0.value))
        else:
            return 0
    
    @property
    def capacity_factor_EL(self):
        model, time = self.model, self.data.time
        if model.K_E_EL.value != 0:
            return sum([model.P_E_EL[t].value for t in time]) / (len(time) * model.K_E_EL.value)
        else:
            return 0
    
    @property
    def capacity_factor_FC(self):
        model, time = self.model, self.data.time
        if model.K_E_FC.value != 0:
            return sum([model.P_E_FC[t].value for t in time]) / (len(time) * model.K_E_FC.value)
        else:
            return 0
            
    @property
    def capacity_factor_MT(self):
        model, time = self.model, self.data.time
        if model.K_CH4_MT.value != 0:
            return sum([model.P_CH4_MT[t].value for t in time]) / (len(time) * model.K_CH4_MT.value)
        else:
            return 0

    @property
    def capacity_factor_H2S(self):
        model, time = self.model, self.data.time
        if model.S_H2.value != 0:
            return sum([model.E_H2[t].value for t in time]) / (len(time) * model.S_H2.value)
        else:
            return 0
        
    @property
    def capacity_factor_NG_OCGT(self):
        model, time = self.model, self.data.time
        if model.K_E_OCGT.value != 0:
            return sum([model.P_OCGT[t].value for t in time]) / (len(time) * model.K_E_OCGT.value)
        else:
            return 0
        
    @property
    def capacity_factor_NG_CCGT(self):
        model, time = self.model, self.data.time
        if model.K_E_CCGT.value != 0:
            return sum([model.P_CCGT[t].value for t in time]) / (len(time) * model.K_E_CCGT.value)
        else:
            return 0
        
    @property
    def capacity_factor_BM(self):
        model, time = self.model, self.data.time
        if model.kappa_BM.value != 0:
            return sum([model.P_BM[t].value for t in time]) / (len(time) * model.kappa_BM.value)
        else:
            return 0
        
    @property
    def capacity_factor_CHP(self):
        model, time = self.model, self.data.time
        if model.kappa_CHP.value != 0:
            return sum([model.P_CHP[t].value for t in time]) / (len(time) * model.kappa_CHP.value)
        else:
            return 0
        
    @property
    def capacity_factor_WS(self):
        model, time = self.model, self.data.time
        if model.kappa_WS.value != 0:
            return sum([model.P_WS[t].value for t in time]) / (len(time) * model.kappa_WS.value)
        else:
            return 0

    @property
    def capacity_factor_NK(self):
        model, time = self.model, self.data.time
        if model.kappa_NK.value != 0:
            return sum([model.P_E_NK[t].value for t in time]) / (len(time) * model.kappa_NK.value)
        else:
            return 0
        
    @property
    def capacity_factor_SMR(self):
        model, time = self.model, self.data.time
        if model.K_H2_SMR.value != 0:
            return sum([model.P_H2_SMR[t].value for t in time]) / (len(time) * model.K_H2_SMR.value)
        else:
            return 0
        
    @property
    def capacity_factor_imports_electricity(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_I[t].value for t in time]) / (len(time) * model.kappa_IE.value)
    
    @property
    def capacity_factor_exports_electricity(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_E[t].value for t in time]) / (len(time) * model.kappa_IE.value)
    
    @property
    def curtailment_solar_PV(self):
        model, time = self.model, self.data.time
        P_S, kappa_0, K_S, gamma_S = [model.P_S[t].value for t in time], self.data.kappa_S_0, model.K_S.value, [model.gamma_S[t] for t in time]
        return sum([(kappa_0 + K_S) * gamma_S[t] for t in time]) - sum([P_S[t] for t in time])
    
    @property
    def curtailment_wind_onshore(self):
        model, time = self.model, self.data.time
        P, kappa_0, K, gamma = [model.P_W_on[t].value for t in time], self.data.kappa_W_on_0, model.K_W_on.value, [model.gamma_W_on[t] for t in time]
        return sum([(kappa_0 + K) * gamma[t] for t in time]) - sum([P[t] for t in time])
    
    @property
    def curtailment_wind_offshore(self):
        model, time = self.model, self.data.time
        P, kappa_0, K, gamma = [model.P_W_off[t].value for t in time], self.data.kappa_W_off_0, model.K_W_off.value, [model.gamma_W_off[t] for t in time]
        return sum([(kappa_0 + K) * gamma[t] for t in time]) - sum([P[t] for t in time])
    
    @property
    def total_curtailment(self):
        return self.curtailment_solar_PV + self.curtailment_wind_onshore + self.curtailment_wind_offshore
    
    @property
    def curtailment_percentage(self):
        return (self.total_curtailment / self.electricity_production_RES) * 100
     
    @property
    def curtailment_solar_PV_hours(self):
        model, time = self.model, self.data.time
        dual_prod_S = [model.dual[model.solar_PV_power_output_definition[t]] for t in time]
        return sum([1 for t in time if dual_prod_S[t] == 0.]) / self.data.n_y
    
    @property
    def curtailment_wind_onshore_hours(self):
        model, time = self.model, self.data.time
        dual_prod_W_on = [model.dual[model.wind_onshore_power_output_definition[t]] for t in time]
        return sum([1 for t in time if dual_prod_W_on[t] == 0.]) / self.data.n_y
    
    @property
    def curtailment_wind_offshore_hours(self):
        model, time = self.model, self.data.time
        dual_prod_W_off = [model.dual[model.wind_offshore_power_output_definition[t]] for t in time]
        return sum([1 for t in time if dual_prod_W_off[t] == 0.]) / self.data.n_y
        

    @property
    def n_cycles_H2S(self):
        model, time = self.model, self.data.time
        if model.S_H2.value != 0:
            return sum([model.P_H2_EL[t].value for t in time]) / model.S_H2.value
        else:
            return 0
        
    @property
    def n_cycles_B(self):
        model, time = self.model, self.data.time
        if model.S_B.value != 0:
            return sum([model.P_PtB[t].value for t in time]) / model.S_B.value
        else:
            return 0

    @property
    def n_cycles_NGS(self):
        model, time = self.model, self.data.time
        if model.xi_NGS.value != 0:
            return sum([model.P_NGtNGS[t].value for t in time]) / model.xi_NGS.value
        else:
            return 0

    @property
    def electricity_production_total(self):
        model, time = self.model, self.data.time
        return sum([model.P_S[t].value + model.P_W_on[t].value + model.P_W_off[t].value + model.P_E_FC[t].value + \
                    model.P_E_BM[t].value + model.P_E_WS[t].value + model.P_E_CHP[t].value + model.P_E_OCGT[t].value + \
                    model.P_E_CCGT[t].value + model.P_E_NK[t].value for t in time]) * model.delta_t.value
    
    @property
    def electricity_production_S(self):
        model, time = self.model, self.data.time
        return sum([model.P_S[t].value for t in time]) * model.delta_t.value
    
    @property
    def electricity_production_W_on(self):
        model, time = self.model, self.data.time
        return sum([model.P_W_on[t].value for t in time]) * model.delta_t.value
    
    @property
    def electricity_production_W_off(self):
        model, time = self.model, self.data.time
        return sum([model.P_W_off[t].value for t in time]) * model.delta_t.value
    
    @property
    def electricity_production_RES(self):
        return self.electricity_production_S + self.electricity_production_W_on + self.electricity_production_W_off
    
    @property
    def electricity_production_surplus_RES(self):
        model, time = self.model, self.data.time
        mismatch = [model.P_S[t].value + model.P_W_on[t].value + model.P_W_off[t].value - model.lambda_E[t] - model.lambda_E_HT[t] - model.L_E_TR[t].value for t in time]
        surplus = sum([value for value in mismatch if value > 0])
        return surplus

    @property
    def electricity_production_FC(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_FC[t].value for t in time]) * model.delta_t.value

    @property
    def net_electricity_production_BM(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_BM[t].value for t in time]) * model.delta_t.value

    @property
    def net_electricity_production_WS(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_WS[t].value for t in time]) * model.delta_t.value

    @property
    def net_electricity_production_CHP(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_CHP[t].value for t in time]) * model.delta_t.value

    @property
    def net_electricity_production_NG_OCGT(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_OCGT[t].value for t in time]) * model.delta_t.value

    @property
    def net_electricity_production_NG_CCGT(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_CCGT[t].value for t in time]) * model.delta_t.value

    @property
    def electricity_consumption_EL(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_EL[t].value for t in time]) * model.delta_t.value
    
    @property
    def water_consumption_EL(self):
        model, time = self.model, self.data.time
        return sum(model.m_H2O[t].value for t in time) * model.delta_t.value
    
    @property
    def oxygen_production_EL(self):
        model, time = self.model, self.data.time
        return sum(model.m_O2[t].value for t in time) * model.delta_t.value

    @property
    def hydrogen_production_EL(self):
        model, time = self.model, self.data.time
        return sum([model.P_H2_EL[t].value for t in time]) * model.delta_t.value
    
    @property
    def hydrogen_consumption_MT(self):
        model, time = self.model, self.data.time
        return sum([model.P_H2_MT[t].value for t in time]) * model.delta_t.value
    
    @property
    def methane_production_MT(self):
        model, time = self.model, self.data.time
        return sum([model.P_CH4_MT[t].value for t in time]) * model.delta_t.value
    
    @property
    def electricity_consumption_CCS_BM(self):
        model, time = self.model, self.data.time
        return sum([model.P_BM_CCS[t].value for t in time]) * model.delta_t.value
    
    @property
    def electricity_consumption_CCS_WS(self):
        model, time = self.model, self.data.time
        return sum([model.P_WS_CCS[t].value for t in time]) * model.delta_t.value
    
    @property
    def electricity_consumption_CCS_CHP(self):
        model, time = self.model, self.data.time
        return sum([model.P_CHP_CCS[t].value for t in time]) * model.delta_t.value

    @property
    def electricity_consumption_CCS_OCGT(self):
        model, time = self.model, self.data.time
        return sum([model.P_OCGT_CCS[t].value for t in time]) * model.delta_t.value
    
    @property
    def electricity_consumption_CCS_CCGT(self):
        model, time = self.model, self.data.time
        return sum([model.P_CCGT_CCS[t].value for t in time]) * model.delta_t.value

    @property
    def electricity_consumption_CCS_SMR(self):
        model, time = self.model, self.data.time
        return sum([model.P_SMR_CCS[t].value for t in time]) * model.delta_t.value

    @property
    def electricity_consumption_ACC(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_ACC[t].value for t in time]) * model.delta_t.value

    @property
    def NG_consumption_ACC(self):
        model, time = self.model, self.data.time
        return sum([model.P_NG_ACC[t].value for t in time]) * model.delta_t.value

    @property
    def electricity_annual_load(self):
        model = self.model
        return sum(model.lambda_E[:]) * model.delta_t.value

    @property
    def electricity_annual_load_HT(self):
        model = self.model
        return sum(model.lambda_E_HT[:]) * model.delta_t.value

    @property
    def electricity_annual_load_TR(self):
        model = self.model
        return sum(model.L_E_TR[:].value) * model.delta_t.value
    
    @property
    def electricity_imports(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_I[t].value for t in time]) * model.delta_t.value
    
    @property
    def electricity_exports(self):
        model, time = self.model, self.data.time
        return sum([model.P_E_E[t].value for t in time]) * model.delta_t.value

    @property
    def hydrogen_annual_load_ID(self):
        model = self.model
        return sum(model.lambda_H2_ID[:]) * model.delta_t.value

    @property
    def hydrogen_annual_load_TR(self):
        model = self.model
        return sum(model.lambda_H2_TR[:]) * model.delta_t.value

    @property
    def hydrogen_imports(self):
        model, time = self.model, self.data.time
        return sum([model.P_H2_I[t].value for t in time]) * model.delta_t.value
    
    @property
    def NG_imports(self):
        model, time = self.model, self.data.time
        return sum([model.P_NG_I[t].value for t in time]) * model.delta_t.value

    @property
    def NG_annual_load_TR(self):
        model, time = self.model, self.data.time
        return sum([model.lambda_NG_TR[t] for t in time]) * model.delta_t.value

    @property
    def NG_annual_load_HT(self):
        model, time = self.model, self.data.time
        return sum([model.lambda_NG_HT[t] for t in time]) * model.delta_t.value

    @property
    def NG_annual_load_ID(self):
        model, time = self.model, self.data.time
        return sum([model.lambda_NG_ID[t] for t in time]) * model.delta_t.value

    @property
    def NG_annual_load_SRM(self):
        model, time = self.model, self.data.time
        return sum([model.lambda_NGtH2[t] for t in time]) * model.delta_t.value

    @property
    def ENS_E(self):
        model, time = self.model, self.data.time
        return sum([model.L_E_ENS[t].value for t in time]) * model.delta_t.value

    @property
    def ENS_NG(self):
        model, time = self.model, self.data.time
        return sum([model.L_NG_ENS[t].value for t in time]) * model.delta_t.value

    @property
    def ENS_H2(self):
        model, time = self.model, self.data.time
        return sum([model.L_H2_ENS[t].value for t in time]) * model.delta_t.value

    @property
    def CAPEX_W_on(self):
        model = self.model
        return model.zeta_W_on.value * model.K_W_on.value

    @property
    def FOM_W_on(self):
        model = self.model
        return model.theta_W_on_f.value * model.n_y.value * model.K_W_on.value

    @property
    def VOM_W_on(self):
        model, time = self.model, self.data.time
        return model.theta_W_on_v.value * sum([model.P_W_on[t].value for t in time]) * model.delta_t.value

    @property
    def total_cost_W_on(self):
        return self.CAPEX_W_on + self.FOM_W_on + self.VOM_W_on

    @property
    def CAPEX_W_off(self):
        model = self.model
        return model.zeta_W_off.value * model.K_W_off.value

    @property
    def FOM_W_off(self):
        model = self.model
        return model.theta_W_off_f.value * model.n_y.value * model.K_W_off.value

    @property
    def VOM_W_off(self):
        model, time = self.model, self.data.time
        return model.theta_W_off_v.value * sum([model.P_W_off[t].value for t in time]) * model.delta_t.value

    @property
    def total_cost_W_off(self):
        return self.CAPEX_W_off + self.FOM_W_off + self.VOM_W_off

    @property
    def CAPEX_S(self):
        model = self.model
        return model.zeta_S.value * model.K_S.value

    @property
    def FOM_S(self):
        model = self.model
        return model.theta_S_f.value * model.n_y.value * model.K_S.value

    @property
    def VOM_S(self):
        model, time = self.model, self.data.time
        return model.theta_S_v.value * sum([model.P_S[t].value for t in time]) * model.delta_t.value

    @property
    def total_cost_S(self):
        return self.CAPEX_S + self.FOM_S + self.VOM_S

    @property
    def CAPEX_EL(self):
        model = self.model
        return model.zeta_EL.value * model.K_E_EL.value

    @property
    def FOM_EL(self):
        model = self.model
        return model.theta_EL_f.value * model.K_E_EL.value * model.n_y.value

    @property
    def total_cost_EL(self):
        return self.CAPEX_EL + self.FOM_EL

    @property
    def CAPEX_FC(self):
        model = self.model
        return model.zeta_FC.value * model.K_E_FC.value

    @property
    def FOM_FC(self):
        model = self.model
        return model.theta_FC_f.value * model.K_E_FC.value * model.n_y.value

    @property
    def VOM_FC(self):
        model, time = self.model, self.data.time
        return model.theta_FC_v.value * sum(model.P_E_FC[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_FC(self):
        return self.CAPEX_FC + self.FOM_FC + self.VOM_FC

    @property
    def CAPEX_H2S(self):
        model = self.model
        return model.zeta_H2S.value * model.S_H2.value

    @property
    def FOM_H2S(self):
        model = self.model
        return model.theta_H2S_f.value * model.S_H2.value * model.n_y.value

    @property
    def total_cost_H2S(self):
        return self.CAPEX_H2S + self.FOM_H2S

    @property
    def FOM_NGS(self):
        model = self.model
        return model.theta_NGS_f.value * model.xi_NGS.value * model.n_y.value

    @property
    def total_cost_NGS(self):
        return self.FOM_NGS

    @property
    def CAPEX_MT(self):
        model = self.model
        return model.zeta_MT.value * model.K_CH4_MT.value

    @property
    def FOM_MT(self):
        model = self.model
        return model.theta_MT_f.value * model.K_CH4_MT.value * model.n_y.value

    @property
    def total_cost_MT(self):
        return self.CAPEX_MT + self.FOM_MT

    @property
    def CAPEX_B_e(self):
        model = self.model
        return model.zeta_B_E.value * model.S_B.value

    @property
    def FOM_B_e(self):
        model = self.model
        return model.theta_B_E_f.value * model.S_B.value * model.n_y.value

    @property
    def total_cost_B_e(self):
        return self.CAPEX_B_e + self.FOM_B_e

    @property
    def CAPEX_B_p(self):
        model = self.model
        return model.zeta_B_P.value * model.K_B.value

    @property
    def FOM_B_p(self):
        model = self.model
        return model.theta_B_P_f.value * model.K_B.value * model.n_y.value

    @property
    def total_cost_B_p(self):
        return self.CAPEX_B_p + self.FOM_B_p

    @property
    def total_cost_B(self):
        return self.total_cost_B_e + self.total_cost_B_p

    @property
    def CAPEX_OCGT(self):
        model = self.model
        return  model.zeta_OCGT.value * model.K_E_OCGT.value

    @property
    def FOM_OCGT(self):
        model = self.model
        return model.theta_OCGT_f.value * model.K_E_OCGT.value * model.n_y.value

    @property
    def VOM_OCGT(self):
        model, time = self.model, self.data.time
        return model.theta_OCGT_v.value * sum([model.P_OCGT[t].value for t in time]) * model.delta_t.value

    @property
    def estimated_fuel_costs_NG_OCGT(self):
        model, time = self.model, self.data.time
        return sum(model.P_NG_OCGT[t].value * model.theta_NG_I[t] for t in time) * model.delta_t.value

    @property
    def CO2_cost_OCGT(self):
        model, time = self.model, self.data.time
        return model.theta_CO2.value * sum(model.Q_CO2_OCGT[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_OCGT(self):
        return self.CAPEX_OCGT + self.FOM_OCGT + self.VOM_OCGT + self.CO2_cost_OCGT

    @property
    def CAPEX_CCS_OCGT(self):
        model = self.model
        return model.zeta_NG_CCS.value * model.K_CO2_OCGT_CCS.value

    @property
    def FOM_CCS_OCGT(self):
        model = self.model
        return model.theta_NG_CCS_f.value * model.n_y.value * model.K_CO2_OCGT_CCS.value

    @property
    def VOM_CCS_OCGT(self):
        model, time = self.model, self.data.time
        return model.theta_NG_CCS_v.value * sum(model.Q_CO2_OCGT_CCS[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_OCGT_CCS(self):
        return self.CAPEX_CCS_OCGT + self.FOM_CCS_OCGT + self.VOM_CCS_OCGT

    @property
    def total_cost_OCGT_with_CCS(self):
        return self.total_cost_OCGT + self.total_cost_OCGT_CCS

    @property
    def CAPEX_CCGT(self):
        model = self.model
        return  model.zeta_CCGT.value * model.K_E_CCGT.value

    @property
    def FOM_CCGT(self):
        model = self.model
        return model.theta_CCGT_f.value * model.K_E_CCGT.value * model.n_y.value

    @property
    def VOM_CCGT(self):
        model, time = self.model, self.data.time
        return model.theta_CCGT_v.value * sum([model.P_CCGT[t].value for t in time]) * model.delta_t.value

    @property
    def estimated_fuel_costs_NG_CCGT(self):
        model, time = self.model, self.data.time
        return sum(model.P_NG_CCGT[t].value * model.theta_NG_I[t] for t in time) * model.delta_t.value

    @property
    def CO2_cost_CCGT(self):
        model, time = self.model, self.data.time
        return model.theta_CO2.value * sum(model.Q_CO2_CCGT[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_CCGT(self):
        return self.CAPEX_CCGT + self.FOM_CCGT + self.VOM_CCGT + self.CO2_cost_CCGT

    @property
    def CAPEX_CCS_CCGT(self):
        model = self.model
        return model.zeta_NG_CCS.value * model.K_CO2_CCGT_CCS.value

    @property
    def FOM_CCS_CCGT(self):
        model = self.model
        return model.theta_NG_CCS_f.value * model.n_y.value * model.K_CO2_CCGT_CCS.value

    @property
    def VOM_CCS_CCGT(self):
        model, time = self.model, self.data.time
        return model.theta_NG_CCS_v.value * sum(model.Q_CO2_CCGT_CCS[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_CCGT_CCS(self):
        return self.CAPEX_CCS_CCGT + self.FOM_CCS_CCGT + self.VOM_CCS_CCGT

    @property
    def total_cost_CCGT_with_CCS(self):
        return self.total_cost_CCGT + self.total_cost_CCGT_CCS

    @property
    def FOM_CHP(self):
        model, data = self.model, self.data
        return data.theta_CHP_f * model.kappa_CHP.value * model.n_y.value

    @property
    def VOM_CHP(self):
        model, time = self.model, self.data.time
        return model.theta_CHP_v.value * sum([model.P_CHP[t].value for t in time]) * model.delta_t.value

    @property
    def estimated_fuel_costs_CHP(self):
        model, time = self.model, self.data.time
        return sum(model.P_NG_CHP[t].value * model.theta_NG_I[t] for t in time) * model.delta_t.value

    @property
    def CO2_cost_CHP(self):
        model, time = self.model, self.data.time
        return model.theta_CO2.value * sum(model.Q_CO2_CHP[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_CHP(self):
        return self.FOM_CHP + self.VOM_CHP + self.CO2_cost_CHP

    @property
    def CAPEX_CCS_CHP(self):
        model = self.model
        return model.zeta_CHP_CCS.value * model.K_CO2_CHP_CCS.value

    @property
    def FOM_CCS_CHP(self):
        model = self.model
        return model.theta_CHP_CCS_f.value * model.n_y.value * model.K_CO2_CHP_CCS.value

    @property
    def VOM_CCS_CHP(self):
        model, time = self.model, self.data.time
        return model.theta_CHP_CCS_v.value * sum(model.Q_CO2_CHP_CCS[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_CHP_CCS(self):
        return self.CAPEX_CCS_CHP + self.FOM_CCS_CHP + self.VOM_CCS_CHP

    @property
    def total_cost_CHP_with_CCS(self):
        return self.total_cost_CHP + self.total_cost_CHP_CCS

    @property
    def FOM_BM(self):
        model, data = self.model, self.data
        return data.theta_BM_f * model.kappa_BM.value * model.n_y.value

    @property
    def VOM_BM(self):
        model, time = self.model, self.data.time
        return model.theta_BM_v.value * sum([model.P_BM[t].value for t in time]) * model.delta_t.value

    @property
    def fuel_costs_BM(self):
        model, time = self.model, self.data.time
        return model.theta_BM_fuel.value * sum(model.P_BM[t].value for t in time) * model.delta_t / model.eta_BM.value

    @property
    def CO2_cost_BM(self):
        model, time = self.model, self.data.time
        return model.theta_CO2.value * sum(model.Q_CO2_BM[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_BM(self):
        return self.FOM_BM + self.VOM_BM + self.fuel_costs_BM + self.CO2_cost_BM

    @property
    def CAPEX_CCS_BM(self):
        model = self.model
        return model.zeta_BM_CCS.value * model.K_CO2_BM_CCS.value

    @property
    def FOM_CCS_BM(self):
        model = self.model
        return model.theta_BM_CCS_f.value * model.n_y.value * model.K_CO2_BM_CCS.value

    @property
    def VOM_CCS_BM(self):
        model, time = self.model, self.data.time
        return model.theta_BM_CCS_v.value * sum(model.Q_CO2_BM_CCS[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_BM_CCS(self):
        return self.CAPEX_CCS_BM + self.FOM_CCS_BM + self.VOM_CCS_BM

    @property
    def total_cost_BM_with_CCS(self):
        return self.total_cost_BM + self.total_cost_BM_CCS

    @property
    def FOM_WS(self):
        model, data = self.model, self.data
        return data.theta_WS_f * model.kappa_WS.value * model.n_y.value

    @property
    def VOM_WS(self):
        model, time = self.model, self.data.time
        return model.theta_WS_v.value * sum([model.P_WS[t].value for t in time]) * model.delta_t.value

    @property
    def fuel_costs_WS(self):
        model, time = self.model, self.data.time
        return model.theta_WS_fuel.value * sum(model.P_WS[t].value for t in time) * model.delta_t / model.eta_WS.value

    @property
    def CO2_cost_WS(self):
        model, time = self.model, self.data.time
        return model.theta_CO2.value * sum(model.Q_CO2_WS[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_WS(self):
        return self.FOM_WS + self.VOM_WS + self.fuel_costs_WS + self.CO2_cost_WS

    @property
    def CAPEX_CCS_WS(self):
        model = self.model
        return model.zeta_WS_CCS.value * model.K_CO2_WS_CCS.value

    @property
    def FOM_CCS_WS(self):
        model = self.model
        return model.theta_WS_CCS_f.value * model.n_y.value * model.K_CO2_WS_CCS.value

    @property
    def VOM_CCS_WS(self):
        model, time = self.model, self.data.time
        return model.theta_WS_CCS_v.value * sum(model.Q_CO2_WS_CCS[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_WS_CCS(self):
        return self.CAPEX_CCS_WS + self.FOM_CCS_WS + self.VOM_CCS_WS

    @property
    def total_cost_WS_with_CCS(self):
        return self.total_cost_WS + self.total_cost_WS_CCS

    @property
    def CAPEX_SMR(self):
        model = self.model
        return  model.zeta_SMR.value * model.K_H2_SMR.value

    @property
    def FOM_SMR(self):
        model = self.model
        return model.theta_SMR_f.value * model.K_H2_SMR.value * model.n_y.value

    @property
    def VOM_SMR(self):
        model, time = self.model, self.data.time
        return model.theta_SMR_v.value * sum(model.P_H2_SMR[t].value for t in time) * model.delta_t

    @property
    def estimated_fuel_costs_SMR(self):
        model, time = self.model, self.data.time
        return sum(model.P_NG_SMR[t].value * model.theta_NG_I[t] for t in time) * model.delta_t.value

#    @property
#    def CO2_cost_SMR(self):
#        model, time = self.model, self.data.time
#        return sum(model.theta_CO2.value * model.Q_CO2_SMR[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_SMR(self):
        return self.CAPEX_SMR + self.FOM_SMR + self.VOM_SMR

    @property
    def CAPEX_CCS_SMR(self):
        model = self.model
        return model.zeta_SMR_CCS.value * model.K_CO2_SMR_CCS.value

    @property
    def FOM_CCS_SMR(self):
        model = self.model
        return model.theta_SMR_CCS_f.value * model.n_y.value * model.K_CO2_SMR_CCS.value

    @property
    def VOM_CCS_SMR(self):
        model, time = self.model, self.data.time
        return model.theta_SMR_CCS_v.value * sum(model.Q_CO2_SMR_CCS[t].value for t in time) * model.delta_t.value

    @property
    def total_cost_SMR_CCS(self):
        return self.CAPEX_CCS_SMR + self.FOM_CCS_SMR + self.VOM_CCS_SMR

    @property
    def total_cost_SMR_with_CCS(self):
        return self.total_cost_SMR + self.total_cost_SMR_CCS

    @property
    def total_cost_CCS(self):
        return self.total_cost_OCGT_CCS + self.total_cost_CCGT_CCS + self.total_cost_CHP_CCS + self.total_cost_BM_CCS + self.total_cost_WS_CCS + self.total_cost_SMR_CCS

    @property
    def FOM_NK(self):
        model, theta_NK_f = self.model, self.data.theta_NK_f
        return theta_NK_f * model.kappa_NK.value * model.n_y.value
    
    @property
    def VOM_NK(self):
        model, time = self.model, self.data.time
        return model.theta_NK_v.value * sum([model.P_E_NK[t].value for t in time]) * model.delta_t.value

    @property
    def fuel_costs_NK(self):
        model, time = self.model, self.data.time
        return model.theta_NK_fuel.value * sum(model.P_E_NK[t].value for t in time) * model.delta_t / model.eta_NK.value

    @property
    def total_cost_NK(self):
        return self.FOM_NK + self.VOM_NK + self.fuel_costs_NK

    @property
    def total_cost_NG_imports(self):
        model, time = self.model, self.data.time
        return sum([model.theta_NG_I[t] * model.P_NG_I[t].value for t in time])

    @property
    def total_cost_E_imports(self):
        model, time = self.model, self.data.time
        return sum([model.theta_IE[t] * model.P_E_I[t].value for t in time])

    @property
    def total_cost_H2_imports(self):
        model, time = self.model, self.data.time
        return sum([model.theta_H2_IE[t] * model.P_H2_I[t].value for t in time])

    @property
    def CAPEX_ACC(self):
        model = self.model
        return model.zeta_ACC.value * model.K_CO2_ACC.value

    @property
    def FOM_ACC(self):
        model = self.model
        return model.theta_ACC_f.value * model.K_CO2_ACC.value * model.n_y.value

    @property
    def estimated_fuel_costs_ACC(self):
        model, time = self.model, self.data.time
        return sum(model.P_NG_ACC[t].value * model.theta_NG_I[t] for t in time)

    @property
    def total_cost_ACC(self):
        return self.CAPEX_ACC + self.FOM_ACC

    @property
    def VOM_PH(self):
        model, time = self.model, self.data.time
        return model.theta_PH_v.value * sum([model.P_PtPH[t].value for t in time]) * model.delta_t.value

    @property
    def FOM_PH(self):
        model, theta_PH_f = self.model, self.data.theta_PH_f
        return theta_PH_f * model.xi_PH.value * model.n_y.value

    @property
    def total_cost_PH(self):
        return self.FOM_PH + self.VOM_PH

    @property
    def OPEX_ENS_E(self):
        model, time = self.model, self.data.time
        return model.varsigma_ENS_E.value * sum([model.L_E_ENS[t].value for t in time]) * model.delta_t.value

    @property
    def OPEX_ENS_NG(self):
        model, time = self.model, self.data.time
        return model.varsigma_ENS_NG.value * sum([model.L_NG_ENS[t].value for t in time]) * model.delta_t.value

    @property
    def OPEX_ENS_H2(self):
        model, time = self.model, self.data.time
        return model.varsigma_ENS_H2.value * sum([model.L_H2_ENS[t].value for t in time]) * model.delta_t.value

    @property
    def CAPEX_CO2_storage(self):
        model = self.model
        return model.zeta_CO2S.value * model.S_CO2.value

    @property
    def FOM_CO2_storage(self):
        model = self.model
        return model.theta_CO2S_f.value * model.S_CO2.value * model.n_y.value

    @property
    def total_cost_CO2S(self):
        return self.CAPEX_CO2_storage + self.FOM_CO2_storage

    @property
    def CO2_exports(self):
        model, time = self.model, self.data.time
        return sum(model.Q_CO2_E[t].value for t in time) * model.delta_t.value

    @property
    def CO2_captured(self):
        model, time = self.model, self.data.time
        return sum(model.Q_CO2_OCGT_CCS[t].value + model.Q_CO2_CCGT_CCS[t].value + model.Q_CO2_BM_CCS[t].value + \
                   model.Q_CO2_WS_CCS[t].value + model.Q_CO2_CHP_CCS[t].value + model.Q_CO2_SMR_CCS[t].value + model.Q_CO2_ACC[t].value for t in time)

    @property
    def electricity_system_cost(self):
        model, time = self.model, self.data.time
        return (self.total_cost_W_on + self.total_cost_W_off + self.total_cost_S + self.total_cost_OCGT + self.total_cost_CCGT + self.total_cost_BM + \
                self.total_cost_WS + self.total_cost_CHP + self.total_cost_B + self.total_cost_PH + self.total_cost_FC + self.total_cost_NK + self.total_cost_E_imports) / sum(model.lambda_E[t] + model.lambda_E_HT[t] + model.L_E_TR[t].value + model.P_E_EL[t].value + model.P_E_SMR[t].value + model.P_OCGT_CCS[t].value + model.P_CCGT_CCS[t].value + model.P_CHP_CCS[t].value + model.P_BM_CCS[t].value + model.P_WS_CCS[t].value + model.P_SMR_CCS[t].value + model.P_E_ACC[t].value - model.L_E_ENS[t].value for t in time)

    @property
    def hydrogen_cost(self):
        model, time = self.model, self.data.time
        return (self.total_cost_SMR + self.total_cost_EL + self.total_cost_H2S + self.total_cost_H2_imports) / sum(model.P_H2_MT[t].value + model.P_E_FC[t].value + model.lambda_H2_ID[t] + model.lambda_H2_TR[t] - model.L_H2_ENS[t].value for t in time)

    @property
    def hydrogen_system_cost(self):
        return self.total_cost_SMR + self.total_cost_EL + self.total_cost_H2S + self.total_cost_H2_imports

    @property
    def natural_gas_cost(self):
        model, time = self.model, self.data.time
        return (self.total_cost_NG_imports + self.total_cost_MT + self.total_cost_NGS) / sum(model.P_NG_ACC[t].value + model.P_NG_SMR[t].value + model.P_NG_CCGT[t].value + model.P_NG_OCGT[t].value + model.P_NG_CHP[t].value + model.lambda_NG_HT[t] + model.lambda_NG_TR[t] + model.lambda_NG_ID[t] - model.lambda_NGtH2[t] - model.L_NG_ENS[t].value for t in time)

    @property
    def natural_gas_system_cost(self):
        return self.total_cost_NG_imports + self.total_cost_MT + self.total_cost_NGS

    @property
    def CO2_system_cost(self):
        return self.total_cost_OCGT_CCS + self.total_cost_CCGT_CCS + self.total_cost_BM_CCS + self.total_cost_WS_CCS + self.total_cost_CHP_CCS + self.total_cost_SMR_CCS + self.total_cost_ACC + self.total_cost_CO2S + self.total_cost_CO2_E

    @property
    def total_cost_ENS_E(self):
        model, time = self.model, self.data.time
        return model.varsigma_ENS_E.value * sum(model.L_E_ENS[t].value for t in time)

    @property
    def total_cost_ENS_H2(self):
        model, time = self.model, self.data.time
        return model.varsigma_ENS_H2.value * sum(model.L_H2_ENS[t].value for t in time)

    @property
    def total_cost_ENS_NG(self):
        model, time = self.model, self.data.time
        return model.varsigma_ENS_NG.value * sum(model.L_NG_ENS[t].value for t in time)

    @property
    def total_cost_CO2_E(self):
        model, time = self.model, self.data.time
        return model.theta_CO2_E.value * sum(model.Q_CO2_E[t].value for t in time)

    @property
    def system_cost(self):
        return(self.total_cost_WS_with_CCS + self.total_cost_BM_with_CCS + self.total_cost_CHP_with_CCS + self.total_cost_OCGT_with_CCS + self.total_cost_CCGT_with_CCS + self.total_cost_SMR_with_CCS + self.total_cost_S + self.total_cost_W_off + self.total_cost_W_on + \
               self.total_cost_CO2S + self.total_cost_H2S + self.total_cost_NGS + self.total_cost_B + self.total_cost_ACC + self.total_cost_NG_imports + self.total_cost_H2_imports + self.total_cost_E_imports + self.total_cost_EL + self.total_cost_MT + self.total_cost_FC + \
               self.total_cost_PH + self.total_cost_NK + self.total_cost_ENS_E + self.total_cost_ENS_H2 + self.total_cost_ENS_NG + self.total_cost_CO2_E)
