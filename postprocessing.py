import numpy as np
import pandas as pd


class PostProcess:

    def __init__(self, path):
        self.path = path
        self.rs = pd.read_csv(path)
        self.T = len(self.rs['global.demand_el_base'])

    @property
    def solar_PV_generation(self):
        return sum(self.rs['SOLAR_PV_PLANTS.electricity'])

    @property
    def onshore_wind_generation(self):
        return sum(self.rs['ONSHORE_WIND_PLANTS.electricity'])

    @property
    def offshore_wind_generation(self):
        return sum(self.rs['OFFSHORE_WIND_PLANTS.electricity'])

    @property
    def res_generation(self):
        return self.solar_PV_generation + self.onshore_wind_generation + self.offshore_wind_generation

    @property
    def solar_PV_curtailment(self):
        max_gen = sum(self.rs['SOLAR_PV_PLANTS.capacity_factor']*(self.rs['SOLAR_PV_PLANTS.pre_installed_capacity'][0]+self.rs['SOLAR_PV_PLANTS.capacity'][0]))
        return (max_gen - self.solar_PV_generation)


    @property
    def onshore_wind_curtailment(self):
        max_gen = sum(self.rs['ONSHORE_WIND_PLANTS.capacity_factor']*(self.rs['ONSHORE_WIND_PLANTS.pre_installed_capacity'][0]+self.rs['ONSHORE_WIND_PLANTS.capacity'][0]))
        return (max_gen - self.onshore_wind_generation)

    @property
    def offshore_wind_curtailment(self):
        max_gen = sum(self.rs['OFFSHORE_WIND_PLANTS.capacity_factor']*(self.rs['OFFSHORE_WIND_PLANTS.pre_installed_capacity'][0]+self.rs['OFFSHORE_WIND_PLANTS.capacity'][0]))
        return (max_gen - self.offshore_wind_generation)

    @property
    def res_curtailment(self):
        return self.solar_PV_curtailment + self.onshore_wind_curtailment + self.offshore_wind_curtailment

    @property
    def smr_cf(self):
        return sum(self.rs['SMR_w_PCCC.SMR.hydrogen'])/(self.rs['SMR_w_PCCC.SMR.capacity'][0]*self.T)

    @property
    def exported_co2(self):
        return sum(self.rs['CARBON_DIOXIDE_EXPORTS.carbon_dioxide'])

    @property
    def ocgt_cf(self):
        return sum(self.rs['OCGT_w_PCCC.OCGT.electricity'])/(self.rs['OCGT_w_PCCC.OCGT.capacity'][0]*self.T)

    @property
    def ccgt_cf(self):
        return sum(self.rs['CCGT_w_PCCC.CCGT.electricity'])/(self.rs['CCGT_w_PCCC.CCGT.capacity'][0]*self.T)

    @property
    def ccgt_net_cf(self):
        return sum(self.rs['CCGT_w_PCCC.electricity'])/(self.rs['CCGT_w_PCCC.CCGT.capacity'][0]*self.T)

    @property
    def solar_pv_cost(self):
        inv_cost = self.rs['SOLAR_PV_PLANTS.investment_cost'][0]
        fom_cost = self.rs['SOLAR_PV_PLANTS.fom_cost'][0]
        vom_cost = self.rs['SOLAR_PV_PLANTS.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def onshore_wind_cost(self):
        inv_cost = self.rs['ONSHORE_WIND_PLANTS.investment_cost'][0]
        fom_cost = self.rs['ONSHORE_WIND_PLANTS.fom_cost'][0]
        vom_cost = self.rs['ONSHORE_WIND_PLANTS.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def offshore_wind_cost(self):
        inv_cost = self.rs['OFFSHORE_WIND_PLANTS.investment_cost'][0]
        fom_cost = self.rs['OFFSHORE_WIND_PLANTS.fom_cost'][0]
        vom_cost = self.rs['OFFSHORE_WIND_PLANTS.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def wpp_cost(self):
        fom_cost = self.rs['WASTE_POWER_PLANTS_w_PCCC.WASTE_POWER_PLANTS.fom_cost'][0]
        fuel_cost = self.rs['WASTE_POWER_PLANTS_w_PCCC.WASTE_POWER_PLANTS.fuel_cost'][0]
        non_fuel_vom_cost = self.rs['WASTE_POWER_PLANTS_w_PCCC.WASTE_POWER_PLANTS.non_fuel_vom_cost'][0]
        return fom_cost + fuel_cost + non_fuel_vom_cost

    @property
    def pccc_wpp_cost(self):
        inv_cost = self.rs['WASTE_POWER_PLANTS_w_PCCC.PCCC.investment_cost'][0]
        fom_cost = self.rs['WASTE_POWER_PLANTS_w_PCCC.PCCC.fom_cost'][0]
        vom_cost = self.rs['WASTE_POWER_PLANTS_w_PCCC.PCCC.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def wpp_w_pccc_cost(self):
        emissions_cost = self.rs['WASTE_POWER_PLANTS_w_PCCC.emissions_cost'][0]
        return self.wpp_cost + self.pccc_wpp_cost + emissions_cost

    @property
    def bmpp_cost(self):
        fom_cost = self.rs['BIOMASS_POWER_PLANTS_w_PCCC.BIOMASS_POWER_PLANTS.fom_cost'][0]
        fuel_cost = self.rs['BIOMASS_POWER_PLANTS_w_PCCC.BIOMASS_POWER_PLANTS.fuel_cost'][0]
        non_fuel_vom_cost = self.rs['BIOMASS_POWER_PLANTS_w_PCCC.BIOMASS_POWER_PLANTS.non_fuel_vom_cost'][0]
        return fom_cost + fuel_cost + non_fuel_vom_cost

    @property
    def pccc_bmpp_cost(self):
        inv_cost = self.rs['BIOMASS_POWER_PLANTS_w_PCCC.PCCC.investment_cost'][0]
        fom_cost = self.rs['BIOMASS_POWER_PLANTS_w_PCCC.PCCC.fom_cost'][0]
        vom_cost = self.rs['BIOMASS_POWER_PLANTS_w_PCCC.PCCC.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def bmpp_w_pccc_cost(self):
        emissions_cost = self.rs['BIOMASS_POWER_PLANTS_w_PCCC.emissions_cost'][0]
        return self.bmpp_cost + self.pccc_bmpp_cost + emissions_cost

    @property
    def ncpp_cost(self):
        fom_cost = self.rs['NUCLEAR_POWER_PLANTS.fom_cost'][0]
        fuel_cost = self.rs['NUCLEAR_POWER_PLANTS.fuel_cost']
        non_fuel_vom_cost = self.rs['NUCLEAR_POWER_PLANTS.non_fuel_vom_cost'][0]
        return fom_cost + fuel_cost + non_fuel_vom_cost

    @property
    def ng_imports_cost(self):
        return self.rs['NATURAL_GAS_IMPORTS.imports_cost'][0]

    @property
    def h2_imports_cost(self):
        return self.rs['HYDROGEN_IMPORTS.imports_cost'][0]

    @property
    def el_imports_cost(self):
        return self.rs['ELECTRICITY_INTERCONNECTION.imports_cost'][0]

    @property
    def co2_exports_cost(self):
        return self.rs['CARBON_DIOXIDE_EXPORTS.exports_cost'][0]

    @property
    def electrolysis_cost(self):
        inv_cost = self.rs['ELECTROLYSIS_PLANTS.investment_cost'][0]
        fom_cost = self.rs['ELECTROLYSIS_PLANTS.fom_cost'][0]
        vom_cost = self.rs['ELECTROLYSIS_PLANTS.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def h2_fc_cost(self):
        inv_cost = self.rs['HYDROGEN_FUEL_CELLS.investment_cost'][0]
        fom_cost = self.rs['HYDROGEN_FUEL_CELLS.fom_cost'][0]
        vom_cost = self.rs['HYDROGEN_FUEL_CELLS.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def chp_cost(self):
        fom_cost = self.rs['CHP_PLANTS_w_PCCC.CHP_PLANTS.fom_cost'][0]
        non_fuel_vom_cost = self.rs['CHP_PLANTS_w_PCCC.CHP_PLANTS.non_fuel_vom_cost'][0]
        return inv_cost + non_fuel_vom_cost

    @property
    def pccc_chp_cost(self):
        inv_cost = self.rs['CHP_PLANTS_w_PCCC.PCCC.investment_cost'][0]
        fom_cost = self.rs['CHP_PLANTS_w_PCCC.PCCC.fom_cost'][0]
        vom_cost = self.rs['CHP_PLANTS_w_PCCC.PCCC.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def chp_w_pccc_cost(self):
        emissions_cost = self.rs['CHP_PLANTS_w_PCCC.emissions_cost'][0]
        return self.chp_cost + self.pccc_chp_cost + emissions_cost

    @property
    def ocgt_cost(self):
        inv_cost = self.rs['OCGT_w_PCCC.OCGT.investment_cost'][0]
        fom_cost = self.rs['OCGT_w_PCCC.OCGT.fom_cost'][0]
        non_fuel_vom_cost = self.rs['OCGT_w_PCCC.OCGT.non_fuel_vom_cost'][0]
        return inv_cost + fom_cost + non_fuel_vom_cost

    @property
    def pccc_ocgt_cost(self):
        inv_cost = self.rs['OCGT_w_PCCC.PCCC.investment_cost'][0]
        fom_cost = self.rs['OCGT_w_PCCC.PCCC.fom_cost'][0]
        vom_cost = self.rs['OCGT_w_PCCC.PCCC.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def ocgt_w_pccc_cost(self):
        emissions_cost = self.rs['OCGT_w_PCCC.emissions_cost'][0]
        return self.ocgt_cost + self.pccc_ocgt_cost + emissions_cost

    @property
    def ccgt_cost(self):
        inv_cost = self.rs['CCGT_w_PCCC.CCGT.investment_cost'][0]
        fom_cost = self.rs['CCGT_w_PCCC.CCGT.fom_cost'][0]
        non_fuel_vom_cost = self.rs['CCGT_w_PCCC.CCGT.non_fuel_vom_cost'][0]
        return inv_cost + fom_cost + non_fuel_vom_cost

    @property
    def pccc_ccgt_cost(self):
        inv_cost = self.rs['CCGT_w_PCCC.PCCC.investment_cost'][0]
        fom_cost = self.rs['CCGT_w_PCCC.PCCC.fom_cost'][0]
        vom_cost = self.rs['CCGT_w_PCCC.PCCC.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def ccgt_w_pccc_cost(self):
        emissions_cost = self.rs['CCGT_w_PCCC.emissions_cost'][0]
        return self.ccgt_cost + self.pccc_ccgt_cost + emissions_cost

    @property
    def methanation_cost(self):
        inv_cost = self.rs['METHANATION_PLANTS.investment_cost'][0]
        fom_cost = self.rs['METHANATION_PLANTS.fom_cost'][0]
        vom_cost = self.rs['METHANATION_PLANTS.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def smr_cost(self):
        inv_cost = self.rs['SMR_w_PCCC.SMR.investment_cost'][0]
        fom_cost = self.rs['SMR_w_PCCC.SMR.fom_cost'][0]
        vom_cost = self.rs['SMR_w_PCCC.SMR.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def pccc_smr_cost(self):
        inv_cost = self.rs['SMR_w_PCCC.PCCC.investment_cost'][0]
        fom_cost = self.rs['SMR_w_PCCC.PCCC.fom_cost'][0]
        vom_cost = self.rs['SMR_w_PCCC.PCCC.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def smr_w_pccc_cost(self):
        emissions_cost = self.rs['SMR_w_PCCC.emissions_cost'][0]
        return self.smr_cost + self.pccc_smr_cost + emissions_cost

    @property
    def dac_cost(self):
        inv_cost = self.rs['DIRECT_AIR_CAPTURE_PLANTS.investment_cost'][0]
        fom_cost = self.rs['DIRECT_AIR_CAPTURE_PLANTS.fom_cost'][0]
        vom_cost = self.rs['DIRECT_AIR_CAPTURE_PLANTS.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def bs_cost(self):
        inv_cost = self.rs['BATTERY_STORAGE.investment_cost'][0]
        fom_cost = self.rs['BATTERY_STORAGE.fom_cost'][0]
        vom_cost = self.rs['BATTERY_STORAGE.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def phs_cost(self):
        fom_cost = self.rs['PUMPED_HYDRO_STORAGE.fom_cost'][0]
        vom_cost = self.rs['PUMPED_HYDRO_STORAGE.vom_cost'][0]
        return fom_cost + vom_cost

    @property
    def ngs_cost(self):
        fom_cost = self.rs['NATURAL_GAS_STORAGE.fom_cost'][0]
        vom_cost = self.rs['NATURAL_GAS_STORAGE.vom_cost'][0]
        return fom_cost + vom_cost

    @property
    def h2s_cost(self):
        inv_cost = self.rs['HYDROGEN_STORAGE.investment_cost'][0]
        fom_cost = self.rs['HYDROGEN_STORAGE.fom_cost'][0]
        vom_cost = self.rs['HYDROGEN_STORAGE.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def co2s_cost(self):
        inv_cost = self.rs['CARBON_DIOXIDE_STORAGE.investment_cost'][0]
        fom_cost = self.rs['CARBON_DIOXIDE_STORAGE.fom_cost'][0]
        vom_cost = self.rs['CARBON_DIOXIDE_STORAGE.vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def ens_el_cost(self):
        return self.rs['ELECTRICAL_DEMAND_RESPONSE.load_shedding_cost'][0]

    @property
    def ens_ng_cost(self):
        return self.rs['NATURAL_GAS_DEMAND_RESPONSE.load_shedding_cost'][0]

    @property
    def ens_h2_cost(self):
        return self.rs['HYDROGEN_DEMAND_RESPONSE.load_shedding_cost'][0]

    @property
    def demand_el_base(self):
        return sum(self.rs['global.demand_el_base'])

    @property
    def demand_el_ht(self):
        return sum(self.rs['global.demand_el_ht'])

    @property
    def demand_el_tr(self):
        return sum(self.rs['global.demand_el_tr'])

    @property
    def demand_el(self):
        return self.demand_el_base + self.demand_el_ht + self.demand_el_tr

    @property
    def demand_h2_id(self):
        return sum(self.rs['links.industry_demand_hydrogen'])

    @property
    def demand_h2_tr(self):
        return sum(self.rs['links.transport_demand_hydrogen'])

    @property
    def demand_h2(self):
        return self.demand_h2_id + self.demand_h2_tr

    @property
    def demand_ng_ht(self):
        return sum(self.rs['heat_demand_natural_gas'])

    @property
    def demand_ng_id(self):
        return sum(self.rs['industry_demand_natural_gas'])

    @property
    def demand_ng_tr(self):
        return sum(self.rs['transport_demand_natural_gas'])

    @property
    def legacy_demand_ng_smr(self):
        return sum(self.rs['legacy_SMR_demand_natural_gas'])

    @property
    def demand_ng(self):
        return self.demand_ng_ht + self.demand_ng_id + self.demand_ng_tr - self.legacy_demand_ng_smr

    @property
    def ens_el(self):
        return sum(self.rs['ELECTRICAL_DEMAND_RESPONSE.load_shedding'])

    @property
    def ens_h2(self):
        return sum(self.rs['HYDROGEN_DEMAND_RESPONSE.load_shedding'])

    @property
    def ens_ng(self):
        return sum(self.rs['NATURAL_GAS_DEMAND_RESPONSE.load_shedding'])

    @property
    def el_cost(self):
        res_cost = self.solar_pv_cost + self.onshore_wind_cost + self.offshore_wind_cost
        gt_cost = self.ocgt_w_pccc_cost + self.ccgt_w_pccc_cost
        pol_disp_cost = self.wpp_w_pccc_cost + self.bmpp_w_pccc_cost + self.chp_w_pccc_cost
        clean_disp_cost = self.h2_fc_cost + self.ncpp_cost
        stor_cost = self.battery_cost + self.phs_cost
        imports_cost = self.el_imports_cost
        cost = res_cost + gt_cost + pol_disp_cost + clean_disp_cost + stor_cost + imports_cost
        return cost / (demand_el - ens_el)

    @property
    def h2_cost(self):
        prod_cost = self.SMR_w_PCCC_cost + self.electrolysis_cost
        stor_cost = self.h2s_cost
        imports_cost = self.h2_imports_cost
        cost = prod_cost + stor_cost + imports_cost
        return cost / (demand_h2 - ens_h2)

    @property
    def ng_cost(self):
        prod_cost = self.methanation_cost
        stor_cost = self.ngs_cost
        imports_cost = self.ng_imports_cost
        cost = prod_cost + stor_cost + imports_cost
        return cost / (demand_ng - ens_ng)
