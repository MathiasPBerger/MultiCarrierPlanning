import numpy as np
import pandas as pd


class PostProcess:

    def __init__(self, path):
        self.path = path
        self.rs = pd.read_csv(path)
        self.T = len(self.rs['global.demand_el_base'])

    @property
    def solar_PV_generation(self):
        node_id = 'SOLAR_PV_PLANTS.'
        return sum(self.rs[node_id+'electricity'])

    @property
    def onshore_wind_generation(self):
        node_id = 'ONSHORE_WIND_PLANTS.'
        return sum(self.rs[node_id+'electricity'])

    @property
    def offshore_wind_generation(self):
        node_id = 'OFFSHORE_WIND_PLANTS.'
        return sum(self.rs[node_id+'electricity'])

    @property
    def res_generation(self):
        return (self.solar_PV_generation + self.onshore_wind_generation
                + self.offshore_wind_generation)

    @property
    def solar_PV_curtailment(self):
        node_id = 'SOLAR_PV_PLANTS.'
        max_gen = (sum(self.rs[node_id+'capacity_factor']
                   * (self.rs[node_id+'pre_installed_capacity'][0]
                      + self.rs[node_id+'capacity'][0])))
        return (max_gen - self.solar_PV_generation)


    @property
    def onshore_wind_curtailment(self):
        node_id = 'ONSHORE_WIND_PLANTS.'
        max_gen = (sum(self.rs[node_id+'capacity_factor']
                   * (self.rs[node_id+'pre_installed_capacity'][0]
                      + self.rs[node_id+'capacity'][0])))
        return (max_gen - self.onshore_wind_generation)

    @property
    def offshore_wind_curtailment(self):
        node_id = 'OFFSHORE_WIND_PLANTS.'
        max_gen = (sum(self.rs[node_id+'capacity_factor']
                   * (self.rs[node_id+'pre_installed_capacity'][0]
                      + self.rs[node_id+'capacity'][0])))
        return (max_gen - self.offshore_wind_generation)

    @property
    def res_curtailment(self):
        return (self.solar_PV_curtailment + self.onshore_wind_curtailment
                + self.offshore_wind_curtailment)

    @property
    def smr_cf(self):
        node_id = 'SMR_w_PCCC.SMR.'
        return 100 * (sum(self.rs[node_id+'hydrogen'])
                / (self.rs[node_id+'capacity'][0] * self.T))

    @property
    def ng_smr_cons(self):
        node_id = 'SMR_w_PCCC.'
        return sum(self.rs[node_id+'natural_gas'])

    @property
    def ocgt_cf(self):
        node_id = 'OCGT_w_PCCC.OCGT.'
        if self.rs[node_id+'capacity'][0] > 1e-3:
            return (100 * sum(self.rs[node_id+'electricity'])
                    / (self.rs[node_id+'capacity'][0] * self.T))
        else:
            return 0.0

    @property
    def ccgt_cf(self):
        node_id = 'CCGT_w_PCCC.CCGT.'
        if self.rs[node_id+'capacity'][0] > 1e-3:
            return (100 * sum(self.rs[node_id+'electricity'])
                    / (self.rs[node_id+'capacity'][0] * self.T))
        else:
            return 0.0

    @property
    def ccgt_net_cf(self):
        node_id = 'CCGT_w_PCCC.'
        if self.rs[node_id+'CCGT.capacity'][0] > 1e-3:
            return (100 * sum(self.rs[node_id+'electricity'])
                    / (self.rs[node_id+'CCGT.capacity'][0] * self.T))
        else:
            return 0.0

    @property
    def solar_pv_cost(self):
        node_id = 'SOLAR_PV_PLANTS.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def onshore_wind_cost(self):
        node_id = 'ONSHORE_WIND_PLANTS.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def offshore_wind_cost(self):
        node_id = 'OFFSHORE_WIND_PLANTS.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def wpp_cost(self):
        node_id = 'WASTE_POWER_PLANTS_w_PCCC.WASTE_POWER_PLANTS.'
        fom_cost = self.rs[node_id+'fom_cost'][0]
        fuel_cost = self.rs[node_id+'fuel_cost'][0]
        non_fuel_vom_cost = self.rs[node_id+'non_fuel_vom_cost'][0]
        return fom_cost + fuel_cost + non_fuel_vom_cost

    @property
    def pccc_wpp_cost(self):
        node_id = 'WASTE_POWER_PLANTS_w_PCCC.PCCC.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def wpp_w_pccc_cost(self):
        node_id = 'WASTE_POWER_PLANTS_w_PCCC.'
        emissions_cost = self.rs[node_id+'emissions_cost'][0]
        return self.wpp_cost + self.pccc_wpp_cost + emissions_cost

    @property
    def bmpp_cost(self):
        node_id = 'BIOMASS_POWER_PLANTS_w_PCCC.BIOMASS_POWER_PLANTS.'
        fom_cost = self.rs[node_id+'fom_cost'][0]
        fuel_cost = self.rs[node_id+'fuel_cost'][0]
        non_fuel_vom_cost = self.rs[node_id+'non_fuel_vom_cost'][0]
        return fom_cost + fuel_cost + non_fuel_vom_cost

    @property
    def pccc_bmpp_cost(self):
        node_id = 'BIOMASS_POWER_PLANTS_w_PCCC.PCCC.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def bmpp_w_pccc_cost(self):
        node_id = 'BIOMASS_POWER_PLANTS_w_PCCC.'
        emissions_cost = self.rs[node_id+'emissions_cost'][0]
        return self.bmpp_cost + self.pccc_bmpp_cost + emissions_cost

    @property
    def ncpp_cost(self):
        node_id = 'NUCLEAR_POWER_PLANTS.'
        fom_cost = self.rs[node_id+'fom_cost'][0]
        fuel_cost = self.rs[node_id+'fuel_cost'][0]
        non_fuel_vom_cost = self.rs[node_id+'non_fuel_vom_cost'][0]
        if rs[node_id+'pre_installed_capacity'][0] > 1e-3:
            return fom_cost + fuel_cost + non_fuel_vom_cost
        else:
            return 0.0

    @property
    def ng_imports_cost(self):
        node_id = 'NATURAL_GAS_IMPORTS.'
        return self.rs[node_id+'imports_cost'][0]

    @property
    def h2_imports_cost(self):
        node_id = 'HYDROGEN_IMPORTS.'
        return self.rs[node_id+'imports_cost'][0]

    @property
    def el_imports_cost(self):
        node_id = 'ELECTRICITY_INTERCONNECTION.'
        return self.rs[node_id+'imports_cost'][0]

    @property
    def co2_exports_cost(self):
        node_id = 'CARBON_DIOXIDE_EXPORTS.'
        return self.rs[node_id+'exports_cost'][0]

    @property
    def electrolysis_cost(self):
        node_id = 'ELECTROLYSIS_PLANTS.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def h2_fc_cost(self):
        node_id = 'HYDROGEN_FUEL_CELLS.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def chp_cost(self):
        node_id = 'CHP_PLANTS_w_PCCC.CHP_PLANTS.'
        fom_cost = self.rs[node_id+'fom_cost'][0]
        non_fuel_vom_cost = self.rs[node_id+'non_fuel_vom_cost'][0]
        return fom_cost + non_fuel_vom_cost

    @property
    def pccc_chp_cost(self):
        node_id = 'CHP_PLANTS_w_PCCC.PCCC.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def chp_w_pccc_cost(self):
        node_id = 'CHP_PLANTS_w_PCCC.'
        emissions_cost = self.rs[node_id+'emissions_cost'][0]
        return self.chp_cost + self.pccc_chp_cost + emissions_cost

    @property
    def ocgt_cost(self):
        node_id = 'OCGT_w_PCCC.OCGT.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        non_fuel_vom_cost = self.rs[node_id+'non_fuel_vom_cost'][0]
        return inv_cost + fom_cost + non_fuel_vom_cost

    @property
    def pccc_ocgt_cost(self):
        node_id = 'OCGT_w_PCCC.PCCC.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def ocgt_w_pccc_cost(self):
        node_id = 'OCGT_w_PCCC.'
        emissions_cost = self.rs[node_id+'emissions_cost'][0]
        return self.ocgt_cost + self.pccc_ocgt_cost + emissions_cost

    @property
    def ccgt_cost(self):
        node_id = 'CCGT_w_PCCC.CCGT.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        non_fuel_vom_cost = self.rs[node_id+'non_fuel_vom_cost'][0]
        return inv_cost + fom_cost + non_fuel_vom_cost

    @property
    def pccc_ccgt_cost(self):
        node_id = 'CCGT_w_PCCC.PCCC.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def ccgt_w_pccc_cost(self):
        node_id = 'CCGT_w_PCCC.'
        emissions_cost = self.rs[node_id+'emissions_cost'][0]
        return self.ccgt_cost + self.pccc_ccgt_cost + emissions_cost

    @property
    def methanation_cost(self):
        node_id = 'METHANATION_PLANTS.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def smr_cost(self):
        node_id = 'SMR_w_PCCC.SMR.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def pccc_smr_cost(self):
        node_id = 'SMR_w_PCCC.PCCC.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def smr_w_pccc_cost(self):
        node_id = 'SMR_w_PCCC.'
        emissions_cost = self.rs[node_id+'emissions_cost'][0]
        return self.smr_cost + self.pccc_smr_cost + emissions_cost

    @property
    def dac_cost(self):
        node_id = 'DIRECT_AIR_CAPTURE_PLANTS.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def bs_cost(self):
        node_id = 'BATTERY_STORAGE.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def phs_cost(self):
        node_id = 'PUMPED_HYDRO_STORAGE.'
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return fom_cost + vom_cost

    @property
    def ngs_cost(self):
        node_id = 'NATURAL_GAS_STORAGE.'
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return fom_cost + vom_cost

    @property
    def h2s_cost(self):
        node_id = 'HYDROGEN_STORAGE.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def co2s_cost(self):
        node_id = 'CARBON_DIOXIDE_STORAGE.'
        inv_cost = self.rs[node_id+'investment_cost'][0]
        fom_cost = self.rs[node_id+'fom_cost'][0]
        vom_cost = self.rs[node_id+'vom_cost'][0]
        return inv_cost + fom_cost + vom_cost

    @property
    def ens_el_cost(self):
        node_id = 'ELECTRICAL_DEMAND_RESPONSE.'
        return self.rs[node_id+'load_shedding_cost'][0]

    @property
    def ens_ng_cost(self):
        node_id = 'NATURAL_GAS_DEMAND_RESPONSE.'
        return self.rs[node_id+'load_shedding_cost'][0]

    @property
    def ens_h2_cost(self):
        node_id = 'HYDROGEN_DEMAND_RESPONSE.'
        return self.rs[node_id+'load_shedding_cost'][0]

    @property
    def demand_el_base(self):
        return sum(self.rs['global.demand_el_base'])

    @property
    def demand_el_ht(self):
        return sum(self.rs['global.demand_el_ht'])

    @property
    def demand_el_tr(self):
        return sum(self.rs['global.demand_el_ev'][0:int(self.T/24)])

    @property
    def demand_el_exo(self):
        return self.demand_el_base + self.demand_el_ht + self.demand_el_tr

    @property
    def demand_h2_id(self):
        return sum(self.rs['links.industry_demand_hydrogen'])

    @property
    def demand_h2_tr(self):
        return sum(self.rs['links.transport_demand_hydrogen'])

    @property
    def demand_h2_exo(self):
        return self.demand_h2_id + self.demand_h2_tr

    @property
    def demand_ng_ht(self):
        return sum(self.rs['links.heat_demand_natural_gas'])

    @property
    def demand_ng_id(self):
        return sum(self.rs['links.industry_demand_natural_gas'])

    @property
    def demand_ng_tr(self):
        return sum(self.rs['links.transport_demand_natural_gas'])

    @property
    def legacy_demand_ng_smr(self):
        return sum(self.rs['links.legacy_SMR_demand_natural_gas'])

    @property
    def demand_ng_exo(self):
        return (self.demand_ng_ht + self.demand_ng_id + self.demand_ng_tr
                - self.legacy_demand_ng_smr)

    @property
    def ens_el(self):
        node_id = 'ELECTRICAL_DEMAND_RESPONSE.'
        return sum(self.rs[node_id+'load_shedding'])

    @property
    def ens_h2(self):
        node_id = 'HYDROGEN_DEMAND_RESPONSE.'
        return sum(self.rs[node_id+'load_shedding'])

    @property
    def ens_ng(self):
        node_id = 'NATURAL_GAS_DEMAND_RESPONSE.'
        return sum(self.rs[node_id+'load_shedding'])

    @property
    def av_ng_import_cost(self):
        node_id = 'NATURAL_GAS_IMPORTS.'
        return 1000 * self.rs[node_id+'import_cost'].mean()

    @property
    def el_imported(self):
        node_id = 'ELECTRICITY_INTERCONNECTION.'
        return sum(self.rs[node_id+'electricity_imported'])

    @property
    def el_exported(self):
        node_id = 'ELECTRICITY_INTERCONNECTION.'
        return sum(self.rs[node_id+'electricity_exported'])

    @property
    def h2_imported(self):
        node_id = 'HYDROGEN_IMPORTS.'
        return sum(self.rs[node_id+'hydrogen'])

    @property
    def ng_imported(self):
        node_id = 'NATURAL_GAS_IMPORTS.'
        return sum(self.rs[node_id+'natural_gas'])

    @property
    def co2_exported(self):
        node_id = 'CARBON_DIOXIDE_EXPORTS.'
        return sum(self.rs[node_id+'carbon_dioxide'])

    @property
    def el_cost(self):
        res_cost = (self.solar_pv_cost + self.onshore_wind_cost
                    + self.offshore_wind_cost)
        gt_cost = self.ocgt_w_pccc_cost + self.ccgt_w_pccc_cost
        pol_disp_cost = (self.wpp_w_pccc_cost + self.bmpp_w_pccc_cost
                         + self.chp_w_pccc_cost)
        clean_disp_cost = self.h2_fc_cost + self.ncpp_cost
        stor_cost = self.bs_cost + self.phs_cost
        imports_cost = self.el_imports_cost
        cost = (res_cost + gt_cost + pol_disp_cost + clean_disp_cost
                + stor_cost + imports_cost)
        demand_dac = sum(self.rs['DIRECT_AIR_CAPTURE_PLANTS.electricity'])
        demand_electrolysis = sum(self.rs['ELECTROLYSIS_PLANTS.electricity'])
        demand_smr = sum(self.rs['SMR_w_PCCC.electricity'])
        demand = (self.demand_el_exo + demand_dac + demand_electrolysis
                  + demand_smr)
        return 1000 * (cost / (demand - self.ens_el))

    @property
    def h2_cost(self):
        tech_cost = self.smr_w_pccc_cost + self.electrolysis_cost
        stor_cost = self.h2s_cost
        imports_cost = self.h2_imports_cost
        cost = tech_cost + stor_cost + imports_cost
        demand_mt = sum(self.rs['METHANATION_PLANTS.hydrogen'])
        demand_h2_fc = sum(self.rs['HYDROGEN_FUEL_CELLS.hydrogen'])
        demand = self.demand_h2_exo + demand_mt + demand_h2_fc
        return 1000 * (cost / (self.demand_h2_exo - self.ens_h2))

    @property
    def ng_cost(self):
        prod_cost = self.methanation_cost
        stor_cost = self.ngs_cost
        imports_cost = self.ng_imports_cost
        cost = prod_cost + stor_cost + imports_cost
        demand_smr = sum(self.rs['SMR_w_PCCC.natural_gas'])
        demand_dac = sum(self.rs['DIRECT_AIR_CAPTURE_PLANTS.natural_gas'])
        demand_chp = sum(self.rs['CHP_PLANTS_w_PCCC.natural_gas'])
        demand_ocgt = sum(self.rs['OCGT_w_PCCC.natural_gas'])
        demand_ccgt = sum(self.rs['CCGT_w_PCCC.natural_gas'])
        demand = (self.demand_ng_exo + demand_smr + demand_dac + demand_chp
                  + demand_ocgt + demand_ccgt)
        return 1000 * (cost / (demand - self.ens_ng))

    @property
    def co2_cost(self):
        pccc_cost = (self.pccc_wpp_cost + self.pccc_bmpp_cost
                     + self.pccc_chp_cost + self.pccc_ocgt_cost
                     + self.pccc_ccgt_cost)
        dac_cost = self.dac_cost
        stor_cost = self.co2s_cost
        exports_cost = self.co2_exports_cost
        cost = pccc_cost + dac_cost + stor_cost + exports_cost
        if self.co2_exported > 1e-3:
            return 1000 * (cost / self.co2_exported)
        else:
            return 0.0

    @property
    def el_price_mean_dual(self):
        return 1000 * self.rs['links.power_balance.dual'].mean()

    @property
    def el_price_median_dual(self):
        return 1000 * self.rs['links.power_balance.dual'].median()

    @property
    def h2_price_mean_dual(self):
        return 1000 * self.rs['links.hydrogen_balance.dual'].mean()

    @property
    def h2_price_median_dual(self):
        return 1000 * self.rs['links.hydrogen_balance.dual'].median()

    @property
    def ng_price_mean_dual(self):
        return 1000 * self.rs['links.natural_gas_balance.dual'].mean()

    @property
    def ng_price_median_dual(self):
        return 1000 * self.rs['links.natural_gas_balance.dual'].median()

    @property
    def h2_electrolysis_cost(self):
        tech_cost = self.electrolysis_cost
        feed_cost = (sum(self.rs['ELECTROLYSIS_PLANTS.electricity']
                     * self.rs['links.power_balance.dual']))
        produced_volume = sum(self.rs['ELECTROLYSIS_PLANTS.hydrogen'])
        return (tech_cost + feed_cost) / produced_volume

    @property
    def h2_smr_cost(self):
        tech_cost = self.smr_w_pccc_cost
        feed_cost = (sum(self.rs['SMR_w_PCCC.natural_gas']
                     * self.rs['links.natural_gas_balance.dual']))
        produced_volume = sum(self.rs['SMR_w_PCCC.hydrogen'])
        if self.rs['SMR_w_PCCC.SMR.capacity'][0] > 1e-3:
            return (tech_cost + feed_cost) / produced_volume
        else:
            return 0.0
