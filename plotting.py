from matplotlib import pyplot as plt
import pandas as pd
#import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter
from matplotlib.dates import MonthLocator, DateFormatter, DayLocator
from datetime import datetime, timedelta

## TIME INDICES

n_h = 43800
n_d = 1825

## PLOT 1: daily electricity, natural gas and hydrogen demands

# d_index = [d for d in range(365)]
# d_st, d_nd = d_index[0], d_index[-1]
# d_dates = list()
# d_dates.append(datetime.strptime('2014-01-01', '%Y-%m-%d'))
# for d in d_index[1:]:
#     d_dates.append(d_dates[0] + timedelta(days = d))
#
# path_scen1 = 'results/results_scenario1.csv'
# rs1 = pd.read_csv(path_scen1)
#
# demand_el_base = rs1['global.demand_el_base']
# demand_el_ht = rs1['global.demand_el_ht']
# demand_el_tr = rs1['global.demand_el_ev'][0:n_d]
#
# d_demand_el_base = [sum(demand_el_base[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
# d_demand_el_ht = [sum(demand_el_ht[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
# d_demand_el = [(d_demand_el_base[d] + d_demand_el_ht[d] + demand_el_tr[d]) for d in range(n_d)]
#
# demand_ng_ht = rs1['links.heat_demand_natural_gas']
# demand_ng_id = rs1['links.industry_demand_natural_gas']
# demand_ng_tr = rs1['links.transport_demand_natural_gas']
# demand_ng_leg_smr = rs1['links.legacy_SMR_demand_natural_gas']
#
# d_demand_ng_ht = [sum(demand_ng_ht[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
# d_demand_ng_id = [sum(demand_ng_id[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
# d_demand_ng_tr = [sum(demand_ng_tr[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
# d_demand_leg_smr = [sum(demand_ng_leg_smr[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
# d_demand_ng = [(d_demand_ng_ht[d] + d_demand_ng_id[d] + d_demand_ng_tr[d] - d_demand_leg_smr[d]) for d in range(n_d)]
#
# demand_h2_id = rs1['links.industry_demand_hydrogen']
# demand_h2_tr = rs1['links.transport_demand_hydrogen']
#
# d_demand_h2_id = [sum(demand_h2_id[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
# d_demand_h2_tr = [sum(demand_h2_tr[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
# d_demand_h2 = [(d_demand_h2_id[d] + d_demand_h2_tr[d]) for d in range(n_d)]
#
# fig1 = plt.figure()
# ax = plt.gca()
# labels = ['Electricity', 'Natural Gas', 'Hydrogen']
# colors = ['gold', 'salmon', 'lightblue']
# ax.fill_between(d_dates, d_demand_el[d_st:d_nd+1], color=colors[0], label=labels[0], alpha=.75)
# ax.fill_between(d_dates, d_demand_ng[d_st:d_nd+1], color=colors[1], label=labels[1], alpha=.55)
# ax.fill_between(d_dates, d_demand_h2[d_st:d_nd+1], color=colors[2], label=labels[2], alpha=.95)
# ax.xaxis.set_major_locator(MonthLocator())
# ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
# ax.xaxis.set_major_formatter(NullFormatter())
# ax.xaxis.set_minor_formatter(DateFormatter('%b'))
# ax.legend(loc='upper center', ncol=3, fontsize=10)
# ax.set_xlim(d_dates[0], d_dates[-1])
# plt.ylabel('Energy Demand [GWh/day]')
# plt.show()

## PLOT 2: solar PV capacity factors

path = 'results/results_scenario1.csv'
rs = pd.read_csv(path)

solar_pv_cf = rs['SOLAR_PV_PLANTS.capacity_factor']
onshore_wind_cf = rs['ONSHORE_WIND_PLANTS.capacity_factor']
offshore_wind_cf = rs['OFFSHORE_WIND_PLANTS.capacity_factor']

# h_index = [h for h in range(720)]
# h_st, h_nd = h_index[0], h_index[-1]
# h_dates = list()
# h_dates.append(datetime.strptime('2014-08-01', '%Y-%m-%d'))
# for h in h_index[1:]:
#     h_dates.append(h_dates[0] + timedelta(hours = h))
#
# fig2 = plt.figure()
# labels = ['Solar PV', 'Onshore Wind', 'Offshore Wind']
# colors = ['gold', 'seagreen', 'royalblue']
# ax = plt.gca()
# ax.plot(h_dates, solar_pv_cf[h_st:h_nd+1], color=colors[0], label=labels[0])
# ax.plot(h_dates, onshore_wind_cf[h_st:h_nd+1], color=colors[1], label=labels[1])
# ax.plot(h_dates, offshore_wind_cf[h_st:h_nd+1], color=colors[2], label=labels[2])
# ax.xaxis.set_major_locator(MonthLocator())
# ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
# ax.xaxis.set_major_formatter(NullFormatter())
# ax.xaxis.set_minor_formatter(DateFormatter('%b'))
# ax.legend(loc='best', fontsize=10)
# ax.set_xlim(h_dates[0], h_dates[-1])
# plt.title('Capacity Factors of Renewable Resources')
# plt.ylabel('Capacity Factor [-]')
# plt.show()

# h_index = [h for h in range(43800)]
# h_dates = list()
# h_dates.append(datetime.strptime('2014-01-01', '%Y-%m-%d'))
# for h in h_index[1:]:
#     h_dates.append(h_dates[0] + timedelta(hours = h))
#
# labels = ['Solar PV', 'Onshore Wind', 'Offshore Wind']
# colors = ['gold', 'seagreen', 'royalblue']
#
# fig2,(ax1,ax2,ax3,ax4) = plt.subplots(1,4,sharey=True, facecolor='w')
# ax1.plot(h_dates, solar_pv_cf, color=colors[0], label=labels[0])
# ax1.plot(h_dates, onshore_wind_cf, color=colors[1], label=labels[1])
# ax1.plot(h_dates, offshore_wind_cf, color=colors[2], label=labels[2])
# ax1.xaxis.set_minor_locator(DayLocator(interval=2))
# ax1.xaxis.set_major_formatter(DateFormatter('%d/%m'))
# ax1.tick_params(axis='x', rotation=90)
# ax1.set_xlim(h_dates[216], h_dates[384])
# ax1.set_ylabel('Capacity Factors [-]', fontsize=15)
# ax1.spines['right'].set_visible(False)
# ax1.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
# ax2.plot(h_dates, solar_pv_cf, color=colors[0], label=labels[0])
# ax2.plot(h_dates, onshore_wind_cf, color=colors[1], label=labels[1])
# ax2.plot(h_dates, offshore_wind_cf, color=colors[2], label=labels[2])
# ax2.xaxis.set_minor_locator(DayLocator(interval=2))
# ax2.xaxis.set_major_formatter(DateFormatter('%d/%m'))
# ax2.tick_params(axis='x', rotation=90)
# ax2.set_xlim(h_dates[2376], h_dates[2544])
# ax2.spines['left'].set_visible(False)
# ax2.spines['right'].set_visible(False)
# ax2.tick_params(axis='y', color='w')
# ax3.plot(h_dates, solar_pv_cf, color=colors[0], label=labels[0])
# ax3.plot(h_dates, onshore_wind_cf, color=colors[1], label=labels[1])
# ax3.plot(h_dates, offshore_wind_cf, color=colors[2], label=labels[2])
# ax3.xaxis.set_minor_locator(DayLocator(interval=2))
# ax3.xaxis.set_major_formatter(DateFormatter('%d/%m'))
# ax3.tick_params(axis='x', rotation=90)
# ax3.tick_params(axis='y', color='w')
# ax3.set_xlim(h_dates[4560], h_dates[4728])
# ax3.spines['left'].set_visible(False)
# ax3.spines['right'].set_visible(False)
# ax3.legend(loc='upper right', fontsize=10)
# ax4.plot(h_dates, solar_pv_cf, color=colors[0], label=labels[0])
# ax4.plot(h_dates, onshore_wind_cf, color=colors[1], label=labels[1])
# ax4.plot(h_dates, offshore_wind_cf, color=colors[2], label=labels[2])
# ax4.xaxis.set_minor_locator(DayLocator(interval=2))
# ax4.xaxis.set_major_formatter(DateFormatter('%d/%m'))
# ax4.tick_params(axis='x', rotation=90)
# ax4.tick_params(axis='y', color='w')
# ax4.set_xlim(h_dates[6768], h_dates[6936])
# ax4.spines['left'].set_visible(False)
# plt.show()

## PLOT 3: hydrogen imports capacity

path = 'results/results_scenario1.csv'
rs = pd.read_csv(path)

h2_import_capacity = rs['HYDROGEN_IMPORTS.capacity']

h_index = [h for h in range(43800)]
h_dates = list()
h_dates.append(datetime.strptime('2014-01-01', '%Y-%m-%d'))
for h in h_index[1:]:
    h_dates.append(h_dates[0] + timedelta(hours = h))

colors = ['grey']

fig2,(ax1,ax2,ax3,ax4) = plt.subplots(1,4,sharey=True, facecolor='w')
ax1.plot(h_dates, h2_import_capacity, color=colors[0])
ax1.xaxis.set_minor_locator(DayLocator(interval=2))
ax1.xaxis.set_major_formatter(DateFormatter('%d/%m'))
ax1.tick_params(axis='x', rotation=90)
ax1.set_xlim(h_dates[216], h_dates[384])
ax1.set_ylabel('Import Capacity [GWh/h]', fontsize=15)
ax1.spines['right'].set_visible(False)
ax1.set_yticks([0.0, 1.75, 3.5, 5.25, 7.0])
ax2.plot(h_dates, h2_import_capacity, color=colors[0])
ax2.xaxis.set_minor_locator(DayLocator(interval=2))
ax2.xaxis.set_major_formatter(DateFormatter('%d/%m'))
ax2.tick_params(axis='x', rotation=90)
ax2.set_xlim(h_dates[2376], h_dates[2544])
ax2.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.tick_params(axis='y', color='w')
ax3.plot(h_dates, h2_import_capacity, color=colors[0])
ax3.xaxis.set_minor_locator(DayLocator(interval=2))
ax3.xaxis.set_major_formatter(DateFormatter('%d/%m'))
ax3.tick_params(axis='x', rotation=90)
ax3.tick_params(axis='y', color='w')
ax3.set_xlim(h_dates[4560], h_dates[4728])
#ax3.set_title('Hydrogen Tanker Schedule')
ax3.spines['left'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax4.plot(h_dates, h2_import_capacity, color=colors[0])
ax4.xaxis.set_minor_locator(DayLocator(interval=2))
ax4.xaxis.set_major_formatter(DateFormatter('%d/%m'))
ax4.tick_params(axis='x', rotation=90)
ax4.tick_params(axis='y', color='w')
ax4.set_xlim(h_dates[6768], h_dates[6936])
ax4.spines['left'].set_visible(False)
plt.show()

## PLOT 4: plant capacities


## PLOT 5: storage operation

# h_index = [h for h in range(8760)]
# h_st, h_nd = h_index[0], h_index[-1]
# h_dates = list()
# h_dates.append(datetime.strptime('2014-01-01', '%Y-%m-%d'))
# for h in h_index[1:]:
#     h_dates.append(h_dates[0] + timedelta(hours = h))
#
# path_scen2 = 'results/results_scenario2.csv'
# rs2 = pd.read_csv(path_scen2)
#
# bs_soc = rs2['BATTERY_STORAGE.electricity_stored']
# h2s_soc = rs2['HYDROGEN_STORAGE.hydrogen_stored']
# ng_soc = rs2['NATURAL_GAS_STORAGE.natural_gas_stored']
#
# ttl_fnt = 10
#
# fig7 = plt.figure()
# colors = ['gold', 'lightblue', 'salmon']
# ax1 = plt.subplot(311)
# ax1.fill_between(h_dates, bs_soc[h_st:h_nd+1], color=colors[0], alpha=.55)
# ax1.set_title('Battery Storage', fontsize=ttl_fnt)
# ax1.set_yticks([0, 20, 40])
# ax1.axes.get_xaxis().set_visible(False)
# ax2 = plt.subplot(312, sharex=ax1)
# ax2.fill_between(h_dates, h2s_soc[h_st:h_nd+1], color=colors[1], alpha=.65)
# ax2.set_title('Hydrogen Storage', fontsize=ttl_fnt)
# ax2.axes.get_xaxis().set_visible(False)
# ax2.set_ylabel('State of Charge [GWh]')
# ax2.set_yticks([0, 100, 200, 300])
# ax3 = plt.subplot(313, sharex=ax1)
# ax3.fill_between(h_dates, ng_soc[h_st:h_nd+1], color=colors[2], alpha=.75)
# ax3.set_title('Natural Gas Storage', fontsize=ttl_fnt)
# ax3.set_yticks([0, 2500, 5000, 7500])
# ax3.xaxis.set_major_locator(MonthLocator())
# ax3.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
# ax3.xaxis.set_major_formatter(NullFormatter())
# ax3.xaxis.set_minor_formatter(DateFormatter('%b'))
# ax3.set_xlim(h_dates[0], h_dates[-1])
# plt.show()
