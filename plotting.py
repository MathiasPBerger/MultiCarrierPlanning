from matplotlib import pyplot as plt
import pandas as pd
#import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter
from matplotlib.dates import MonthLocator, DateFormatter
from datetime import datetime, timedelta


## PLOT 1: daily electricity, natural gas and hydrogen demands

path_scen1 = 'results/results_scenario1.csv'
rs1 = pd.read_csv(path_scen1)

n_h = 43800
n_d = 1825

demand_el_base = rs1['global.demand_el_base']
demand_el_ht = rs1['global.demand_el_ht']
demand_el_tr = rs1['global.demand_el_ev'][0:n_d]

d_demand_el_base = [sum(demand_el_base[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_el_ht = [sum(demand_el_ht[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_el = [(d_demand_el_base[d] + d_demand_el_ht[d] + demand_el_tr[d]) for d in range(n_d)]

demand_ng_ht = rs1['links.heat_demand_natural_gas']
demand_ng_id = rs1['links.industry_demand_natural_gas']
demand_ng_tr = rs1['links.transport_demand_natural_gas']
demand_ng_leg_smr = rs1['links.legacy_SMR_demand_natural_gas']

d_demand_ng_ht = [sum(demand_ng_ht[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_ng_id = [sum(demand_ng_id[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_ng_tr = [sum(demand_ng_tr[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_leg_smr = [sum(demand_ng_leg_smr[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_ng = [(d_demand_ng_ht[d] + d_demand_ng_id[d] + d_demand_ng_tr[d] - d_demand_leg_smr[d]) for d in range(n_d)]

demand_h2_id = rs1['links.industry_demand_hydrogen']
demand_h2_tr = rs1['links.transport_demand_hydrogen']

d_demand_h2_id = [sum(demand_h2_id[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_h2_tr = [sum(demand_h2_tr[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_h2 = [(d_demand_h2_id[d] + d_demand_h2_tr[d]) for d in range(n_d)]

d_index = [d for d in range(365)]
st, nd = d_index[0], d_index[-1]
d_dates = list()
d_dates.append(datetime.strptime('2014-01-01', '%Y-%m-%d'))
for d in d_index[1:]:
    d_dates.append(d_dates[0] + timedelta(days = d))

fig1 = plt.figure()
ax = plt.gca()
labels = ['Electricity', 'Natural Gas', 'Hydrogen']
colors = ['gold', 'salmon', 'skyblue']
ax.fill_between(d_dates, d_demand_el[st:nd+1], color=colors[0], label=labels[0], alpha=.75)
ax.fill_between(d_dates, d_demand_ng[st:nd+1], color=colors[1], label=labels[1], alpha=.55)
ax.fill_between(d_dates, d_demand_h2[st:nd+1], color=colors[2], label=labels[2], alpha=.95)
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
ax.xaxis.set_major_formatter(NullFormatter())
ax.xaxis.set_minor_formatter(DateFormatter('%b'))
ax.legend(loc='upper center', ncol=3, fontsize=6)
ax.set_xlim(d_dates[0], d_dates[-1])
plt.ylabel('Energy Demand [GWh/day]')
plt.show()

## PLOT 2: solar PV capacity factors


## PLOT 3: onshore wind capacity factors


## PLOT 4: offshore wind capacity factors


## PLOT 5: plant capacities


## PLOT 6: storage operation

path_scen2 = 'results/results_scenario2.csv'
rs2 = pd.read_csv(path_scen2)

bs_soc = rs2['BATTERY_STORAGE.electricity_stored']
h2s_soc = rs2['HYDROGEN_STORAGE.hydrogen_stored']
ng_soc = rs2['NATURAL_GAS_STORAGE.natural_gas_stored']

h_index = [h for h in range(8760)]
st, nd = h_index[0], h_index[-1]
h_dates = list()
h_dates.append(datetime.strptime('2014-01-01', '%Y-%m-%d'))
for h in h_index[1:]:
    h_dates.append(h_dates[0] + timedelta(hours = h))

ttl_fnt = 10

fig1 = plt.figure()
colors = ['gold', 'skyblue', 'salmon']
ax1 = plt.subplot(311)
ax1.fill_between(h_dates, bs_soc[st:nd+1], color=colors[0], alpha=.55)
ax1.set_title('Battery Storage', fontsize=ttl_fnt)
ax1.axes.get_xaxis().set_visible(False)
ax2 = plt.subplot(312, sharex=ax1)
ax2.fill_between(h_dates, h2s_soc[st:nd+1], color=colors[1], alpha=.65)
ax2.set_title('Hydrogen Storage', fontsize=ttl_fnt)
ax2.axes.get_xaxis().set_visible(False)
ax2.set_ylabel('State of Charge [GWh]')
ax3 = plt.subplot(313, sharex=ax1)
ax3.fill_between(h_dates, ng_soc[st:nd+1], color=colors[2], alpha=.75)
ax3.set_title('Natural Gas Storage', fontsize=ttl_fnt)
ax3.xaxis.set_major_locator(MonthLocator())
ax3.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
ax3.xaxis.set_major_formatter(NullFormatter())
ax3.xaxis.set_minor_formatter(DateFormatter('%b'))
ax3.set_xlim(h_dates[0], h_dates[-1])
plt.show()
