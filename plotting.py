from matplotlib import pyplot as plt
import pandas as pd
#import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter
from matplotlib.dates import MonthLocator, DateFormatter
from datetime import datetime, timedelta

path = 'results/results_scenario1.csv'
rs = pd.read_csv(path)

## PLOT 1: daily electricity, natural gas and hydrogen demands

n_h = 43800
n_d = 1825

demand_el_base = rs['global.demand_el_base']
demand_el_ht = rs['global.demand_el_ht']
demand_el_tr = rs['global.demand_el_ev'][0:n_d]

d_demand_el_base = [sum(demand_el_base[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_el_ht = [sum(demand_el_ht[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_el = [(d_demand_el_base[d] + d_demand_el_ht[d] + demand_el_tr[d]) for d in range(n_d)]

demand_ng_ht = rs['links.heat_demand_natural_gas']
demand_ng_id = rs['links.industry_demand_natural_gas']
demand_ng_tr = rs['links.transport_demand_natural_gas']
demand_ng_leg_smr = rs['links.legacy_SMR_demand_natural_gas']

d_demand_ng_ht = [sum(demand_ng_ht[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_ng_id = [sum(demand_ng_id[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_ng_tr = [sum(demand_ng_tr[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_leg_smr = [sum(demand_ng_leg_smr[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_ng = [(d_demand_ng_ht[d] + d_demand_ng_id[d] + d_demand_ng_tr[d] - d_demand_leg_smr[d]) for d in range(n_d)]

demand_h2_id = rs['links.industry_demand_hydrogen']
demand_h2_tr = rs['links.transport_demand_hydrogen']

d_demand_h2_id = [sum(demand_h2_id[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_h2_tr = [sum(demand_h2_tr[d+t] for t in range(24)) for d in range(n_h) if not d % 24]
d_demand_h2 = [(d_demand_h2_id[d] + d_demand_h2_tr[d]) for d in range(n_d)]

index = [d for d in range(365)]
st, nd = index[0], index[-1]
dates = list()
dates.append(datetime.strptime('2014-01-01', '%Y-%m-%d'))
for d in index[1:]:
    dates.append(dates[0] + timedelta(days = d))

fig1 = plt.figure()
ax = plt.gca()
labels = ['Electricity', 'Natural Gas', 'Hydrogen']
colors = ['royalblue', 'firebrick', 'seagreen']
ax.fill_between(dates, d_demand_el[st:nd+1], color=colors[0], label=labels[0], alpha=.55)
ax.fill_between(dates, d_demand_ng[st:nd+1], color=colors[1], label=labels[1], alpha=.65)
ax.fill_between(dates, d_demand_h2[st:nd+1], color=colors[2], label=labels[2], alpha=.75)
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
ax.xaxis.set_major_formatter(NullFormatter())
ax.xaxis.set_minor_formatter(DateFormatter('%b'))
ax.legend(loc='upper center', ncol=3, fontsize=6)
ax.set_xlim(dates[0], dates[-1])
plt.ylabel('Energy Demand [GWh/day]')
plt.show()

## PLOT 2: solar PV capacity factors



## PLOT 3: onshore wind capacity factors


## PLOT 4: offshore wind capacity factors


## PLOT 5: plant capacities


## PLOT 6: storage operation
