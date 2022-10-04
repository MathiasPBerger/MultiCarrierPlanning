from matplotlib import pyplot as plt
import pandas as pd
from matplotlib.ticker import NullFormatter
from matplotlib.dates import MonthLocator, DateFormatter, DayLocator
from datetime import datetime, timedelta
import numpy as np

## Plot selection

show_demand = False
show_cf = False
show_h2_imp = False
show_stor_op = False
show_ng_bal = False
show_gen_capas = True
show_stor_capas = True
show_gas_prod_capas = True
show_ev_charge_dist = True

## PLOT 1: daily electricity, natural gas and hydrogen demands

if show_demand:
    n_h = 43800
    n_d = int(43800/24)

    d_index = [d for d in range(365)]
    d_st, d_nd = d_index[0], d_index[-1]
    d_dates = list()
    d_dates.append(datetime.strptime('2014-01-01', '%Y-%m-%d'))
    for d in d_index[1:]:
        d_dates.append(d_dates[0] + timedelta(days = d))

    path_scen1 = 'results/results_scenario1_cplex.csv'
    rs1 = pd.read_csv(path_scen1)

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

    fig1 = plt.figure()
    ax = plt.gca()
    labels = ['Electricity', 'Natural Gas', 'Hydrogen']
    colors = ['gold', 'salmon', 'lightblue']
    ax.fill_between(d_dates, d_demand_el[d_st:d_nd+1], color=colors[0], label=labels[0], alpha=.75)
    ax.fill_between(d_dates, d_demand_ng[d_st:d_nd+1], color=colors[1], label=labels[1], alpha=.55)
    ax.fill_between(d_dates, d_demand_h2[d_st:d_nd+1], color=colors[2], label=labels[2], alpha=.95)
    ax.xaxis.set_major_locator(MonthLocator())
    ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.xaxis.set_minor_formatter(DateFormatter('%b'))
    ax.legend(loc='upper center', ncol=3, fontsize=10)
    ax.set_xlim(d_dates[0], d_dates[-1])
    plt.ylabel('Energy Demand [GWh/day]')
    plt.show()

## PLOT 2: capacity factors

if show_cf:
    path = 'results/results_scenario1_cplex.csv'
    rs = pd.read_csv(path)

    solar_pv_cf = rs['SOLAR_PV_PLANTS.capacity_factor']
    onshore_wind_cf = rs['ONSHORE_WIND_PLANTS.capacity_factor']
    offshore_wind_cf = rs['OFFSHORE_WIND_PLANTS.capacity_factor']

    h_index = [h for h in range(43800)]
    h_dates = list()
    h_dates.append(datetime.strptime('2014-01-01', '%Y-%m-%d'))
    for h in h_index[1:]:
        h_dates.append(h_dates[0] + timedelta(hours = h))

    labels = ['Solar PV', 'Onshore Wind', 'Offshore Wind']
    colors = ['gold', 'seagreen', 'royalblue']

    fig2,(ax1,ax2,ax3,ax4) = plt.subplots(1,4,sharey=True, facecolor='w')
    ax1.plot(h_dates, solar_pv_cf, color=colors[0], label=labels[0])
    ax1.plot(h_dates, onshore_wind_cf, color=colors[1], label=labels[1])
    ax1.plot(h_dates, offshore_wind_cf, color=colors[2], label=labels[2])
    ax1.xaxis.set_minor_locator(DayLocator(interval=2))
    ax1.xaxis.set_major_formatter(DateFormatter('%d/%m'))
    ax1.tick_params(axis='x', rotation=90)
    ax1.set_xlim(h_dates[216], h_dates[384])
    ax1.set_ylabel('Capacity Factors [-]', fontsize=15)
    ax1.spines['right'].set_visible(False)
    ax1.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.plot(h_dates, solar_pv_cf, color=colors[0], label=labels[0])
    ax2.plot(h_dates, onshore_wind_cf, color=colors[1], label=labels[1])
    ax2.plot(h_dates, offshore_wind_cf, color=colors[2], label=labels[2])
    ax2.xaxis.set_minor_locator(DayLocator(interval=2))
    ax2.xaxis.set_major_formatter(DateFormatter('%d/%m'))
    ax2.tick_params(axis='x', rotation=90)
    ax2.set_xlim(h_dates[2376], h_dates[2544])
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.tick_params(axis='y', color='w')
    ax3.plot(h_dates, solar_pv_cf, color=colors[0], label=labels[0])
    ax3.plot(h_dates, onshore_wind_cf, color=colors[1], label=labels[1])
    ax3.plot(h_dates, offshore_wind_cf, color=colors[2], label=labels[2])
    ax3.xaxis.set_minor_locator(DayLocator(interval=2))
    ax3.xaxis.set_major_formatter(DateFormatter('%d/%m'))
    ax3.tick_params(axis='x', rotation=90)
    ax3.tick_params(axis='y', color='w')
    ax3.set_xlim(h_dates[4560], h_dates[4728])
    ax3.spines['left'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.legend(loc='upper right', fontsize=10)
    ax4.plot(h_dates, solar_pv_cf, color=colors[0], label=labels[0])
    ax4.plot(h_dates, onshore_wind_cf, color=colors[1], label=labels[1])
    ax4.plot(h_dates, offshore_wind_cf, color=colors[2], label=labels[2])
    ax4.xaxis.set_minor_locator(DayLocator(interval=2))
    ax4.xaxis.set_major_formatter(DateFormatter('%d/%m'))
    ax4.tick_params(axis='x', rotation=90)
    ax4.tick_params(axis='y', color='w')
    ax4.set_xlim(h_dates[6768], h_dates[6936])
    ax4.spines['left'].set_visible(False)
    plt.show()

## PLOT 3: hydrogen imports capacity

if show_h2_imp:
    path = 'results/results_scenario1_cplex.csv'
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

## PLOT 4: storage operation

if show_stor_op:
    h_index = [h for h in range(8760)]
    h_st, h_nd = h_index[0], h_index[-1]
    h_dates = list()
    h_dates.append(datetime.strptime('2014-01-01', '%Y-%m-%d'))
    for h in h_index[1:]:
        h_dates.append(h_dates[0] + timedelta(hours = h))

    path_scen2 = 'results/results_scenario2_cplex.csv'
    rs2 = pd.read_csv(path_scen2)

    bs_soc = rs2['BATTERY_STORAGE.electricity_stored']
    h2s_soc = rs2['HYDROGEN_STORAGE.hydrogen_stored']
    ng_soc = rs2['NATURAL_GAS_STORAGE.natural_gas_stored']

    ttl_fnt = 10

    fig5 = plt.figure()
    colors = ['gold', 'lightblue', 'salmon']
    ax1 = plt.subplot(311)
    ax1.fill_between(h_dates, bs_soc[h_st:h_nd+1], color=colors[0], alpha=.55)
    ax1.set_title('Battery Storage', fontsize=ttl_fnt)
    ax1.set_yticks([0, 20, 40])
    ax1.axes.get_xaxis().set_visible(False)
    ax2 = plt.subplot(312, sharex=ax1)
    ax2.fill_between(h_dates, h2s_soc[h_st:h_nd+1], color=colors[1], alpha=.65)
    ax2.set_title('Hydrogen Storage', fontsize=ttl_fnt)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.set_ylabel('State of Charge [GWh]')
    ax2.set_yticks([0, 100, 200, 300])
    ax3 = plt.subplot(313, sharex=ax1)
    ax3.fill_between(h_dates, ng_soc[h_st:h_nd+1], color=colors[2], alpha=.75)
    ax3.set_title('Natural Gas Storage', fontsize=ttl_fnt)
    ax3.set_yticks([0, 2500, 5000, 7500])
    ax3.xaxis.set_major_locator(MonthLocator())
    ax3.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
    ax3.xaxis.set_major_formatter(NullFormatter())
    ax3.xaxis.set_minor_formatter(DateFormatter('%b'))
    ax3.set_xlim(h_dates[0], h_dates[-1])
    plt.show()

# PLOT 5: natural gas balance

if show_ng_bal:
    path = 'results_scenario3_cplex.csv'
    rs = pd.read_csv(path)
    tmp = list(rs['NATURAL_GAS_STORAGE.natural_gas'])
    ngsp = np.zeros(len(tmp))
    ngsm = np.zeros(len(tmp))
    for i, v in enumerate(tmp):
        if v > 0.0:
            ngsp[i] = v
        else:
            ngsm[i] = v
    ngsp, ngsm = pd.Series(ngsp), pd.Series(ngsm)

    plt.fill_between(rs.index, 0, rs['METHANATION_PLANTS.methane'])
    plt.fill_between(rs.index, rs['METHANATION_PLANTS.methane'], rs['METHANATION_PLANTS.methane']+rs['NATURAL_GAS_IMPORTS.natural_gas'])
    plt.fill_between(rs.index, rs['METHANATION_PLANTS.methane']+rs['NATURAL_GAS_IMPORTS.natural_gas'], rs['METHANATION_PLANTS.methane']+rs['NATURAL_GAS_IMPORTS.natural_gas']+rs['NATURAL_GAS_DEMAND_RESPONSE.load_shedding'])
    plt.fill_between(rs.index, rs['METHANATION_PLANTS.methane']+rs['NATURAL_GAS_IMPORTS.natural_gas']+rs['NATURAL_GAS_DEMAND_RESPONSE.load_shedding'], rs['METHANATION_PLANTS.methane']+rs['NATURAL_GAS_IMPORTS.natural_gas']+rs['NATURAL_GAS_DEMAND_RESPONSE.load_shedding']+ngsp)
    plt.fill_between(rs.index, -rs['links.heat_demand_natural_gas']-rs['links.industry_demand_natural_gas']+rs['links.legacy_SMR_demand_natural_gas']-rs['links.transport_demand_natural_gas']-rs['SMR_w_PCCC.natural_gas']-rs['CHP_PLANTS_w_PCCC.natural_gas']-rs['CCGT_w_PCCC.natural_gas']-rs['OCGT_w_PCCC.natural_gas']-ngsm)
    plt.show()

## PLOT 6: power generation capacity

if show_gen_capas:
    name = 'results/results_cplex_scenario'
    nr = [str(i+1) for i in range(5)]
    ext = '.csv'
    files = [name+n+ext for n in nr]
    rss = [pd.read_csv(f) for f in files]

    ccgt, ocgt, fc, nk, bm, wst, chp, bt, phs = [], [], [], [], [], [], [], [], []
    bwst, bchp, bccgt, bocgt, bfc, bnk, bbt, bphs = [], [], [], [], [], [], [], []
    for rs in rss:
        bm.append(rs['BIOMASS_POWER_PLANTS_w_PCCC.BIOMASS_POWER_PLANTS.pre_installed_capacity'][0])
        wst.append(rs['WASTE_POWER_PLANTS_w_PCCC.WASTE_POWER_PLANTS.pre_installed_capacity'][0])
        bwst.append(bm[-1])
        chp.append(rs['CHP_PLANTS_w_PCCC.CHP_PLANTS.pre_installed_capacity'][0])
        bchp.append(bwst[-1]+wst[-1])
        nk.append(rs['NUCLEAR_POWER_PLANTS.pre_installed_capacity'][0])
        bnk.append(bchp[-1]+chp[-1])
        ccgt.append(rs['CCGT_w_PCCC.CCGT.capacity'][0])
        bccgt.append(bnk[-1]+nk[-1])
        ocgt.append(rs['OCGT_w_PCCC.OCGT.capacity'][0])
        bocgt.append(bccgt[-1]+ccgt[-1])
        fc.append(rs['HYDROGEN_FUEL_CELLS.capacity'][0])
        bfc.append(bocgt[-1]+ocgt[-1])
        bt.append(rs['BATTERY_STORAGE.capacity_flow'][0])
        bbt.append(bfc[-1]+fc[-1])
        phs.append(rs['PUMPED_HYDRO_STORAGE.pre_installed_capacity_flow'][0])
        bphs.append(bbt[-1]+bt[-1])


    labels = ['Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5']
    width = 0.35   # the width of the bars: can also be len(x) sequence

    fig, ax = plt.subplots()

    ax.bar(labels, bm, width, label='Biomass', color='yellowgreen')
    ax.bar(labels, wst, width, bottom=bwst, label='Waste', color='grey')
    ax.bar(labels, chp, width, bottom=bchp, label='CHP', color='pink')
    ax.bar(labels, nk, width, bottom=bnk, label='Nuclear', color='khaki')
    ax.bar(labels, ccgt, width, bottom=bccgt, label='CCGT', color='salmon')
    ax.bar(labels, ocgt, width, bottom=bocgt, label='OCGT', color='lightsalmon')
    ax.bar(labels, fc, width, bottom=bfc, label='Fuel Cells', color='lightblue')
    ax.bar(labels, bt, width, bottom=bbt, label='Battery', color='gold')
    ax.bar(labels, phs, width, bottom=bphs, label='Pumped-Hydro', color='royalblue')

    ax.set_ylabel('Capacity [GW]', fontsize=10)
    ax.set_ylim(bottom=0.0, top=32.5)
    ax.legend(ncol=3, fontsize=9)

    plt.show()

## PLOT 7: storage capacity

if show_stor_capas:
    name = 'results/results_cplex_scenario'
    nr = [str(i+1) for i in range(5)]
    ext = '.csv'
    files = [name+n+ext for n in nr]
    rss = [pd.read_csv(f) for f in files]

    ch4s, h2s, bt, phs, co2s = [],[],[],[],[]
    bch4s, bh2s, bbt, bphs = [],[],[],[]
    for rs in rss:
        ch4s.append(rs['NATURAL_GAS_STORAGE.pre_installed_capacity_stock'][0])
        h2s.append(rs['HYDROGEN_STORAGE.capacity_stock'][0])
        bh2s.append(ch4s[-1])
        bt.append(rs['BATTERY_STORAGE.capacity_stock'][0])
        bbt.append(bh2s[-1]+h2s[-1])
        phs.append(rs['PUMPED_HYDRO_STORAGE.pre_installed_capacity_stock'][0])
        bphs.append(bbt[-1]+bt[-1])

    labels = ['Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5']
    width = 0.35   # the width of the bars: can also be len(x) sequence

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [7, 1]})
    ax1.bar(labels, ch4s, width, label='Methane Storage', color='salmon')
    ax1.bar(labels, h2s, width, bottom=bh2s, label='Hydrogen Storage', color='lightblue')
    ax1.bar(labels, bt, width, bottom=bbt, label='Battery Storage', color='gold')
    ax1.bar(labels, phs, width, bottom=bphs, label='Pumped-Hydro Storage', color='royalblue')
    ax1.set_ylim(bottom=7925.0, top=8400.0)
    ax2.bar(labels, ch4s, width, label='Methane Storage', color='salmon')
    ax2.bar(labels, h2s, width, bottom=bh2s, label='Hydrogen Storage', color='lightblue')
    ax2.bar(labels, bt, width, bottom=bbt, label='Battery Storage', color='gold')
    ax2.bar(labels, phs, width, bottom=bphs, label='Pumped-Hydro Storage', color='royalblue')
    ax2.set_ylim(bottom=0.0, top=600.0)
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.tick_params(bottom=False, right=True, labeltop=False, labelright=True)  # don't put tick labels at the top
    ax2.tick_params(bottom=True, right=True, labeltop=False, labelright=True)

    #
    d = .005  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (0, 0), **kwargs)        # top-left diagonal
    ax1.plot((1 - d, 1 + d), (0, 0), **kwargs)  # top-right diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1, 1), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1, 1), **kwargs)  # bottom-right diagonal

    # y label and legend
    ax1.set_ylabel('Capacity [GWh]', fontsize=10)
    ax1.yaxis.set_label_coords(-0.1, 0.35)
    ax1.legend(ncol=1, fontsize=9)

    plt.show()

## PLOT 8: gas production capacity

if show_gas_prod_capas:
    name = 'results/results_cplex_scenario'
    nr = [str(i+1) for i in range(5)]
    ext = '.csv'
    files = [name+n+ext for n in nr]
    rss = [pd.read_csv(f) for f in files]

    el, smr, mt = [],[],[]
    bsmr, bmt = [],[]
    for rs in rss:
        el.append(rs['ELECTROLYSIS_PLANTS.capacity'][0])
        smr.append(rs['SMR_w_PCCC.SMR.capacity'][0])
        bsmr.append(el[-1])
        mt.append(rs['METHANATION_PLANTS.capacity'][0])
        bmt.append(bsmr[-1]+smr[-1])

    labels = ['Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4', 'Scenario 5']
    width = 0.35   # the width of the bars: can also be len(x) sequence

    fig, ax = plt.subplots()
    ax.bar(labels, el, width, label='Electrolysis', color='lightblue')
    ax.bar(labels, smr, width, bottom=bsmr, label='SMR', color='steelblue')
    ax.bar(labels, mt, width, bottom=bmt, label='Methanation', color='salmon')
    ax.set_ylim(bottom=0.0, top=4.0)

    # y label and legend
    ax.set_ylabel('Capacity [GWh/h]', fontsize=10)
    ax.legend(ncol=3, fontsize=9)

    plt.show()

## PLOT 9: ev charging distribution

if show_ev_charge_dist:

    path = 'results/results_cplex_scenario1.csv'
    rs = pd.read_csv(path)

    fig, (ax1,ax2) = plt.subplots(2,1,sharex=True, gridspec_kw={'height_ratios':[7,1]})

    ax1.hist(rs['ELECTRICAL_DEMAND_RESPONSE.ev_charge'], 1000, histtype='step', cumulative=True, density=True)
    ax1.set_ylim(bottom=0.625, top=1.05)
    ax1.tick_params(bottom=False)
    ax2.hist(rs['ELECTRICAL_DEMAND_RESPONSE.ev_charge'], 1000, histtype='step', cumulative=True, density=True)
    ax2.set_ylim(bottom=0.0, top=0.125)
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)

    d = .005  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (0, 0), **kwargs)        # top-left diagonal
    ax1.plot((1 - d, 1 + d), (0, 0), **kwargs)  # top-right diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1, 1), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1, 1), **kwargs)  # bottom-right diagonal

    ax1.set_ylabel('Empirical Cumulative Probability [-]', fontsize=10)
    ax1.yaxis.set_label_coords(-0.1, 0.35)
    ax2.set_xlabel('EV Charge [GW]', fontsize=10)

    plt.show()
