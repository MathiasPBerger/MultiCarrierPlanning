import matplotlib
import matplotlib.pyplot as plt
from os import makedirs
from numpy import asarray
from time import strftime
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

class Plot(object):
    def __init__(self, model, data, save_fig=True):
        
        self.model = model
        self.data = data
        self.save_fig = save_fig
           
    def plot(self, folder, start=0, end=-1, start_zoom=4000, end_zoom=4500, 
             figs=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20,21,22,23,24]):
        #figs=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

        font = {'family' : 'serif',
                    'weight' : 'normal',
                    'size'   : 35}
        lines={'linewidth':1}
        xticks={'labelsize':25}
        yticks={'labelsize':25}
        figsizes={'figsize':(40,20)}
        
        legends={'loc':"upper right",'fontsize': 25}
        matplotlib.rc('font', **font)
        matplotlib.rc('lines', **lines) 
        matplotlib.rc('legend', **legends) 
        matplotlib.rc('xtick', **xticks) 
        matplotlib.rc('ytick', **yticks)
        matplotlib.rc('figure', **figsizes)
        
        time = self.data.time[start:end]
        load = list(asarray(list(self.data.lambda_E.values())))[start:end]
        solar_prod = [self.model.P_S[t].value for t in time]
        w_on_prod = [self.model.P_W_on[t].value for t in time]
        solar_onshore_wind_prod = [solar_prod[t]+w_on_prod[t] for t in time]
        w_off_prod = [self.model.P_W_off[t].value for t in time]
        renewable_prod = [solar_prod[t]+w_on_prod[t]+w_off_prod[t] for t in time]
        mismatch = [renewable_prod[t]-load[t] for t in time]
        mismatch_p = []
        for v in mismatch:
            if v>0:
                mismatch_p.append(v)
            else:
                mismatch_p.append(0)
        nuclear_prod = [self.model.P_NK[t].value for t in time]
        PH_prod_pos = []
        PH_prod_neg = []
        PH_prod = [self.model.P_PH[t].value for t in time]
        for t in time:
            if PH_prod[t]>=0.0:
                PH_prod_pos.append(PH_prod[t])
                PH_prod_neg.append(0.0)
            else:
                PH_prod_pos.append(0.0)
                PH_prod_neg.append(-PH_prod[t])
        disp_prod = [self.model.P_BM_e[t].value + self.model.P_WS_e[t].value + self.model.P_CHP_e[t].value for t in time]
        H2_prod = [self.model.P_H2[t].value for t in time]
        NG_prod = [self.model.P_NG_CCGT[t].value for t in time]
        trs = [self.model.P_IE[t].value for t in time]
        trs_im = [self.model.P_I[t].value for t in time]
        trs_ex = [self.model.P_E[t].value for t in time]
        ENS = [self.model.L_E[t].value for t in time]
        
        E_H2 = [self.model.E_H2[t].value for t in time]
        E_PH = [self.model.E_PH[t].value for t in time]
#        E_CH4 = [self.model.E_CH4[t].value for t in time]
        
        P_H2_TRM = [self.model.P_H2_TRM[t].value for t in time]
        P_H2_IDM = [self.model.P_H2_IDM[t].value for t in time]
        P_NG_TRM = [self.model.P_NG_TRM[t].value for t in time]
        
        P_PtG = [self.model.P_PtG[t].value for t in time]
        P_PtH2 = [self.model.P_PtH2[t].value for t in time]
        P_H2 = [self.model.P_H2[t].value for t in time]
        P_H2tCH4 = [self.model.P_H2tCH4[t].value for t in time]
#        P_CH4 = [self.model.P_CH4[t].value for t in time]
        P_CH4tNG = [self.model.P_CH4tNG[t].value for t in time]
        P_CH4tNG_m = [-self.model.P_CH4tNG[t].value for t in time]
        P_PH = [self.model.P_PH[t].value for t in time]
        P_PtPH = [self.model.P_PtPH[t].value for t in time]
#        P_PHtP = [self.model.eta_PHtP * self.model.P_PHtP[t].value for t in time]
        P_PHtP_neg = [-self.model.eta_PHtP * self.model.P_PHtP[t].value for t in time]
        PH_prod_m = [-self.model.P_PH[t].value for t in time]
        P_H2_out_m = [-self.model.P_H2_out[t].value for t in time]
        P_H2tP_m = [-self.model.P_H2tP[t].value for t in time]
        P_H2tP_TRM  = [(-self.model.P_H2_TRM[t].value + P_H2tP_m[t]) for t in time]
        P_H2tP_TRM_IDM = [(P_H2tP_TRM[t] - self.model.P_H2_IDM[t].value) for t in time]
        #P_B = [self.model.P_B[t].value for t in time]
        
        P_PtB = [self.model.P_PtB[t].value*self.model.eta_PtB for t in time]
        P_BtP = [self.model.P_BtP[t].value for t in time]
        P_BtP_neg = [-self.model.P_BtP[t].value for t in time]
        E_B = [self.model.E_B[t].value for t in time]
        S_B = [self.model.S_B.value for t in time]
        K_B = [self.model.K_B.value*self.model.rho_B.value for t in time]
        K_B_dis = [-self.model.K_B.value for t in time]
        
        P_NGtNGS = [self.model.P_NGtNGS[t].value*self.model.eta_NGS for t in time]
#        P_NGStNG = [self.model.P_NGStNG[t].value for t in time]
        P_NGStNG_neg = [-self.model.P_NGStNG[t].value for t in time]
        E_NGS = [self.model.E_NGS[t].value for t in time]
        xi_NGS = [self.model.xi_NGS.value for t in time]
        kappa_NGS = [self.model.kappa_NGS.value*self.model.rho_NGS.value for t in time]
        kappa_NGS_dis = [-self.model.kappa_NGS.value for t in time]
        
        total_curtail = [self.model.P_C[t].value for t in time]
        curtail = [self.model.P_C[t].value for t in time]
        #curtail = [self.model.P_C_p[t].value for t in time]
        #negative_curtail = [self.model.P_C_n[t].value for t in time]
        negative_curtail = [0 for t in time]
        PtG = [self.model.P_PtG[t].value for t in time]
        load_curtail = [load[t]+curtail[t] for t in time]
        load_curtail_PtG = [load[t]+curtail[t]+PtG[t] for t in time]
        load_curtail_PtG_PtPH = [load_curtail_PtG[t] + P_PtPH[t] for t in time]
        load_curtail_PtG_PtPH_PtB = [load_curtail_PtG_PtPH[t] + P_PtB[t] for t in time]
        total_load = [load_curtail_PtG_PtPH_PtB[t] + trs_ex[t] for t in time]
        total_prod = [nuclear_prod[t]+solar_prod[t]+w_on_prod[t]+w_off_prod[t]+PH_prod_pos[t]+disp_prod[t]+H2_prod[t]+NG_prod[t]+ENS[t]+trs_im[t]+P_BtP[t] for t in time]
        
        heat = [self.model.lambda_NG_HIDM[t] for t in time]
        heat_power = [self.model.lambda_NG_HIDM[t] + self.model.P_NGtP_CCGT[t].value for t in time]
        heat_power_NGS = [self.model.lambda_NG_HIDM[t] + self.model.P_NGtP_CCGT[t].value + self.model.P_NGtNGS[t].value for t in time]
        CH4tNG = [-self.model.P_CH4tNG[t].value for t in time]
        CH4tNG_im = [-self.model.P_CH4tNG[t].value -self.model.P_NG_I[t].value for t in time]
        CH4tNG_im_NGS = [-self.model.P_NGStNG[t].value -self.model.P_CH4tNG[t].value -self.model.P_NG_I[t].value for t in time]
        CH4tNG_im_NGS_ENS = [-self.model.P_NGStNG[t].value -self.model.P_CH4tNG[t].value -self.model.P_NG_I[t].value - self.model.L_NG[t].value for t in time]
        
        inst_ren_capa = self.data.kappa_S_0 + self.data.kappa_W_on_0 + self.data.kappa_W_off_0 +\
                        self.model.K_S.value+self.model.K_W_on.value+self.model.K_W_off.value
        inst_ren_capa_vec = [inst_ren_capa for t in time]
        kappa_PtPH = -self.data.kappa_PtPH
        kappa_PtPH_vec = [kappa_PtPH for t in time]
        kappa_PHtP = self.data.kappa_PHtP*self.data.eta_PHtP
        kappa_PHtP_vec = [kappa_PHtP for t in time]
        kappa_disp_max = self.data.kappa_CHP+self.data.kappa_BM+self.data.kappa_WS
        kappa_disp_max_vec = [kappa_disp_max for t in time]
        kappa_disp_min = self.data.mu_CHP*self.data.kappa_CHP+self.data.mu_BM*self.data.kappa_BM+self.data.mu_WS*self.data.kappa_WS
        kappa_disp_min_vec = [kappa_disp_min for t in time]
        kappa_NG_0 = 0
        kappa_NG = self.model.K_NG_CCGT.value
        kappa_NG_vec = [kappa_NG_0+kappa_NG for t in time]
        kappa_NGNet = self.data.kappa_NGNet
        kappa_NGNet_vec = [kappa_NGNet for t in time]   
        kappa_NK = self.data.kappa_NK
        kappa_NK_vec = [kappa_NK[t] for t in time]
        kappa_H2_TRM = self.data.kappa_H2_TRM; kappa_H2_TRM_vec = [kappa_H2_TRM[t] for t in time]
        kappa_H2_IDM = self.data.kappa_H2_IDM; kappa_H2_IDM_vec = [kappa_H2_IDM[t] for t in time]
        kappa_NG_TRM = self.data.kappa_NG_TRM; kappa_NG_TRM_vec = [kappa_NG_TRM[t] for t in time]
        
        el_price= [1000*self.model.theta_IE[t] for t in time]
                            
        if self.save_fig:
            makedirs(folder+"/plots")
            dir_plots = folder+"/plots"
                            
        if 1 in figs:
            
            fig1=plt.figure()
            fig1.suptitle('ENS', fontweight='bold')
            ax1=plt.subplot(211)
            ax1.set_ylabel('ENS (MWh)',fontweight='bold')
            ax1.plot(time, ENS,'b',label="ENS")
            ax1.fill_between(time, 0, ENS, alpha=0.3, label="ENS")
            ax1.legend(loc='upper right',fontsize=10)
            ax2=plt.subplot(2,1,2)
            ax2.set_ylabel('Prod - Cons Mismatch (MWh)',fontweight='bold')
            ax2.set_xlabel('Time',fontweight='bold')
            ax2.plot(time, mismatch,label="Load - Renewable Generation Mismatch")
            ax2.fill_between(time, 0, mismatch, alpha=0.3, label="Prod - Cons")
            ax2.legend(loc='upper right',fontsize=10)
            
            if self.save_fig:
                    fig1.savefig(dir_plots+"/ENS_plot.pdf")
                    plt.show()
                
        if 2 in figs:
            
            fig2 = plt.figure()
            fig2.suptitle("Renewable Production", fontweight="bold")
            plt.ylabel("Renewable Production [MWh/h]", fontweight="bold")
            plt.xlabel("Time",fontweight="bold")
            plt.plot(time,renewable_prod, "w", label="Renewable Production")
            plt.plot(time,inst_ren_capa_vec, "r", label="Total Installed Renewable Capacity")
            plt.fill_between(time, 0, solar_prod, alpha=0.3,label="Solar")
            plt.fill_between(time, solar_prod, solar_onshore_wind_prod, alpha=0.3, label="Onshore Wind")
            plt.fill_between(time, solar_onshore_wind_prod, renewable_prod, alpha=0.3, label="Offshore Wind")
            fig2.legend(loc='upper right',fontsize= 10)
            
            if self.save_fig:
                fig2.savefig(dir_plots+"/renewable_prod_plot.pdf")
                plt.show()
                
        if 3 in figs:        
            fig3 = plt.figure()
            fig3.suptitle("Curtailed Production", fontweight="bold")
            ax1 = plt.subplot(2, 1, 1)
            ax1.set_ylabel("Curtailed Production [MWh/h]", fontweight="bold")
            ax1.plot(time, curtail, "k", label="Curtailed Production")
            ax1.fill_between(time, 0, curtail, alpha=0.3,label="Curtailed Production")
            #ax1.fill_between(time, 0, curtail_S, alpha=0.3,label="Curtailed Solar")
            #ax1.fill_between(time, curtail_S, curtail_S_W_on, alpha=0.3, label="Curtailed Onshore Wind")
            #ax1.fill_between(time, curtail_S_W_on, curtail, alpha=0.3, label="Curtailed Offshore Wind")
            ax1.legend(loc='upper right',fontsize= 10)
            ax1.xaxis.set_visible(False)
            ax2 = plt.subplot(2, 1, 2, sharex=ax1)
            ax2.plot(time, negative_curtail, label="Negative Curtailment (slack)")
            ax2.set_xlabel("Time",fontweight="bold")
            ax2.legend(loc='upper right',fontsize= 10)
            
            if self.save_fig:
                fig3.savefig(dir_plots+"/test_curtailment_plot.pdf")
                plt.show()
        
        if 4 in figs:
                 
            fig4 = plt.figure()
            fig4.suptitle("Energy Balance", fontweight="bold")
            ax1 = plt.subplot(2,1,1)
            ax1.set_ylabel("Energy Balance [MWh/h]", fontweight="bold")
            ax1.plot(time,total_load, "w", label="Load + Curtailment + PtG + PtPH")
            ax1.fill_between(time, 0, load, alpha=0.3,label="Load")
            ax1.fill_between(time, load, load_curtail, alpha=0.3, label="Curtailment")
            ax1.fill_between(time, load_curtail, load_curtail_PtG, alpha=0.3, label="PtG")
            ax1.fill_between(time, load_curtail_PtG, load_curtail_PtG_PtPH, alpha=0.3, label="PtPH")
            ax1.fill_between(time, load_curtail_PtG_PtPH, load_curtail_PtG_PtPH_PtB, alpha=0.3, label="PtB")
            ax1.fill_between(time, load_curtail_PtG_PtPH_PtB, total_load, alpha=0.3, label="Exports")
            ax1.legend(loc='upper right',fontsize=7)
            ax1.set_ylim([0, 1.1*max(total_load)])
            ax2=plt.subplot(2,1,2, sharex=ax1)
            ax2.set_xlabel('Time',fontweight='bold')
            ax2.plot(time,total_prod,'w',label="Solar + Wind + PHtP + Disp + H2 + NG + ENS")
            ax2.fill_between(time, 0, nuclear_prod, alpha=0.3, label="Nuclear")
            ax2.fill_between(time, nuclear_prod, [nuclear_prod[t]+solar_prod[t] for t in time], alpha=0.3,label="Solar")
            ax2.fill_between(time, [nuclear_prod[t]+solar_prod[t] for t in time], [nuclear_prod[t]+solar_prod[t]+w_on_prod[t] for t in time], alpha=0.3,label="Onshore Wind")
            ax2.fill_between(time, [nuclear_prod[t]+solar_prod[t]+w_on_prod[t] for t in time], [nuclear_prod[t]+renewable_prod[t] for t in time], alpha=0.3,label="Offshore Wind")
            ax2.fill_between(time, [nuclear_prod[t]+renewable_prod[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t] for t in time], alpha=0.3, label="PHtP")
            ax2.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t] for t in time], alpha=0.3,label="NG")
            ax2.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t] for t in time], alpha=0.3, label="H2")
            ax2.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t] for t in time], alpha=0.3, label="Disp")
            ax2.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t]+trs_im[t] for t in time], alpha=0.3, label="Imports")
            ax2.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t]+trs_im[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t]+trs_im[t]+P_BtP[t] for t in time], alpha=0.3, label="BtP")
            ax2.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t]+trs_im[t]+P_BtP[t] for t in time], total_prod, alpha=0.3, label="ENS")
            ax2.legend(loc='upper right', fontsize=5)
            ax2.set_ylim([0, 1.1*max(total_prod)])
            
            if self.save_fig:
                fig4.savefig(dir_plots+"/energy_balance_plot.pdf")
                plt.show()
         
        if 5 in figs:
            fig5 = plt.figure()
            fig5.suptitle("Energy Storage SOCs")
            ax1 = plt.subplot(2, 1, 1)
            ax1.plot(time, E_H2, 'b',label="H2 Storage")
            ax1.fill_between(time, 0, E_H2, alpha=0.3, label="H2 Storage SOC")
            ax1.legend(loc='upper right',fontsize=10)
            ax2 = plt.subplot(2, 1, 2)
            ax2.plot(time,E_PH, 'b',label="PH Storage")
            ax2.fill_between(time, 0, E_PH, alpha=0.3, label="PH Storage SOC")
            ax2.legend(loc='upper right', fontsize=10)
            
            if self.save_fig:
                fig5.savefig(dir_plots+"/SOC_plot.pdf")
                plt.show()
                
        if 6 in figs:
            fig6 = plt.figure()
            fig6.suptitle("PtG System Powers")
            ax1 = plt.subplot(4, 1, 1)
            ax1.plot(time, P_PtG, 'b',label="P_PtG Curve")
            ax1.fill_between(time, 0, P_PtG, alpha=0.3,label="P_PtG")
            ax1.xaxis.set_visible(False)
            ax1.legend(loc='upper right',fontsize=10)
            ax2 = plt.subplot(4, 1, 2, sharex=ax1)
            ax2.plot(time, P_H2, 'b',label="P_H2 Curve")
            ax2.fill_between(time, 0, P_H2, alpha=0.3,label="P_H2")
            ax2.xaxis.set_visible(False)
            ax2.legend(loc='upper right', fontsize=10)
            ax3 = plt.subplot(4, 1, 3, sharex=ax1)
            ax3.plot(time, P_H2tCH4, 'b',label="P_H2tCH4 Curve")
            ax3.fill_between(time, 0, P_H2tCH4, alpha=0.3,label="P_H2tCH4")
            ax3.xaxis.set_visible(False)
            ax3.legend(loc='upper right', fontsize=10)
            ax4 = plt.subplot(4, 1, 4, sharex=ax1)
            ax4.plot(time, P_CH4tNG, 'b',label="P_CH4tNG Curve")
            ax4.fill_between(time, 0, P_CH4tNG, alpha=0.3,label="P_CH4tNG")
            ax4.set_xlabel("Time")
            ax4.legend(loc='upper right', fontsize=10)
            
            if self.save_fig:
                fig6.savefig(dir_plots+"/PtG_powers_plot.pdf")
                plt.show()
                
        if 7 in figs:
            
            fig7 = plt.figure()
            fig7.suptitle("Load - Renewable Generation Mismatch")
            ax1 = plt.subplot(3, 1, 1)
            ax1.plot(time, load, 'b',label="Load Curve")
            ax1.fill_between(time, 0, load, alpha=0.3,label="Load")
            ax1.legend(loc='upper right',fontsize=10)
            ax2 = plt.subplot(3, 1, 2, sharex=ax1)
            ax2.plot(time, renewable_prod, 'b',label="Renewable Generation Curve")
            ax2.fill_between(time, 0, renewable_prod, alpha=0.3,label="Renewable Generation")
            ax2.legend(loc='upper right', fontsize=10)
            ax3 = plt.subplot(3, 1, 3, sharex=ax1)
            ax3.plot(time, mismatch, 'b',label="Load - Renewable Generation Mismatch Curve")
            ax3.fill_between(time, 0, mismatch, alpha=0.3,label="Load - Renewable Generation Mismatch")
            ax3.set_xlabel("Time")
            ax3.legend(loc='upper right', fontsize=10)
            
            if self.save_fig:
                fig7.savefig(dir_plots+"/load_renewable_mismatch_plot.pdf")
                plt.show()
                
        if 8 in figs:
            fig8 = plt.figure()
            fig8.suptitle("PH, Disp, NG Powers")
            ax1 = plt.subplot(3, 1, 1)
            ax1.plot(time, P_PH, 'b',label="P_PH Curve")
            ax1.plot(time, kappa_PtPH_vec, "r", label="PtPH Capacity")
            ax1.plot(time, kappa_PHtP_vec, "r", label="PHtP Capacity")
            ax1.fill_between(time, 0, P_PH, alpha=0.3,label="P_PH")
            ax1.legend(loc='upper right',fontsize=10)
            ax2 = plt.subplot(3, 1, 2, sharex=ax1)
            ax2.plot(time,disp_prod, 'b',label="P_disp Curve")
            ax2.plot(time, kappa_disp_max_vec, "r", label="Disp Capacity")
            ax2.plot(time, kappa_disp_min_vec, "r", label="Disp Min Prod Level")
            ax2.fill_between(time, 0, disp_prod, alpha=0.3,label="P_disp")
            ax2.legend(loc='upper right', fontsize=10)
            ax3 = plt.subplot(3, 1, 3, sharex=ax1)
            ax3.plot(time, NG_prod, 'b',label="P_NG Curve")
            ax3.plot(time, kappa_NG_vec, "r", label="NG Capacity")
            ax3.fill_between(time, 0, NG_prod, alpha=0.3,label="P_NG")
            ax3.set_xlabel("Time")
            ax3.legend(loc='upper right', fontsize=10)
            
            if self.save_fig:
                fig8.savefig(dir_plots+"/PH_disp_NG_powers_plot.pdf")
                plt.show()
                
        if 9 in figs:
                 
            fig9 = plt.figure()
            fig9.suptitle("Energy Balance", fontweight="bold")
            ax1 = plt.subplot(3, 1, 1)
            ax1.plot(time,total_prod,'w',label="Total Production")
            ax1.fill_between(time, 0, nuclear_prod, alpha=0.3, label="Nuclear")
            ax1.fill_between(time, nuclear_prod, [nuclear_prod[t]+solar_prod[t] for t in time], alpha=0.3,label="Solar")
            ax1.fill_between(time, [nuclear_prod[t]+solar_prod[t] for t in time], [nuclear_prod[t]+solar_prod[t]+w_on_prod[t] for t in time], alpha=0.3,label="Onshore Wind")
            ax1.fill_between(time, [nuclear_prod[t]+solar_prod[t]+w_on_prod[t] for t in time], [nuclear_prod[t]+renewable_prod[t] for t in time], alpha=0.3,label="Offshore Wind")
            ax1.fill_between(time, [nuclear_prod[t]+renewable_prod[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t] for t in time], alpha=0.3, label="PHtP")
            ax1.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t] for t in time], alpha=0.3,label="NG")
            ax1.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t] for t in time], alpha=0.3, label="H2")
            ax1.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t] for t in time], alpha=0.3, label="Disp")
            ax1.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t]+trs_im[t] for t in time], alpha=0.3, label="Imports")
            ax1.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t]+trs_im[t] for t in time], [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t]+trs_im[t]+P_BtP[t] for t in time], alpha=0.3, label="BtP")
            ax1.fill_between(time, [nuclear_prod[t]+renewable_prod[t]+PH_prod_pos[t]+NG_prod[t]+H2_prod[t]+disp_prod[t]+trs_im[t]+P_BtP[t] for t in time], total_prod, alpha=0.3, label="ENS")
            ax1.legend(loc='upper right', fontsize=5)
            ax1.set_ylim([0, 1.1*max(total_prod)])
            ax2 = plt.subplot(3, 1, 2, sharex=ax1)
            ax2.plot(time, E_H2, 'b',label="H2 Storage")
            ax2.fill_between(time, 0, E_H2, alpha=0.3,label="H2 Storage SOC")
            ax2.legend(loc='upper right',fontsize=10)
            ax3 = plt.subplot(3, 1, 3, sharex=ax1)
            ax3.plot(time, E_PH, 'b',label="PH Storage")
            ax3.fill_between(time, 0, E_PH, alpha=0.3,label="PH Storage SOC")
            ax3.legend(loc='upper right', fontsize=10)
            ax3.set_xlabel('Time',fontweight='bold')
            
            if self.save_fig:
                fig9.savefig(dir_plots+"/PH_disp_NG_powers_plot.pdf")
                plt.show()
                
        if 10 in figs:
            
            fig10 = plt.figure()
            fig10.suptitle("H2 Storage", fontweight="bold")
            ax1=plt.subplot(2,1,1)
            ax1.plot(time, E_H2 , 'b',label="H2 Storage")
            ax1.fill_between(time, 0, E_H2, alpha=0.3,label="H2 Storage SOC")
            ax1.legend(loc='upper right', fontsize=10)
            ax2 = plt.subplot(2,1,2, sharex=ax1)
            ax2.plot(time, P_PtH2, 'g',label="P_PtH2")
            ax2.plot(time, P_H2_out_m, 'r',label="P_H2_out")
            ax2.fill_between(time, P_H2_out_m, P_H2tP_TRM_IDM, alpha=0.3,label="P_H2tCH4")
            ax2.fill_between(time, P_H2tP_TRM_IDM, P_H2tP_TRM, alpha=0.3,label="P_H2_IDM")
            ax2.fill_between(time, P_H2tP_TRM, P_H2tP_m, alpha=0.3,label="P_H2_TRM")
            ax2.fill_between(time, P_H2tP_m, 0, alpha=0.3,label="P_H2tP")
            ax2.fill_between(time, 0, P_PtH2, alpha=0.3,label="P_PtH2")
            ax2.legend(loc='upper right',fontsize=10)
            ax2.set_xlabel('Time',fontweight='bold')
            
            if self.save_fig:
                fig10.savefig(dir_plots+"/H2_storage_plot.pdf")
                plt.show()
                
#        if 11 in figs:
#            
#            fig11 = plt.figure()
#            fig11.suptitle("CH4 Storage", fontweight="bold")
#            ax1=plt.subplot(2,1,1)
#            ax1.plot(time,E_CH4, 'b',label="CH4 Storage")
#            ax1.fill_between(time, 0, E_CH4, alpha=0.3,label="CH4 Storage SOC")
#            ax1.legend(loc='upper right', fontsize=10)
#            ax2 = plt.subplot(2,1,2, sharex=ax1)
#            ax2.plot(time, P_CH4, 'r',label="P_CH4")
#            ax2.plot(time, P_CH4tNG_m, 'b',label="P_CH4tNG")
#            ax2.fill_between(time, P_CH4tNG_m, 0, alpha=0.3,label="P_CH4tNG")
#            ax2.fill_between(time, 0, P_CH4, alpha=0.3,label="P_CH4")
#            ax2.legend(loc='upper right',fontsize=10)
#            ax2.set_xlabel('Time',fontweight='bold')
#            
#            if self.save_fig:
#                fig11.savefig(dir_plots+"/CH4_storage_plot.pdf")
#                plt.show()
                
        if 12 in figs:
            
            fig12 = plt.figure()
            fig12.suptitle("NG Linepack", fontweight="bold")
            plt.plot(time, kappa_NGNet_vec, "r", label="Linepack Capacity")
            plt.fill_between(time, 0, heat, label="Heating")
            plt.fill_between(time, heat, heat_power, label="Power")
            plt.fill_between(time, heat_power, heat_power_NGS, label="NGtNGS")
            plt.fill_between(time, CH4tNG_im_NGS_ENS, CH4tNG_im_NGS, label="ENS")
            plt.fill_between(time, CH4tNG_im_NGS, CH4tNG_im, label="NGStNG")
            plt.fill_between(time, CH4tNG_im, CH4tNG, label="Imports")
            plt.fill_between(time, CH4tNG, 0, label="CH4tNG")
            plt.legend(loc='upper right',fontsize=10)
            
            if self.save_fig:
                fig12.savefig(dir_plots+"/NG_network_plot.pdf")
                plt.show()
                
        if 13 in figs:
            
            fig13 = plt.figure()
            fig13.suptitle("PH Storage", fontweight="bold")
            ax1=plt.subplot(3, 1, 1)
            ax1.plot(time, E_PH, 'b',label="PH Storage SOC")
            ax1.fill_between(time, 0, E_PH, alpha=0.3,label="PH Storage SOC")
            ax1.legend(loc='upper right', fontsize=10)
            ax2 = plt.subplot(3, 1, 2, sharex=ax1)
            ax2.plot(time, PH_prod_m, 'b',label="P_PH")
            ax2.fill_between(time, PH_prod_m, 0, alpha=0.3, color="blue",label="P_PH")
            ax2.fill_between(time, 0, PH_prod_m, alpha=0.3, color="blue",)
            ax2.legend(loc='upper right',fontsize=10)
            ax3 = plt.subplot(3, 1, 3)
            ax3.set_xlabel('Time',fontweight='bold')
            ax3.plot(time, P_PtPH, 'b',label="P_PtPH")
            ax3.plot(time, P_PHtP_neg, 'b',label="P_PHtP")
            ax3.fill_between(time, P_PHtP_neg, 0, alpha=0.3, color="blue",label="P_PHtP")
            ax3.fill_between(time, 0, P_PtPH, alpha=0.3, color="blue",label="P_PtPH")
            ax3.legend(loc='upper right',fontsize=10)
            
            if self.save_fig:
                fig13.savefig(dir_plots+"/PH_storage_plot.pdf")
                plt.show()
                
        if 14 in figs:
            
            fig14 = plt.figure()
            fig14.suptitle("PH Power Analysis", fontweight="bold")
            ax1=plt.subplot(3, 1, 1)
            ax1.plot(time, mismatch, 'b',label="Load - Renewable Generation Mismatch Curve")
            ax1.fill_between(time, 0, mismatch, alpha=0.3,label="Load - Renewable Generation Mismatch")
            ax1.legend(loc="upper right", fontsize=10)
            ax1.xaxis.set_visible(False)
            ax2=plt.subplot(3, 1, 2, sharex=ax1)
            ax2.plot(time, P_PtPH, 'b',label="P_PtPH")
            ax2.plot(time, P_PHtP_neg, 'b',label="P_PHtP")
            ax2.fill_between(time, P_PHtP_neg, 0, alpha=0.3, color="blue",label="P_PHtP")
            ax2.fill_between(time, 0, P_PtPH, alpha=0.3, color="blue",label="P_PtPH")
            ax2.legend(loc='upper right',fontsize=10)
            ax2.xaxis.set_visible(False)
            ax3=plt.subplot(3, 1, 3, sharex=ax1)
            ax3.plot(time,curtail, "b", label="Curtailed Production")
            #ax3.plot(time, mismatch_p, "r")
            ax3.fill_between(time, 0, curtail, alpha=0.3,label="Curtailed Power")
            ax3.legend(loc='upper right',fontsize=10)
            ax3.set_xlabel('Time',fontweight='bold')
            
            if self.save_fig:
                fig14.savefig(dir_plots+"/renewables_PH_output_plot.pdf")
                plt.show()
                
        if 15 in figs:
            
            fig15 = plt.figure()
            fig15.suptitle("Interconnection Use", fontweight="bold")
            plt.plot(time, trs, 'b',label="Interconnection Power")
            plt.fill_between(time, 0, trs, alpha=0.3,label="Interconnection Power")
            plt.xlabel('Time',fontweight='bold')
            plt.ylabel('Interconnection Power')
            plt.legend(loc='upper right',fontsize=10)
            
            if self.save_fig:
                fig15.savefig(dir_plots+"/interconnection_plot.pdf")
                plt.show()
                
        if 16 in figs:
            
            fig16 = plt.figure()
            fig16.suptitle("Curtailment Analysis", fontweight="bold")
            ax1=plt.subplot(3, 1, 1)
            ax1.plot(time,mismatch, 'b',label="Load - Renewable Generation Mismatch Curve")
            ax1.fill_between(time, 0, mismatch, alpha=0.3,label="Load - Renewable Generation Mismatch")
            ax2=plt.subplot(3, 1, 2, sharex=ax1)
            ax2.plot(time,trs, "w", label="Interconnection")
            ax2.fill_between(time, 0, trs, alpha=0.3,label="Interconnection")
            ax2.legend(loc='upper right',fontsize=10)
            ax3=plt.subplot(3, 1, 3, sharex=ax1)
            ax3.plot(time, curtail, "b", label="Curtailed Production")
            ax3.plot(time, mismatch_p, "r", label="Renewable Production Surplus")
            ax3.fill_between(time, 0, curtail, alpha=0.3,label="Curtailed Power")
            ax3.legend(loc='upper right',fontsize=10)
            ax3.set_xlabel('Time',fontweight='bold')
            
            if self.save_fig:
                fig16.savefig(dir_plots+"/curtailment_analysis_plot.pdf")
                plt.show()
                
        if 17 in figs:
            
            fig17 = plt.figure()
            fig17.suptitle("PH Power Analysis", fontweight="bold")
            ax1=plt.subplot(4, 1, 1)
            ax1.plot(time, mismatch, 'b',label="Load - Renewable Generation Mismatch Curve")
            ax1.fill_between(time, 0, mismatch, alpha=0.3,label="Load - Renewable Generation Mismatch")
            ax1.legend(loc='upper right', fontsize=10)
            ax2=plt.subplot(4, 1, 2, sharex=ax1)
            ax2.plot(time, P_PtPH, 'b', label="P_PtPH")
            ax2.plot(time, P_PHtP_neg, 'b', label="P_PHtP")
            ax2.fill_between(time, P_PHtP_neg, 0, alpha=0.3, color="blue",label="P_PHtP")
            ax2.fill_between(time, 0, P_PtPH, alpha=0.3, color="blue",label="P_PtPH")
            ax2.legend(loc='upper right',fontsize=10)
            ax3=plt.subplot(4, 1, 3, sharex=ax1)
            ax3.plot(time, curtail, "b", label="Curtailed Production")
            ax3.plot(time, mismatch_p, "r", label="Renewable Production Surplus")
            ax3.fill_between(time, 0, curtail, alpha=0.3,label="Curtailed Power")
            ax3.legend(loc='upper right',fontsize=10)
            ax3.set_xlabel('Time',fontweight='bold')
            ax4=plt.subplot(4, 1, 4, sharex=ax1)
            ax4.plot(time, disp_prod, "b", label="Dispatchable")
            ax4.fill_between(time, 0, disp_prod, alpha=0.3,label="Dispatchable")
            ax4.legend(loc='upper right',fontsize=10)
            ax4.set_xlabel('Time',fontweight='bold')
            
            if self.save_fig:
                fig17.savefig(dir_plots+"/curtailment_analysis_plot_2.pdf")
                plt.show()
                
        if 18 in figs:
            
            fig18 = plt.figure()
            fig18.suptitle("Nuclear Power", fontweight="bold")
            plt.plot(time, nuclear_prod, label="Nuclear")
            plt.plot(time, kappa_NK_vec, label="Nuclear Capacity")
            plt.fill_between(time, nuclear_prod, alpha=0.3, label="Nuclear")
            plt.xlabel("Time", fontweight="bold")
            plt.ylabel("Nuclear Production [MW]")
            plt.legend(loc='upper right',fontsize=10)
            
            if self.save_fig:
                fig18.savefig(dir_plots+"/nuclear_plot.pdf")
                plt.show()
                
        if 19 in figs:
            
            fig19 = plt.figure()
            fig19.suptitle("Curtailment Analysis II", fontweight="bold")
            ax1=plt.subplot(3, 1, 1)
            ax1.plot(time, total_curtail, 'b',label="P_C")
            ax1.fill_between(time, 0, total_curtail, alpha=0.3, label="P_C")
            ax2=plt.subplot(3, 1, 2, sharex=ax1)
            ax2.plot(time, curtail, "b", label="P_C_p")
            ax2.fill_between(time, 0, curtail, alpha=0.3,label="P_C_p")
            ax2.legend(loc='upper right',fontsize=10)
            ax3=plt.subplot(3, 1, 3, sharex=ax1)
            ax3.plot(time, negative_curtail, "b", label="P_C_n")
            ax3.fill_between(time, 0, negative_curtail, alpha=0.3,label="P_C_n")
            ax3.legend(loc='upper right',fontsize=10)
            ax3.set_xlabel('Time',fontweight='bold')
            
            if self.save_fig:
                fig19.savefig(dir_plots+"/curtailment_analysis_II_plot.pdf")
                plt.show()
    
        if 20 in figs:
            
            fig20, ax = plt.subplots()
            fig20.suptitle("Renewable Production", fontweight="bold")
            ax.set_ylabel("Renewable Production [MWh/h]", fontweight="bold")
            ax.set_xlabel("Time",fontweight="bold")
            ax.plot(time,renewable_prod, "w", label="Renewable Production")
            #ax.plot(time,inst_ren_capa_vec, "r", label="Total Installed Renewable Capacity")
            ax.fill_between(time, 0, solar_prod, alpha=0.3,label="Solar")
            ax.fill_between(time, solar_prod, solar_onshore_wind_prod, alpha=0.3, label="Onshore Wind")
            ax.fill_between(time, solar_onshore_wind_prod, renewable_prod, alpha=0.3, label="Offshore Wind")
            fig20.legend(loc='upper left',fontsize= 10)
            
            axins = inset_axes(ax, width="25%", height="30%", loc=1)
            x1, x2, y1, y2 = start_zoom, end_zoom, 0, max(renewable_prod[start_zoom:end_zoom])
            axins.set_xlim(x1, x2)
            axins.set_ylim(y1, y2)
            ax.plot(time,renewable_prod, "w", label="Renewable Production")
            axins.fill_between(time, 0, solar_prod, alpha=0.3)
            axins.fill_between(time, solar_prod, solar_onshore_wind_prod, alpha=0.3)
            axins.fill_between(time, solar_onshore_wind_prod, renewable_prod, alpha=0.3)
            plt.yticks(visible=False)
            plt.xticks(visible=False)
            mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0")
            
            if self.save_fig:
                fig20.savefig(dir_plots+"/renewable_prod_plot.pdf")
                plt.show()
                
        if 21 in figs:
            
            fig21 = plt.figure()
            fig21.suptitle("Battery Storage", fontweight="bold")
            ax1=plt.subplot(2, 1, 1)
            ax1.plot(time, E_B, 'b',label="Battery Storage SOC")
            ax1.plot(time, S_B, color='red', label="Battery Energy Capacity")
            ax1.fill_between(time, 0, E_B, alpha=0.3,label="Battery Storage SOC")
            ax1.legend(loc='upper right', fontsize=10)
            ax2 = plt.subplot(2, 1, 2, sharex=ax1)
            ax2.set_xlabel('Time',fontweight='bold')
            ax2.plot(time, P_PtB, 'b',label="P_PtB")
            ax2.plot(time, P_BtP_neg, 'b',label="P_BtP")
            ax2.plot(time, K_B_dis, color='red', label="Battery Dis. Power Capacity")
            ax2.plot(time, K_B, color='red', label="Battery Ch. Power Capacity")
            ax2.fill_between(time, P_BtP_neg, 0, alpha=0.3, color="blue",label="P_BtP")
            ax2.fill_between(time, 0, P_PtB, alpha=0.3, color="blue",label="P_PtB")
            ax2.legend(loc='upper right',fontsize=10)
            
            if self.save_fig:
                fig21.savefig(dir_plots+"/battery_plot.pdf")
                plt.show()
                
        if 22 in figs:
            
            fig22 = plt.figure()
            fig22.suptitle("NG Storage", fontweight="bold")
            ax1=plt.subplot(2, 1, 1)
            ax1.plot(time, E_NGS, 'b',label="NG Storage SOC")
            ax1.plot(time, xi_NGS, color='red', label="NG Energy Capacity")
            ax1.fill_between(time, 0, E_NGS, alpha=0.3,label="NG Storage SOC")
            ax1.legend(loc='upper right', fontsize=10)
            ax2 = plt.subplot(2, 1, 2, sharex=ax1)
            ax2.set_xlabel('Time',fontweight='bold')
            ax2.plot(time, P_NGtNGS, 'b',label="P_NGtNGS")
            ax2.plot(time, P_NGStNG_neg, 'b',label="P_NGStNG")
            ax2.plot(time, kappa_NGS_dis, color='red', label="NGS Dis. Power Capacity")
            ax2.plot(time, kappa_NGS, color='red', label="NGS Ch. Power Capacity")
            ax2.fill_between(time, P_NGStNG_neg, 0, alpha=0.3, color="blue",label="P_NGStNG")
            ax2.fill_between(time, 0, P_NGtNGS, alpha=0.3, color="blue",label="P_NGtNGS")
            ax2.legend(loc='upper right',fontsize=10)
            
            if self.save_fig:
                fig22.savefig(dir_plots+"/NGS_plot.pdf")
                plt.show()
                
        if 23 in figs:
            
            fig23 = plt.figure()
            fig23.suptitle("RES - Wholesale Electricity Prices", fontweight="bold")
            ax1=plt.subplot(4, 1, 1)
            ax1.fill_between(time, 0, solar_prod, alpha=0.3,color="yellow", label="Solar")
            ax1.legend(loc='upper right', fontsize=10)
            ax2 = plt.subplot(4, 1, 2, sharex=ax1)
            ax2.fill_between(time, 0, w_on_prod, alpha=0.3, color="green",label="Won")
            ax2.legend(loc='upper right',fontsize=10)
            ax3 = plt.subplot(4, 1, 3, sharex=ax1)
            ax3.fill_between(time, 0, w_off_prod, alpha=0.3, color="red",label="Woff")
            ax3.legend(loc='upper right',fontsize=10)
            ax4 = plt.subplot(4, 1, 4, sharex=ax1)
            ax4.fill_between(time, 0, el_price, alpha=0.3, color="blue",label="Elprice")
            ax4.legend(loc='upper right',fontsize=10)
            
            if self.save_fig:
                fig23.savefig(dir_plots+"/elprice_plot.pdf")
                plt.show()
                
                
        if 24 in figs:
            
            fig24 = plt.figure()
            fig24.suptitle("Other Gas Markets", fontweight="bold")
            ax1=plt.subplot(3, 1, 1)
            ax1.plot(time, P_H2_TRM, 'g', label="H2 TR demand")
            ax1.plot(time, kappa_H2_TRM_vec, color='r', label="H2 TR cap")
            ax1.fill_between(time, 0, P_H2_TRM, alpha=0.3, color='g', label="H2 TR demand")
            ax1.legend(loc='upper right', fontsize=10)
            ax2 = plt.subplot(3, 1, 2, sharex=ax1)
            ax2.plot(time, P_H2_IDM, 'b',label="H2 IN demand")
            ax2.plot(time, kappa_H2_IDM_vec, 'r',label="H2 IN cap")
            ax2.fill_between(time, 0, P_H2_IDM, alpha=0.3, color='b', label="H2 ID demand")            
            ax2.legend(loc='upper right',fontsize=10)
            ax3 = plt.subplot(3, 1, 3, sharex=ax1)
            ax3.plot(time, P_NG_TRM, 'olive',label="NG TR demand")
            ax3.plot(time, kappa_NG_TRM_vec, 'r',label="NG TR cap")
            ax3.fill_between(time, 0, P_NG_TRM, alpha=0.3, label="H2 ID demand")            
            ax3.legend(loc='upper right',fontsize=10)            
            ax3.set_xlabel('Time',fontweight='bold')
            if self.save_fig:
                fig24.savefig(dir_plots+"/other_markets.pdf")
                plt.show()
                
        return None
            