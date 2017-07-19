import sys
sys.path.append("../bin/")
import seaborn as sns
cl = sns.color_palette()
from kinetics.ala_kinetics import *
from kinetics.plt_kinetics import *
import os
from kinetics.sorted_lifetimes import *
#from mkdir_p import mkdir_p
import scipy as sp
from kinetics.arrhenius import *
#from kinetics.plot_rna import *


recalc=True
dt=1
delta_t = str(5*dt)
run_name="remd_st1/rb_remd2-1_25jan17_q0.7_dt{}".format(delta_t)
run_name2 = "remd_st2/rb_remd2-2_25jan17_q0.7_dt{}".format(delta_t)

#for st in ["remd_st1", "remd_st2"]:
#   mkdir_p(st)


rep_ar = np.genfromtxt("remd_st1_q/replica_temp.xvg")
tmp_l = np.genfromtxt("remd/tmp_dir_list_r2_n-rb.txt")

# /home/tb/lustelzl/from-bio
q_path = "remd_st1_q/_demux_rbn2-remd2-154ns_rna_f*"



#q_d = make_q_dict(q_path, split_char= "_f")
q_d_tr1 = make_q_dict(q_path, split_char= "_f", last_index=30000, step=dt,
                      verbose=True)
s= [[0.1, 0.7]]

o_remd_st1 = get_rates_pop_multi_temp(q_d_tr1, rep_ar[:30000:dt,:], run_name, state_def=s,
                                      return_raw_events=True, remd=True,
                                      recalc=recalc, split_state_tp=True)

state_df_st1, all_tp_temp_df_st1, pt_st1, rates_st1, events_st1 = o_remd_st1
# reweight
tp_rew, rew_r = calc_ptpq_all_temp(events_st1, q_d_tr1, (0,1), rep_ar[:30000:dt,:], pt_st1)
pd.to_pickle(rew_r,run_name+"_rew_{}.pickle".format(delta_t))
pd.to_pickle(tp_rew, run_name+"_tp-rew_{}.pickle".format(delta_t))

# symmetrize rates
total_time =  rep_ar[:30000:dt,:][-1, 0] - rep_ar[:30000:dt,:][0,0]
print total_time
trans_names=[(0,1), (1,0)]
#rew_sym_r_st1 = sym_counts_calc_rate(rew_r, pt_st1, rew_r.type.unique(), total_time,
#                                     verbose=False)
#pd.to_pickle(rew_sym_r_st1, run_name+"_rew_rs_{}.pickle".format(delta_t)) 

# st2
rep_ar2 = np.genfromtxt("remd_st2_q/replica_temp.xvg")
q_r2="remd_st2_q/_demux_rbn2-remd2-2_rna_f*"

q_d_tr2 = make_q_dict(q_r2, split_char= "_f", last_index=-1,step=dt)


o_remd_st2 = get_rates_pop_multi_temp(q_d_tr2, rep_ar2[::dt], run_name2, state_def=s,
                                      return_raw_events=True, remd=True,
                                      recalc=recalc, split_state_tp=True)

state_df_st2, all_tp_temp_df_st2, pt_st2, rates_st2, events_st2 = o_remd_st2

tp_rew2, rew_r2 = calc_ptpq_all_temp(events_st2, q_d_tr2, (0,1), rep_ar2[::dt,:], pt_st2)

pd.to_pickle(rew_r2, run_name+"_rew_{}.pickle".format(delta_t))
pd.to_pickle(tp_rew2, run_name2+"_tp-rew_{}.pickle".format(delta_t))

# symmetrize counts
#total_time2 = rep_ar2[::dt,:][-1,0] - rep_ar2[::dt,:][0,0]
#rew_sym_r_st2 = sym_counts_calc_rate(rew_r2, pt_st2, rew_r2.type.unique(), total_time2, verbose=True)
#average the two stages
#sym_rew_r_df = av_sym_rates_stages([rew_sym_r_st1 , rew_sym_r_st2 ])


# combine the two REMD stages

total_time_st1 = rep_ar[-1,0] - rep_ar[0,0]
total_time_st2 = rep_ar2[-1,0] - rep_ar2[0,0]
print total_time, total_time_st2
remd_comb_time = total_time + total_time_st2    
trans_l = [(0,1), (1,0)]

r_c_q07, pt_c_q07, tp_c_q07 = combine_rates_stages_simple([pt_st1, pt_st2],
                                                          [tp_rew, tp_rew2],
                                                          remd_comb_time,2 , trans_l,
                                                          return_all=True)
                                                          
                                                      
                                                          
r_cs_q07 = sym_counts_calc_rate(r_c_q07, pt_c_q07, trans_names,remd_comb_time)
r_cs_q07_err = err_log_rate(r_cs_q07, weight_name='sym_weight', diff_from_est="True")

pd.to_pickle(r_cs_q07_err, "rb_remd_r_sym_st1-2_dt{}.pickle".format(delta_t))


# MD

path_md = "md_q"
#fn= "rbn_run2_md_con_*.00"
fn="rbn_run2_md_complete_cf_q*.00"

run_name="rb_md_19dec16_0.7_dt5"
q_dict_pickle=run_name+"q_dict.pickle"

if os.path.exists(q_dict_pickle):
   with open(q_dict_pickle, "rb") as p:
        q_dict_rtr = pickle.load(p)
else:     
     q_dict = make_q_dict(path_md + fn, "_q", first_index=0, last_index=-1)
     test_q_dict_len_d = {}
     for q, val in q_dict.items():
         test_q_dict_len_d[q] = val[-1,0]
     
     print "shortest traj ends at: {}".format(np.min(test_q_dict_len_d.values()))
     q_dict_tr = make_q_dict(path_md + fn, "_q", first_index=0, last_index=998601)

     q_dict_rtr = {}
     for ti, t in enumerate(tmp_l):
         q_dict_rtr[ti] = q_dict_tr[t][::dt*5]
       
     with open(q_dict_pickle,"wb") as w:
          pickle.dump(q_dict_rtr, w)   
             
rep_fn=run_name+"_rep.txt"
if os.path.exists(rep_fn):
   rep_ar = np.genfromtxt(rep_fn)
else:             
     t_pts = q_dict_rtr[10][:,0]
     rep_ar = gen_rep_ar_md(tmp_l, sim_length=len(t_pts), time_pts= t_pts )
     np.savetxt(rep_fn, rep_ar)
     
     
_r_md = get_rates_pop_multi_temp(q_dict_rtr, rep_ar, run_name, state_def=[[0.1, 0.7]], remd=False,
                                 return_raw_events=True, split_state_tp=True, recalc=recalc)

md_state_df, md_all_tp_temp_df, md_pt, md_rates, transition_events = _r_md

md_rs_q07 = sym_counts_calc_rate(md_rates, md_pt, trans_l, rep_ar[-1,0] - rep_ar[0,0])
md_rs_q07_err = err_log_rate(md_rs_q07, weight_name='sym_weight', diff_from_est="True")

# Activation energies


R=8.3144598
Tm=98.0

## Activation energy REMD

rln_sremd_f_c = np.log(r_cs_q07_err[r_cs_q07_err.type==(0,1)].rate.values.astype(np.float64))
rln_sremd_uf_c =np.log(r_cs_q07_err[r_cs_q07_err.type==(1,0)].rate.values.astype(np.float64))

opt_uf13_r = curve_fit(arrhenius, tmp_l[5:-6], rln_sremd_uf_c[5:-5], maxfev=5000)
print opt_uf13_r[0][0] / (Tm*R)
Ea_rremd_uf13c, A_rremd_uf13c  = opt_uf13_r[0]

opt_f13_r = curve_fit(arrhenius, tmp_l[5:-6], rln_sremd_f_c[5:-6], maxfev=2000)
print opt_f13_r[0][0] / (Tm*R)
Ea_rremd_f13c, A_rremd_f13c  = opt_f13_r[0]

remd_xl = [Tm / _t for _t in gen_x_axis_t_rna(r_cs_q07_err, tmp_l)]

## Activation energy MD


# Plot comparison


md_xl= [ Tm / _t for _t in gen_x_axis_t_rna(md_rs_q07_err, tmp_l )]
t_ticks = np.arange(0.94, 1.04, step=0.02)

fig, ax = plt.subplots(figsize=(3.5,3))
sns.set_style('ticks')
#plt.style.use('seaborn-ticks')
md_cl = [cl[5], cl[4]]
md_label_l = ["MD", "MD"]
ax, md_leg = plot_rate_est_err(ax, md_rs_q07_err , trans_l, md_label_l, md_xl,
                               md_cl, return_legend_handles=True, symbol='o',
                               fill_symbol=True, capsize=3.5, elinewidth=1.5)
                               
remd_label_l = ["replica exchange", "replica exchange"]
remd_cl = [cl[0], cl[2]]

ax, remd_leg = plot_rate_est_err(ax, r_cs_q07_err,
                                 trans_l, remd_label_l, remd_xl,
                                 remd_cl, return_legend_handles=True, symbol="s",
                                 capsize=3.5, elinewidth=1.5)
                                 
                                 
l1 = ax.legend([remd_leg[1], md_leg[1]], [remd_label_l[1],md_label_l[1]],
loc=[0.14,0.74], handletextpad=-0.48, fontsize=13,labelspacing=0.2)
l2 = ax.legend([md_leg[0], remd_leg[0]], [md_label_l[0],remd_label_l[0]],
loc=[0.14,-0.02], handletextpad=-0.48, fontsize=13, labelspacing=0.15)

ax.add_artist(l1)
ax.semilogy()
ax.set_xlim(0.937, 1.047)
ax.set_ylim(10**-3.3, 10**1.7)
ax.tick_params(axis='both', which="major", labelsize=13, right=True)
ax.tick_params(axis='both', which="minor", right=True)
ax.set_xlabel(r"$\mathregular{T_m/T}$", fontsize=14)
ax.set_ylabel(r"$\mathregular{k}$ [$\mathregular{1/10^6 \, MD \, steps}$]",
fontsize=14)
#ax.set_ylabel("Rate coefficient", fontsize=14)
#ax.set_xlabel("Inverse temperature", fontsize=14)
ax.set_xticks(t_ticks)
ax.yaxis.set_label_coords(-0.18,0.5)
ax.xaxis.set_label_coords(0.5, -0.18)
ax.tick_params(which='both', direction='out', pad=1.8,
               labelsize=12)

ax.plot(1.0/ (tmp_l[5:-6]/ Tm), A_rremd_uf13c* np.exp( - Ea_rremd_uf13c / (R* tmp_l[5:-6]) ), c=cl[2],lw=1.5)
ax.plot(1.0/ (tmp_l[5:-6] /Tm), A_rremd_f13c * np.exp(-Ea_rremd_f13c / (R* tmp_l[5:-6]) ), '-', c=cl[0], lw=1.5)   
fig.tight_layout()

fig.savefig("rate_plot/test.png")
