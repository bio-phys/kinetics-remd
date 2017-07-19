import seaborn as sns
cl = sns.color_palette()
from kinetics.ala_kinetics import *
from kinetics.kinetics import make_q_dict
import os
from kinetics.sorted_lifetimes import *
import matplotlib
import glob
from mkdir_p import mkdir_p
import argparse

parser = argparse.ArgumentParser() 
parser.add_argument("dt", help="observation interval")
parser.add_argument("run_name", help="name of refinement run")
args = parser.parse_args()

dt=int(args.dt)

print dt

# ## REMD 1
tmp_l = np.genfromtxt(
    "../temp_list_formatted")
#q_path="/home/tb/lustelzl/Projects/kinetics/alanine-dipeptide/remd-1/out-st1/demux_1ps/process/phi_psi/*_phi_psi_p.txt"
q_path="../../ala2/remd/*_phi_psi_p.txt"
q_dict = make_q_dict(q_path, "_", split_pos=0, last_index=-1, step=dt)
print q_dict[0].shape
#rep_ar=np.genfromtxt("/home/tb/lustelzl/Projects/kinetics/alanine-dipeptide/remd-1/out-st1/demux_1ps/ad_demux/replica_temp_full.txt")
rep_ar = np.genfromtxt("../../ala2/remd/replica_temp_full.txt")

print rep_ar[:3]


run_name ="remd_dt{}/ala_remd_st1_dt{}ps_{}".format(dt,dt, args.run_name)
remd_dir_name="remd_dt{}".format(dt)
mkdir_p(remd_dir_name)


s=[[-100, 50], [0, 150]]

o_remd_st1 = get_rates_pop_multi_temp(q_dict, rep_ar[::dt], run_name, state_def=s,
                                      return_raw_events=True, remd=True,
                                      recalc=True, split_state_tp=True,
                                      stop_tp_at_next_end=True, verbose=False)



remd_state_df, remd_all_tp_temp_df, remd_pt, remd_rates, remd_events = o_remd_st1


remd_total_time = rep_ar[-1,0] - rep_ar[0,0]
print remd_total_time


r_s = sym_counts_calc_rate(remd_rates, remd_pt, remd_rates.type.unique(), remd_total_time)


r_s_ln = err_log_rate(r_s, weight_name="sym_weight", diff_from_est=True)
pd.to_pickle(r_s_ln, "{}_rates_sym_ln_{}.pickle".format(run_name, args.run_name))


# ### MD

md1_d = {}
for i, t in enumerate(tmp_l):
    md1_d[i] = np.genfromtxt(
    "../../ala2/md_1/ala1_md_st1_{:.2f}_f.xtc_phi_psi_p.txt".format(
            t))[1::dt]


md2_d = {}
for i, t in enumerate(tmp_l):
    md2_d[i] = np.genfromtxt(
        "../../ala2/md_2/ala1_md_st1_{:.2f}_f.xtc_phi_psi_p.txt".format(
            t))[1::dt]

md3_d = {}
for i, t in enumerate(tmp_l):
    md3_d[i] = np.genfromtxt(
        "../../ala2/md_3/ala1_md_st3_{:.2f}_f.xtc_phi_psi_p.txt".format(
            t))[1::dt]


d = {'md_st1' : md1_d, 'md_st2' : md2_d, 'md_st3' : md3_d}


md_run_name="md_dt{}_{}".format(dt, args.run_name)

md_dir_name="md_dt{}_st".format(dt)
for i in range(1,4):
    mkdir_p(md_dir_name+str(i))

_c = combine_stages_df(md_run_name, s, d, None, remd=False,tmp_l=tmp_l,
                      dir_name=md_dir_name, return_comb_pt=True, recalc=True,
                      n_blocks=None, stop_tp_at_next_end=True)


r_df_md, [st_st_c, tp_st_c, p_st_c, r_st_c], pt_c_md = _c


total_time = 0
for k, v in d.items():
    print v[0].__len__()
    #total_time += v[0].__len__()
    total_time +=  v[0][-1,0] - v[0][0,0]
print total_time

md_dir_c = "c_md_dt{}".format(dt)
mkdir_p(md_dir_c)

rep_md = gen_rep_ar_md(tmp_l, sim_length=total_time)
np.savetxt(md_dir_c+"/rep_ar_md_st1-3_{}.txt".format(args.run_name), rep_md)


pd.to_pickle(r_df_md, md_dir_c+"/rates_md_st1-3_{}.txt".format(args.run_name))

pd.to_pickle(pt_c_md, md_dir_c+"/pt_md_st1-3_{}.txt".format(args.run_name))


# ### symmetrise counts
md_r_s = sym_counts_calc_rate(r_df_md, pt_c_md, r_df_md.type.unique(), total_time)
md_r_s_ln = err_log_rate(md_r_s, weight_name="sym_weight", diff_from_est=True)
pd.to_pickle(md_r_s_ln, md_dir_c+"/rates_md_st1-3_sym_ln_{}.pickle".format(args.run_name))

# analyse lifetimes
# for now just REMD remd_state_df, remd_all_tp_temp_df
trans_from_h = [t for t in possible_transitions(4) if t[:2] == (0,0,)]
trans_from_c = [t for t in possible_transitions(4) if t[:2] == (0,1,)]

_tp0 = remd_all_tp_temp_df[remd_all_tp_temp_df.temperature==0]
_tba0 = remd_state_df[remd_state_df.temperature==0]
    
_dw_h = loop_dwell_trans_temp(_tp0, _tba0, trans_from_h, tu=dt)
_dw_c = loop_dwell_trans_temp(_tp0, _tba0, trans_from_c, tu=dt)
pd.to_pickle(_dw_h, "{}_dw_h_{}.pickle".format(run_name, args.run_name))
pd.to_pickle(_dw_c, "{}_dw_c_{}.pickle".format(run_name, args.run_name))

