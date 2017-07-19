
# coding: utf-8

# # Table of Contents
#  <p><div class="lev1 toc-item"><a data-toc-modified-id="REMD-1" href="#REMD"><span class="toc-item-num">1&nbsp;&nbsp;</span>REMD</a></div><div class="lev1 toc-item"><a data-toc-modified-id="MD-2" href="#MD"><span class="toc-item-num">2&nbsp;&nbsp;</span>MD</a></div>

# In[82]:

import seaborn as sns
cl = sns.color_palette()
from kinetics.ala_kinetics import *
from kinetics.kinetics import make_q_dict
import os
from kinetics.sorted_lifetimes import *
import matplotlib
import glob


# In[83]:

recalc=True


# In[84]:

# "../../../alanine-dipeptide/remd-1/temp_list_formatted"
tmp_l = np.genfromtxt("tem")


# # REMD

# In[85]:

#q_path="/home/tb/lustelzl/Projects/kinetics/alanine-dipeptide/remd-1/out-st1/demux_1ps/process/phi_psi/*_phi_psi_p.txt"
q_path = "./remd/*_phi_psi_p.txt"


# In[86]:

get_ipython().system(u' ls -lrt /home/tb/lustelzl/Projects/kinetics/alanine-dipeptide/remd-1/out-st1/demux_1ps/process/phi_psi/*_phi_psi_p.txt')


# In[87]:

get_ipython().system(u' ls remd')


# In[ ]:




# In[88]:

q_dict = make_q_dict(q_path, "_", split_pos=0, last_index=-1)


# In[89]:

rep_ar = np.genfromtxt('remd/replica_temp_full.txt')
run_name ="remd-out/ala_remd_st1_dt1ps_23may17"
s=[[-100, 50], [0, 150]]


# In[90]:

o_remd_st1 = get_rates_pop_multi_temp(q_dict, rep_ar, run_name, state_def=s,
                                      return_raw_events=True, remd=True,
                                      recalc=recalc, split_state_tp=True,
                                      stop_tp_at_next_end=True, verbose=False)


# In[91]:

state_df, all_tp_temp_df, pt, rates, events = o_remd_st1


# # MD

# In[92]:

# "/home/tb/lustelzl/Projects/kinetics/alanine-dipeptide/run1/out-st1-multi-temp/phi_psi/ala1_md_st1_{:.2f}_f.xtc_phi_psi_p.txt"
# /home/tb/lustelzl/Projects/kinetics/alanine-dipeptide/run1/out-st2-multi-temp/phi_psi/ala1_md_st1_{:.2f}_f.xtc_phi_psi_p.txt
# /home/tb/lustelzl/Projects/kinetics/alanine-dipeptide/run1/out-st3-multi-temp/phi_psi/ala1_md_st3_{:.2f}_f.xtc_phi_psi_p.txt


# In[93]:

get_ipython().system(u' ls /home/tb/lustelzl/Projects/kinetics/alanine-dipeptide/run1/out-st2-multi-temp/phi_psi/*_p.txt | tail')


# In[94]:

md1_d = {}
for i, t in enumerate(tmp_l):
    md1_d[i] = np.genfromtxt("md_1/ala1_md_st1_{:.2f}_f.xtc_phi_psi_p.txt".format(
            t))

md2_d = {}
for i, t in enumerate(tmp_l):
    md2_d[i] = np.genfromtxt(
        "md_2/ala1_md_st1_{:.2f}_f.xtc_phi_psi_p.txt".format(
            t))

md3_d = {}
for i, t in enumerate(tmp_l):
    md3_d[i] = np.genfromtxt(
        "md_3/ala1_md_st3_{:.2f}_f.xtc_phi_psi_p.txt".format(
            t))


# In[95]:

d = {'md_st1' : md1_d, 'md_st2' : md2_d, 'md_st3' : md3_d}


# In[96]:

st_l =  ['md_st1', 'md_st2', 'md_st3']


# In[97]:

corr_out_d = {}

bl_l = []

for i, st_name in enumerate(st_l):
    print i
    md_d = d[st_name]
    _rep_fn ="md-out/st{}/rep_md_".format(i+1)+st_name+'_corr.npy'
    if os.path.exists(_rep_fn):
        rep_md = np.load(_rep_fn)
    else:
        rep_md = gen_rep_ar_md(tmp_l, time_pts= md_d[0][:,0])
        np.save(_rep_fn,rep_md)
    _run_name = "md-out/st{}/".format(i+1)+st_name+'_corr'
    print _run_name
    md_st_df, md_tp_df, md_p, md_r =  get_rates_pop_multi_temp(md_d, rep_md,
                                                               _run_name,
                                                               state_def=s, remd=False,
                                                               recalc=False,
                                                               split_state_tp=True,
                                                               stop_tp_at_next_end=True, verbose=False)
    
   
    md_p_5bl, md_r_5b =  block_rates_pop_m(md_st_df, rep_md, md_tp_df, s, 5)
    
    md_r_5b.dropna(inplace=True)
    md_r_5b["rate"] = md_r_5b["rate"].astype(float)   
    # don't average here concat the blocks
    bl_l.append(md_r_5b)
    
    corr_out_d[st_name] =  md_st_df, md_tp_df, md_p, md_r
    
md_df_corr = pd.concat(bl_l)


# need to find combined rate for the three set of MD simulations

# In[98]:

corr_out_d.keys()


# In[99]:

corr_out_d['md_st1'].__len__() # md_st_df, md_tp_df, md_p


# In[100]:

pt_l_md = [ corr_out_d[st][2] for st  in st_l]


# In[101]:

tp_l_md = [ corr_out_d[st][1] for st in st_l]


# In[102]:

total_time = 0
for k, v in d.items():
    print v[0].__len__()
    total_time += v[0].__len__()
    print  v[0][-1,0] - v[0][0,0]
total_time


# In[103]:

ala_ttime = total_time - 3


# In[104]:

trans_l = possible_transitions(4)


# In[105]:

mdc_r, mdc_pt, mdc_tp = combine_rates_stages_simple(pt_l_md, tp_l_md, ala_ttime, 4, trans_l, return_all=True)


# In[106]:

# now symmetrise and calculate erros


# In[107]:

# quick comparison
mdc_r.head()


# In[108]:

rates.head()


# In[ ]:



