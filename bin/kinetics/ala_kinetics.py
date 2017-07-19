__author__ = 'lustelzl'

import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import defaultdict
import pandas as pd
import itertools
import string
from collections import Counter
from numpy.testing import assert_array_equal
import cPickle as pickle
from waiting_time import get_bin_centers
#from waiting_time import *



def make_q_dict(path, split_char, first_index=0, last_index=None, step=1, split_pos=-1):
    q_l = glob.glob(path)
    q_dict = {}
    for q_fn in q_l:
        trj = os.path.basename(q_fn).split(split_char)[split_pos] # "_c"
        q_dict[int(float(trj))] = np.genfromtxt(q_fn)[first_index:last_index:step,:]
    return q_dict



def multi_tp(rc_ar, state_def, verbose=False, stop_tp_at_next_end=False):
    """

    :param rc_ar:
    :param state_def:
    :return:
    """
    assert len(rc_ar[0, 1:]) == len(state_def)

    o = np.zeros(rc_ar.shape)
    o[:,0] = rc_ar[:,0]


    for ci, col in enumerate(rc_ar[0, 1:]):
        # assiging to closest core state
        if np.abs(col - state_def[ci][0]) < np.abs(col - state_def[ci][1]):
            o[0, ci+1] = 0.0
        else:
            o[0, ci+1] = 1.0

    if verbose:
        print o[0, :]

    #o, events_d = find_multi_tp(rc_ar, state_def, o, verbose=verbose)
    o = assign_multi_states(rc_ar, o, state_def, verbose=verbose)
    events_d = find_events_multi_rc(o, rc_ar, state_def, verbose=verbose, stop_tp_at_next_end=stop_tp_at_next_end)
    return o, events_d
       
       
def assign_multi_states(rc_ar, out_ar, state_def, verbose=False):
    out_ar[:,0] = rc_ar[:,0]
    prev_row = out_ar[0, 1:]
    for ri, row in enumerate(rc_ar[1:, :]):
        out_l = state_trans_row_rc(row, prev_row, state_def, verbose=verbose)
        out_ar[ri + 1, 1:] = out_l
        prev_row = out_l
    return out_ar  
    
    
def state_trans_row_rc(rc_row, prev_row, state_def, verbose=False):
    """
    Is any RC indicating a transition?
    """
    out_l = []   
    #if verbose:
    #    print rc_row, prev_row, state_def
    
    for ei, e in enumerate(rc_row[1:]):
        if e <= state_def[ei][0]:
            out_l.append(0)
        elif e >= state_def[ei][1]:
            out_l.append(1)
        else:
            out_l.append(prev_row[ei])
    return out_l
    
    
def find_events_multi_rc(out_ar, rc_ar, state_def, verbose=False, stop_tp_at_next_end=False):
    """
    Parameters:
    -----------
    out_ar: array_like
    rc_ar: array_like
    state_def: array_like
    
    Returns:
    --------
    events_time_d: dictionary
    
    """
    events_end_d = find_end_points(out_ar)
    events_time_d = find_start_points(out_ar, rc_ar, events_end_d, state_def,
                                      verbose=verbose, stop_tp_at_next_end=stop_tp_at_next_end)
    return events_time_d
    
    
def find_end_points(out_ar):
    """
    Parameters:
    -----------
    out_ar: array_like
    
    Returns:
    --------
    events_end_d: dictionary
    """ 
    #events_end_d = defaultdict(list)
    events_end_d = {}
    prev_row = out_ar[0, 1:]
    for ri, row in enumerate(out_ar[1:,:]):
        #print row, prev_row
        config = row[1:]
        if not np.allclose(config, prev_row):
            #events_end_d[tuple(prev_row) + tuple(config)].append(row[0])
            events_end_d[row[0]] = tuple(prev_row) + tuple(config)
        prev_row = config    
    return events_end_d


def find_start_points(out_ar, rc_ar, events_end_time_d, state_def, verbose=False,
                      stop_tp_at_next_end=False):
    """
    Loop backward in the RC array starting from the end points of transitions.
    
    The function should be able to deal correctly with trajectories that start on a transition path,
    but an extra modification might still be required. At the moment a transition is only counted
    if two stable states are connected. 
    
    I could assign the start time of the simulation to the (few) events for which no stable starting
    structure is detected.     
    """
    events_time_d = defaultdict(list)
    
    reversed_rc_ar = rc_ar[::-1]
    reversed_out_ar = out_ar[::-1]
    # alternatively use state_trans_row_rc on end_config
    
    for event_end_time, event in events_end_time_d.items():
        event_start = None
        if verbose:
           # print event, event_end_time
           pass
           
        end_index = np.where(reversed_rc_ar[:,0] == event_end_time)       
        end_state = reversed_out_ar[end_index][0][1:]
        
        for row in reversed_rc_ar[end_index[0][0]:,: ]:
            out_l = state_trans_row_rc(row, end_state, state_def)                 
            if not np.allclose(end_state, out_l):
                #stable starting configuration found and we can start to analyse next event
                event_start = row[0]
               # print out_l, event_start
                events_time_d[event].append((event_start, event_end_time))
                break
            elif (stop_tp_at_next_end == True) and row[0] in events_end_time_d.keys():
                 if event[:2] == events_end_time_d[row[0]][2:]:
                 # could check that value of events_end_time_d is first half of event
                 #event_start_i = np.where(np.array(events_end_time_d.keys()) == row[0])[0]
                    #event_start = events_end_time_d.keys()[event_start_i]
                    event_start = row[0]
                    events_time_d[event].append((event_start, event_end_time))
                    if verbose:
                       print event, event_start, event_end_time
                    break  
                 
        if event_start == None:
            if verbose:
               print out_l, rc_ar[0,0]
               print event, event_end_time
            events_time_d[event].append((rc_ar[0,0], event_end_time))
                
    return events_time_d


def state_row(rc_row, prev_row, state_def, on_tp, events_start_time): #Remove
    """
    Is any RC indicating a transition?
    """
    on_tp_l = []
    out_l = []
    for ei, e in enumerate(rc_row[1:]):
        if e <= state_def[ei][0]:
            out_l.append(0)
            on_tp_l.append(False)
        elif e >= state_def[ei][1]:
            out_l.append(1)
            on_tp_l.append(False)
        else:
            out_l.append(prev_row[ei])
            if not on_tp:
                on_tp_l.append(True)
     
    if not events_start_time: 
        if np.any(on_tp_l):         
           events_start_time = rc_row[0]
   # print on_tp_l
    return out_l, np.any(on_tp_l), events_start_time



def transition_path_data_pd(data_ar, start_time, end_time):
    df = pd.DataFrame(data_ar)
    return df[(df.ix[:, 0] >= start_time) & (df.ix[:, 0] <= end_time)]


def multi_pop_pd(state_data_f, state_def, total_frames, state_names=None):
    """
    Will work for two reaction coordinates with
    :param state_data_f:
    :param state_def:
    :return:
    """
    if not state_names:
       state_names = generate_state_names(state_def)
    state_pop = {}
    col_str = string.lowercase[0:len(state_def)]
    
    for i, st in enumerate(state_names):
        #mask = np.zeros(state_data_f.shape, dtype=bool)
        #mask_l = []
        for i, column in enumerate(col_str):
            #print i, column
            #mask[:,i] = np.logical_and(state_data_f[column] == st[i], mask)    
            #mask_l.append() 
            if i == 0:
               mask = state_data_f[column] == st[i]
            else:
                 mask = np.logical_and(mask, state_data_f[column] == st[i])   
        state_pop[st] = np.float64(len(state_data_f[mask])) / total_frames
    return state_pop


def ala4_temp_pop(state_data_f, state_def, rep_ar, total_frames=None, state_names=None):

    if total_frames is None:
       total_frames = np.float64(len(rep_ar[:,1]))

    state_pop_temp_l = []
    for t in rep_ar[0,1:]:
        _dict = multi_pop_pd(state_data_f[state_data_f.temperature == t], state_def, total_frames,
                             state_names=state_names)
        _dict["temperature"] = t
        state_pop_temp_l.append(_dict)
    return pd.DataFrame(state_pop_temp_l)
    

def generate_state_names(state_def):
    return list(itertools.product([0,1], repeat=len(state_def)))
    

def possible_transitions(number_states):
    '''
    number_states: number of states
    
    Enumerate possible transitions. For n=2 we have 4 possible transitions
    '''
    trans_names = list(itertools.product([0,1], repeat=number_states)) 
    trans_names = [x for x in trans_names if not x[:2] == x[2:]]
    return  trans_names
    

def count_events_at_temp(all_tp_temp_df, number_states, trans_names=None):
    if not trans_names:
       trans_names = possible_transitions(number_states)

    col = ["temperature", "type", "events", "sum_weight"]
    tra_df = pd.DataFrame(data=np.zeros((0,len(col))), columns=col )

    for tp in trans_names:
       _transition =  all_tp_temp_df[all_tp_temp_df.type == tp ]
       t = _transition.temperature.mean()
       # should assert that all temperatures in the selection are equal!
       tra_df = tra_df.append({'temperature' : t, "type" : tp, "events": _transition.fraction.size,
                               "sum_weight" : _transition.fraction.sum()}, ignore_index=True)
    return tra_df


def ala4st_temp_rates(state_pop_temp, all_tp_temp_df):
    col = ["temperature", "type","rate", "events", "sum_weight"]
    r_df = pd.DataFrame(data=np.zeros((0,len(col))), columns=col )

    for t in sorted(all_tp_temp_df["temperature"].unique()):
        tra_df = count_events_at_temp(all_tp_temp_df[all_tp_temp_df.temperature == t], 4)
        r_df = pd.concat([r_df, tra_df], axis=0)

    return r_df
    

def temp_frac_tp(temp_col):
    event_temp_prob = {}
    tp_length = len(temp_col)
    for temp, c in Counter(temp_col).items():
        event_temp_prob[temp] = np.float64(c) / np.float64(tp_length)
    return event_temp_prob


def event_weights_multi(trj_event_d, trj_i, rep_ar, df, remd=True, verbose=False):
    for key, event_list in trj_event_d.items():
        #print key
        for _tuple in event_list:
            _s = transition_path_data_pd(rep_ar, _tuple[0], _tuple[1]).values[:, trj_i + 1]

            if remd:
               frac = temp_frac_tp(_s)
               
               if verbose:
                  print frac
               
               for ti, w in frac.items():
                   df = df.append({'temperature': ti, 'type': key, 'traj' : trj_i, 'start': _tuple[0],
                                   'stop': _tuple[1], 'fraction': w}, ignore_index=True)
            else:
            # if we don not re weight,
                df = df.append({'temperature': trj_i, 'type': key, 'traj' : trj_i, 'start': _tuple[0],
                                'stop': _tuple[1], 'fraction': 1.0}, ignore_index=True)
    return df


def run_multi_tp(trj_rc_dict, rep_ar, state_def, return_raw_events=False,
                 remd=True, stop_tp_at_next_end=False, verbose=False):
    """
    Find transitions and can deal with more than a single reaction 
    coordinate.

    Parameters:
    -----------
    trj_rc_dict: dictionary
    rep_ar: array_like
    state_def: array_like    

    Returns:
    --------
    state_df: pandas dataframe
    tp_temp_df: pandas dataframe
    raw_events (optionally): dictionary

    """
    raw_events = {}

    col = ['temperature', 'type', 'traj', 'start', 'stop', 'fraction']
    tp_temp_df = pd.DataFrame(data=np.zeros((0, len(col))), columns=col)

    state_df = setup_state_table(np.zeros((0, len(state_def) + 3) ), state_def)

    for trj_i, rc in trj_rc_dict.items():
        print trj_i

        trj_bin, _events = multi_tp(rc, state_def, stop_tp_at_next_end=stop_tp_at_next_end,
                                    verbose=verbose)

        assert_array_equal(rep_ar[:,0], trj_bin[:,0])

        dat = np.column_stack(([trj_i]*len(trj_bin), trj_bin[:,0],
                               rep_ar[:, trj_i + 1], trj_bin[:, 1:]))

        _df = setup_state_table(dat, state_def)

        state_df = pd.concat([state_df, _df], axis=0)

        raw_events[trj_i] = _events
        tp_temp_df = event_weights_multi(_events, trj_i, rep_ar, tp_temp_df, remd=remd)

    if return_raw_events:
        return state_df, tp_temp_df, raw_events
    else:
         return state_df, tp_temp_df


def setup_state_table(ar, state_def):
    s_df = pd.DataFrame(ar, columns=['trj_i', 'time','temperature']
                        + [string.lowercase[i] for i in range(len(state_def))])
    s_df.set_index(['trj_i', 'time'], inplace='True')
    return s_df
    
    

def rates_ala4_multi_temp(all_tp_temp_df, pt, total_time, number_states,
                          time_unit_factor=1000.0, trans_names=None ):
    col = ["temperature", "type", "events", "sum_weight", "rate"]
    r_df = pd.DataFrame(data=np.zeros((0,len(col))), columns=col )
    
    # set rates to 0 in r_df

    for t in sorted(all_tp_temp_df["temperature"].unique()):
        tra_df = count_events_at_temp(all_tp_temp_df[all_tp_temp_df.temperature == t], number_states,
                                      trans_names=trans_names)
        tra_df["rate"] = None

        # loop through all types of transitions and find starting state
        for typ in tra_df["type"].unique():
           _pt = pt[pt.temperature == t]
           _pt = _pt[typ[:len(typ) /2]]
           _tra =  tra_df[tra_df.type == typ]
           _w_ev =  _tra.events * (np.float64(_tra.sum_weight) / _tra.events)
           r  = _w_ev / (np.float64(_pt) * total_time / time_unit_factor)
       # print _w_ev, r
       # print r.values
           tra_df = tra_df.set_value(np.where(tra_df.type == typ)[0][0], 'rate', r.values[0])
        #tra_df[np.where(tra_df.type==typ)[0][0]]["rate"] = r
        r_df = pd.concat([r_df, tra_df], axis=0)
    return r_df


def get_rates_pop_multi_temp(q_dict, replica_ar,run_name, state_def=[[-100, 50], [0, 150]], return_raw_events=False,
                             remd=True, recalc=False, time_unit_factor=1000.0, split_state_tp=False, verbose=False,
                             stop_tp_at_next_end=False): #verbose not used atm 
    """
    returns:
    state_df, all_tp_temp_df, pt, rates
    events_d optional

    save/ load previous analysis by pickling all key objects
    
    time_unit_factor is set to 10^3 by default to give out rates in ns. (Gromacs uses ps)
    
    Set split_state_tp=True to correctly assign TPs to the reactent and product states.
    """
    tba_fn = run_name+"_tba.pickle"
    tp_fn = run_name+"_tp.pickle"
    if return_raw_events:
       events_fn = run_name+"events.pickle"
       fn_l = [tba_fn, tp_fn, events_fn]
    else:
         fn_l = [tba_fn, tp_fn]

    if all([os.path.exists(fn) for fn in fn_l]) and not recalc:
       state_df = pd.read_pickle(tba_fn)
       all_tp_temp_df = pd.read_pickle(tp_fn)
       if return_raw_events:
          with open(events_fn, "r") as wp:
               events_d = pickle.load(wp)             
               
    else:
         _r_multi_tp = run_multi_tp(q_dict, replica_ar, state_def, return_raw_events=return_raw_events, remd=remd,
                                    stop_tp_at_next_end=stop_tp_at_next_end, verbose=verbose)
         if return_raw_events:
            state_df, all_tp_temp_df, events_d = _r_multi_tp
            with open(events_fn,"w") as wp:
                 pickle.dump(events_d, wp)
         else:
              state_df, all_tp_temp_df = _r_multi_tp
              # correct TBA assign half of TPs to reactents the other half to products.
         if split_state_tp:
            state_df = loop_corr_tba(state_df, all_tp_temp_df)
         state_df.to_pickle(tba_fn)
         all_tp_temp_df.to_pickle(tp_fn)

    pt_fn = run_name + "_pop.pickle"
    if os.path.exists(pt_fn) and not recalc:
       pt = pd.read_pickle(pt_fn)
    else:
         pt = ala4_temp_pop(state_df, state_def, replica_ar)
         pt.to_pickle(pt_fn)

    rate_fn = run_name + "_rates.pickle"
    if os.path.exists(rate_fn) and not recalc:
       print "prev calculated rates {}".format(rate_fn)
       rates = pd.read_pickle(rate_fn)
    else:  
         total_time = replica_ar[-1, 0] - replica_ar[0, 0]
         number_states = len(generate_state_names(state_def))
         rates = rates_ala4_multi_temp(all_tp_temp_df, pt, total_time, number_states, time_unit_factor=time_unit_factor)
         rates.to_pickle(rate_fn)

    if return_raw_events:
       return state_df, all_tp_temp_df, pt, rates, events_d
    else:
         return state_df, all_tp_temp_df, pt, rates   


def block_rates_pop_m(state_df, replica_ar, all_tp_temp_df, state_defin, number_blocks, verbose=False, time_step=1, time_unit_factor=1000.0):
    """
    """
    block_size = replica_ar[:,0].__len__() * time_step  / number_blocks
    blks = np.arange(replica_ar[0,0], replica_ar[-1,0] + block_size, step=block_size) 
    print block_size  
    pt_l = []
    r_l = []
    
    print blks
       
    for bi, block in enumerate(blks[1:]):     
        prev_bl = blks[bi]
        if verbose: 
           print prev_bl, block
        _rep_ar = replica_ar[(replica_ar[:,0] >= prev_bl ) & (replica_ar[:,0] < block)]
        
        #print _rep_ar
        
        #NB we dont give a proper temperature list to ala4_temp_pop; could just columns in rep_ar
        # The function will work correctly as it will take first line from rep_ar to determine a temperature list
        # the temperature won't be in a nice order- as in the full rep_ar, but each entry exists only once and the correct
        # beheviour will occur
        _pt = ala4_temp_pop(state_df.query("time >= {} and time < {}".format(prev_bl,block)),state_defin, _rep_ar)
        pt_l.append(_pt)     
        number_states = len(generate_state_names(state_defin))
        _r =  rates_ala4_multi_temp(all_tp_temp_df[(all_tp_temp_df.start >= prev_bl ) & (all_tp_temp_df.stop < block)],
                                    _pt, _rep_ar[-1,0] - _rep_ar[0,0], number_states, time_unit_factor=time_unit_factor)
       # _r.loc[_r['rate'].isnull(), 'rate'] = 0.0 #  sala.loc[sala['N10'].isnull(), 'N10'] = 'vcv'
       # _r["rate"].replace('None', 0.0, inplace=True) 
        _r = _r.set_value(_r['rate'].isnull(), 'rate', 0.0)                          
        _r.dropna(inplace=True)
        _r["rate"] = _r["rate"].astype(np.float64)
        r_l.append(_r)
        
        if verbose:
           print bi
    
    pt_concat = pd.concat(pt_l)
    pt_gr = pt_concat.groupby("temperature")  
    rt_concat = pd.concat(r_l)   
    rt_concat = rt_concat.dropna()    
    return pt_concat, rt_concat



def gen_rep_ar_md(temp_l, sim_length=150*1000 + 1, time_pts=None):
    """
    
    Generate an array for the replica MD simulations. The columns do not change as
    no exchanges occur during the simulations
    """
    num_rep = len(temp_l)
    fake_rep_ar = [0] +  [x  for x in range(num_rep)] 
    fake_rep_ar = np.array(fake_rep_ar).reshape((1, num_rep + 1))
    rep_md = np.zeros((sim_length, num_rep + 1))
    rep_md[0,:] = fake_rep_ar
    for i, t in enumerate(rep_md[0,:]): 
        rep_md[:, i] = t

    # http://stackoverflow.com/questions/11295609/how-can-i-check-whether-the-numpy-array-is-empty-or-not
    if time_pts is None:
       print "generate time axis" 
       time_pts = np.arange(sim_length)
    rep_md[:,0] = time_pts # rep_md[:,0] = md_dict[0][:,0]
    return rep_md


def state_bin2num(x, state_num_dict=None):
    """
    Simple function to assign a single state index based on serval columns in a dataframe.
    Perhaps I should use something like this to do transition based assignments?

    """
    if state_num_dict is None:
       state_num_dict = {(0,0) : 0, (0,1) : 1, (1,0) : 2, (1,1) : 3 }
    return state_num_dict[(x[0], x[1])]


def av_blks_df(rt_bl, pad_number_blocks=None, return_event_weight=None, verbose=False):
    """
    Input: rt_bl
    output: av_gr_df, errors
    """
    temp_l = rt_bl["temperature"].unique()
    trans_l = rt_bl["type"].unique()
    
    if pad_number_blocks:
       for temp in temp_l:
           for trans in trans_l:
               if verbose:
                  print temp, trans 
               _s_bl = rt_bl[(rt_bl["temperature"] == temp) & (rt_bl['type'] == trans) ]
               if len(_s_bl) < pad_number_blocks:
                  add_d = {'temperature' : temp, 'type' : trans, 'events' : 0.0, 'sum_weight' : 0.0, 'rate' : 0.0  }
                  rt_bl = rt_bl.append(add_d, ignore_index=True)             
    
    av_gr_df = rt_bl["rate"].groupby([rt_bl["temperature"], rt_bl["type"]]).mean().unstack()
    errors = rt_bl["rate"].groupby([rt_bl["temperature"], rt_bl["type"]]).std().unstack()
   
    if return_event_weight:
       sum_weight = rt_bl["sum_weight"].groupby([rt_bl["temperature"], rt_bl["type"]]).sum().unstack() 
       return av_gr_df, errors, sum_weight
    else:
         return av_gr_df, errors
         
         
def delta_T_tps(trans_type, events_d, rep_ar, trj_l=np.arange(12), return_tp_temperatures=False):
    """
    Are trajectories switching temperature during the transition paths?
    Based on _delta_T() in temperature_tp.ipynb in:
    /home/lustelzl/Projects/kinetics/alanine-dipeptide/remd-1/out-st1/demux_1ps/transition_paths/
    """
    delta_t_d = {}
    tp_temp_d = {}
    for trj in trj_l:
        for trans, trans_times in events_d[trj].items():
            if trans == trans_type:
               for ev_times in trans_times:
                   _temp = transition_path_data_pd(rep_ar, ev_times[0], ev_times[1]).values[:, trj +1]
                   
                   if return_tp_temperatures:
                      tp_temp_d[(trj, ev_times)] = _temp
                   
                   delta_t_d[(trj, ev_times)] = _temp[-1] - _temp[0]
    if return_tp_temperatures:
       return delta_t_d, tp_temp_d
    else:
      return delta_t_d


def analyse_num_rex_tp(transition_events,rep_ar, transition_types=[(0,1), (1,0)], trj_l=[0,1], return_dict=False):
    """
    Get number of replica exchanges on forward and backward transition paths.

    Parameters
    ----------
    transition_events : pandas dataframe
    rep_ar : array-like
    transitions: dictionary of transition events
    trj_l : list of trajectories
    return_dict : to return the dictionaries for the reaction

    Returns
    -------
    av_n_rex_tp : average number of replica exchanges on the transition paths
    (rex_f_d, rex_uf_d) : optional dictionaries for forward and backward reaction

    """
    deltaT_f, deltaT_f_detail = delta_T_tps(transition_types[0], transition_events,
                                            rep_ar, return_tp_temperatures=True, trj_l=trj_l)
    deltaT_uf, deltaT_uf_detail = delta_T_tps(transition_types[1], transition_events,
                                              rep_ar, return_tp_temperatures=True, trj_l=trj_l)

    rex_f_d, rex_f_stretch = analyse_rex_on_tp(deltaT_f_detail)
    rex_uf_d, rex_uf_stretch = analyse_rex_on_tp(deltaT_uf_detail)
    rex_tp_val = rex_f_d.values() + rex_uf_d.values()
    av_n_rex_tp = np.mean(rex_tp_val)
    if return_dict:
        return av_n_rex_tp, (rex_f_d, rex_uf_d)
    else:
        return av_n_rex_tp
    

def analyse_rex_on_tp(deltaT_detail):
    rex_d = {}
    stretch_btw_rex_d  = {}
    for k, v in deltaT_detail.items():
        _stretch = [sum(1 for i in g) for j,g in itertools.groupby(v)]
        _stretch_len = len(_stretch)
        # if there are more than 1 stretes we need to subract one from the 
        # number of stretches. The first stretch does not come from REX
        #if _stretch_len > 1:
        _stretch_len = _stretch_len -1          
        rex_d[k] =_stretch_len
        stretch_btw_rex_d[k] = _stretch
    return rex_d, stretch_btw_rex_d        


def temperature_rc_bin(tp_rc_d, tp_temp_d, step=0.025, lower=0.1, upper=0.6, row_index=1):
    """
    tp_rc_temp_d: dictionary with RC values during the TPs
    tp_temp_d: dictionary with temperature for each
               point on TP (also called deltaT_f_st_detail)
    """
    tp_bins = np.arange(lower, upper + step, step=step)
    rc_hist_temp_d = defaultdict(list)

    assert tp_rc_d.keys().__len__() == tp_temp_d.keys().__len__()

    for key, value in tp_rc_d.items():
        for ri, row in enumerate(value):
            # np.argmin returns index of smallest element
            _rc_bin = tp_bins.flat[np.abs(tp_bins - row[row_index]).argmin()]
            _temp_ri = tp_temp_d[key][ri]
            rc_hist_temp_d[_rc_bin].append(_temp_ri)
    return rc_hist_temp_d


def loop_events_tp_rc(trans_type, events_d, rc_dict, trj_l, verbose=False):
    """
    Loop through events and extract the reaction coordinate values for
    them. 
    
    trans_type: tuple specifying the transition
    events_d: dictionary containing the raw events for each trajectory
    rd_dict: dictionary with reaction coordinate calculated for each trajectory
    trj_l: list of trajectories
    verbose: verbosity
    
    returns
    tp_rc_d: dictionary containing the RC during TPs 
    
    """
    # not tested yet!
    event_counter = 0
    tp_rc_d = {}
    for trj in trj_l:
        for trans, trans_times in events_d[trj].items():
            if trans == trans_type: 
                for ev_times in trans_times:
                    event_counter = event_counter + 1
                    tp_rc_d[(trj, ev_times)] = transition_path_data_pd(
                        rc_dict[trj], ev_times[0], ev_times[1]).values
    assert event_counter  ==  tp_rc_d.keys().__len__()
    if verbose:
        print event_counter, tp_rc_d.keys().__len__()
    return tp_rc_d

    
def _time_multi_index_gr_seq(df, time_1, time_2):
    """
    Extract frames with time index greater time_1 and smaller or equal time_2
    """
    _s = (df.index.get_level_values('time') > time_1) & (df.index.get_level_values('time') <= time_2)
    return _s


def split_react_prod_tp(tp_df, tba_state, verbose=False):
    """
    Go through a list of events. Half of the transition paths belongs to the reactent state
    the other half to the product state.
    
    The correction probably has to be applied on a per-trajectory basis.
    """
    # copy state data frame. Only the copy will be modified
    tba = tba_state.copy()
    len_type = np.unique(tp_df.type.values)
    #print len_type
    len_type = len_type[0]
    #print len_type
    
    for ri, row in tp_df.iterrows() :
        _tp_length = int(row.stop - row.start)
        _2nd_half = _tp_length /2
        _2nd_half = row.start + _2nd_half       
        if verbose:
            print _tp_length, _2nd_half, row.start, row.stop   
            
        #_s = (tba.index.get_level_values('time') > _2nd_half) & (tba.index.get_level_values('time') <= row.stop)
        _s = _time_multi_index_gr_seq(tba, _2nd_half, row.stop)
        
        if len(len_type) == 4: # not tested; apply to Ala2
            tba.loc[_s,'a'] = row['type'][2]
            tba.loc[_s,'b'] = row['type'][3]
        else:
            tba.loc[_s, 'a'] = row['type'][1]
        
    return tba
    
    
def loop_corr_tba(sim_tba, tp_df, verbose=False):
    """
    Loop through the trajectories and the apply a correction to the transition based assignment. Half of the 
    transition paths should be assigned to the reactent state, the other half to the product state. 
    """  
    trj_l = sim_tba.index.get_level_values('trj_i').unique()
    #print trj_l
    c = 0
    for traj in trj_l:
        if  traj in tp_df.traj.unique():
            _st =sim_tba[sim_tba.index.get_level_values('trj_i') == traj]    
            if c == 0:
                r_st_tba =  split_react_prod_tp(tp_df[tp_df.traj == traj], _st, verbose=verbose)
            else:
                _st_tba = split_react_prod_tp(tp_df[tp_df.traj == traj], _st, verbose=verbose)
                r_st_tba = r_st_tba.append(_st_tba, ignore_index=False)
            c = c + 1
    return r_st_tba
    
    
# error analysis

def count_err(N):
    """
    \exp(\pm 1 / sqrt(N)
    """
    return np.exp(- 1.0 / np.sqrt(N)), np.exp(1.0 / np.sqrt(N))
    
 
def diff_from_rate_est(rate, uncertainity):
    """
    difference between est +/- uncertainity - for plotting
    """
    return np.absolute(np.float64(rate) - uncertainity) 
    

def err_log_rate(r, weight_name='sum_weight', diff_from_est=False):
    """
    sigma = k^ / sqrt(N)
    """
    c = r.copy()
    # what is dropna doing exactly?
    c.dropna(subset=['rate'],inplace=True)
    #c.fillna(value=0)
    #r['err_p'] = r.apply(lambda row: (np.float64(row['rate'])*np.exp(1.0/row['sum_weight'])) ,
    #    axis=1)
    #r['err_m'] = r.apply(lambda row: (np.float64(row['rate'])*np.exp(-1.0/row['sum_weight'])) ,
    #    axis=1)
    #  np.absolute(np.float64(x['rate']) - error  
    # http://stackoverflow.com/questions/26614465/python-pandas-apply-function-if-a-column-value-is-not-null
    
    # np.exp(1.0/ np.sqrt(x[weight_name]))
    # np.exp(-1.0/ np.sqrt(x[weight_name]))
    c['std_p'] = c.apply(lambda x: np.float64(x['rate'])* count_err(x[weight_name])[1] if pd.notnull(x['rate'])  else 0, axis=1)
    c['std_m'] = c.apply(lambda x: np.float64(x['rate'])* count_err(x[weight_name])[0] if pd.notnull(x['rate'])  else 0, axis=1)
    
    # New Flage diff_from_est ? 
    if diff_from_est:
       # np.absolute(md_r_s_ln[md_r_s_ln.type==(1,0)].rate - md_r_s_ln[md_r_s_ln.type==(1,0)].std_m
       c['err_m'] = c.apply(lambda x: diff_from_rate_est(x['rate'], x['std_m']) if pd.notnull(x['rate']) else 0, axis=1)
       c['err_p'] = c.apply(lambda x: diff_from_rate_est(x['rate'], x['std_p']) if pd.notnull(x['rate']) else 0, axis=1)
    
    return c


def err_maxL_rate(r):
    c = r.copy()
    _err = c.apply(lambda row: np.float64(row['rate']) / row['sum_weight'] if pd.notnull(row['rate']) else 0, axis=1)
    c['err_p'] = _err
    c['err_m'] = _err
    return c
    
    
def read_av_pickle(fn):
    with open(fn, "rb") as p:
        r, err, w = pickle.load(p)
    return r, err, w
    
# combining stages



def combine_stages_df(name_suffix, s, rc_d, rep_ar_d, stage_l=['md_st1', 'md_st2', 'md_st3'], n_blocks=5,
                      recalc=True, remd=True, tmp_l=None, time_step=1.0, verbose=False, time_unit_factor=1000.0,
                      dir_name='st1', return_comb_pt=False, stop_tp_at_next_end=False): 
    """
    Combining Ala2 MD stages, and CG-RNA REMD calculation stages to collect aggregate statistics.
    The function also allows me to repeat the error estimation with different block sizes. 
    
    Parameters
    
    name_suffix: generic name of the stages
    s: state
    md_d: dictionary with simulation time series
    rep_ar_d: dictionary of replica index files
    n_blocks: Number of blocks per calculation stage
    recalc: Set to False to load previous calculation
    remd: bool
    tmp_l: If replica index array need to be generated a list of temperature
           or indeces has to be supplied  
          
    """
   # combined the stages to get aggregated statistics
    st_st_c = []
    tp_st_c = []
    p_st_c = []
    r_st_c = []
    
    # list for block averging
    bl_l = []
    total_blocks = 0
    total_time = 0
    
    for i, st_name in enumerate(stage_l):
        print i
        md_d = rc_d[st_name]
        if rep_ar_d:
           rep_ar = rep_ar_d[st_name]
        else:   
             rep_ar = gen_rep_ar_md(tmp_l, sim_length=len(md_d[0]), time_pts= md_d[0][:,0])
        _run_name = "{}{}/{}_bl{}".format(dir_name,i+1, st_name+name_suffix, n_blocks) #'_29mar16'
        print _run_name
        md_st_df, md_tp_df, md_p, md_r =  get_rates_pop_multi_temp(md_d, rep_ar,
                                                               _run_name, 
                                                               state_def=s, remd=remd,
                                                               recalc=recalc,
                                                               split_state_tp=True,
                                                               verbose=verbose,
                                                               time_unit_factor=time_unit_factor,
                                                               stop_tp_at_next_end=stop_tp_at_next_end)
        # combining statistics across stages 
        st_st_c.append(md_st_df)
        tp_st_c.append(md_tp_df)
        p_st_c.append(md_p)
        r_st_c.append(md_r) 
        
        # total time in block
        total_time = total_time + (md_d[0][-1,0] - md_d[0][0,0] )
                                                                  
        # for block averaging
        if n_blocks:
        
           total_blocks = total_blocks + n_blocks
           md_p_bl, md_r_bl =  block_rates_pop_m(md_st_df, rep_ar, md_tp_df, s, n_blocks, time_step=time_step,
                                              time_unit_factor=time_unit_factor)
           md_r_bl.dropna(inplace=True)
           md_r_bl["rate"] = md_r_bl["rate"].astype(np.float64)
           # don't average here concat the blocks
           bl_l.append(md_r_bl)
    
    if n_blocks:    
       r_bl_df = pd.concat(bl_l)
       avr_r_bl, err_bl, w_bl = av_blks_df(r_bl_df, pad_number_blocks=total_blocks, return_event_weight=True)
        
    # generate/return list of tp_dfs, p_dfs and state_dfs that can then be combined.
    # duplicate frames?
    # rate calculation could be included here
       
    # for rate calculation I have to find the mean population for each temperature 

    _pt = pd.concat(p_st_c).groupby("temperature").mean().reset_index()
    n_states = _pt.ix[:, _pt.columns != 'temperature'].columns.__len__()
    r_df = rates_ala4_multi_temp(pd.concat(tp_st_c), _pt, total_time, n_states, time_unit_factor=time_unit_factor)

    # r_df and err_bl are probably what I will use on a regular basis
    if not n_blocks:
       return r_df, [st_st_c, tp_st_c, p_st_c, r_st_c], _pt 
    elif return_comb_pt:
       return r_df, err_bl,  [avr_r_bl, w_bl], [st_st_c, tp_st_c, p_st_c, r_st_c], _pt
    else:
         return r_df, err_bl,  [avr_r_bl, w_bl], [st_st_c, tp_st_c, p_st_c, r_st_c]

         
       
def _float_with_NaN(val):
    try:
        return np.float(val)
    except:
        return np.nan  


def name_rev_transition(trans_name, _len_half_trans_name):
    """Generate the name of the reverse transition"""
    return trans_name[_len_half_trans_name:] +  trans_name[:_len_half_trans_name]


def sym_counts_calc_rate(r_df, p_df, trans_names, total_time, time_unit_factor=1000, verbose=False):
    """
    for each transition take forward and reverse counts. We loop though through all transitions,
    at all temperatures. 
    """
    col = ["temperature", "type", "rate", "events", "sum_weight", "rev_events", "rev_sum_weight", "sym_weight"]
    r_symc_df = pd.DataFrame(data=np.zeros((0,len(col))), columns=col )
    _temp_ar = p_df.temperature.unique()
    
    for _typ in trans_names:
        _len_half_name = len(_typ) /2
        rev_name = name_rev_transition(_typ, _len_half_name)
        for t in _temp_ar:
            # get transition counts for forward and backward rates
            _rf_t = r_df[(r_df.temperature== t) & (r_df.type == _typ)]
            _rb_t = r_df[(r_df.temperature== t) & (r_df.type == rev_name)] 
            _p_t = p_df[(p_df.temperature == t) ][_typ[:_len_half_name]]
            #print _rb_t.sum_weight.values, _rf_t.sum_weight.values
            #print np.nansum([_float_with_NaN(_rf_t.sum_weight.values),_float_with_NaN(_rb_t.sum_weight.values)]) 
            _sym_w_counts = np.nansum([_float_with_NaN(_rf_t.sum_weight.values),
                                       _float_with_NaN(_rb_t.sum_weight.values)]) 
            _sym_rate = _sym_w_counts / (2.0* np.float64(_p_t)*total_time / time_unit_factor)
            #if not _sym_rate:
            #   _sym_rate = 0
            
            if verbose:
                print _typ
                print _rf_t.sum_weight.values
                print _rb_t.sum_weight.values
                print _sym_w_counts, _sym_rate
            # write out new data frame with additional colums for sym counts
            r_symc_df =  r_symc_df.append({"temperature" : t, "type" : _typ, "rate": _sym_rate,
                             "events" : _float_with_NaN(_rf_t.events.values),
                             "sum_weight" :  _float_with_NaN(_rf_t.sum_weight.values),
                             "rev_events" :  _float_with_NaN(_rb_t.events.values),
                                           "rev_sum_weight" : _float_with_NaN(_rb_t.sum_weight.values),
                             "sym_weight" : _sym_w_counts}, ignore_index=True)
            
    r_symc_df.rev_events.fillna(0, inplace=True)
    r_symc_df.rev_sum_weight.fillna(0, inplace=True)
    return r_symc_df
    
    
def weight_count_tp_prior(temp_col, rc_ar, prior, bins_prior,
                          prior_t=False, verbose=False, return_tp_p = False,
                         neglect_small_weights=1.0/ 200.0):
    """
    The analysis of REMD is complicated if transition paths are long and contain
    intra-well relaxation.     

    For each frame on the TP we estimate how likely such a structure
    is to be found on a transition path. 
    """ 
    event_temp_prob = {}
    
    if prior_t:
        p_tp = prior[int(temp)]
    else:
         p_tp = prior
        
    w_l = []
    t_l = []
    
    # find the p_tp_q for a give point on the transition path
    for i, frame_rc in enumerate(rc_ar):
        prior_index = [np.abs(bins_prior - frame_rc[1]).argmin()]
        _p_tp_at_q = p_tp[prior_index]
        t_frame  = temp_col[i]
        w_l.append(_p_tp_at_q)
        t_l.append(t_frame)
    
    w_sum = np.sum(w_l)
    w_ar = np.array(w_l) / w_sum
    
    t_ar = np.array(t_l)
    for t_i in t_ar:
        event_temp_prob[t_i] = np.sum(w_ar[t_ar == t_i]) / w_sum
        
    norm_prob = {}
    
    if neglect_small_weights:
        sum_w = 0.0
        for t_i, t_p in event_temp_prob.items():
            if t_p > neglect_small_weights:
                sum_w += t_p
    else:
          sum_w = sum(event_temp_prob.values())
    
    for t_i, t_p in event_temp_prob.items():
        if neglect_small_weights:
            if t_p > neglect_small_weights:
                norm_prob[t_i] = t_p / sum_w
        else:
            norm_prob[t_i] = t_p / sum_w


    np.testing.assert_almost_equal(sum(norm_prob.values()), 1.0)

    return norm_prob
    
    
def loop_reweight_events(events_st1, rep_ar, tp_rc_d, p_tp_q_all_not, tp_bc,
                        neglect_small_weights=10**-2, verbose=False, t_l=np.arange(24)):
    col = ['temperature', 'type', 'traj', 'start', 'stop', 'fraction']
    df = pd.DataFrame(data=np.zeros((0, len(col))), columns=col)

    event_counter = 0
    _c = 0

    tp_prob_d = {}

    for trj_i in t_l:
        _ev_trj_i = events_st1[trj_i]
        for key, event_list in _ev_trj_i.items():
            #print key, event_list
            _c = _c +  len(event_list)
            for _tuple in event_list:
                #print trj_i, _tuple
                event_counter = event_counter + 1
                _s = transition_path_data_pd(rep_ar, _tuple[0], _tuple[1]).values[:, trj_i + 1]
            #frac = _temp_frac_tp_prior(_s, p_t_tp)
                rc_for_event = tp_rc_d[(trj_i,_tuple)]
                frac = weight_count_tp_prior(_s, rc_for_event,
                                         p_tp_q_all_not, tp_bc, verbose=verbose,
                                             neglect_small_weights=neglect_small_weights)
                tp_prob_d[(trj_i,_tuple)] = frac
            
                for ti, w in frac.items():
                    if neglect_small_weights:
                        if not w ==  np.inf:
                            df = df.append({'temperature': ti, 'type': key, 'traj' : trj_i, 'start': _tuple[0],
                                'stop': _tuple[1], 'fraction': w}, ignore_index=True)
                    else:
                        df = df.append({'temperature': ti, 'type': key, 'traj' : trj_i, 'start': _tuple[0],
                                'stop': _tuple[1], 'fraction': w}, ignore_index=True)
    return df, tp_prob_d
    
    
def calculate_transition_path_length(tp_df, temperature, trans_type, REMD=False):
    """
    Calculate the length of the transition paths from MD and REMD for a particular
    transition at a given temperature. 
    """
    tp_trans_temp_df = tp_df[(tp_df.temperature == temperature) & (tp_df.type ==trans_type)]
    tp_length = tp_trans_temp_df.stop.values - tp_trans_temp_df.start.values
    if REMD:
       return tp_length, tp_trans_temp_df.fraction
    else:
         return tp_length
         

def calc_ptpq_all_temp(raw_events_d, rc_dict, trans_type, rep_ar, tba_df,
                       trj_l=np.arange(0,24), time_unit_factor=1000.0):
    '''
    Calculate the probability to be on a transition path given the value of the reaction coordinate.
    '''
    trans_rev = trans_type[::-1]
    _, temperature_forward_d = delta_T_tps(trans_type, raw_events_d, rep_ar, trj_l=trj_l, return_tp_temperatures=True)
    _, temperature_backward_d = delta_T_tps(trans_rev, raw_events_d, rep_ar, trj_l=trj_l, return_tp_temperatures=True)
    
    forward_tp_rc_d = loop_events_tp_rc(trans_type, raw_events_d, rc_dict, trj_l)
    backward_tp_rc_d = loop_events_tp_rc(trans_rev, raw_events_d, rc_dict, trj_l)  
    tp_rc_d = forward_tp_rc_d.copy()
    tp_rc_d.update(backward_tp_rc_d)
    # probability distribution of RC given that we are on a transition path
    p_q_tp = np.vstack((tp_rc_d.values()))
    
    # equilibrium distribution at ALL temperatures
    p_eq_all = []
    for key, values in rc_dict.items():
        p_eq_all.append(np.vstack(values[:,1]))
    p_eq_all_ar = np.vstack(p_eq_all)
    
    # binned equilibrium distribution 
    bins = np.arange(0,1, step=0.02)
    tp_c, tp_b = np.histogram(p_q_tp[:,1], bins=bins)
    tp_bc =  get_bin_centers(tp_b)
    # get_bin_centres
    peq_c, peq_b = np.histogram(p_eq_all_ar, bins=bins)
    
    p_tp_q_all_not =  np.float64(tp_c) / np.float64(peq_c) 
    
    # reweighting
    tp_rew, tp_prob_d = loop_reweight_events(raw_events_d, rep_ar, tp_rc_d, p_tp_q_all_not,
                              tp_bc, neglect_small_weights=False, verbose=True, t_l=trj_l)
                              
    # total time; if one analyses only a segment, the selection has be to 
    total_time = rep_ar[-1, 0] - rep_ar[0, 0]
    rew_r = rates_ala4_multi_temp(tp_rew, tba_df, total_time, 2, time_unit_factor=time_unit_factor)
                          
    # return new weights and the new rates
    return tp_rew, rew_r
    
   
def av_sym_rates_stages(list_sym_df_r):
    col=['temperature', 'type', 'rate', 'events',
         'sum_weight', 'rev_events','rev_sum_weights', 'sym_weight']
    sym_rew_r_df = pd.DataFrame(data=np.zeros((0,len(col))), columns=col)
    
    for idx, row in list_sym_df_r[0].iterrows():
    #print row
        _rate_l = []
        _events_l = []
        _sum_weight_l = []
        _rev_events_l = []
        _rev_sum_weights_l = []
        _sym_weight_l = []
        
        for _stage in list_sym_df_r:
            _sel_trans =  _stage[(_stage.temperature == row.temperature) & (_stage.type == row.type)]
            _sel_trans[~np.isfinite(_sel_trans.rate)] == np.nan
            if all(_sel_trans.rate.notnull()) and all(np.isfinite(_sel_trans.rate.values)): 
               _rate_l.extend(_sel_trans.rate.values.astype(np.float))
               _events_l.extend(_sel_trans.events.values.astype(np.float))
               _sum_weight_l.extend(_sel_trans.sum_weight.values.astype(np.float))
               _rev_events_l.extend(_sel_trans.rev_events.values.astype(np.float))
               _rev_sum_weights_l.extend(_sel_trans.rev_sum_weight.values.astype(np.float))
               _sym_weight_l.extend(_sel_trans.sym_weight.values.astype(np.float))

        print _rate_l
        sym_rew_r_df = sym_rew_r_df.append({'temperature' : row.temperature, 'type' : row.type,
                                       'rate' : np.mean(_rate_l) ,
                                       'events': np.sum(_events_l),
                                       'sum_weight' : np.sum(_sum_weight_l),
                                       'rev_events' : np.sum(_rev_events_l),
                                       'rev_sum_weights' : np.sum(_rev_sum_weights_l),
                                       'sym_weight' : np.sum(_sym_weight_l)},
                                       ignore_index=True)
    return sym_rew_r_df # _rev_sum_weight_l
    
    
    
def combine_rates_stages_simple(pt_l, tp_df_l, combined_total_time, n_states, trans_names,
                         return_all=False):
    """
    Parameters
    ----------
    pt_l:
    tp_df_l:
    combined_total_time:
    n_states:
    trans_names:
    return_all:
    
    Returns:
    --------
    r_c
    pt (optional)
    tp_comb_df (optional)
    
    """                     
    _pt = pd.concat(pt_l).groupby("temperature").mean().reset_index()
    tp_comb_df = pd.concat(tp_df_l)
    r_c = rates_ala4_multi_temp(tp_comb_df, _pt, combined_total_time, n_states, trans_names=trans_names)
    if return_all:
        return r_c, _pt, tp_comb_df
    else:
        return r_c    
