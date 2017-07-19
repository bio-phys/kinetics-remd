#from kinetics.waiting_time import *
#from kinetics.waiting_time import *
from waiting_time import *
from ala_kinetics import *


def loop_dwell_trans_temp(tp_df, states_df, transition, tu=1.0, verbose=False):
    """
  
    Loop through alll the temperatures. For each temperature find the unique starting
    points of transitions.
    
    """
    col = ['temperature', 'type', 'traj', 'start', 'stop',
           'weight','wait_T', 'wait', 'prev_state']
    df = pd.DataFrame(data=np.zeros((0, len(col))), columns=col)
    

    #t1 = remd_state[remd_state.temperature == 10].sort_index(level='time')
    for temp_i in np.unique(tp_df.temperature.values):
        _s_t = states_df[states_df.temperature == temp_i].sort_index(level='time')
        # order transitions
        gr0 = tp_df[tp_df.temperature==temp_i].groupby(['start'])
        index = [gp_keys[0] for gp_keys in gr0.groups.values()]
        unique_df = tp_df.reindex(index)
        s = unique_df.sort_values(by='start')
        #s_type = s[s.type==_transition]
        s_type = s[s.type.isin(transition)]
        if s_type[s_type.temperature == temp_i].values.size > 0:
           df = dwell_trans_temp(df, s_type, _s_t, tu=tu, verbose=verbose)
    return df


def dwell_trans_temp(df, s_type, s_t,tu=1, verbose=False):
    
    #gr0 = tp_df[tp_df.temperature==t].groupby(['start'])
    #index = [gp_keys[0] for gp_keys in gr0.groups.values()]
    #unique_df = tp_df.reindex(index)
    #s = unique_df.sort_values(by='start')
    _prev = 0
    #s_type = s[s.type==(0,0,0,1)]
    len_type = np.unique(s_type.type.values)
    #print len_type[0]
    if verbose:
       print len_type
    len_type = len(len_type[0])
    if verbose:
       print len_type

    for ri, row in s_type.iterrows():
        end_dw = row.start
        start_dw = _prev
            
        #_tp = tp_df[tp_df.start ==end_dw]
        _tp = s_type[s_type.start == end_dw]
        _tp_type = np.unique(_tp.type)
        _tp_temp = np.unique(s_t.temperature.values)[0]
        w = _tp[_tp.temperature == _tp_temp].fraction.values
        _tp_traj = _tp.traj.unique()[0].astype(int)
        ar = s_t.ix[(s_t.index.get_level_values('time') >= start_dw) &
                    (s_t.index.get_level_values('time') <= end_dw)]
        if len_type == 4: 
           _dwell_reactent_state = ar[(ar.a == _tp_type[0][0]) & (ar.b == _tp_type[0][1])]
        else:
             _dwell_reactent_state = ar[(ar.a == _tp_type[0][0])]
                   
        total_wait = np.float64(end_dw) - np.float64(start_dw)
        
        #print w[0]
        
        _row_dict = {'temperature': _tp_temp, 'type': _tp_type[0], 'traj' : _tp_traj,
                     'start' : start_dw, 'stop' : end_dw, 'weight': w[0],
                     'wait_T': len(_dwell_reactent_state)*tu, 'wait': total_wait}
        df = df.append(_row_dict, ignore_index=True)
        _prev = row.stop
    return df
    
    
def weight_dwell_trans_temp(remd_dw_df, temperature):
    _remd_dw_df = remd_dw_df[remd_dw_df.temperature == temperature]
    return _remd_dw_df.wait_T / _remd_dw_df.weight, _remd_dw_df.weight
