__author__ = 'lustelzl'

from ala_kinetics import *


def determine_block_size(trans, state_df, rep_ar, all_tp_temp_df, s, block_n, verbose=False, time_step=1):
    """
    Vary block sizes in the rate calculation to determine the optimal size of the blocks. 
    """
    bl_m ={}
    bl_err = {}
    for n in block_n:
        _pt, _rt = block_rates_pop_m(state_df, rep_ar, all_tp_temp_df, s, n, verbose=verbose, time_step=time_step )
        _bl = av_blks_df(_rt[_rt.type == trans], pad_number_blocks=n, return_event_weight=False, verbose=verbose)
        bl_m[n] = _bl[0]
        bl_err[n] = _bl[1]
    return bl_m, bl_err

    
def get_block_mean_err(df_l, key_l, loc, length=150.0*1000):
    return np.array([(length / n, df_l[0][n].loc[loc], df_l[1][n].loc[loc] ) for n in key_l])
