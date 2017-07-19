import numpy as np

def mean_tau_w(tau_a, w_a):
    return np.average( tau_a / w_a, weights=w_a)


def var_tau_w(tau_a, w_a, av):
    return np.sum( w_a*(tau_a/w_a - av)**2) / (np.sum(w_a) )


def check_moments_w(tau_a, w_a):
    _av = mean_tau_w(tau_a, w_a)
    _var = var_tau_w(tau_a, w_a, _av)
    return _av, _var


def loop_check_moments_w(dt_l, dw_dict, exclude_empty_waits=True):
    var_l = []
    av_l = []
    x_l = []
    for dt in sorted(dt_l):
        dw_df = dw_dict[dt]
        if np.any(dw_df.wait_T.values > 0):
            if exclude_empty_waits:
               m = dw_df.weight > 0
               _dw = dw_df[m]
            else:
                 _dw = dw_df  
        
            _exp = check_moments_w(_dw.wait_T,
                                   _dw.weight)
            var_l.append(_exp[1])
            av_l.append(_exp[0])
            x_l.append(dt)
    return np.column_stack((x_l, av_l, var_l))



