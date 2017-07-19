from ala_kinetics import *
from scipy.optimize import curve_fit
import scipy.stats as stats


cl = sns.color_palette()


def get_bin_centers(b):
    return (b[:-1] + b[1:]) *0.5
    
    
def cumulative_hist(data, bins=100, weights=None, density=False):
    #print weights
    if weights is None:
        h, b = np.histogram(data, bins=bins,
                   density=density)
    else:
        h, b = np.histogram(data, bins=bins, weights=weights,
                           density=density)

    bc = get_bin_centers(b)
    cum = np.cumsum(h)
    if density:
        cum_n = cum 
    else:
       # normalization
       cum_n = cum /  np.float64(cum.max())
    
    return bc, cum_n
  

def exp_dist_func(x, l):
    return l*np.exp(-l*x)


def cumulative_dist_func(x,l):
    return 1.0 - np.exp(-l*x)
    

def biexp_cdf(t, w, tau1, tau2):
    return 1 - w*np.exp(-t/tau1) - (1-w)*np.exp(-t/tau2)
    

def fit_waiting_cdf(data, bins=100, weights=None,
                   return_noncumulative=False, return_bin_edges=False,
                   fit_cumulative_cdf=True, density=False):
    """
    fit_func: cumlative_dist_func or exp_dist_func
    """
    
    bc, cum_n  = cumulative_hist(data, bins=bins, weights=weights,
                                 density=density)
    
    if fit_cumulative_cdf:
       popt = curve_fit(cumulative_dist_func, bc, cum_n, p0=0.002)
    else:
         popt = curve_fit(exp_dist_func, bc, cum_n, p0=0.002)   
    
    return_l = [bc, cum_n, popt[0]]
    
    if return_noncumulative:
        return_l.append(h)
    if return_bin_edges:
        return_l.append(b)
    
    return return_l

    
def fit_plot_cdf(ax, wait_T, dist_label='CDF', fit_label="_nolegend_", bins=1000, weights=None,
                 return_h=False, plot_fit=True, color=None, fit_color=None, density=False):
                 
    b, h, l = fit_waiting_cdf(wait_T, bins=bins, weights=weights,
                              density=density)
    if color:
       ax.plot(b, h, lw=2, label=dist_label, c=color)
    else:   
        ax.plot(b, h, lw=2, label=dist_label)
    if plot_fit:
       if fit_color:
          ax.plot(b, cumulative_dist_func(b, l), c=fit_color, label=fit_label)
       else:
            ax.plot(b, cumulative_dist_func(b, l), label=fit_label)

    if return_h:
       return ax, b, l, h
    else:   
         return ax, b, l


def format_cdf(fig, ax, time_factor=1.0, time_units_label="x 1000 MD steps", fontsize=14):
    ax.set_xlabel('Time '+time_units_label, fontsize=fontsize)
    ax.set_xlim(10**1*time_factor, 10**5.5*time_factor)
    ax.set_ylabel("CDF", fontsize=14)
    fig.tight_layout()
    return fig, ax


def plot_cdf_multi_temp(dw, ax, temp_l, bins=100, temp_label_l=None, dist_label_l=None,
                        rate_type=None, rate_factor=1):
                        
    
    for ti, temp in enumerate(temp_l):
        _dw = dw[dw.temperature == temp]
        if temp_label_l:
           t_label=' T={:.0f} '.format(temp_label_l[ti])
        else:
             t_label=' T={:.0f}'.format(temp)
        if dist_label_l:
           dist_label = dist_label_l[ti] + t_label
        else:
             dist_label = t_label   
        ax[ti], b, l = fit_plot_cdf(ax[ti], _dw.wait_T / _dw.weight, bins=bins, dist_label=dist_label)
        ax[ti].legend(loc=8)
        
        if rate_type is not None:
           _r = rate_type[rate_type.temperature == temp].rate.values[0] / rate_factor
           ax[ti].plot(b, cumlative_dist_func(b, _r), '--') # , label='rate coeff'
    return ax


def plot_cdf_tp(dw, ax, temp_l, bins=100, temp_label_l=None, dist_label_l=None,
                rate_type=None, rate_factor=1, plot_fit=True):
    
    for ti, temp in enumerate(temp_l):
        _dw = dw[dw.temperature == temp]
        _tp = _dw.stop - _dw.start
        #print _tp.head()
        if temp_label_l:
           t_label=' T={:.0f} '.format(temp_label_l[ti])
        else:
             t_label=' T={:.0f}'.format(temp)
        if dist_label_l:
           dist_label = dist_label_l[ti] + t_label
        else:
             dist_label = t_label   
        ax[ti], b, l = fit_plot_cdf(ax[ti], _tp, weights=_dw.fraction.values,  bins=bins, dist_label=dist_label,
                                    plot_fit=plot_fit)
        ax[ti].legend(loc=8)
        
        if rate_type is not None:
           _r = rate_type[rate_type.temperature == temp].rate.values[0] / rate_factor
           ax[ti].plot(b, cumlative_dist_func(b, _r), '--')
    return ax


    
def sum_inv_lifetimes_two_state_lagtime(t,k):
    """
    The sum of the inverse lifetimes in a two state systems depends on the
    sum of the rate coefficients. 
    """
    return (1 - np.exp(-k*t)) /t
    
 
def inv_lifetime_two_state_lagtime(t,p,k): 
    return p*(1 - np.exp(-k*t)) /t


def inv_lifetime_three_state_lagtime(t, p, k1, k2, w):
    """
    expected apparent rate coefficient for a three-state kinetic model
    """
    return  p*(1 - w*np.exp(-k1*t) - (1.0 - w)*np.exp(-k2*t)) / t



def plot_compare_lifetime_distributions(ax, bins, dist1, dist2, dist1_label, dist2_label, dist1_cl, dist2_cl,
                                        hw = 2 * 10, hl = 0.04 * 2, arrow=None, arrow_sign=1, verbose=False,
                                        arrow_start=0.22, exclude_empty_waits=True, plot_fit=False,
                                        return_lambda=False, ks_test=False, weights1=None, weights2=None,
                                        arw1_off=0, arw2_off=0):

    """

    :param ax:
    :param bins:
    :param dist1:
    :param dist2:
    :param dist1_label:
    :param dist2_label:
    :param dist1_cl:
    :param dist2_cl:
    :param _hw:
    :param _hl:
    :param arrow: None, mean, median, fit
    :param arrow_sign:
    :param verbose:
    :param plot_fit: False,
    :param return_lambda: False,
    :param ks_test: False
    :return:
    """
    ax.semilogx()
    
    #exclude empty waits here
    if exclude_empty_waits:
       
        m1 = dist1 > 0
        m2 = dist2 > 0
        dist1_ = dist1[m1]
        dist2_ = dist2[m2]
        #weights1 = weights1[m1]
        if np.any(weights2):
           weights2 = weights2[m2]  
    
    else:
        dist1_ = dist1
        dist2_ = dist2

    ax, b_md, l_md, hist_md_dw = fit_plot_cdf(ax, dist1_, plot_fit=plot_fit, color=dist1_cl,
                                             return_h=True, bins=bins, dist_label=dist1_label,
                                             fit_color=dist1_cl, weights=weights1)

    ax, b_remd, l_remd, hist_remd_dw = fit_plot_cdf(ax, dist2_, color=dist2_cl,
                                                   return_h=True, plot_fit=plot_fit, bins=bins,
                                                   dist_label=dist2_label, fit_color=dist2_cl,
                                                   weights=weights2)
    
    if arrow == "mean":
        arrow_md = np.average(dist1_, weights=weights1)
        arrow_remd = np.average(dist2_, weights=weights2)
    elif arrow == "median":
        arrow_md = np.median(dist1_)
        arrow_remd = np.median(dist2_)
    elif arrow == "fit":
        arrow_md = l_md
        arrow_remd = l_remd

    if arrow and verbose:
        print (1.0/arrow_md, 1.0/arrow_remd)

    if arrow == "fit":
        ax.arrow(1.0 / arrow_md, 1 - np.exp(-1), 0, arrow_sign*arrow_start, color=dist1_cl,
                 head_width=hw, head_length=hl,
                 lw=2)
        ax.arrow(1.0 / arrow_remd, 1 - np.exp(-1), 0, arrow_sign*arrow_start, color=dist2_cl,
                 head_width=hw, head_length=hl,lw=2)
    
    elif arrow == "mean":
         ax.arrow(arrow_md, 1 - np.exp(-1) - arw1_off, 0, arrow_sign*arrow_start, color=dist1_cl, head_width=hw, head_length=hl,
                 lw=2, zorder=50)
         ax.arrow(arrow_remd, 1 - np.exp(-1) - arw2_off, 0, arrow_sign*arrow_start, color=dist2_cl, head_width=hw, head_length=hl,
                  lw=2, zorder=50)
    elif arrow == "median":
         ax.arrow(arrow_md, 0.5, 0, arrow_sign*0.22, color=dist1_cl, head_width=hw, head_length=hl, lw=2)
         ax.arrow(arrow_remd, 0.5, 0, arrow_sign*0.22, color=dist2_cl, head_width=hw, head_length=hl, lw=2)
         
    if return_lambda and arrow == "fit":
       return ax, l_md, l_remd 
    elif ks_test:
        l_fit = [l_md, l_remd] 
        for i, d in enumerate([dist1, dist2]):
            _ks = stats.kstest(d, cumulative_dist_func, args=l_fit[i])
            print (i, _ks)
        return ax, _ks
    else:
        return ax

