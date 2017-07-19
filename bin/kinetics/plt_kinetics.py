__author__ = 'lustelzl'

from ala_kinetics import *


def plot_rama(phi_psi_ar):

    g = sns.jointplot(x=phi_psi_ar[:,1], y= phi_psi_ar[:,2],
                      kind="hexbin", xlim=(-180,180),  ylim=(-180,180))
    g.set_axis_labels("$\phi \;(deg)$", "$\psi \;(deg)$")
    ax0 = g.fig.axes[0]

    ax0.plot([-180.0, 180], [75]*2, "--", c="black")
    ax0.plot([-180.0, 180], [-100]*2, "--", c="black")
    ax0.plot([0,0], [-180, 180], "--", c="black")
    ax0.plot([-180,-180], [-180, 180], "--", c="black")
    f = g.fig
    f.tight_layout()
    return f, ax0


def plot_all_rama(trj_d):
    for trj_nu, phi_psi_ar in trj_d.items():
        f, ax = plot_rama(phi_psi_ar)
        f.savefig("psi_phi_{}.png".format(trj_nu))


def ck_bl_av_trans(av_gr_df, errors, rates, xtick_labels, fig_ax=np.array([]),
                   label="5x30 ns {}", fig_name=None, number_states=4):
    if np.size(fig_ax) == 0:
       fig ,ax = plt.subplots(4,3, figsize = (15, 12))
    else:
         fig, ax = fig_ax
    a = ax.flat

    poss_trans = possible_transitions(number_states)

    for tr_i, trans in enumerate(poss_trans):
        if trans in av_gr_df.keys():   
           _r = rates[rates.type == trans]
           if np.all(_r.events > 0):    
              a[tr_i].errorbar(av_gr_df[trans].index, av_gr_df[trans].values,yerr=errors[trans], fmt="o",
                               label=label)    
              a[tr_i].plot(_r.temperature, _r.rate, "s", label="150 ns {}".format(trans) )
              a[tr_i].legend(loc=2)
              a[tr_i].set_xticklabels(xtick_labels)
              a[tr_i].set_ylabel("Rate ($\mathrm{s^{-1}}$)")
    
    if fig_name: 
       fig.savefig(fig_name)
    else:
         return fig, ax 
         

def _plot_prev_calc(rate_df, rate_ar, fig_ax=np.array([])):
    """
    comparing new (Pandas) and old NumPy calculations
    """
    if np.size(fig_ax) == 0:
       fig ,ax = plt.subplots(figsize = (6, 6))
    else:
         fig, ax = fig_ax
    cl = sns.color_palette()    
       
    plt.plot(rate_df[rate_df['type'] == (0,1)].temperature,
             rate_df[rate_df['type'] == (0,1)].rate/ 1000.0, "o")

    for ri, row in enumerate(rate_ar):
        if row[1] > 0:
           plt.plot(ri, row[1], '.', c=cl[1])
    return fig, ax       
           
           
def _corr_plot_prev_curr_log(rate_df, rate_ar, tmp_dir_l, verbose=False):
    fig, ax = plt.subplots(figsize=(5,5))
    ax.loglog()
    plt.plot([10**-6, 10**-1],[10**-6, 10**-1] )
    diff_l = []
     
    for t in rate_df[rate_df['type']== (0,1)].temperature:
        _t = rate_df[(rate_df['type']== (0,1)) & (rate_df['temperature']==t) ]
        t = tmp_dir_l[_t.temperature.values[0]]
        _md = rate_ar[rate_ar[:,0] == t]
        delta = _md[0][1] - _t.rate.values/ 1000.0
        if verbose:
           print t, delta
        diff_l.append(delta)   
        ax.plot(_md[0][1], _t.rate.values/ 1000.0, "o", c="blue", alpha=0.5)
    return fig, ax, diff_l    
    
    
def convert_temp_indeces(tmp_dir_l, t_values):
    """
    Generate list of temperature for rate plots
    """
    t_l = [tmp_dir_l[int(i)] if not np.isnan(i) else np.nan  for i in t_values]
    return np.array(t_l)
    
    
def plt_av_rate(av_gr_df,errors, fig_ax, xtick_labels):
    fig, ax = fig_ax
    a = ax.flat
    poss_trans = possible_transitions(4) 
    
    for tr_i, trans in enumerate(poss_trans):
        if trans in av_gr_df.keys():   
 
              a[tr_i].errorbar(av_gr_df[trans].index, av_gr_df[trans].values,yerr=errors[trans], fmt="o",
                               label="{} 150 ns".format(trans))    
              a[tr_i].legend(loc=2)
              a[tr_i].set_xticklabels(xtick_labels)
              a[tr_i].set_ylabel("Rate ($\mathrm{s^{-1}}$)")
    return fig, ax
    


def plot_transition(ax, t_type, rates_df, label=None):
    _r = rates_df[rates_df.type == t_type]
    if not label:
        label=t_type
    ax.plot(_r.temperature, _r.rate, "o", label=label)
    return ax

def plot_twelve_trans(rates_md, rates_remd, tmp_l, fig_obj=None, fig_fn=None, label1="MD", label2="REMD"):
    if fig_obj:
       fig, ax = fig_obj
    else:
         fig, ax = plt.subplots(4,3, figsize = (15, 12))
    poss_trans = possible_transitions(4)
    for it, trans in enumerate(poss_trans):
        if it < 3:
           plot_transition(ax[0,it], trans, rates_md, label="{} {}".format(label1, trans))    
           plot_transition(ax[0,it], trans, rates_remd, label="{} {}".format(label2, trans))
           ax[0,it].legend(loc=2)
        if it >= 3 and it < 6:
           plot_transition(ax[1, it -3], trans, rates_md, label="{} {}".format(label1, trans))
           plot_transition(ax[1, it -3 ], trans, rates_remd, label="{} {}".format(label2, trans))
           ax[1,it-3].legend(loc=2)
        if it >= 6 and it < 9:
           plot_transition(ax[2, it -6], trans, rates_md, label="{} {}".format(label1, trans))
           plot_transition(ax[2, it -6 ], trans, rates_remd, label="{} {}".format(label2, trans))
           ax[2,it-6].legend(loc=2)
        if it >= 9 and it < 12:
           plot_transition(ax[3, it -9], trans, rates_md, label="{} {}".format(label1, trans))
           plot_transition(ax[3, it -9 ], trans, rates_remd, label="{} {}".format(label2, trans))
           ax[3,it-9].legend(loc=2 )

    for a in ax.flat:
        a.set_ylabel("$\mathrm{Rate}$ ($\mathrm{ns^{-1}}$)")
        a.set_xticks(np.arange(tmp_l.__len__())[::2])
        a.set_xticklabels( tmp_l[::2])
        a.set_xlabel("Temperature (K)")  
    fig.tight_layout() 
    if fig_fn:
       fig.savefig(fig_fn+".png")
       fig.savefig(fig_fn+".pdf")
    return fig, ax
    
    
def plot_kji(r, md_rt, fig_ax = None, clip_err=True, err_is_err_bar=True,
             highlight_trans=None, highlight_colour='red', verbose=False,
             highlight_alpha=1.0, other_pts_alpha=0.7, highlight_only=False):
    # can be simplified a lot by just passing the right selections         
    cl = sns.color_palette()
    count_plotted_pts = 0
    
    if not fig_ax:
        fig, ax = plt.subplots(figsize=(4,4))
    else:
        fig, ax = fig_ax
    
    #sns.set_style("ticks",  {"xtick.major.size": 12, "ytick.major.size": 12})
    
    if clip_err:
        ax.set_xscale('log', nonposx='clip' )
        ax.set_yscale('log', nonposy='clip')
    else:
        ax.set_xscale('log')
        ax.set_yscale('log')
        
    pts = [0, 10**3]        
    plt.plot(pts, pts, '--', c='grey')
    
    poss_trans = possible_transitions(4)
    _high_c = 0
    for tr_i, trans in enumerate(poss_trans):
        if verbose:
           print trans
        
        if highlight_trans and trans in highlight_trans:
           c= highlight_colour
           a = highlight_alpha
        else:
             c=cl[0]
             a =other_pts_alpha
    
        for temp in np.arange(0,12):
            if highlight_trans and trans in highlight_trans and verbose:
               _high_c += 1
            _m = md_rt[(md_rt.temperature == temp) & (md_rt.type == trans)]
            _r = r[(r.temperature== temp) & (r.type == trans)]
            if _m.size> 0 and (_r.size > 0 ) and (np.all(_m.rate > 0)) and (np.all(_r.rate > 0)):
                count_plotted_pts += 1
                if highlight_only and not trans in highlight_trans:
                   pass
                else:
                     ax.plot(_m.rate, _r.rate,  '.', c=c)
                     if err_is_err_bar:
                        ax.errorbar(_m.rate, _r.rate, xerr=[_m['std_m'], _m['std_p']],
                               yerr=[_r['std_m'], _r['std_p']], c=c, alpha=a)
                     else:
                          ax.errorbar(_m.rate, _r.rate, xerr=[np.abs(_m.rate - _m['std_m']), np.abs(_m.rate - _m['std_p'])],
                                      yerr=[np.abs(_r.rate - _r['std_m']), np.abs(_r.rate - _r['std_p'])], c=c, alpha=a)

    ax.set_xlim(10**-2.2, 10**2.2)        
    ax.set_ylim(10**-2.2, 10**2.2)
    
    if verbose:
       print count_plotted_pts, _high_c
    return fig, ax
    
    
def plot_temp_rc(rc_hist_l, fig_obj=None, colour=None, fmt="o"):
    """
    Plot temperature along the reaction coordinate
    """
    if fig_obj:
       fig, ax = fig_obj
    else:
         fig, ax = plt.subplots(figsize = (4.2, 2.5))
    
    merged_dict = {}
    for key, val in rc_hist_l[0].items():
        _temp_data = []
        for data in rc_hist_l:
            _temp_data.append(data[key])
        merged_dict[key] = np.hstack(_temp_data)
        
    for key, value in merged_dict.items():
        #print value
        if colour:
           ax.errorbar(key, np.mean(value), yerr=np.std(value), fmt=fmt, c=colour)
        else:
             ax.errorbar(key, np.mean(value), yerr=np.std(value), fmt=fmt)
    fig.tight_layout()
    return fig, ax
       

def format_hc_rate_plot(fig, ax, temperature_list, font_size=14, rotation=0):
    ax.legend(loc=2)
    ax.set_xticks(np.arange(12))
    ax.set_xticklabels(temperature_list.astype(int), rotation=rotation)
    ax.set_ylabel(r"$\mathregular{k}_{h/c}}$ [$\mathregular{ns^{-1}}$]", fontsize=font_size)
    ax.set_xlabel("$\mathregular{T \, [K]}$", fontsize=font_size)
    ax.set_xlim(-0.5,11.5)

    fig.tight_layout()
    return fig, ax
    

def plot_rate_est_err(ax,rate_df, tra_l, label_l, x_l, cl_l, return_legend_handles=False, fill_symbol=False, symbol="o",
                      **kwargs):
    """
    useful primarily to plot n-rb rates vs (inverse) temperature Tm/T
    """
    legend_l = []
    for ci, tra in enumerate(tra_l):
       _remd_tra = rate_df[rate_df.type==tra]
       if fill_symbol:
          _fill = cl_l[ci]
       else:
            _fill = 'None'
            
       l, = ax.plot(x_l[ci], _remd_tra.rate,symbol,c=cl_l[ci], label=label_l[ci],
                mfc=_fill, mec=cl_l[ci], mew=2)
                
       ax.errorbar(x_l[ci], _remd_tra.rate, yerr=[_remd_tra.err_m,_remd_tra.err_p],
           fmt=None, c=cl_l[ci], capthick=2, ecolor=cl_l[ci], label="_nolegend_",
           **kwargs)
           
       legend_l.append(l)
       
    if return_legend_handles:
       return ax,  legend_l
    else:       
         return ax
         
         
def gen_x_axis_t_rna(rate_df, tmp_dir_l):
    """
    generate x axis (temperature) for rate plots for two state folding.
    """
    xr_uf=np.array(
                  [tmp_dir_l[int(t)]  if pd.notnull(t) else 0 for t in rate_df[rate_df.type==(1,0)].temperature])
    xr_f=np.array(
                 [tmp_dir_l[int(t)]  if pd.notnull(t) else 0 for t in rate_df[rate_df.type==(0,1)].temperature])
    return xr_f, xr_uf
