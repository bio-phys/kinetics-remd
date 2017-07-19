from kinetics import *
import seaborn as sns
cl = sns.color_palette()
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def calc_dwell_time(rc_ar, plot_fold_unfold=False):
    out_ar, f_events, u_events = find_transition_paths(rc_ar)

    out_l = []

    dwell_time()

    return out_l


def dwell_time(events_0, events_1):
    """
    Dwell time in state 0: First time in state 1 - first time in state 0

    :param events_0:
    :param events_1:
    :return:
    """
    dw_l = []
    for event_index, event in events_0:
        dw_l.append(events_1[event_index][1] - event[1])
    return np.array(dw_l)


def cumulative_exp(x, k):
    return 1 - np.exp(-k*x)


def fit_cum_exp(dw_times):
    h, bins, patches = plt.hist(dw_times, cumulative=True, bins=100, normed=1, alpha=0.4)
    bin_c = bin_centres(bins)
    popt = curve_fit(cumulative_exp, bin_c, h, p0=0.0001)
    return bin_c, h, popt[0]


def bin_centres(bins):
    return (bins[:-1] + bins[1:]) / 2.0


def plot_dw_times(bin_c, h, rate):

   # plt.plot(bin_c, h, alpha)
    plt.plot(bin_c, cumulative_exp(bin_c, rate))
