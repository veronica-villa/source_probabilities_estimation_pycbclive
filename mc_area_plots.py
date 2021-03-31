# Functions to make plots related to the source probability calculation.
# By V. Villa-Ortega, March 2021

import numpy as np
from matplotlib import use; use('Agg')
from matplotlib import pyplot as plt
from pycbc.conversions import mass2_from_mchirp_mass1 as m2mcm1
from pycbc import mchirp_area

def insert_args(parser):
    massplot_group = parser.add_argument_group("Arguments for plotting "
                                               "chirp mass contour on m1m2 "
                                               "plane.")
    massplot_group.add_argument('--x-limits', nargs=2, default=[0.0, 45.0])
    massplot_group.add_argument('--y-limits', nargs=2, default=[0.0, 20.0])

def from_cli(args):
    return {'x_limits': args.x_limits, 'y_limits': args.y_limits}

def contour_mass_plot(trig_mc_det, mass_limits, mass_bdary, z, plot_args):
    trig_mc = mchirp_area.src_mass_from_z_det_mass(z["central"], z["delta"],
                                  trig_mc_det["central"], trig_mc_det["delta"])
    mcb = trig_mc[0] + trig_mc[1]
    mcs = trig_mc[0] - trig_mc[1]
    m2_min = mass_limits["min_m2"]
    m1_max = mass_limits["max_m1"]
    ns_max = mass_bdary["ns_max"]
    gap_max = mass_bdary["gap_max"]
    mib = (2.**0.2) * mcb
    mis = (2.**0.2) * mcs

    lim_m1b = min(m1_max, m2mcm1(mcb, m2_min))
    m1b = np.linspace(mib, lim_m1b, num=100)
    m2b = m2mcm1(mcb, m1b)

    lim_m1s = min(m1_max, m2mcm1(mcs, m2_min))
    m1s = np.linspace(mis, lim_m1s, num=100)
    m2s = m2mcm1(mcs, m1s)
   
    fig = plt.figure()

#plot contour
    if mib > m1_max:
        plt.plot((m1_max, m1_max), (m2mcm1(mcs, lim_m1s), m1_max), "b")
    else:
        plt.plot(m1b, m2b, "b")
        plt.plot((m1_max, m1_max), (m2mcm1(mcs, lim_m1s),
                                    m2mcm1(mcb, lim_m1b)),"b")
    if mis >= m2_min:
        plt.plot(m1s, m2s, "b")
        plt.plot((lim_m1s, lim_m1b), (m2_min, m2_min), "b")
    else:
        plt.plot((m2_min, lim_m1b), (m2_min, m2_min), "b")
#plot limits
    plt.plot((m2_min, m1_max), (m2_min, m1_max), "k--")
    plt.plot((ns_max, ns_max), (m2_min, ns_max), "k:")
    plt.plot((gap_max, gap_max), (m2_min, gap_max), "k:")
    plt.plot((ns_max, m1_max), (ns_max, ns_max), "k:")
    plt.plot((gap_max, m1_max), (gap_max, gap_max), "k:")
#colour plot
    plt.fill_between(np.arange(0.0, ns_max-0.01, 0.01), gap_max,
                     m1_max, color=source_color('NSBH'), alpha=0.5)
    plt.fill_between(np.arange(gap_max, m1_max, 0.01), 0.0,
                     ns_max, color=source_color('NSBH'))
    plt.fill_between(np.arange(ns_max, gap_max, 0.01),
                     np.arange(ns_max, gap_max, 0.01), m1_max,
                     color=source_color('Mass Gap'), alpha=0.5)
    plt.fill_between(np.arange(0.0, ns_max, 0.01), ns_max,
                     gap_max, color=source_color('Mass Gap'), alpha=0.5)
    plt.fill_between(np.arange(gap_max, m1_max, 0.01),
                     np.arange(gap_max, m1_max, 0.01), m1_max,
                     color=source_color('BBH'), alpha=0.5)
    plt.fill_between(np.arange(gap_max, m1_max, 0.01),
                     np.arange(gap_max, m1_max, 0.01), gap_max,
                     color=source_color('BBH'))
    plt.fill_between(np.arange(0.0, ns_max, 0.01), 0.0,
                     np.arange(0.0, ns_max, 0.01),
                     color=source_color('BNS'))
    plt.fill_between(np.arange(0.0, ns_max, 0.01), ns_max,
                     np.arange(0.0, ns_max, 0.01),
                     color=source_color('BNS'), alpha=0.5)
    plt.fill_between(np.arange(ns_max, gap_max, 0.01),
                     np.arange(ns_max, gap_max, 0.01),
                     color=source_color('Mass Gap'))
    plt.fill_between(np.arange(gap_max, m1_max, 0.01), ns_max,
                     gap_max, color=source_color('Mass Gap'))
#colour contour
    x1 = np.arange(mis, mib+0.01, 0.01)
    plt.fill_between(x1, x1, m2mcm1(mcs, x1),
                     facecolor=(1,1,1,0.5), edgecolor=(0,0,0,0))
    x2 = np.arange(mib, lim_m1b, 0.01)
    plt.fill_between(x2, m2mcm1(mcb, x2), m2mcm1(mcs, x2),
                     facecolor=(1,1,1,0.5), edgecolor=(0,0,0,0))
#plot_details
    plt.xlim(left=plot_args['x_limits'][0], right=plot_args['x_limits'][1])
    plt.ylim(bottom=plot_args['y_limits'][0], top=plot_args['y_limits'][1])
    plt.xlabel("M1")
    plt.ylabel("M2")

    return fig


def probabilities_plot(probabilities):
    prob_plot = {k: v for (k, v) in probabilities.items() if v != 0.0}
    labels, sizes = zip(*prob_plot.items())
    colors = [source_color(label) for label in labels]
    fig, ax = plt.subplots()
    ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
           textprops={'fontsize': 15})
    ax.axis('equal')
    return fig

_source_color_map = {
    'BNS': '#A2C8F5',   # light blue
    'NSBH': '#FFB482',  # light orange
    'BBH': '#FE9F9B',   # light red
    'Mass Gap': '#8EE5A1',  # light green
    'GNS': '#98D6CB',   # turquoise
    'GG': '#79BB87',    # green
    'BHG': '#C6C29E'    # dark khaki
}

def source_color(source):
    return _source_color_map[source]

