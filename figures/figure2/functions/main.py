# -*- coding: utf-8 -*-
"""
Helper functions for generating Python summaries of MFAST results.

Author: Conor Bacon
Date: 18 October 2021
"""

# --- Import libraries and setup styles ---
import pathlib
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import obspy
import pandas as pd

try:
    plt.style.use("mfast_summary")
except OSError:
    warnings.warn("You have not added the 'mfast_summary.mplstyle' stylesheet to your"
                  " matplotlib config - continuing with matplotlib defaults.")
mpl.rcParams["font.family"] = "Helvetica"


# --- i/o function space ---
def resolve(parent, name):
    """Resolve a wildcard glob and return first result."""
    return list(parent.glob(name))[0]


def read_cluster_analysis(path, event, station, fb=1):
    """
    Read in the files containing cluster analysis information.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the
        origin time).
    station : str
        Name of the station to produce summary plot for.
    fb : int, optional {1, 2, 3}
        Which filter band to use, of the choices determined
        by MFAST.

    Returns
    -------
    results : numpy.ndarray of float, shape(7, n_windows)
        Values of phi/dt for each window.
    clusters : numpy.ndarray of float, shape(2, n_clusters)
        Positions of clusters identified.

    """

    uid = f"{event}.{station}"

    results = pd.read_csv(resolve(path / uid, f"{uid}.*.fb{fb}.clustxy"),
                          delim_whitespace=True, header=None)
    clusters = pd.read_csv(resolve(path / uid, f"{uid}.*.fb{fb}.clusters"),
                           delim_whitespace=True, header=None)    

    return results.values.T, clusters.values.T


def read_energy_grid(path, event, station, fb=1):
    """
    Read in the files containing the energy grid.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the
        origin time).
    station : str
        Name of the station to produce summary plot for.
    fb : int, optional {1, 2, 3}
        Which filter band to use, of the choices determined
        by MFAST.

    Returns
    -------
    energy : numpy.ndarray of float, shape(3, n_phi*n_dt)
        Energy for each (phi, dt) pair.

    """

    uid = f"{event}.{station}"

    return pd.read_csv(resolve(path / uid, f"{uid}.*.fb{fb}.error"),
                       delim_whitespace=True, header=None).values.T


def read_hodograms(path, event, station, fb=1):
    """
    Read in the files containing particle motion information.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the
        origin time).
    station : str
        Name of the station to produce summary plot for.
    fb : int, optional {1, 2, 3}
        Which filter band to use, of the choices determined
        by MFAST.

    Returns
    -------
    original_pm : numpy.ndarray of float, shape(2, n_samples)
        Uncorrected particle motion.
    corrected_pm : numpy.ndarray of float, shape(2, n_samples)
        Corrected particle motion.

    """

    uid = f"{event}.{station}"

    original_pm = pd.read_csv(resolve(path / uid, f"{uid}.*.fb{fb}.pm"),
                              delim_whitespace=True, header=None)
    corrected_pm = pd.read_csv(resolve(path / uid, f"{uid}.*.fb{fb}.pmc"),
                               delim_whitespace=True, header=None)

    return original_pm.values.T, corrected_pm.values.T


def read_waveforms(path, event, station, fb=1):
    """
    Read in the waveform files.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the
        origin time).
    station : str
        Name of the station to produce summary plot for.
    fb : int, optional {1, 2, 3}
        Which filter band to use, of the choices determined
        by MFAST.

    Returns
    -------
    st : obspy.Stream object
       Waveform data in the form of a Stream of Traces. 

    """

    uid = f"{event}.{station}"

    st = obspy.read(str(resolve(path / uid, f"{uid}.*.fb1.e")))
    st += obspy.read(str(resolve(path / uid, f"{uid}.*.fb1.n")))
    st += obspy.read(str(resolve(path / uid, f"{uid}.*.fb1.z")))
    
    return st


def read_rotated_waveforms(path, event, station, sampling_delta, fb=1):
    """
    Read in the files containing rotated waveform information.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the
        origin time).
    station : str
        Name of the station to produce summary plot for.
    fb : int, optional {1, 2, 3}
        Which filter band to use, of the choices determined
        by MFAST.

    Returns
    -------
    times : numpy.ndarray of float, shape(n_samples)
        Timestamps (in seconds relative to 0) for waveform data.
    tc : numpy.ndarray of float, shape(n_samples)
        Corrected transverse component.
    rc : numpy.ndarray of float, shape(n_samples)
        Corrected radial component.
    t : numpy.ndarray of float, shape(n_samples)
        Original transverse component.
    r : numpy.ndarray of float, shape(n_samples)
        Original radial component.

    """

    uid = f"{event}.{station}"
    start_idx, end_idx = int(12.7 / sampling_delta), int(17.3 / sampling_delta) + 1
    # Uncorrected radial/transverse
    r = pd.read_csv(resolve(path / uid, f"{uid}.*.fb1.sro"),
                    delim_whitespace=True, header=None)
    times, r = r.values[start_idx:end_idx, :].T
    t = pd.read_csv(resolve(path / uid, f"{uid}.*.fb1.sto"),
                    delim_whitespace=True, header=None)
    _, t = t.values[start_idx:end_idx, :].T

    # Corrected radial/transverse
    rc = pd.read_csv(resolve(path / uid, f"{uid}.*.fb1.src"),
                     delim_whitespace=True, header=None)
    _, rc = rc.values[start_idx:end_idx, :].T
    tc = pd.read_csv(resolve(path / uid, f"{uid}.*.fb1.stc"),
                     delim_whitespace=True, header=None)
    _, tc = tc.values[start_idx:end_idx, :].T

    return times, tc, rc, t, r


def read_fast_slow(path, event, station, fb=1):
    """
    Read in the files containing the fast/slow information.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the
        origin time).
    station : str
        Name of the station to produce summary plot for.
    fb : int, optional {1, 2, 3}
        Which filter band to use, of the choices determined
        by MFAST.

    Returns
    -------
    fast : numpy.ndarray of float, shape(n_samples)
        Fast component.
    slow : numpy.ndarray of float, shape(n_samples)
        Slow component.
    fast_cor : numpy.ndarray of float, shape(n_samples)
        Corrected fast component.
    slow_cor : numpy.ndarray of float, shape(n_samples)
        Corrected slow component.

    """

    uid = f"{event}.{station}"

    # Uncorrected slow/fast waves
    fast = pd.read_csv(resolve(path / uid, f"{uid}.*.fb1.sf"),
                       delim_whitespace=True, header=None).values.T
    slow = pd.read_csv(resolve(path / uid, f"{uid}.*.fb1.ss"),
                       delim_whitespace=True, header=None).values.T

    # Corrected slow/fast waves
    fast_cor = pd.read_csv(resolve(path / uid, f"{uid}.*.fb1.sfc"),
                       delim_whitespace=True, header=None).values.T
    slow_cor = pd.read_csv(resolve(path / uid, f"{uid}.*.fb1.ssc"),
                       delim_whitespace=True, header=None).values.T

    return fast, slow, fast_cor, slow_cor


# --- Plotting function space ---
def plot_summary(path, event, station, results, output_path, fb=1):
    """
    Utility function bringing together all of the plotting
    methods for each panel.

    Parameters
    ----------
    path : pathlib.Path object
        Directory containing results folders.
    event : str
        Unique identifier for the event (constructed from the
        origin time).
    station : str
        Name of the station to produce summary plot for.
    fb : int, optional {1, 2, 3}
        Which filter band to use, of the choices determined
        by MFAST.

    """

    # --- Create figure and add axes ---
    fig = _build_grid()

    # --- Filtered waveforms ---
    st = read_waveforms(path, event, station)
    _filtered_waveforms(fig.axes[0], st, results)
    sampling_delta = st[0].stats.delta
    print(f"This station has a sampling rate of {st[0].stats.sampling_rate}")
    # --- R/T waveforms ---
    times, tc, rc, t, r = read_rotated_waveforms(path, event, station,
                                                 sampling_delta)

    _rotated_waveforms(fig.axes[1], times, tc, rc, t, r, results)

    # --- Fast/Slow waveforms ---
    fast, slow, fast_cor, slow_cor = read_fast_slow(path, event, station)
    _fast_slow(fig.axes[2], fast, slow, results, "c")
    _fast_slow(fig.axes[4], fast_cor, slow_cor, results, "d")
    fig.axes[4].legend(fontsize=8, loc=3)

    # --- Hodograms ---
    opm, cpm = read_hodograms(path, event, station)
    _hodogram(fig.axes[3], opm[0], opm[1], "e")
    _hodogram(fig.axes[5], cpm[0], cpm[1], "f")

    # --- Cluster analysis ---
    cluster_results, clusters = read_cluster_analysis(path, event, station)
    _cluster_analysis(fig.axes[6:9], cluster_results, clusters)

    # --- Energy grid ---
    dts, phis, vals = read_energy_grid(path, event, station)
    _energy_grid(fig.axes[9], dts, phis, vals, results)

    # --- Summary information ---
    _summary_information(fig.axes[-1], results)

    plt.savefig(output_path / "figure2.pdf", dpi=400, bbox_inches="tight")
    plt.close()

def _build_grid():
    """
    Utility function to construct the grid of panels in the figure.
    
    Axes are:
        0: Filtered waveforms
        1: Parallel/Perpendicular waveforms before and after correction
      2-5: Corrected waveforms (zoomed) and hodograms
      6-8: Cluster analysis
        9: Energy grid
       10: Summary information

    Returns
    -------
    fig : matplotlib.pyplot.Figure object
        A fully-specified figure containing a grid of panels.

    """

    fig = plt.figure(figsize=(15, 8.65))
    grid_specs = {"nrows": 10,
                  "ncols": 18,
                  "wspace": 2.5,
                  "hspace": 2.5}

    for i in [0, 5]:
        spec = GridSpec(**grid_specs).new_subplotspec((i, 0), colspan=6, rowspan=5)
        fig.add_subplot(spec)

    for i in range(4):
        spec = GridSpec(**grid_specs).new_subplotspec((0 + 3*(i % 2), 6 + 3*(i // 2)), colspan=3, rowspan=3)
        fig.add_subplot(spec)

    for i in [6, 8]:
        spec = GridSpec(**grid_specs).new_subplotspec((i, 6), colspan=3, rowspan=2)
        fig.add_subplot(spec)

    spec = GridSpec(**grid_specs).new_subplotspec((6, 9), colspan=3, rowspan=4)
    fig.add_subplot(spec)

    for i in [0, 6]:
        spec = GridSpec(**grid_specs).new_subplotspec((i, 12), colspan=6, rowspan=6)
        fig.add_subplot(spec)

    return fig


def _filtered_waveforms(ax, st, results):
    """
    Plot the filtered waveforms as ZNE components.

    Parameters
    ----------
    ax : matplotlib.pyplot.Axes object
        Axes on which to plot waveforms.
    st : obspy.Stream object
        Stream containing the filtered waveforms.

    """

    # Taper the data
    for tr in st:
        tr.taper(0.15)

    norm = 1.5*max([max(abs(tr.data)) for tr in st])
    for i, comp in enumerate("ENZ"):
        tr = st.select(component=comp)[0]
        ax.plot(tr.times(), tr.data / norm + 2*i, linewidth=0.6, zorder=2)
        if comp != "Z":
            ax.fill_between([results["wbeg"], results["wend"]], y1=-.75 + 2*i, y2=0.75 + 2*i,
                            color="gray", alpha=0.15, zorder=1)

    ax.vlines(15, -1.5, 5.5, linewidth=0.8, linestyle="--", color="black", alpha=0.6)

    ax.set_xlim([st[0].times()[0], st[0].times()[-1]])
    ax.set_ylim([-1.1, 5.1])
    ax.set_yticks(range(0, 5, 2))
    ax.set_yticklabels("ZNE")
    ax.set_xlabel("Seconds")
    ax.set_ylabel("Component")
    ax.text(0.04, 0.95, "(a)", va="center", ha="center", fontweight="bold",
            transform=ax.transAxes)


def _rotated_waveforms(ax, times, tc, rc, t, r, results):
    """
    Plot the rotated waveforms in the ray frame before and after correction.

    Parameters
    ----------
    ax : matplotlib.pyplot.Axes object
        Axes on which to plot waveforms.
    times : numpy.ndarray of float, shape(n_samples)
        Timestamps (in seconds relative to 0) for waveform data.
    tc : numpy.ndarray of float, shape(n_samples)
        Corrected transverse component.
    rc : numpy.ndarray of float, shape(n_samples)
        Corrected radial component.
    t : numpy.ndarray of float, shape(n_samples)
        Original transverse component.
    r : numpy.ndarray of float, shape(n_samples)
        Original radial component.
    results : pandas.DataFrame object
        Contains general results information for the MFAST measurement.
    
    """

    norm = 1.2*max([max(abs(tr)) for tr in [tc, rc]])
    ax.plot(times, tc / norm, linewidth=1, color="k", zorder=2)
    ax.plot(times, rc /  norm + 2, linewidth=1, color="k", zorder=2)
    norm = 1.2*max([max(abs(tr)) for tr in [t, r]])
    ax.plot(times, t / norm + 4, linewidth=1, color="k", zorder=2)
    ax.plot(times, r /  norm + 6, linewidth=1, color="k", zorder=2)
    for i in range(2):
        ax.fill_between([results["wbeg"], results["wend"]], y1=-0.95+2*i, y2=0.95+2*i,
                        color="gray", alpha=0.15, zorder=1)
    ax.set_xlim([times[0], times[-1]])

    ax.vlines(15, -1.2, 7.2, linewidth=0.8, color="k", alpha=0.6, linestyle="--")

    ax.set_ylim([-1.1, 7.1])
    ax.set_yticks(range(0, 7, 2))
    ax.set_yticklabels(["Tcor", "Pcor", "T", "P"])
    ax.set_xlabel("Seconds")
    ax.set_ylabel("Component")
    ax.text(0.04, 0.95, "(b)", va="center", ha="center", fontweight="bold",
            transform=ax.transAxes)


def _fast_slow(ax, fast, slow, results, letter):
    """
    Plot the fast and slow components.
    
    Parameters
    ----------
    ax : matplotlib.pyplot.Axes object
        Axes on which to plot fast and slow waveforms.
    fast : numpy.ndarray of float, shape(n_samples)
        Fast component.
    slow : numpy.ndarray of float, shape(n_samples)
        Slow component.
    results : pandas.DataFrame object
        Contains general results information for the MFAST measurement.
    letter : string
        Panel label.
    
    """

    ax.plot(fast[0], fast[1], linewidth=1, zorder=2, color="k", label="Fast")
    ax.plot(slow[0], slow[1], linewidth=1, zorder=2, color="k", linestyle="--",
            label="Slow")
    ax.fill_between([results["wbeg"], results["wend"]], y1=-1, y2=1,
                    color="gray", alpha=0.15, zorder=1)
    ax.set_xlim([fast[0, 0], fast[0, -1]])
    ax.set_xlabel("Seconds")
    ax.set_yticklabels([])
    ax.text(0.03, 0.94, f"({letter})", va="center", ha="left",
            fontweight="bold", transform=ax.transAxes)


def _hodogram(ax, x, y, letter):
    """
    Plot the fast and slow components.
    
    Parameters
    ----------
    ax : matplotlib.pyplot.Axes object
        Axes on which to plot fast and slow waveforms.
    x : numpy.ndarray of float, shape(n_samples)
        Component 1 for particle motion diagram.
    y : numpy.ndarray of float, shape(n_samples)
        Component 2 for particle motion diagram.
    letter : string
        Panel label.

    """
    
    ax.plot(x, y, linewidth=1, color="k", alpha=0.9)
    ax.set_xlim([-1, 1])
    ax.set_xlabel("Slow comp.")
    ax.set_xticklabels([])
    ax.set_ylim([-1, 1])
    ax.set_ylabel("Fast comp.")
    ax.set_yticklabels([])
    ax.text(0.03, 0.94, f"({letter})", va="center", ha="left",
            fontweight="bold", transform=ax.transAxes)


def _cluster_analysis(axes, results, clusters):
    """
    Plot the cluster analysis results.
    
    Parameters
    ----------
    axes : matplotlib.pyplot.Axes object
        Axes on which to plot cluster analysis results.
    results : pandas.DataFrame object
        Contains general results information for the MFAST measurement.
    results : numpy.ndarray of float, shape(7, n_windows)
        Values of phi/dt for each window.

    """

    windows, _, _, phis, phi_err, dts, dt_err = results

    # Fast orientation as a function of window number
    axes[0].errorbar(windows, phis, yerr=phi_err, ms=1.5, fmt=".",
                     linewidth=0.5, capsize=1.2, mew=0.5)
    axes[0].text(0.03, 0.9, "(g)", va="center", ha="left",
                 fontweight="bold", transform=axes[0].transAxes)

    # Delay time as a function of window number
    axes[1].errorbar(windows, dts, yerr=dt_err, ms=1.5, fmt=".",
                     linewidth=0.5, capsize=1.2, mew=0.5)
    axes[1].text(0.03, 0.9, "(h)", va="center", ha="left",
                 fontweight="bold", transform=axes[1].transAxes)

    # Fast orientation against delay time
    axes[2].scatter(dts, phis, s=1, marker="x")
    axes[2].scatter(clusters[0], clusters[1], marker="o",
                    s=80, edgecolors="red", facecolors="none")
    axes[2].text(0.04, 0.94, "(i)", va="center", ha="left",
                 fontweight="bold", transform=axes[2].transAxes)

    # Tweak y limits/labels for phi
    for i in [0, 2]:
        axes[i].set_ylim([-90, 90])
        axes[i].set_yticks(range(-90, 91, 30))
        axes[i].set_yticklabels(range(-90, 91, 30))
        axes[i].set_ylabel("Fast orientation, °")

    # Tweak limits/labels for dt
    axes[1].set_ylim([0., 0.6])
    axes[1].set_ylabel("\u03B4t, s")
    axes[2].set_xlim([0., 0.6])
    axes[2].set_xlabel("\u03B4t, s")

    # Tweak x limits/labels for windows
    axes[1].set_xlim([windows[0], windows[-1]])
    axes[1].set_xlabel("Window")
    
    # Setup axes sharing
    axes[0].get_shared_x_axes().join(*axes[:1])
    axes[0].set_xticklabels([])


def _energy_grid(ax, dts, phis, vals, results):

    # Parse into X, Y, and Z
    Z = vals.reshape(len(set(phis)), len(set(dts))).T
    diffs = [(v[-1] - v[0]) / (len(v) * 2) for v in [dts, phis]]
    X, Y = np.mgrid[dts[0]-diffs[0]:dts[-1]+diffs[0]:Z.shape[0]*1j,
                    phis[0]-diffs[1]:phis[-1]+diffs[1]:Z.shape[1]*1j]

    ax.pcolormesh(X, Y, Z, edgecolors="face",
                  cmap="inferno_r")
    ax.contour(X, Y, Z, 9, colors="#cccccc")
    dt, phi = results["tlag"], results["fast"]
    ax.scatter(dt, phi, marker="+", s=50,
               label="Best-fitting (\u03C6, \u03B4t)")

    ax.set_ylabel("Fast orientation, °", fontsize=8)
    ax.set_xlabel("Delay time, s", fontsize=8)
    ax.set_xlim([0, 0.3])
    ax.set_ylim([-90, 90])
    ax.set_yticks(range(-90, 91, 15))
    ax.set_yticklabels(range(-90, 91, 15))
    ax.text(0.05, 0.95, "(j)", va="center", ha="center", fontweight="bold",
            transform=ax.transAxes)


def _summary_information(ax, results):
    ax.set_axis_off()

    # Event info
    ax.text(0.5, 0.9, "Event information", ha="center", fontsize=13)
    stat = results.stat
    uid = results["1event"].split(f".{stat}.")[0] + "." + stat
    _position_text(ax, 0.45, 0.8, "Event UID:", uid)
    _position_text(ax, 0.45, 0.71, "Depth:", f"{results.depthkm:5.2f} km")
    dist = results["10dist(ev-stat)"]
    _position_text(ax, 0.45, 0.62, "Distance:", f"{dist:5.2f} km")

    # Results info
    ax.text(0.5, 0.45, "Results", ha="center", fontsize=14)
    _position_text(ax, 0.45, 0.35, "Grade:", f"{results.gradeABCNR}")
    _position_text(ax, 0.45, 0.26, "Fast orientation:",
                   f"{results.fast:3.1f} \u00b1 {results.Dfast:3.1f} °")
    _position_text(ax, 0.45, 0.17, "Delay time:",
                   f"{results.tlag:5.3f} \u00b1 {results.Dtlag:5.3f} s")
    dspol = results["15Dspol"]
    _position_text(ax, 0.45, 0.08, "Source polarisation:",
                   f"{results.spol:5.1f} \u00b1 {dspol:5.1f} °")


def _position_text(ax, x, y, text1, text2):
    ax.text(x, y, text1, ha="right", va="center", fontsize=10)
    ax.text(x + 0.02, y, text2, ha="left", va="center", fontsize=10)
