"""
diagbox_plots.py
================
Matplotlib port of the legacy NCAR-graphics plotting layer of the
Interannual Diagnostic Box Model (IDBM).

Fortran source                    ->  Python
---------------------------------------------------------------------
bkgroundcolor.F / simplecolor.F   ->  gscr(), bkgroundcolor(), simplecolor(), COLOR_TABLE
calc_grid.F                       ->  calc_grid()
plot_results.F                    ->  plot_results()
generate_ncar_plots.F             ->  generate_plots()
generate_mc_ncar_plots.F          ->  generate_mc_plots()

NCAR call                         ->  matplotlib equivalent
---------------------------------------------------------------------
opngks / clsgks / frame           ->  GKSWorkstation (PdfPages; one page per frame)
gsclip(0)                         ->  clip_on=False on every artist
gscr(1, idx, r, g, b)             ->  COLOR_TABLE[idx] = (r, g, b)
gsplci / gstxci                   ->  color= kwarg
gsln(1..5)                        ->  linestyle= (DASH_PATTERNS)
getset / set(vpl,vpr,vpb,vpt,..)  ->  fig.add_axes([left, bottom, width, height]) + set_xlim/ylim
tick4(-12,-8,-12,-8)              ->  tick_params(direction="out")
gridal(mx,nx,my,ny,...)           ->  box + FixedLocator major/minor ticks
curved                            ->  ax.plot
plchhq(x,y,text,size,angle,cntr)  ->  fig.text / ax.text
pcsetc('FC','$') function codes   ->  mathtext  ($S$..$N$ -> ^{..},  $G$ -> Greek)
getdate                           ->  datetime.now()

Array conventions (same as the f2py return values of run_simulation):
all process arrays are (process, time), i.e. D_C[0:7, :nstep], fluxes[0:10, :nstep].
Only the first nstep columns are used; the rest is padding.
"""

from __future__ import annotations

import datetime as _dt
import math
import re
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

# ======================================================================
# Constants from the Fortran parameter statements
# ======================================================================
SEC_DAY = 24.0 * 3600.0
DAY_YEAR = 365.0
G_MOL = 12.01

# Viewport of the data window inside the frame (NCAR "set" call, NDC units)
XPOSL, XPOSR, YPOSB, YPOST = 0.15, 0.95, 0.2, 0.8

# NCAR character heights are fractions of the frame width. FONT_SCALE is a
# single knob to make them comfortable on screen (1.0 = literal NCAR size).
FIG_SIZE = (8.0, 8.0)
FONT_SCALE = 1.4

MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# ======================================================================
# GKS colour table  (bkgroundcolor.F + simplecolor.F)
# ======================================================================
COLOR_TABLE: dict[int, tuple[float, float, float]] = {}


def gscr(idx: int, r: float, g: float, b: float) -> None:
    """NCAR GKS 'set colour representation' (workstation id dropped)."""
    COLOR_TABLE[idx] = (r, g, b)


def bkgroundcolor() -> None:
    """Background white, foreground 1 = black, 2 = grey (bkgroundcolor.F)."""
    gscr(0, 1.0, 1.0, 1.0)
    gscr(1, 0.0, 0.0, 0.0)
    gscr(2, 0.7, 0.7, 0.7)


def simplecolor() -> None:
    """Colour table for indices 2..19 (simplecolor.F). Overrides index 2."""
    gscr(2, 1.0, 0.0, 0.0)      # red
    gscr(3, 0.0, 0.0, 1.0)      # blue
    gscr(4, 0.6, 0.0, 0.8)      # darkviolet
    gscr(5, 1.0, 0.0, 1.0)      # purple
    gscr(6, 0.0, 1.0, 0.0)      # green
    gscr(7, 1.0, 1.0, 0.0)      # yellow
    gscr(8, 0.7, 1.0, 0.2)      # greenyellow
    gscr(9, 0.5, 1.0, 0.0)      # chartreuse
    gscr(10, 0.2, 1.0, 0.5)     # celeste
    gscr(11, 1.0, 0.0, 0.2)     # orangered
    gscr(12, 0.0, 0.9, 1.0)     # aqua
    gscr(13, 0.0, 0.75, 1.0)    # deepskyblue
    gscr(14, 0.25, 0.45, 0.95)  # royalblue
    gscr(15, 0.4, 0.35, 0.8)    # slateblue
    gscr(16, 1.0, 0.85, 0.0)    # gold
    gscr(17, 0.8, 0.8, 1.0)     # lavender
    gscr(18, 0.85, 0.45, 0.8)   # orchid
    gscr(19, 1.0, 0.65, 0.0)    # orange


# GKS line types 1..5 -> matplotlib linestyles
DASH_PATTERNS = {
    1: "-",                       # solid
    2: (0, (6, 3)),               # dashed
    3: (0, (1, 2)),               # dotted
    4: (0, (6, 2, 1, 2)),         # dash-dot
    5: (0, (6, 2, 1, 2, 1, 2)),   # dash-dot-dot
}


# ======================================================================
# Small helpers
# ======================================================================
def _fortran_fmt(fmt: str):
    """'(f7.2)' -> a callable producing '%.2f'-style strings (width ignored)."""
    m = re.fullmatch(r"\(\s*[fF](\d+)\.(\d+)\s*\)", fmt.strip())
    if not m:
        return lambda v: f"{v:g}"
    nd = int(m.group(2))
    return lambda v: f"{v:.{nd}f}"


def _nint(v: float) -> int:
    """Fortran NINT: round half away from zero."""
    return int(math.floor(v + 0.5)) if v >= 0 else -int(math.floor(-v + 0.5))


def _pt(size: float, fig) -> float:
    """Convert an NCAR character height (fraction of frame width) to points."""
    return size * fig.get_figwidth() * 72.0 * FONT_SCALE


def nstep_per_year(ts: float) -> int:
    # Fortran: int(1.d0/ts). round() protects against 399.9999999 -> 399.
    return int(round(1.0 / ts))


# ======================================================================
# GKS workstation: opngks / frame / clsgks
# ======================================================================
class GKSWorkstation:
    """Stands in for opngks/frame/clsgks.  Every frame() call writes one
    page of a multi-page PDF (like the NCAR metafile had one frame per
    plot).  Set show=True to also open interactive windows at close()."""

    def __init__(self, filename: Optional[str] = "diagbox_plots.pdf", show: bool = False):
        self.show = show
        self._pdf = PdfPages(filename) if filename else None
        self._figs = []
        bkgroundcolor()

    def new_frame(self):
        fig = plt.figure(figsize=FIG_SIZE, facecolor=COLOR_TABLE[0])
        return fig

    def frame(self, fig) -> None:
        if self._pdf is not None:
            self._pdf.savefig(fig)
        if self.show:
            self._figs.append(fig)
        else:
            plt.close(fig)

    def close(self) -> None:          # clsgks
        if self._pdf is not None:
            self._pdf.close()
        if self.show:
            plt.show()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


# ======================================================================
# calc_grid.F
# ======================================================================
def calc_grid(x):
    """Nice axis limits and tick counts for automatically scaled axes.

    Returns (axismin, axismax, nmajticks, nminticks).  nmajticks is the
    number of major *intervals*; nminticks the number of minor intervals
    per major interval (as gridal interprets them).
    """
    nmajdef = 10
    x = np.asarray(x, dtype=float).ravel()
    x = x[np.isfinite(x)]
    nvalid = x.size

    if nvalid <= 1:
        if nvalid == 1:
            delta = 0.5 * max(1.0, abs(x[0]))
            return x[0] - delta, x[0] + delta, 5, 4
        return 0.0, 1.0, 5, 4

    xmin, xmax = float(x.min()), float(x.max())

    delta = xmax - xmin
    if delta <= 0.0:
        delta = 0.05 * max(1.0, abs(xmin))
        if delta == 0.0:
            delta = 1.0
        xmin -= delta
        xmax += delta
        delta = xmax - xmin

    dx = delta / nmajdef

    if abs(dx) < 1:
        iz = int(abs(math.log10(abs(dx))) + 1)
    else:
        iz = -int(math.log10(abs(dx)))

    fac = 10.0 ** iz
    dx = float(int(dx * fac + math.copysign(0.5, dx))) / fac   # round to 1 sig. digit
    if abs(dx) < 1.0e-30:
        dx = 1.0

    xl = float(int(xmin / dx)) * dx     # int() truncates toward zero, as in Fortran
    xr = float(int(xmax / dx)) * dx
    if xl > xmin:
        xl -= dx
    if xr < xmax:
        xr += dx

    nmajticks = _nint((xr - xl) / dx)
    nminticks = 5
    return xl, xr, nmajticks, nminticks


# ======================================================================
# plot_results.F
# ======================================================================
def _line_style(ns, colored, auto_style, line_dashpat, line_color):
    """Colour index / dash pattern selection (ns is 0-based here)."""
    if colored:
        colindx = ns + 2 if auto_style else line_color[ns]
    else:
        colindx = 1
    dashpat = (ns % 5) + 1 if auto_style else line_dashpat[ns]
    return COLOR_TABLE[colindx], DASH_PATTERNS[dashpat]


def plot_results(fig, plottype, xsets, ysets,
                 xautoscale, yautoscale, xplot, yplot, nticks, tickformat,
                 legend, legendpos, legendlabel,
                 auto_style, line_dashpat, line_color, colored):
    """Draw one frame's worth of lines, ticks, tick labels and legend.

    xsets, ysets : lists of 1-D arrays (one per set; lengths may differ,
                   replacing xvalplot/yvalplot/nvalplot)
    xplot, yplot : (min, max) of the axes (used if not autoscaled)
    nticks       : [[xmajor, xminor], [ymajor, yminor]]
    tickformat   : (x_fortran_format, y_fortran_format), e.g. ('(f7.0)', '(f7.2)')

    Returns (ax, xplot, yplot, nticks) with autoscaled values filled in.
    """
    nsets = len(xsets)
    xplot = list(xplot) if xplot is not None else [0.0, 1.0]
    yplot = list(yplot) if yplot is not None else [0.0, 1.0]
    nticks = [list(nticks[0]), list(nticks[1])] if nticks is not None else [[5, 4], [5, 4]]

    # ---- autoscale --------------------------------------------------
    if xautoscale:
        xplot[0], xplot[1], nticks[0][0], nticks[0][1] = calc_grid(
            np.concatenate([np.asarray(a, float).ravel() for a in xsets]))
    if yautoscale:
        yplot[0], yplot[1], nticks[1][0], nticks[1][1] = calc_grid(
            np.concatenate([np.asarray(a, float).ravel() for a in ysets]))

    xmin, xmax = xplot
    ymin, ymax = yplot
    xrange_, yrange_ = xmax - xmin, ymax - ymin

    if colored:
        simplecolor()
    black = COLOR_TABLE[1]

    # ---- viewport / user coordinates  (getset + set) ------------------
    ax = fig.add_axes([XPOSL, YPOSB, XPOSR - XPOSL, YPOST - YPOSB])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_facecolor(COLOR_TABLE[0])

    # ---- tick marks (tick4 + gridal): outside ticks, box, no grid ------
    xmaj = np.linspace(xmin, xmax, max(1, nticks[0][0]) + 1)
    ymaj = np.linspace(ymin, ymax, max(1, nticks[1][0]) + 1)

    def minors(maj, nminor):
        if nminor < 2:
            return []
        out = []
        for a, b in zip(maj[:-1], maj[1:]):
            out.extend(np.linspace(a, b, nminor + 1)[1:-1])
        return out

    ax.xaxis.set_major_locator(FixedLocator(xmaj))
    ax.yaxis.set_major_locator(FixedLocator(ymaj))
    ax.xaxis.set_minor_locator(FixedLocator(minors(xmaj, nticks[0][1])))
    ax.yaxis.set_minor_locator(FixedLocator(minors(ymaj, nticks[1][1])))
    ax.tick_params(which="major", direction="out", length=6, colors=black)
    ax.tick_params(which="minor", direction="out", length=3, colors=black)
    for s in ax.spines.values():
        s.set_color(black)

    # ---- tick labels -------------------------------------------------
    xfmt = _fortran_fmt(tickformat[0])
    yfmt = _fortran_fmt(tickformat[1])
    ax.tick_params(axis="both", labelsize=_pt(0.011, fig))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: yfmt(v)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    if plottype == "ts":
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: xfmt(v)))
    elif plottype == "sa":
        # month names centred in each major interval, no numeric labels
        ax.xaxis.set_major_formatter(NullFormatter())
        dist = xrange_ / max(1, nticks[0][0])
        for i in range(min(nticks[0][0], 12)):
            ax.text(xmin + (i + 0.5) * dist, -0.03, MONTH_LABELS[i],
                    transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=_pt(0.011, fig), color=black, clip_on=False)

    # ---- the lines (curved) ---------------------------------------------
    for ns in range(nsets):
        x = np.asarray(xsets[ns], float)
        y = np.asarray(ysets[ns], float)
        good = np.isfinite(x) & np.isfinite(y)          # NaN / Inf filter
        if good.sum() >= 2:                              # need >= 2 valid points
            color, ls = _line_style(ns, colored, auto_style, line_dashpat, line_color)
            ax.plot(x[good], y[good], color=color, linestyle=ls, clip_on=False)

    # ---- legend: drawn by hand in user coordinates, like the original ----
    if legend:
        lx = xmin + legendpos[0] * xrange_
        ly = ymin + legendpos[1] * yrange_
        for ns in range(nsets):
            if ns < len(legendlabel) and legendlabel[ns] != "":
                x0, y0 = lx, ly - ns * 0.05 * yrange_
                x1 = lx + 0.05 * xrange_
                color, ls = _line_style(ns, colored, auto_style, line_dashpat, line_color)
                ax.plot([x0, x1], [y0, y0], color=color, linestyle=ls, clip_on=False)
                ax.text(x1 + 0.05 * xrange_, y0, legendlabel[ns], ha="left", va="center",
                        fontsize=_pt(0.009, fig), color=black, clip_on=False)

    return ax, xplot, yplot, nticks


# ======================================================================
# Frame decoration + "frame" (the block that is repeated ~25x in Fortran)
# ======================================================================
def _annotate(fig, title, subtitle, notes, ylabel, author_tag):
    """Title, subtitle, notes, y-axis label and date stamp in NDC coordinates."""
    black = COLOR_TABLE[1]
    kw = dict(transform=fig.transFigure, color=black)

    fig.text(XPOSL + (XPOSR - XPOSL) / 2.0, YPOST + (1.0 - YPOST) / 1.5, title,
             ha="center", va="center", fontsize=_pt(0.017, fig), **kw)
    xs = XPOSL + 0.25 * (XPOSR - XPOSL)
    fig.text(xs, YPOST + (1.0 - YPOST) / 2.5, subtitle,
             ha="left", va="center", fontsize=_pt(0.012, fig), **kw)
    if notes:
        fig.text(xs, YPOST + (1.0 - YPOST) / 4.0, notes,
                 ha="left", va="center", fontsize=_pt(0.012, fig), **kw)
    fig.text(0.2 * XPOSL, YPOSB + (YPOST - YPOSB) / 2.0, ylabel,
             rotation=90, ha="center", va="center", fontsize=_pt(0.012, fig), **kw)
    stamp = f"{author_tag} {_dt.datetime.now():%d-%b-%Y %H:%M}".strip()
    fig.text(0.1, 0.03, stamp, ha="left", va="center", fontsize=_pt(0.007, fig), **kw)


def _draw_frame(gks, sets, *, plottype, title, subtitle, ylabel, notes,
                legendlabel, legendpos=(0.05, 0.95),
                xautoscale=True, xplot=None, nticks=None,
                tickformat=("(f7.0)", "(f7.2)"),
                auto_style=True, line_dashpat=None, line_color=None,
                colored=True, author_tag="(ng)"):
    """One complete NCAR 'frame': set up -> plot_results -> titles -> frame."""
    fig = gks.new_frame()
    xsets = [s[0] for s in sets]
    ysets = [s[1] for s in sets]
    plot_results(fig, plottype, xsets, ysets,
                 xautoscale, True, xplot, None, nticks, tickformat,
                 True, legendpos, legendlabel,
                 auto_style, line_dashpat, line_color, colored)
    _annotate(fig, title, subtitle, notes, ylabel, author_tag)
    gks.frame(fig)


# ---- axis labels: NCAR function codes translated to mathtext -----------
# Original: 'dC/dt [$G$l$R$mol kg$S$-1$N$ d$S$-1$N$]'
# ($G$ = Greek font, $R$ = Roman, $S$ = superscript, $N$ = normal).
# The units are micromol, so mu is used here.
YL_RATE = r"dC/dt [$\mu$mol kg$^{-1}$ d$^{-1}$]"
YL_RATE_ANOM = r"Anomaly dC/dt [$\mu$mol kg$^{-1}$ d$^{-1}$]"
YL_FLUX = r"fluxes [mol m$^{-2}$ y$^{-1}$]"
YL_FLUX_ANOM = r"Anomaly fluxes [mol m$^{-2}$ y$^{-1}$]"
YL_VEL = r"u [m s$^{-1}$]"
YL_C = r"C [$\mu$mol kg$^{-1}$]"
YL_FLUX_GC = r"fluxes [gC m$^{-2}$ d$^{-1}$]"


def _default_names(prefix, n):
    return [f"{prefix} {i + 1}" for i in range(n)]


# ======================================================================
# generate_ncar_plots.F
# ======================================================================
def generate_plots(tstart, tend, ts, nstep, D_C, fluxes, c_sim, c_obs, u_var,
                   exp_name="", notes="", *,
                   D_C_avg=None, fluxes_avg=None,
                   D_C_anom=None, fluxes_anom=None,
                   D_C_runavg=None, fluxes_runavg=None,
                   D_C_a_runavg=None, fluxes_a_runavg=None,
                   proc_rates_name=None, proc_fluxes_name=None,
                   outfile="diagbox_plots.pdf", show=False, author_tag="(ng)"):
    """Port of generate_ncar_plots.  Any of the optional derived arrays that
    is None is simply skipped (the f2py run_simulation wrapper does not return
    the seasonal-mean / anomaly / running-average arrays).

    All arrays are (process, time); D_C needs >= 7 rows, fluxes >= 5 rows
    (plus D_C_avg etc. with >= 5 rows).  proc_*_name are the names that
    lived in processes.h (7 rate names, 5 flux names).
    """
    D_C = np.asarray(D_C, float)
    fluxes = np.asarray(fluxes, float)
    rn = list(proc_rates_name) if proc_rates_name is not None else _default_names("rate", 7)
    fn = list(proc_fluxes_name) if proc_fluxes_name is not None else _default_names("flux", 5)

    npy = nstep_per_year(ts)
    nhalf = npy // 2
    nstart, nstop = nhalf + 1, nstep - nhalf          # 1-based, as in Fortran
    rsl = slice(nstart - 1, nstop)                    # running-average window

    t = tstart + ts * np.arange(nstep)                # interannual time axis
    t_run = tstart + nhalf * ts + ts * np.arange(nstop - nstart + 1)
    t_sa = (tstart - int(tstart)) * 365.0 + ts * 365.0 * np.arange(npy)

    title_ia = "Interannual Diagn. Box Model (IDBM)"
    title_sa = title_ia + ": Seasonal Mean Cycle"
    title_an = title_ia + ": Anomalies"
    title_ra = title_ia + ": Running avg."
    sub = "Experiment : " + exp_name
    rate = 1.0 / (ts * 365.0)           # per-step change -> per day
    flux = SEC_DAY * DAY_YEAR           # mol m-2 s-1    -> mol m-2 yr-1

    def five(arr, scale, xs, sl=slice(None)):
        return [(xs, np.asarray(arr, float)[i, sl] * scale) for i in range(5)]

    common = dict(notes=notes, author_tag=author_tag)
    sa = dict(plottype="sa", xautoscale=False, xplot=(0.0, 365.0), nticks=[[12, 3], [0, 0]])

    with GKSWorkstation(outfile, show) as gks:
        # ---- seasonal mean rates / fluxes -------------------------------
        if D_C_avg is not None:
            _draw_frame(gks, five(np.asarray(D_C_avg)[:, :npy], rate, t_sa), title=title_sa, subtitle=sub,
                        ylabel=YL_RATE, legendlabel=rn[:5], **sa, **common)
        if fluxes_avg is not None:
            _draw_frame(gks, five(np.asarray(fluxes_avg)[:, :npy], flux, t_sa), title=title_sa, subtitle=sub,
                        ylabel=YL_FLUX, legendlabel=fn[:5], **sa, **common)

        ts_kw = dict(plottype="ts", **common)

        # ---- interannual rates, fluxes, velocity ------------------------
        _draw_frame(gks, five(D_C[:, :nstep], rate, t), title=title_ia, subtitle=sub,
                    ylabel=YL_RATE, legendlabel=rn[:5], **ts_kw)
        _draw_frame(gks, five(fluxes[:, :nstep], flux, t), title=title_ia, subtitle=sub,
                    ylabel=YL_FLUX, legendlabel=fn[:5], **ts_kw)
        _draw_frame(gks, [(t, np.asarray(u_var, float)[:nstep])], title=title_ia, subtitle=sub,
                    ylabel=YL_VEL, legendlabel=["horizontal velocity"], **ts_kw)

        # ---- running averages (defined on the central part only) --------
        if D_C_runavg is not None:
            _draw_frame(gks, five(D_C_runavg, rate, t_run, rsl), title=title_ra, subtitle=sub,
                        ylabel=YL_RATE, legendlabel=rn[:5], **ts_kw)
        if fluxes_runavg is not None:
            _draw_frame(gks, five(fluxes_runavg, flux, t_run, rsl), title=title_ra, subtitle=sub,
                        ylabel=YL_FLUX, legendlabel=fn[:5], **ts_kw)

        # ---- anomalies ---------------------------------------------------
        if D_C_anom is not None:
            _draw_frame(gks, five(D_C_anom, rate, t, slice(0, nstep)), title=title_an,
                        subtitle=sub, ylabel=YL_RATE, legendlabel=rn[:5], **ts_kw)
        if fluxes_anom is not None:
            _draw_frame(gks, five(fluxes_anom, flux, t, slice(0, nstep)), title=title_an,
                        subtitle=sub, ylabel=YL_FLUX, legendlabel=fn[:5], **ts_kw)

        # ---- running averages of the anomalies ---------------------------
        if D_C_a_runavg is not None:
            _draw_frame(gks, five(D_C_a_runavg, rate, t_run, rsl), title=title_ra, subtitle=sub,
                        ylabel=YL_RATE_ANOM, legendlabel=rn[:5], **ts_kw)
        if fluxes_a_runavg is not None:
            _draw_frame(gks, five(fluxes_a_runavg, flux, t_run, rsl), title=title_ra,
                        subtitle=sub, ylabel=YL_FLUX_ANOM, legendlabel=fn[:5], **ts_kw)

        # ---- observed vs simulated rate (D_C rows 6 and 7) -----------------
        _draw_frame(gks, [(t, D_C[5, :nstep] * rate), (t, D_C[6, :nstep] * rate)],
                    title=title_ia, subtitle=sub, ylabel=YL_RATE,
                    legendlabel=[rn[5], rn[6]], **ts_kw)

        # ---- observed vs simulated C -------------------------------------
        _draw_frame(gks, [(t, np.asarray(c_sim, float)[:nstep]),
                          (t, np.asarray(c_obs, float)[:nstep])],
                    title=title_ia, subtitle=sub, ylabel=YL_C,
                    legendlabel=["simulated", "observed"], legendpos=(0.75, 0.95),
                    tickformat=("(f7.0)", "(f7.1)"), **ts_kw)


# ======================================================================
# generate_mc_ncar_plots.F
# ======================================================================
def generate_mc_plots(nrun, tstart, tend, ts, nstep,
                      c_sim_std, c_obs,
                      mean_D_C, mean_fluxes, mean_c_sim,
                      wmean_D_C, wmean_fluxes, wmean_c_sim,
                      sg_D_C, sg_fluxes, sg_c_sim,
                      wsg_D_C, wsg_fluxes, wsg_c_sim,
                      nstep2, rsmpl_step, exp_name="", notes="", *,
                      proc_rates_name=None, proc_fluxes_name=None,
                      outfile="diagbox_mc_plots.pdf", show=False, author_tag="(ng)"):
    """Port of generate_mc_ncar_plots (Monte-Carlo mean +/- 1 sigma).

    Every mean_*/sg_* array is (process, sample) with nstep2 samples taken
    every rsmpl_step model steps.  c_sim_std / c_obs have nstep entries;
    mean_c_sim / sg_c_sim are 1-D.  Frames are produced twice: once with
    the unweighted and once with the weighted statistics.  (nrun, tend and
    the D_C_std / fluxes_std arguments of the Fortran routine are unused
    there and therefore omitted/ignored here.)

    Faithfully preserved quirks of the original, worth a look:
      * rates use mean/ts (per year) but the axis label says d-1;
      * fluxes are multiplied by day_year (per year) but labelled gC m-2 d-1;
      * the C plot adds sg_c_sim[nstep2-1] (the *last* sample) at every time.
    """
    rn = list(proc_rates_name) if proc_rates_name is not None else _default_names("rate", 7)
    fn = list(proc_fluxes_name) if proc_fluxes_name is not None else _default_names("flux", 5)

    t_mc = tstart + rsmpl_step * ts * np.arange(nstep2)
    t = tstart + ts * np.arange(nstep)
    sub_u = "unweighted uncertainties; Exp. : " + exp_name
    sub_w = "weighted uncertainties; Exp. : " + exp_name
    fl = SEC_DAY * DAY_YEAR * G_MOL

    def band(mean, sg, scale, n=5):
        """15 sets: mean, mean+sg, mean-sg for the first n processes."""
        mean = np.asarray(mean, float)[:, :nstep2]
        sg = np.asarray(sg, float)[:, :nstep2]
        return ([mean[i] * scale for i in range(n)]
                + [(mean[i] + sg[i]) * scale for i in range(n)]
                + [(mean[i] - sg[i]) * scale for i in range(n)])

    # auto_style = .false. styles used by the MC plots
    dash15 = [1, 2, 3, 4, 1] + [2] * 10
    col5 = [2, 3, 4, 5, 6]
    col15 = col5 + [col5[(i - 1) % 5] for i in range(6, 16)]

    common = dict(plottype="ts", notes=notes, author_tag=author_tag, auto_style=False)

    with GKSWorkstation(outfile, show) as gks:
        for weighted in (False, True):
            sub = sub_w if weighted else sub_u
            mD, sD = (wmean_D_C, wsg_D_C) if weighted else (mean_D_C, sg_D_C)
            mF, sF = (wmean_fluxes, wsg_fluxes) if weighted else (mean_fluxes, sg_fluxes)

            # ---- rates with uncertainty band -----------------------------
            ys = band(mD, sD, 1.0 / ts)
            _draw_frame(gks, [(t_mc, y) for y in ys],
                        title="Interannual Diagnostic Box Model (IDBM)", subtitle=sub,
                        ylabel=YL_RATE, legendlabel=rn[:5] + [""] * 10,
                        line_dashpat=dash15, line_color=col15, **common)

        for weighted in (False, True):
            sub = sub_w if weighted else sub_u
            mF, sF = (wmean_fluxes, wsg_fluxes) if weighted else (mean_fluxes, sg_fluxes)

            # ---- fluxes with uncertainty band -------------------------------
            ys = band(mF, sF, fl)
            _draw_frame(gks, [(t_mc, y) for y in ys],
                        title="Seasonal Diagnostic Box Model (SDBM)", subtitle=sub,
                        ylabel=YL_FLUX_GC, legendlabel=fn[:5] + [""] * 10,
                        line_dashpat=dash15, line_color=col15, **common)

        for weighted in (False, True):
            sub = sub_w if weighted else sub_u
            mD, sD = (wmean_D_C, wsg_D_C) if weighted else (mean_D_C, sg_D_C)
            mD = np.asarray(mD, float)[:, :nstep2]
            sD = np.asarray(sD, float)[:, :nstep2]
            # ---- "observed" (row 6, with band) and "simulated" (row 7) rate --
            ys = [mD[5] / ts, np.asarray(mean_D_C, float)[6, :nstep2] / ts,
                  (mD[5] + sD[5]) / ts, (mD[5] - sD[5]) / ts]
            _draw_frame(gks, [(t_mc, y) for y in ys],
                        title="Seasonal Diagnostic Box Model (SDBM)", subtitle=sub,
                        ylabel=YL_RATE, legendlabel=[rn[5], rn[6], "", ""],
                        line_dashpat=[1, 1, 2, 2], line_color=[2, 3, 2, 2], **common)

        c_sim_std = np.asarray(c_sim_std, float)
        for weighted in (False, True):
            sub = sub_w if weighted else sub_u
            m = np.asarray(wmean_c_sim if weighted else mean_c_sim, float)[:nstep2]
            sg = np.asarray(wsg_c_sim if weighted else sg_c_sim, float)
            k2 = np.arange(nstep2) * rsmpl_step               # (k-1)*rsmpl_step + 1, 0-based
            up = c_sim_std[k2] + sg[nstep2 - 1]
            lo = c_sim_std[k2] - sg[nstep2 - 1]
            sets = [(t, c_sim_std[:nstep]), (t, np.asarray(c_obs, float)[:nstep]),
                    (t_mc, m), (t_mc, up), (t_mc, lo)]
            _draw_frame(gks, sets,
                        title="Seasonal Diagnostic Box Model (SDBM)", subtitle=sub,
                        ylabel=YL_C,
                        legendlabel=["simulated std", "observed", "simulated mean", "", ""],
                        legendpos=(0.70, 0.95), tickformat=("(f7.0)", "(f7.1)"),
                        line_dashpat=[1, 4, 1, 2, 2], line_color=[3, 2, 4, 4, 4], **common)


# ======================================================================
# Glue for the f2py wrapper
# ======================================================================
RUN_SIMULATION_RETURNS = (
    "u_var d_c d_c13 fluxes c_obs dc13_obs h_obs temp_obs sal_obs salk_obs "
    "fco2_o_obs pco2_a_obs fco2_a_obs pco2_a_dry dc13_a_obs ws_obs dcdz_obs "
    "kz_obs c_sim int_rates int_fluxes nyears chisq rsq closure amplitude "
    "amplitude_obs c_min c_max c_min_obs c_max_obs costfn avg_fract "
    "dc13_prc avg_kex").split()


def unpack_run_simulation(out) -> dict:
    """Name the 35 return values of diagbox.run_simulation (order from __doc__)."""
    if len(out) != len(RUN_SIMULATION_RETURNS):
        raise ValueError(f"expected {len(RUN_SIMULATION_RETURNS)} return values, got {len(out)}")
    return dict(zip(RUN_SIMULATION_RETURNS, out))


def derive_seasonal_products(a, ts, nstep):
    """OPTIONAL stand-in for the Fortran seasonal-mean / anomaly / running-mean
    routines (anomalies.h, runavg.h), which were not part of the upload.

    ASSUMED definitions - check them against the Fortran before trusting them:
      avg    : mean over all complete years of each within-year step
      anom   : a - avg (climatology repeated along time)
      runavg : centred moving mean over one year (valid for steps
               nstep_per_halfyr+1 .. nstep-nstep_per_halfyr, NaN elsewhere)
    Input/outputs are (process, time); avg is (process, steps_per_year).
    """
    a = np.asarray(a, float)[:, :nstep]
    npy = nstep_per_year(ts)
    nyears = nstep // npy
    avg = a[:, :nyears * npy].reshape(a.shape[0], nyears, npy).mean(axis=1)
    anom = a - np.tile(avg, (1, nstep // npy + 1))[:, :nstep]
    runavg = np.full_like(a, np.nan)
    half = npy // 2
    kern = np.ones(npy) / npy
    for i in range(a.shape[0]):
        v = np.convolve(a[i], kern, mode="valid")          # length nstep-npy+1
        runavg[i, half:half + v.size] = v
    return avg, anom, runavg


def plot_run_simulation(out, tstart, tend, ts, nstep, exp_name="", notes="",
                        derive=True, **kwargs):
    """Plot the tuple returned by diagbox.run_simulation.

    With derive=True the seasonal-mean / anomaly / running-average plots are
    built with derive_seasonal_products (see its caveats); with derive=False
    only the plots that need nothing but run_simulation's outputs are made.
    """
    r = unpack_run_simulation(out)
    extra = {}
    if derive:
        d_avg, d_anom, d_run = derive_seasonal_products(r["d_c"], ts, nstep)
        f_avg, f_anom, f_run = derive_seasonal_products(r["fluxes"], ts, nstep)
        d_a_avg, _, d_a_run = derive_seasonal_products(d_anom, ts, nstep)
        f_a_avg, _, f_a_run = derive_seasonal_products(f_anom, ts, nstep)
        extra = dict(D_C_avg=d_avg, fluxes_avg=f_avg, D_C_anom=d_anom, fluxes_anom=f_anom,
                     D_C_runavg=d_run, fluxes_runavg=f_run,
                     D_C_a_runavg=d_a_run, fluxes_a_runavg=f_a_run)
    generate_plots(tstart, tend, ts, nstep, r["d_c"], r["fluxes"], r["c_sim"], r["c_obs"],
                   r["u_var"], exp_name, notes, **extra, **kwargs)