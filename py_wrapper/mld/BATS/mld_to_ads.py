#!/usr/bin/env python3
"""
mld_to_ads2.py -- find_mld output  ->  ADS2 file for timeseries_data
                  (MLD written as a pseudo-tracer, two identical levels per station)

usage: mld_to_ads2.py find_mld_output [out.ads2] [options]

  find_mld_output : the full find_mld table (e.g. 1988-2025-MLD), which has
                    lat/lon.  mld_out.set has no position, so use the full file.

Each station gets TWO levels (default 5 m and 10 m) carrying the same MLD, because
timeseries_data wants >= 2 depths per station.  At the timeseries_data prompts use
   shallower boundary : 0     deeper boundary : 11
so the depth average of the two levels is the MLD itself.  The depths are labels,
NOT the mixed layer depth.

Default is one station per DATE (mean MLD of that day's casts), which keeps the
station count (~2000) under the default nstamax of 3000 in timeseries_data.
Use --per-cast for one station per cast (needs nstamax > number of casts).

Layout copied from ctd_csv_to_adsn.py (default 12-tracer build).  The MLD goes in a
spare slot (default 12 = 'phosphate', F6.2, holds up to 999.99).  All other slots are
written as -99.99 with flag -9.
"""
import argparse, os, sys
import numpy as np

BLANK = -99.99
# slot -> (width, decimals) of the default read format in getdata_adsn
SLOTFMT = {1: (6, 2), 2: (7, 3), 3: (7, 2), 4: (7, 2), 5: (7, 2), 6: (7, 2),
           7: (6, 2), 8: (7, 2), 9: (6, 2), 10: (6, 2), 11: (6, 2), 12: (6, 2)}
STOCK_NAMES = ['temperature', 'salinity', 'pressure', 'DIC', 'alkalinity',
               'pCO2 ocean', 'delta 13C', 'oxygen', 'nitrate', 'nitrite',
               'silicate', 'phosphate']
STOCK_UNITS = ['[deg C]', '[psu]', '[dbar]', '[mumol/kg]', '[mueq/kg]', '[ppm]',
               '[per mil]', '[mumol/kg]', '[mumol/kg]', '[mumol/kg]',
               '[mumol/kg]', '[mumol/kg]']


def fmt(v, slot):
    w, d = SLOTFMT[slot]
    s = '%*.*f' % (w, d, v)
    if len(s) > w:
        sys.exit('value %s does not fit slot %d (F%d.%d)' % (v, slot, w, d))
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('infile')
    ap.add_argument('outfile', nargs='?', default=None)
    ap.add_argument('--slot', type=int, default=12, help='tracer slot for MLD (1-12, default 12)')
    ap.add_argument('--depths', type=float, nargs=2, default=[5.0, 10.0],
                    metavar=('Z1', 'Z2'), help='the two label depths (default 5 10)')
    ap.add_argument('--per-cast', action='store_true', help='one station per cast instead of per date')
    ap.add_argument('--bats-box', action='store_true',
                    help='keep only casts inside 31-33N, 63-65.5W (drops Hydrostation S / transects)')
    a = ap.parse_args()
    if not 1 <= a.slot <= 12:
        sys.exit('slot must be 1-12')
    if a.outfile is None:
        a.outfile = os.path.splitext(a.infile)[0] + '.ads2'
    if os.path.abspath(a.outfile) == os.path.abspath(a.infile):
        sys.exit('output would overwrite input')

    # ---- read find_mld table: date decdate julday sta# lat lon cast mld ...
    recs, bad = [], 0
    for ln in open(a.infile):
        p = ln.split()
        if len(p) < 8 or not (p[0].isdigit() and len(p[0]) == 6):
            continue
        try:
            lat, lon, mld = float(p[4]), float(p[5]), float(p[7])
        except ValueError:
            bad += 1
            continue
        if mld <= -99 or not np.isfinite(mld):
            bad += 1
            continue
        recs.append((p[0], lat, lon, mld))
    if a.bats_box:
        n0 = len(recs)
        recs = [r for r in recs if 31 <= r[1] <= 33 and -65.5 <= r[2] <= -63]
        print('box filter dropped %d casts' % (n0 - len(recs)))
    if not recs:
        sys.exit('no usable MLD records found')
    print('read %d casts (%d skipped: failed or unreadable)' % (len(recs), bad))

    # ---- stations (chronological; input is already in date order, but sort to be safe)
    recs.sort(key=lambda r: (('19' if int(r[0][:2]) > 50 else '20') + r[0]))
    if a.per_cast:
        stations = [(d, la, lo, m, 1) for d, la, lo, m in recs]
    else:
        stations, i = [], 0
        while i < len(recs):
            j = i
            while j < len(recs) and recs[j][0] == recs[i][0]:
                j += 1
            blk = recs[i:j]
            stations.append((recs[i][0], np.mean([r[1] for r in blk]),
                             np.mean([r[2] for r in blk]), np.mean([r[3] for r in blk]), len(blk)))
            i = j
    if len(stations) > 9999:
        sys.exit('station number does not fit I4')
    if len(stations) > 3000:
        print('WARNING: %d stations > nstamax (3000 in the default timeseries_data build). '
              'Use the default per-date mode or raise nstamax.' % len(stations))

    # ---- write
    s = a.slot
    names = {k: '' for k in range(1, 13)}
    units = {k: '' for k in range(1, 13)}
    names[s], units[s] = 'MLD', '[m]'
    z1, z2 = sorted(a.depths)
    lines = ['"header", cruise, station number, number of levels, date, latitude, longitude, '
             'bottom depth, minimum depth, maximum depth, ship, data source (find_mld MLD, converted by mld_to_ads2.py)',
             ' bottle, cast, depth,' + ','.join('%-12s' % (names[k] or 'dummy') for k in range(1, 13)) + ', flags',
             ' , ,[m],' + ','.join('%-12s' % (units[k] or '[]') for k in range(1, 13)) + ',']
    for n, (d, la, lo, mld, ncast) in enumerate(stations, 1):
        lines.append('%-9s %-5s %4d %3d %6s %7.2f %7.2f %8.1f %8.1f %8.1f %-5s %-2s' % (
            '"header"', '1', n, 2, d, la, lo, z2, z1, z2, 'MLD', 'bc'))
        for j, z in enumerate((z1, z2), 1):
            vals, flags = [], []
            for k in range(1, 13):
                if k == s:
                    vals.append(fmt(mld, k)); flags.append(2)
                else:
                    vals.append(fmt(BLANK, k)); flags.append(-9)
            lines.append('%3d %3d %5.1f %s %s' % (j, 1, z, ' '.join(vals), ''.join('%2d' % f for f in flags)))
    with open(a.outfile, 'w', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')

    print('wrote %s: %d stations x 2 levels' % (a.outfile, len(stations)))
    print('\nAt the timeseries_data prompts:')
    print('  tracer              : %d' % s)
    print('  flaglevel           : 0')
    print('  set interval        : shallower 0, deeper 11, fileindx 1')
    print('  data-based std dev  : n   (the two levels are identical, so their std is 0)')
    if names[s] != STOCK_NAMES[s - 1]:
        print("\nOptional: in incl/tracer.h rename slot %d to traname '%-12s' traunits '%-10s'"
              % (s, 'MLD', '[m]'))
        print("  (was '%s' / '%s'); only the plot label and output names change." % (STOCK_NAMES[s - 1], STOCK_UNITS[s - 1]))


if __name__ == '__main__':
    main()