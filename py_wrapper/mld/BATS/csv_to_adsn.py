#!/usr/bin/env python3
"""
csv_to_adsn.py  --  BATS CTD profile CSV  ->  ADS file for timeseries_data
                    (default 12-tracer layout: build WITHOUT -Ddericarbquant)
 
usage:  csv_to_adsn.py input.csv [output.ads2] [options]
 
If the output name is omitted it is the input name with extension .ads2
(bats_ctd.csv -> bats_ctd.ads2).
 
Needs: python3, pandas, numpy.
 
WHAT IS CONVERTED
  Each cast (column ID) becomes one ADS station; each depth becomes one level.
  All measured variables are written, each to a tracer slot of the ADS layout
  (menu number in timeseries_data = slot number):
 
      Temperature -> slot 1     Salinity -> slot 2     Pressure -> slot 3
      Oxygen      -> slot 8     BAC      -> slot 4     Flu      -> slot 5
      PAR         -> slot 6
 
  Slots 4,5,6 are DIC/alkalinity/pCO2 in the stock program, so tracer.h must be
  edited to rename them; the script prints the exact lines at the end.
  Change the assignment with  --map NAME=slot  (slot 1-12, or 'none' to drop).
 
QUALITY FLAGS
  timeseries_data keeps a value when  flag >= the flaglevel you type.
  BATS flag -> ADS flag written:   2 (good) -> 2    3 (questionable) -> 1
                                   4 (bad)  -> 0    9 / empty        -> -9 (missing)
  So at the prompt: flaglevel 2 = good only, 1 = good+questionable, 0 = all.
 
STATION FIELDS
  cruise = Cruise_num, station number = Cast, date = UTC deployment date
  (YYMMDD), lat/lon = deployed position, bottom/min/max depth = deepest,
  shallowest and deepest sample of the cast, ship = vessel code,
  source = bc / ba / bb (BATS Core / Bloom A / Bloom B).
 
DEPTH LEVELS
  The ADS format allows at most 999 levels per station and the program's
  nzlevmax is 700, but BATS casts have up to 2201 samples.  Samples are
  therefore averaged into depth bins (default: 2 m down to 500 m, 10 m below).
  Within a bin only values of the best available flag are averaged.
  Use --no-bin to write every sample (fails if a cast has > --nzlevmax levels).
 
DAILY MODE (--daily)
  The program only keeps the date, so casts of the same day share one x value,
  which a smoothing spline cannot handle.  --daily writes one station per
  date (all casts of that date merged and averaged within each depth bin).
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd
 
BLANK = -99.99
# slot -> (width, decimals) of the default ("#else") read format in getdata_adsn
SLOTFMT = {1: (6, 2), 2: (7, 3), 3: (7, 2), 4: (7, 2), 5: (7, 2), 6: (7, 2),
           7: (6, 2), 8: (7, 2), 9: (6, 2), 10: (6, 2), 11: (6, 2), 12: (6, 2)}
STOCK_NAMES = ['temperature', 'salinity', 'pressure', 'DIC', 'alkalinity',
               'pCO2 ocean', 'delta 13C', 'oxygen', 'nitrate', 'nitrite',
               'silicate', 'phosphate']
STOCK_UNITS = ['[deg C]', '[psu]', '[dbar]', '[mumol/kg]', '[mueq/kg]', '[ppm]',
               '[per mil]', '[mumol/kg]', '[mumol/kg]', '[mumol/kg]',
               '[mumol/kg]', '[mumol/kg]']
VARS = {  # csv name: (default slot, tracer name, units)
    'Temperature': (1, 'temperature', '[deg C]'),
    'Salinity':    (2, 'salinity', '[psu]'),
    'Pressure':    (3, 'pressure', '[dbar]'),
    'Oxygen':      (8, 'oxygen', '[mumol/kg]'),
    'BAC':         (4, 'BAC', '[1/m]'),
    'Flu':         (5, 'fluorescence', '[rel]'),
    'PAR':         (6, 'PAR', '[uE/m2/s]'),
}
FLAGMAP = {2: 2, 3: 1, 4: 0, 9: -9}
VESSEL = {'R/V Weatherbird I': 'WB1', 'R/V Weatherbird II': 'WB2',
          'R/V Cape Henlopen': 'CHEN', 'R/V Cape Hatteras': 'CHAT'}
SOURCE = {'BATS Core': 'bc', 'BATS Bloom A': 'ba', 'BATS Bloom B': 'bb'}
 
 
def fmt(v, slot, name):
    w, d = SLOTFMT[slot]
    s = '%*.*f' % (w, d, v)
    if len(s) > w:
        sys.exit("value %s of %s does not fit slot %d (F%d.%d); use --scale or --map"
                 % (v, name, slot, w, d))
    return s
 
 
def parse_kv(items, what, cast):
    out = {}
    for it in items or []:
        if '=' not in it:
            sys.exit("--%s needs NAME=value, got %r" % (what, it))
        k, v = it.split('=', 1)
        if k not in VARS:
            sys.exit("unknown variable %r (known: %s)" % (k, ', '.join(VARS)))
        out[k] = cast(v)
    return out
 
 
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('infile')
    ap.add_argument('outfile', nargs='?', default=None,
                    help='default: input name with extension replaced by .ads2')
    ap.add_argument('--map', action='append', metavar='NAME=slot',
                    help="tracer slot for a variable, e.g. --map PAR=9 or PAR=none")
    ap.add_argument('--scale', action='append', metavar='NAME=factor',
                    help="multiply a variable before writing (keeps precision for "
                         "small values, e.g. --scale BAC=100)")
    ap.add_argument('--daily', action='store_true',
                    help='one station per date (casts of a day merged)')
    ap.add_argument('--no-bin', action='store_true', help='write every sample')
    ap.add_argument('--bin-shallow', type=float, default=2.0, help='bin size (m) above --bin-depth')
    ap.add_argument('--bin-deep', type=float, default=10.0, help='bin size (m) below --bin-depth')
    ap.add_argument('--bin-depth', type=float, default=500.0, help='depth (m) where bins change')
    ap.add_argument('--nzlevmax', type=int, default=700,
                    help='levels per station allowed by your tracer.h (default build: 700)')
    a = ap.parse_args()
 
    if a.outfile is None:
        a.outfile = os.path.splitext(a.infile)[0] + '.ads2'
    if os.path.abspath(a.outfile) == os.path.abspath(a.infile):
        sys.exit('output would overwrite the input file; give a different output name')
 
    # ---- variable -> slot, scale
    slots = {k: v[0] for k, v in VARS.items()}
    for k, v in parse_kv(a.map, 'map', lambda s: None if s.lower() == 'none' else int(s)).items():
        slots[k] = v
    scales = {k: 1.0 for k in VARS}
    scales.update(parse_kv(a.scale, 'scale', float))
    use = {k: s for k, s in slots.items() if s is not None}
    if any(not 1 <= s <= 12 for s in use.values()):
        sys.exit('slots must be 1-12')
    if len(set(use.values())) != len(use):
        sys.exit('two variables map to the same slot: %s' % use)
 
    t0 = time.time()
    need = ['ID', 'ISO_DateTime_UTC_deployed', 'Vessel', 'Latitude_deployed',
            'Longitude_deployed', 'Cruise_type', 'Cruise_num', 'Cast', 'Depth']
    need += list(use) + ['QF_' + k for k in use]
    head = pd.read_csv(a.infile, nrows=0).columns
    miss = [c for c in need if c not in head]
    if miss:
        sys.exit('columns missing in %s: %s' % (a.infile, miss))
    df = pd.read_csv(a.infile, usecols=need, low_memory=False)
    n_in = len(df)
    df = df.dropna(subset=['Depth', 'ISO_DateTime_UTC_deployed'])
    df = df.sort_values(['ISO_DateTime_UTC_deployed', 'ID', 'Depth'], kind='stable').reset_index(drop=True)
    df['date'] = df['ISO_DateTime_UTC_deployed'].str[:10]
    print('read %d rows, %d casts, %d dates (%.1fs)' % (n_in, df['ID'].nunique(), df['date'].nunique(), time.time() - t0))
 
    sid = df['date'] if a.daily else df['ID']          # station key
    df['sid'] = pd.factorize(sid)[0]                    # chronological (data sorted)
 
    # ---- vertical bins
    if a.no_bin:
        df['bin'] = np.arange(len(df))
    else:
        b1, b2, bd = a.bin_shallow, a.bin_deep, a.bin_depth
        n1 = int(np.ceil(bd / b1))
        d = df['Depth'].to_numpy()
        df['bin'] = np.where(d < n1 * b1, np.floor(d / b1), n1 + np.floor((d - n1 * b1) / b2)).astype(np.int64)
    gid = df.groupby(['sid', 'bin'], sort=True).ngroup()
    df['gid'] = gid
    g = df.groupby('gid', sort=True)
 
    lev = pd.DataFrame({'sid': g['sid'].first(), 'z': g['Depth'].mean()})
    unknown_qf = 0
    for k, slot in use.items():
        val = df[k] * scales[k]
        qf = df['QF_' + k]
        fl = qf.map(FLAGMAP)
        unknown_qf += int((fl.isna() & qf.notna()).sum())          # QF values not in FLAGMAP
        fl = fl.fillna(0).astype(int)
        missing = val.isna() | (qf == 9)
        fl = fl.where(~missing, -9)
        val = val.where(~missing)
        tmp = pd.DataFrame({'gid': df['gid'], 'val': val, 'fl': fl})
        best = tmp.groupby('gid')['fl'].transform('max')
        sel = (tmp['fl'] == best) & (tmp['fl'] > -9)
        lev['v%d' % slot] = tmp['val'].where(sel).groupby(tmp['gid']).mean()
        lev['f%d' % slot] = best.groupby(tmp['gid']).max()
    if unknown_qf:
        print('note: %d QF values not in {2,3,4,9} were treated as bad (flag 0)' % unknown_qf)
 
    nlev = lev.groupby('sid').size()
    worst = int(nlev.max())
    limit = min(999, a.nzlevmax)
    if worst > limit:
        sys.exit('a station has %d levels, limit is %d (I3 field = 999, your nzlevmax = %d).\n'
                 'Use coarser bins (--bin-deep / --bin-shallow) or raise nzlevmax in tracer.h '
                 'and pass --nzlevmax.' % (worst, limit, a.nzlevmax))
    if len(nlev) > 3000:
        print('WARNING: %d stations > nstamax (3000 in the default build)' % len(nlev))
 
    # ---- station table
    first = df.groupby('sid').first()
    zmin = df.groupby('sid')['Depth'].min(); zmax = df.groupby('sid')['Depth'].max()   # true sample depths
    cr = first['Cruise_num'].astype(int).astype(str)
    if (cr.str.len() > 5).any():
        sys.exit('Cruise_num longer than 5 characters')
    yrs = first['date'].str[:4].astype(int)
    if (yrs > 2010).any():
        print('WARNING: years after 2010: patch getplotdate.F (year .gt. 10 -> .gt. 50) '
              'or they are read as 19xx.')
 
    # ---- write
    order = sorted(use.items(), key=lambda kv: kv[1])
    colname = {s: '' for s in range(1, 13)}
    colunit = {s: '' for s in range(1, 13)}
    for k, s in use.items():
        nm = VARS[k][1] if scales[k] == 1 else (k + ' x%g' % scales[k])[:12]   # tracer.h names are 12 characters
        colname[s] = nm
        colunit[s] = VARS[k][2]
    lines = ['"header", cruise, station number, number of levels, date, latitude, longitude, '
             'bottom depth, minimum depth, maximum depth, ship, data source (BATS CTD, converted by csv_to_adsn.py)',
             ' bottle, cast, depth,' + ','.join('%-12s' % (colname[s] or 'dummy') for s in range(1, 13)) + ', flags',
             ' , ,[m],' + ','.join('%-12s' % (colunit[s] or '[]') for s in range(1, 13)) + ',']
    levs = lev.sort_values(['sid', 'z'])
    sidv = levs['sid'].to_numpy(); zv = levs['z'].to_numpy()
    colv = {s: levs['v%d' % s].to_numpy() for s in use.values()}
    colf = {s: levs['f%d' % s].to_numpy() for s in use.values()}
    starts = np.flatnonzero(np.r_[True, sidv[1:] != sidv[:-1]])
    ends = np.r_[starts[1:], len(sidv)]
    inv = {s: k for k, s in use.items()}
    for n, (i0, i1) in enumerate(zip(starts, ends)):
        sid_n = int(sidv[i0]); r = first.loc[sid_n]
        station_no = int(r['Cast']) if not a.daily else n + 1
        if station_no > 9999:
            sys.exit('station number %d does not fit I4' % station_no)
        ship = VESSEL.get(r['Vessel'], ''.join(ch for ch in str(r['Vessel']) if ch.isalnum())[:5])
        src = SOURCE.get(r['Cruise_type'], 'xx')
        d = r['date']
        lines.append('%-9s %-5s %4d %3d %6s %7.2f %7.2f %8.1f %8.1f %8.1f %-5s %-2s' % (
            '"header"', cr[sid_n], station_no, i1 - i0, d[2:4] + d[5:7] + d[8:10],
            r['Latitude_deployed'], r['Longitude_deployed'], zmax[sid_n], zmin[sid_n], zmax[sid_n],
            ship, src))
        cast_no = int(r['Cast']) % 1000
        for j in range(i0, i1):
            z = zv[j]
            zs = '%5.1f' % z if z < 1000 else '%5.0f' % z
            if len(zs) > 5:
                sys.exit('depth %s does not fit F5.0' % z)
            vals, flags = [], []
            for s in range(1, 13):
                if s in use.values() and not np.isnan(colv[s][j]):
                    vals.append(fmt(colv[s][j], s, inv[s])); flags.append(int(colf[s][j]))
                else:
                    vals.append(fmt(BLANK, s, 'blank')); flags.append(-9)
            lines.append('%3d %3d %s %s %s' % ((j - i0 + 1) % 1000, cast_no, zs, ' '.join(vals),
                                               ''.join('%2d' % f for f in flags)))
    with open(a.outfile, 'w', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
 
    nd = len(levs)
    print('wrote %s: %d stations, %d levels (max %d per station), %.1fs'
          % (a.outfile, len(starts), nd, worst, time.time() - t0))
    print('\nTracer slots (menu numbers in timeseries_data):')
    for k, s in order:
        w, dec = SLOTFMT[s]
        print('  %2d  %-12s -> %-14s %-11s written with %d decimals%s'
              % (s, k, colname[s], colunit[s], dec,
                 '  (x%g)' % scales[k] if scales[k] != 1 else ''))
    low = [k for k, s in order if k in ('BAC', 'Flu') and scales[k] == 1]
    if low:
        print('\nNote: %s are small numbers and only 2 decimals fit in the ADS columns '
              '(0.01 resolution).\n      Use e.g. --scale BAC=100 --scale Flu=100 if that is too coarse.'
              % ' and '.join(low))
    edits = [(s, k) for k, s in order if colname[s] != STOCK_NAMES[s - 1] or colunit[s] != STOCK_UNITS[s - 1]]
    if edits:
        print('\nEdit incl/tracer.h (data traname / data traunits), keeping the string lengths '
              '(12 and 10 characters):')
        for s, k in edits:
            print("  slot %2d  traname  '%-12s'   traunits '%-10s'   (was '%-12s' / '%-10s')"
                  % (s, colname[s], colunit[s], STOCK_NAMES[s - 1], STOCK_UNITS[s - 1]))
    print('\nAt the prompts: flaglevel 2 = good only, 1 = + questionable, 0 = all data.')
 
 
if __name__ == '__main__':
    main()
