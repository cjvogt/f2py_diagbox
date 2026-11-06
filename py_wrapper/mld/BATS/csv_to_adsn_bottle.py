#!/usr/bin/env python3
"""
csv_to_adsn_bottle.py -- bottle table  ->  ADS file for timeseries_data
Columns: Date LAT LONG Z SAL TEMP d13C AVGDIC AVGALK station Note sDIC sALK
(tab- or comma-separated).  Default 12-tracer layout: build WITHOUT -Ddericarbquant.
 
usage: csv_to_adsn_bottle.py input.csv [output.ads2] [--cruise DIC] [--dayfirst] [--map NAME=slot]
 
Slots (= tracer menu numbers):  TEMP 1, SAL 2, AVGDIC 4, AVGALK 5, d13C 7 match the stock
names; sDIC -> 6 and sALK -> 8 have no slot of their own (the default layout has none, and
only slots 3,4,5,6,8 are wide enough for values near 2000), so tracer.h names for 6 and 8
must be edited -- the script prints the lines.
One station per (station, date); one level per row, sorted by Z.  The program keeps only
the date, so two stations on one day share an x value in time series.
Date: decimal years (1983.75) or calendar dates; either way the program only keeps the day.
Flags: present -> 2 (good), empty -> -9; if the Note column mentions DIC / ALK / Temp the
matching values still present get 1 (questionable).  flaglevel 2 = good only, 1 = + questionable.
If a file mixes sites (station column), run once per site with --station NAME: two stations
on the same date would otherwise share one x value, which a spline cannot handle.
Requires pandas.
"""
import argparse, os, sys
import numpy as np, pandas as pd
 
BLANK = -99.99
FMT = {1: (6, 2), 2: (7, 3), 3: (7, 2), 4: (7, 2), 5: (7, 2), 6: (7, 2), 7: (6, 2), 8: (7, 2),
       9: (6, 2), 10: (6, 2), 11: (6, 2), 12: (6, 2)}
STOCK = {1: ('temperature', '[deg C]'), 2: ('salinity', '[psu]'), 3: ('pressure', '[dbar]'),
         4: ('DIC', '[mumol/kg]'), 5: ('alkalinity', '[mueq/kg]'), 6: ('pCO2 ocean', '[ppm]'),
         7: ('delta 13C', '[per mil]'), 8: ('oxygen', '[mumol/kg]')}
VARS = {'TEMP': (1, 'temperature', '[deg C]'), 'SAL': (2, 'salinity', '[psu]'),
        'AVGDIC': (4, 'DIC', '[mumol/kg]'), 'AVGALK': (5, 'alkalinity', '[mueq/kg]'),
        'd13C': (7, 'delta 13C', '[per mil]'), 'sDIC': (6, 'sDIC', '[mumol/kg]'),
        'sALK': (8, 'sAlk', '[mueq/kg]')}
 
 
def fmt(v, slot, name):
    w, d = FMT[slot]
    s = '%*.*f' % (w, d, v)
    if len(s) > w:
        sys.exit('value %s of %s does not fit slot %d (F%d.%d)' % (v, name, slot, w, d))
    return s
 
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('infile'); ap.add_argument('outfile', nargs='?')
    ap.add_argument('--cruise', default='DIC', help='cruise ID, max 5 characters')
    ap.add_argument('--dayfirst', action='store_true', help='Date is DD/MM/YYYY')
    ap.add_argument('--map', action='append', metavar='NAME=slot')
    ap.add_argument('--station', help='keep only rows whose station column equals this (e.g. BATS or S)')
    a = ap.parse_args()
    out = a.outfile or os.path.splitext(a.infile)[0] + '.ads2'
    if os.path.abspath(out) == os.path.abspath(a.infile):
        sys.exit('output would overwrite the input file')
    slots = {k: v[0] for k, v in VARS.items()}
    for m in a.map or []:
        k, v = m.split('=')
        if k not in VARS: sys.exit('unknown variable %s (known: %s)' % (k, ', '.join(VARS)))
        slots[k] = None if v == 'none' else int(v)
    use = {k: s for k, s in slots.items() if s}
    if len(set(use.values())) != len(use) or any(not 1 <= s <= 12 for s in use.values()):
        sys.exit('bad slot assignment: %s' % use)
 
    df = pd.read_csv(a.infile, sep=None, engine='python')
    df.columns = [c.strip() for c in df.columns]
    if a.station:
        df = df[df['station'].astype(str) == a.station]
        if df.empty: sys.exit('no rows with station == %s' % a.station)
    miss = [c for c in ['Date', 'LAT', 'LONG', 'Z', 'station'] + list(use) if c not in df.columns]
    if miss: sys.exit('columns missing: %s' % miss)
    for c in ['LAT', 'LONG', 'Z'] + list(use):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    dn = pd.to_numeric(df['Date'], errors='coerce')
    if dn.notna().all() and dn.between(1800, 2200).all():      # decimal years, e.g. 1983.75
        yr = np.floor(dn).astype(int); j = np.clip(np.round((dn - yr) * 365), 1, 365).astype(int)
        cum = np.cumsum([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        mo = np.searchsorted(cum, j) + 1; dy = j - np.r_[0, cum][mo - 1]   # the program's 365-day calendar
        df['t'] = dn; df['yr'] = yr
        df['dstr'] = ['%02d%02d%02d' % (y % 100, m, d) for y, m, d in zip(yr, mo, dy)]
        print('Date read as decimal years (rounded to the day, within 0.5 day)')
    else:
        dt = pd.to_datetime(df['Date'], errors='coerce', dayfirst=a.dayfirst)
        df['t'] = dt.dt.year + (dt.dt.dayofyear - 0.5) / 365.0; df['yr'] = dt.dt.year
        df['dstr'] = dt.dt.strftime('%y%m%d')
    bad = df['t'].isna() | df['Z'].isna()
    if bad.any(): print('dropped %d rows with unreadable Date or Z' % bad.sum())
    df = df[~bad].sort_values(['t', 'station', 'Z'], kind='stable').reset_index(drop=True)
    df['sid'] = pd.factorize(df['dstr'] + '|' + df['station'].astype(str))[0]
    num = pd.to_numeric(df['station'], errors='coerce')
    use_num = num.notna().all() and num.max() <= 9999
 
    if (df['yr'] > 2010).any():
        print('WARNING: years after 2010: patch getplotdate.F (year .gt. 10 -> .gt. 50)')
    note = df['Note'].fillna('').astype(str).str.lower() if 'Note' in df.columns else pd.Series('', index=df.index)
    quest = {k: pd.Series(False, index=df.index) for k in use}
    for kw, vs in {'dic': ('AVGDIC', 'sDIC'), 'alk': ('AVGALK', 'sALK'), 'temp': ('TEMP',)}.items():
        for v in vs:
            if v in quest: quest[v] |= note.str.contains(kw)
    names = {s: (VARS[k][1], VARS[k][2]) for k, s in use.items()}
    lines = ['"header", cruise, station number, number of levels, date, latitude, longitude, '
             'bottom depth, minimum depth, maximum depth, ship, data source (converted by csv_to_adsn_bottle.py)',
             ' bottle, cast, depth,' + ','.join('%-12s' % names.get(s, ('dummy',))[0] for s in range(1, 13)) + ', flags',
             ' , ,[m],' + ','.join('%-12s' % names.get(s, ('', '[]'))[1] for s in range(1, 13)) + ',']
    nlev = df.groupby('sid').size().max()
    if nlev > 700: sys.exit('a station has %d levels (nzlevmax 700 / I3 limit 999)' % nlev)
    for sid, g in df.groupby('sid', sort=True):
        r = g.iloc[0]
        stn = int(num[g.index[0]]) if use_num else int(sid) + 1
        lines.append('%-9s %-5s %4d %3d %6s %7.2f %7.2f %8.1f %8.1f %8.1f %-5s %-2s' % (
            '"header"', a.cruise[:5], stn, len(g), r['dstr'],
            np.nan_to_num(g['LAT'].mean()), np.nan_to_num(g['LONG'].mean()),
            g['Z'].max(), g['Z'].min(), g['Z'].max(), 'unkwn', 'bt'))
        for j, (_, x) in enumerate(g.iterrows(), 1):
            z = x['Z']; zs = '%5.1f' % z if z < 1000 else '%5.0f' % z
            if len(zs) > 5: sys.exit('depth %s does not fit F5.0' % z)
            v, f = [], []
            inv = {s: k for k, s in use.items()}
            for s in range(1, 13):
                val = x[inv[s]] if s in inv else np.nan
                ok = pd.notna(val)
                v.append(fmt(val if ok else BLANK, s, inv.get(s, 'blank'))); f.append(-9 if not ok else (1 if s in inv and quest[inv[s]][x.name] else 2))
            lines.append('%3d %3d %s %s %s' % (j % 1000, 1, zs, ' '.join(v), ''.join('%2d' % q for q in f)))
    with open(out, 'w', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('wrote %s: %d stations, %d levels (max %d per station)' % (out, df['sid'].nunique(), len(df), nlev))
    print('tracer slots: ' + ', '.join('%s=%d' % (k, s) for k, s in sorted(use.items(), key=lambda t: t[1])))
    for s, (nm, un) in sorted(names.items()):
        if s not in STOCK or STOCK[s] != (nm, un):
            print("edit incl/tracer.h slot %2d:  traname '%-12s'  traunits '%-10s'" % (s, nm, un))
    print('flaglevel at the prompt: 2 or lower')
 
 
if __name__ == '__main__':
    main()
 
