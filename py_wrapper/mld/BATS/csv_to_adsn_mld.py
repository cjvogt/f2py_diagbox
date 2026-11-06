#!/usr/bin/env python3
"""
Convert an MLD table (tab- or comma-separated, with columns named
decimal_year and mld_m) into the default 12-tracer ADS file layout read by
timeseries_data (build WITHOUT -Ddericarbquant).
 
 * one station per calendar day (all casts of that day are averaged)
 * one depth level per station (depth 0), MLD stored in tracer column 1
 * rows with empty/invalid mld_m are skipped
 
The program rebuilds its time axis from a YYMMDD date as
    year + julday/365      (365-day year, julday = day of year, Jan 1 = 1)
so the date written here is chosen so that x matches decimal_year as
closely as that calendar allows (within ~0.5 day).
 
usage: csv_to_adsn.py input.csv [output.ads2] [--lat 31.67 --lon 64.17 --cruise MLD]
 
If the output name is omitted it is the input name with its extension replaced
by .ads2, in the same folder (mld.csv -> mld.ads2).
"""
import argparse, csv, math, os, sys
from collections import defaultdict
 
DPM = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]   # program's calendar
 
def julday_to_md(j):
    for m, n in enumerate(DPM, 1):
        if j <= n:
            return m, j
        j -= n
    return 12, 31
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile", nargs="?", default=None,
                    help="default: input name with extension replaced by .ads2")
    ap.add_argument("--lat", type=float, default=31.67)
    ap.add_argument("--lon", type=float, default=64.17)
    ap.add_argument("--cruise", default="MLD")
    a = ap.parse_args()
    if a.outfile is None:
        a.outfile = os.path.splitext(a.infile)[0] + ".ads2"
    if os.path.abspath(a.outfile) == os.path.abspath(a.infile):
        sys.exit("output would overwrite the input file; give a different output name")
 
    with open(a.infile, newline="") as fh:
        text = fh.read()
    delim = "\t" if text.split("\n", 1)[0].count("\t") else ","
    rows = list(csv.DictReader(text.splitlines(), delimiter=delim))
    if not rows or "decimal_year" not in rows[0] or "mld_m" not in rows[0]:
        sys.exit("need columns named decimal_year and mld_m; found: %s" % list(rows[0].keys() if rows else []))
 
    days = defaultdict(list)          # (year, julday) -> [mld values]
    nskip = 0; maxerr = 0.0
    for r in rows:
        try:
            dy = float(r["decimal_year"]); v = float(r["mld_m"])
        except (ValueError, TypeError):
            nskip += 1; continue
        if not (math.isfinite(dy) and math.isfinite(v)):
            nskip += 1; continue
        yr = int(math.floor(dy))
        j = min(365, max(1, int(round((dy - yr) * 365.0))))
        maxerr = max(maxerr, abs(yr + j / 365.0 - dy))
        days[(yr, j)].append(v)
 
    keys = sorted(days)
    if len(keys) > 3000:
        print("WARNING: %d stations > nstamax (3000)" % len(keys))
    if any(y > 2010 for y, _ in keys):
        print("WARNING: years after 2010 present. getplotdate.F reads 2-digit years > 10 as 19xx;"
              " patch it (see instructions) or those points land in 19xx.")
 
    out = []
    out.append('"header", cruise, station number, number of levels, date, latitude, longitude, bottom depth, minimum depth, maximum depth, shiptmp, data source')
    out.append(' bottle, cast, depth,MLD         ,dummy       ,dummy       ,dummy       ,dummy       ,dummy       ,dummy       ,dummy       ,dummy       ,dummy       ,dummy       ,dummy       , flags for tracer')
    out.append(' , ,[m],[m]       ,[]        ,[]        ,[]        ,[]        ,[]        ,[]        ,[]        ,[]        ,[]        ,[]        ,[]        ,')
    for n, (yr, j) in enumerate(keys, 1):
        vals = days[(yr, j)]
        mean = sum(vals) / len(vals)
        if abs(mean) >= 1000:
            sys.exit("MLD %.1f does not fit the F6.2 column" % mean)
        mo, d = julday_to_md(j)
        date = "%02d%02d%02d" % (yr % 100, mo, d)
        out.append('%-9s %-5s %4d %3d %6s %7.2f %7.2f %8.1f %8.1f %8.1f %-5s %-2s'
                   % ('"header"', a.cruise[:5], n % 10000, 1, date, a.lat, a.lon,
                      0.0, 0.0, 0.0, "unkwn", "ml"))
        out.append('%3d %3d %5.0f %6.2f %7.3f %7.2f %7.2f %7.2f %7.2f %6.2f %7.2f %6.2f %6.2f %6.2f %6.2f %s'
                   % (1, 1, 0.0, mean, -99.99, -99.99, -99.99, -99.99, -99.99,
                      -99.99, -99.99, -99.99, -99.99, -99.99, -99.99, "".join(" 0" for _ in range(12))))
    with open(a.outfile, "w", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    print("wrote %s" % a.outfile)
    print("rows read: %d, skipped (missing): %d" % (len(rows), nskip))
    print("stations written: %d (casts merged by day)" % len(keys))
    print("largest |x - decimal_year| after date conversion: %.4f yr (%.2f days)" % (maxerr, maxerr * 365))
 
if __name__ == "__main__":
    main()
