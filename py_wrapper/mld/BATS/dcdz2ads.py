#!/usr/bin/env python3
"""
dcdz2ads.py -- repackage find_vergrad.x output (dcdz / vertical gradient)
               into an ADS-format file that timeseries.x can read through
               getdata_adsn.F.

timeseries.x has no "dcdz" tracer, so dcdz is written into an existing
tracer slot (default: slot 2, the F7.3 column -> keeps 3 decimals).
In timeseries.x you then select that tracer (traout) in userinput.

Accepted inputs (auto-detected):
  * the full find_vergrad output  (12 cols: date decdate julday sta# lat lon
                                    mld tramld nobsmld vertgrad sg n)
  * the out.set file               ( 6 cols: decdate julday date mld tramld vertgrad)

Usage:
  python3 dcdz2ads.py INPUT OUTPUT.ads [options]

Options:
  --layout {default,carb}  data-line layout; must match the -D flag your
                           timeseries.x was compiled with:
                             default : no -D flag (12-tracer ADS line)
                             carb    : -Ddericarbquant3/3b/4/5/7, -Ddericarbonate,
                                       -Dderikeelingquant/2  (tra 1,2,3,13-21)
  --slot N                 tracer index to hold dcdz (default 2)
  --depth Z                depth (m) assigned to each value (default 10);
                           choose a timeseries depth interval that contains it
  --flag F                 quality flag for the dcdz value (default 2)
  --missing V              fill value for unused tracers (default -99.99)
  --missflag F             flag for unused tracers (default 9)
  --drop-nomld             drop rows whose tramld is -99.99 (the duplicate
                           "second cast" rows, e.g. -75.334, -44.102)
  --max-abs X              drop rows with |dcdz| > X (outlier cut)
  --min-n N                drop rows with n_vertgrad < N (full file only)
  --start T / --end T      keep only decimal years T_start <= t <= T_end
                           (e.g. --end 2000.9 reproduces the 1988-2000 window)
  --maxdist D              drop stations more than D degrees from BATS
                           (31.67N, 64.17W), e.g. the Oct-2014 survey stations
"""
import argparse, math, sys

# (tracer index, width, decimals) in the order they appear on a data line
LAYOUTS = {
    "default": [(1, 6, 2), (2, 7, 3), (3, 7, 2), (4, 7, 2), (5, 7, 2),
                (6, 7, 2), (7, 6, 2), (8, 7, 2), (9, 6, 2), (10, 6, 2),
                (11, 6, 2), (12, 6, 2)],
    "carb":    [(1, 6, 2), (2, 7, 3), (3, 7, 2)] +
               [(i, 7, 2) for i in range(13, 22)],
}
NFLAGS = 12


def fnum(x):
    try:
        v = float(x)
    except ValueError:
        return float("nan")
    return v


def read_input(path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if len(p) == 12 and p[0].isdigit() and len(p[0]) == 6:
                rows.append(dict(date=p[0], decdate=fnum(p[1]), sta=int(p[3]),
                                 lat=fnum(p[4]), lon=fnum(p[5]), mld=fnum(p[6]),
                                 tramld=fnum(p[7]), dcdz=fnum(p[9]),
                                 sg=fnum(p[10]), n=int(p[11])))
            elif len(p) == 6 and len(p[2]) == 6 and p[2].isdigit():
                rows.append(dict(date=p[2], decdate=fnum(p[0]), sta=0,
                                 lat=31.67, lon=-64.17, mld=fnum(p[3]),
                                 tramld=fnum(p[4]), dcdz=fnum(p[5]),
                                 sg=float("nan"), n=99))
    return rows


def ffmt(v, w, d):
    s = f"{v:{w}.{d}f}"
    if len(s) > w:
        raise ValueError(f"value {v} does not fit F{w}.{d}")
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input"); ap.add_argument("output")
    ap.add_argument("--layout", choices=LAYOUTS, default="default")
    ap.add_argument("--slot", type=int, default=2)
    ap.add_argument("--depth", type=float, default=10.0)
    ap.add_argument("--flag", type=int, default=2)
    ap.add_argument("--missing", type=float, default=-99.99)
    ap.add_argument("--missflag", type=int, default=9)
    ap.add_argument("--drop-nomld", action="store_true")
    ap.add_argument("--max-abs", type=float)
    ap.add_argument("--min-n", type=int)
    ap.add_argument("--start", type=float)
    ap.add_argument("--end", type=float)
    ap.add_argument("--maxdist", type=float)
    a = ap.parse_args()

    layout = LAYOUTS[a.layout]
    slots = [t for t, _, _ in layout]
    if a.slot not in slots:
        sys.exit(f"slot {a.slot} not on a '{a.layout}' data line; choose from {slots}")

    rows = read_input(a.input)
    kept, dropped = [], {}
    def drop(why):
        dropped[why] = dropped.get(why, 0) + 1
    for r in rows:
        v = r["dcdz"]
        if math.isnan(v) or v <= -99.98:
            drop("no dcdz (-99.99/NaN)"); continue
        if a.drop_nomld and r["tramld"] <= -99.98:
            drop("tramld = -99.99"); continue
        if a.max_abs is not None and abs(v) > a.max_abs:
            drop(f"|dcdz| > {a.max_abs}"); continue
        if a.min_n is not None and r["n"] < a.min_n:
            drop(f"n_vertgrad < {a.min_n}"); continue
        yy = int(r["date"][:2]); year = (1900 if yy > 50 else 2000) + yy
        t = year + (r["decdate"] % 1.0)
        if a.start is not None and t < a.start:
            drop("before --start"); continue
        if a.end is not None and t > a.end:
            drop("after --end"); continue
        if a.maxdist is not None and math.hypot(r["lat"] - 31.67,
                                                r["lon"] + 64.17) > a.maxdist:
            drop(f"> {a.maxdist} deg from BATS"); continue
        kept.append(r)

    with open(a.output, "w") as f:
        f.write("Comment         dcdz (vertical gradient below MLD) from find_vergrad.x\n")
        f.write(f"Comment         dcdz stored in tracer slot {a.slot}, depth {a.depth:g} m\n")
        f.write("Type            ADS\n")
        f.write('"header"\n')
        f.write("station header : junk crname sta# nlev yymmdd lat lon zbot zmin zmax ship src\n")
        f.write("data line      : bottle cast z tracers... flags\n")
        for i, r in enumerate(kept, 1):
            sta = r["sta"] if 0 < r["sta"] < 10000 else i
            hdr = (f"{'DCDZ':<9s} {'DCDZ':<5s} {sta:4d} {1:3d} {r['date']:6s} "
                   f"{ffmt(r['lat'],7,2)} {ffmt(r['lon'],7,2)} "
                   f"{ffmt(max(r['mld'],0),8,1)} {ffmt(a.depth,8,1)} {ffmt(a.depth,8,1)} "
                   f"{'VGRAD':<5s} {'FV':<2s}")
            vals, flags = [], []
            for t, w, d in layout:
                if t == a.slot:
                    vals.append(ffmt(r["dcdz"], w, d)); flags.append(a.flag)
                else:
                    vals.append(ffmt(a.missing, w, d)); flags.append(a.missflag)
            flags += [a.missflag] * (NFLAGS - len(flags))
            data = (f"{1:3d} {1:3d} {ffmt(a.depth,5,0)} " + " ".join(vals) + " "
                    + "".join(f"{fl:2d}" for fl in flags[:NFLAGS]))
            f.write(hdr + "\n" + data + "\n")

    print(f"read {len(rows)} rows, wrote {len(kept)} stations to {a.output}")
    for k, v in dropped.items():
        print(f"  dropped {v:4d}  {k}")


if __name__ == "__main__":
    main()
