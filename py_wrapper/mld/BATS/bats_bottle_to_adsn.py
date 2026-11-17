#!/usr/bin/env python3
"""
Convert the BATS discrete-bottle CSV (SO_DateTime_UTC, Bottle_ID, ..., yyyymmdd,
time, decimal_year) into the ADSN format read by getdata_adsn.F
(find_vertgrad and friends).

The fixed-width layout reproduces, column for column, the ADSN files written by
csv_to_adsn_bottle.py / mld_to_ads2.py, e.g.

"header"  BAT      1  24 881021   31.83  -64.17    500.0      1.0    500.0 WB1   bt
  1  81   1.0  -99.99  36.389  -99.99 2008.90  -99.99  -99.99 -99.99  207.00   0.00   0.00   1.04 -99.99  ...flags

Tracer columns follow the order of the program's tracer menu:
   1 temperature   2 salinity   3 pressure   4 DIC        5 alkalinity  6 pCO2
   7 delta 13C     8 oxygen     9 nitrate   10 nitrite   11 silicate   12 phosphate
(pressure, pCO2 and delta 13C are not in the BATS bottle file -> written as missing)

One ADSN "station" = one (cruise, cast).  Bottles within a station are sorted by
depth (calc_vertgrad assumes this).  Values whose QF is not in --keep-flags are
written as missing (-99.99, flag -9).  A side file <outfile>.stations.csv maps
the ADSN station number back to cruise/cast.

usage:
  python3 bats_bottle_to_adsn.py bats_bottle.csv BATS_bottle.ads2
  python3 bats_bottle_to_adsn.py bats_bottle.csv BATS_bottle.ads2 --keep-flags 2,3
  python3 bats_bottle_to_adsn.py bats_bottle.csv BATS_bottle.ads2 --start 1988
"""
import argparse
import csv
import sys
from collections import defaultdict

MISS = -99.99

# (ADSN column label, CSV value column, CSV flag column, fortran width, decimals)
# None for the CSV column -> always missing.
TRACERS = [
    ("temperature", "Temperature",  "QF_Temp",     7, 2),
    ("salinity",    "Salinity",     "QF_Salinity", 8, 3),   # bottle S, CTD S as fallback
    ("pressure",    None,           None,          8, 2),
    ("DIC",         "CO2",          "QF_DIC",      8, 2),
    ("alkalinity",  "Alkalinity",   "QF_Alk",      8, 2),
    ("pCO2 ocean",  None,           None,          8, 2),
    ("delta 13C",   None,           None,          7, 2),
    ("oxygen",      "Oxygen_1",     "QF_Oxygen",   8, 2),
    ("nitrate",     "NO3_plus_NO2", "QF_NO3_NO2",  7, 2),
    ("nitrite",     "NO2",          "QF_NO2",      7, 2),
    ("silicate",    "Silicate",     "QF_Silicate", 7, 2),
    ("phosphate",   "PO4",          "QF_PO4",      7, 2),
]

UNITS = ["[deg C]", "[psu]", "[dbar]", "[mumol/kg]", "[mueq/kg]", "[muatm]",
         "[per mil]", "[mumol/kg]", "[mumol/kg]", "[mumol/kg]", "[mumol/kg]",
         "[mumol/kg]"]

SHIPS = {"r/v weatherbird i": "WB1", "r/v weatherbird ii": "WB2",
         "r/v atlantic explorer": "AE"}


def fnum(x, w, d):
    s = f"{x:{w}.{d}f}"
    if len(s) > w:
        raise ValueError(f"value {x} does not fit F{w}.{d}")
    return s


def get(row, col, flagcol, keep):
    """return (value, flag) or (MISS, -9)"""
    v = row.get(col, "").strip()
    q = row.get(flagcol, "").strip()
    if not v:
        return MISS, -9
    try:
        val = float(v)
        qf = int(float(q)) if q else 2
    except ValueError:
        return MISS, -9
    if val < -90 or qf not in keep:   # -999 style fill value or rejected QF
        return MISS, -9
    return val, qf


def header_line(stanum, cruise, nlev, yymmdd, lat, lon, zbot, zmin, zmax, ship, src):
    return ('"header"  ' + f"{cruise:<5s}" + f"{stanum:5d}" + f"{nlev:4d}" +
            f" {yymmdd:6s}" + f"{lat:8.2f}" + f"{lon:8.2f}" +
            f"{zbot:9.1f}" + f"{zmin:9.1f}" + f"{zmax:9.1f}" +
            f" {ship:<5s}" + f" {src:2s}")


def data_line(bottle, cast, depth, vals, flags):
    s = f"{bottle:3d}{cast:4d}" + fnum(depth, 6, 1)
    for (lab, _, _, w, d), v in zip(TRACERS, vals):
        s += fnum(v, w, d)
    s += " " + "".join(f"{f:2d}" for f in flags)
    return s


def file_header():
    l1 = ('"header", cruise, station number, number of levels, date, latitude, '
          'longitude, bottom depth, minimum depth, maximum depth, ship, '
          'data source (converted by bats_bottle_to_adsn.py)')
    l2 = " bottle, cast, depth," + ",".join(f"{t[0]:<12s}" for t in TRACERS) + ", flags"
    l3 = " , ,[m]," + ",".join(f"{u:<12s}" for u in UNITS) + ","
    return [l1, l2, l3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--keep-flags", default="2",
                    help="comma list of BATS QF values to keep (default: 2 = good)")
    ap.add_argument("--start", type=int, default=None, help="first year to include")
    ap.add_argument("--end", type=int, default=None, help="last year to include")
    ap.add_argument("--min-levels", type=int, default=2,
                    help="drop stations with fewer bottles (default 2)")
    args = ap.parse_args()
    keep = {int(x) for x in args.keep_flags.split(",")}

    casts = defaultdict(list)
    with open(args.infile, newline="") as f:
        for row in csv.DictReader(f):
            try:
                depth = float(row["Depth"])
                ymd = row["yyyymmdd"].strip()
                year = int(ymd[:4])
            except (ValueError, KeyError):
                continue
            if depth < 0 or row.get("QF_Depth", "2").strip() not in ("", "2", "1"):
                continue
            if args.start and year < args.start or args.end and year > args.end:
                continue

            vals, flags = [], []
            for lab, col, qcol, w, d in TRACERS:
                if col is None:
                    v, q = MISS, -9
                else:
                    v, q = get(row, col, qcol, keep)
                    if lab == "salinity" and v == MISS:
                        v, q = get(row, "CTD_Salinity", "QF_CTD_Sal", keep)
                vals.append(v)
                flags.append(q)

            key = (row["Cruise_num"].strip(), int(float(row["Cast_num"])))
            casts[key].append(dict(
                depth=depth, vals=vals, flags=flags, ymd=ymd,
                lat=float(row["Latitude"]), lon=float(row["Longitude"]),
                ship=SHIPS.get(row["Vessel"].strip().lower(), "unkwn")))

    # order stations by date, then cruise, then cast
    keys = sorted(casts, key=lambda k: (casts[k][0]["ymd"], k[0], k[1]))

    out = file_header()
    side = ["adsn_station,cruise,cast,yymmdd,nlev,zmin,zmax"]
    nsta, nlevmax, skipped = 0, 0, 0
    for cruise, cast in keys:
        bot = sorted(casts[(cruise, cast)], key=lambda b: b["depth"])
        # drop bottles with no tracer data at all
        bot = [b for b in bot if any(f != -9 for f in b["flags"])]
        if len(bot) < args.min_levels:
            skipped += 1
            continue
        nsta += 1
        nlev = len(bot)
        nlevmax = max(nlevmax, nlev)
        b0 = bot[0]
        yymmdd = b0["ymd"][2:8]
        zmin, zmax = bot[0]["depth"], bot[-1]["depth"]
        out.append(header_line(nsta, "BAT", nlev, yymmdd, b0["lat"], b0["lon"],
                               zmax, zmin, zmax, b0["ship"], "bt"))
        for i, b in enumerate(bot, 1):
            out.append(data_line(i, cast, b["depth"], b["vals"], b["flags"]))
        side.append(f"{nsta},{cruise},{cast},{yymmdd},{nlev},{zmin},{zmax}")

    with open(args.outfile, "w") as f:
        f.write("\n".join(out) + "\n")
    with open(args.outfile + ".stations.csv", "w") as f:
        f.write("\n".join(side) + "\n")

    print(f"wrote {nsta} stations to {args.outfile} "
          f"({skipped} casts skipped with < {args.min_levels} usable bottles)")
    print(f"max levels per station: {nlevmax}")
    print(f"-> need nstamax >= {nsta} and nzlevmax >= {nlevmax} in tracer.h")


if __name__ == "__main__":
    main()
