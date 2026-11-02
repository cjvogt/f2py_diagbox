"""
Compute mixed layer depth (MLD) for each CTD cast and return it with the cast time.
Works with BCO-DMO BATS CTD files (dataset 3918) and similar long-format CSVs.
 
Columns are auto-detected (case-insensitive). If detection fails, the script prints
your file's headers so you can add the right name to CANDIDATES below.
 
Usage:
    python ctd_mld.py data.csv
    python ctd_mld.py data.csv --threshold 0.125 --out mld.csv --plot
 
Requires: pip install pandas numpy gsw matplotlib
"""
import argparse
import sys
import numpy as np
import pandas as pd
import gsw
 
# Add your own header names here if auto-detection misses them (lowercase).
CANDIDATES = {
    "cast":  ["cast_id", "cast", "cast_number", "profile"],
    "time":  ["time", "iso_datetime_deployed", "isodatetime_deployed", "datetime",
              "date_time", "date"],
    "depth": ["depth", "depth_m", "dep"],
    "temp":  ["temperature", "temp", "temperature_ctd", "t"],
    "sal":   ["salinity", "sal", "salinity_ctd", "practical_salinity"],
}
 
 
def find_columns(df):
    lower = {c.lower().strip(): c for c in df.columns}
    found = {}
    for key, names in CANDIDATES.items():
        found[key] = next((lower[n] for n in names if n in lower), None)
    # Fuzzy fallback for the time column, e.g. "ISO_DateTime_Deployed_UTC"
    if found["time"] is None:
        found["time"] = next((orig for low, orig in lower.items()
                              if "iso" in low and "deploy" in low), None)
    missing = [k for k in ("depth", "temp", "sal") if found[k] is None]
    if missing or (found["time"] is None and "decimal_year" not in lower):
        sys.exit(f"Could not find columns for: {missing or ['time']}\n"
                 f"Your file's headers are:\n  {list(df.columns)}\n"
                 f"Add the right names to CANDIDATES at the top of this script.")
    found["decyear"] = lower.get("decimal_year")
    return found
 
 
def decimal_year_to_datetime(y):
    y = np.asarray(y, float)
    year = np.floor(y).astype(int)
    start = pd.to_datetime(year.astype(str) + "-01-01")
    end = pd.to_datetime((year + 1).astype(str) + "-01-01")
    return start + (end - start) * (y - year)
 
 
def datetime_to_decimal_year(t):
    """Convert datetimes to decimal year (e.g. 1988-07-02 -> ~1988.5), leap-year aware."""
    t = pd.to_datetime(pd.Series(t), utc=True)
    year = t.dt.year
    start = pd.to_datetime(year.astype(str) + "-01-01", utc=True)
    end = pd.to_datetime((year + 1).astype(str) + "-01-01", utc=True)
    return (year + (t - start) / (end - start)).to_numpy()
 
 
def mld_from_cast(depth, temp, sal, threshold=0.03, ref_depth=15.0):
    """MLD = first depth where sigma0 exceeds its value at ref_depth by threshold (kg/m3)."""
    ok = ~(np.isnan(depth) | np.isnan(temp) | np.isnan(sal))
    depth, temp, sal = depth[ok], temp[ok], sal[ok]
    if depth.size < 3 or depth.min() > ref_depth or depth.max() < ref_depth:
        return np.nan
 
    order = np.argsort(depth)
    depth, temp, sal = depth[order], temp[order], sal[order]
 
    sa = gsw.SA_from_SP(sal, depth, 0, 0)       # depth (m) ~ pressure (dbar)
    ct = gsw.CT_from_t(sa, temp, depth)
    sigma0 = gsw.sigma0(sa, ct)
 
    sigma_ref = np.interp(ref_depth, depth, sigma0)
    below = np.where((depth > ref_depth) & (sigma0 >= sigma_ref + threshold))[0]
    if below.size == 0:
        return depth.max()                       # mixed to the bottom of the cast
    i = below[0]
    d0, d1, s0, s1 = depth[i - 1], depth[i], sigma0[i - 1], sigma0[i]
    return d0 + (sigma_ref + threshold - s0) * (d1 - d0) / (s1 - s0)
 
 
def compute_mld_table(df, threshold=0.03):
    c = find_columns(df)
    df = df.copy()
    for k in ("depth", "temp", "sal"):
        df[c[k]] = pd.to_numeric(df[c[k]], errors="coerce").replace(-999, np.nan)
 
    if c["time"] is not None:
        df["_time"] = pd.to_datetime(df[c["time"]], errors="coerce", utc=True)
    else:
        df["_time"] = decimal_year_to_datetime(df[c["decyear"]])
 
    # Split into profiles. A new profile starts whenever the cast label or time
    # changes, or depth goes back up (rows are assumed to be in file order).
    new_prof = (df[c["depth"]].diff() < 0) | (df["_time"] != df["_time"].shift())
    if c["cast"]:
        new_prof |= df[c["cast"]] != df[c["cast"]].shift()
    df["_profile"] = new_prof.cumsum()
 
    rows = []
    for pid, g in df.groupby("_profile"):
        mld = mld_from_cast(g[c["depth"]].to_numpy(float),
                            g[c["temp"]].to_numpy(float),
                            g[c["sal"]].to_numpy(float), threshold)
        label = g[c["cast"]].iloc[0] if c["cast"] else pid
        rows.append({"cast": label, "time": g["_time"].iloc[0], "mld_m": mld})
    out = pd.DataFrame(rows).sort_values("time").reset_index(drop=True)
    out.insert(2, "decimal_year", datetime_to_decimal_year(out["time"]).round(5))
    return out
 
 
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="MLD vs time from CTD data")
    p.add_argument("csv")
    p.add_argument("--threshold", type=float, default=0.03, help="kg/m3 (default 0.03)")
    p.add_argument("--out", default="mld_by_time.csv")
    p.add_argument("--plot", action="store_true")
    args = p.parse_args()
 
    result = compute_mld_table(pd.read_csv(args.csv, low_memory=False), args.threshold)
    result.to_csv(args.out, index=False)
    print(result.head(10).to_string(index=False))
    print("...")
    print(result.tail(5).to_string(index=False))
    print(f"\n{len(result)} casts, {result['time'].min()} to {result['time'].max()}")
    print(f"Full results written to {args.out}")
 
    if args.plot:
        import matplotlib.pyplot as plt
        plt.plot(result["time"], result["mld_m"], ".", ms=3)
        plt.gca().invert_yaxis()
        plt.ylabel("MLD (m)")
        plt.xlabel("Time")
        plt.tight_layout()
        plt.show()
 