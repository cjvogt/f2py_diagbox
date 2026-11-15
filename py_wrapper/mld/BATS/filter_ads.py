#!/usr/bin/env python3
# usage: filter_ads.py in.ads2 out.ads2 [flaglevel]
# used in the MLD workflow, use after ctd_csv_to_adsn.py
import sys
src, dst = sys.argv[1], sys.argv[2]
FLAGLEVEL = int(sys.argv[3]) if len(sys.argv) > 3 else 2
SFC, DEEP, BLANK = 10.0, 50.0, -99.99     # adjust if you like

lines = open(src).read().splitlines()
head, body = lines[:3], lines[3:]
stations, cur = [], None
for ln in body:
    if ln.startswith('"header" '):
        cur = [ln]; stations.append(cur)
    elif cur is not None:
        cur.append(ln)

kept, dropped = 0, 0
with open(dst, 'w') as out:
    out.write('\n'.join(head) + '\n')
    for st in stations:
        good = []
        for ln in st[1:]:
            t = ln.split()
            flags = ln[-24:]                          # 12 flags x I2
            f = [int(flags[2*i:2*i+2]) for i in range(3)]   # T, S, P
            vals = [float(t[3]), float(t[4]), float(t[5])]
            if all(fl >= FLAGLEVEL for fl in f) and all(v != BLANK for v in vals):
                good.append(float(t[2]))
        if good and min(good) <= SFC and max(good) >= DEEP:
            out.write('\n'.join(st) + '\n'); kept += 1
        else:
            dropped += 1
            h = st[0].split()
            print('dropped: cruise %s cast %s date %s (valid %d/%d levels, shallowest %s)'
                  % (h[1], h[2], h[4], len(good), len(st)-1, min(good) if good else 'none'))
print('kept %d, dropped %d' % (kept, dropped))