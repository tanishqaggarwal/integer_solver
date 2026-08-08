#!/usr/bin/env python3
"""h3_weight.py -- hypothesis 3: k has low Hamming weight, or low signed-digit
(NAF-style, digits +-2^i) weight.

Meet in the middle: build a hash table of x(S) for every sum S of <= H distinct
generators, then query x(T - S') for every such S'.  A match means T = S' +- S,
so the whole cross product of size |set|^2 is covered by 2|set| point additions.
With H = 3 unsigned that is every k of Hamming weight <= 6; with H = 4 it is
every k of Hamming weight <= 8.

usage: h3_weight.py <weight|sweight> <maxhalf> <capbits>
"""
import sys, os, subprocess, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from instance import p, n, G, T, PTS, add, mul, neg, sub

SK = os.path.join(HERE, 'sk')
NP = 256

def write_pts(path, target, npts=256):
    with open(path, 'w') as f:
        for i, P in enumerate(PTS[:npts]):
            f.write(f"P{i} {P[0]:064x} {P[1]:064x}\n")
        f.write(f"T {target[0]:064x} {target[1]:064x}\n")

def decode(v, signed):
    """decode the packed generator-index list produced by sk.c"""
    cnt = (v >> 54) & 0x3FF
    terms = []
    for z in range(cnt):
        gi = (v >> (9 * z)) & 0x1FF
        idx = gi % NP
        sgn = 1 if gi < NP else -1
        terms.append(sgn * (1 << idx))
    return terms

def run(mode, maxhalf, capbits, target, tag, npts=256):
    ptsfile = os.path.join(HERE, f'pts_{tag}.txt')
    write_pts(ptsfile, target, npts)
    t0 = time.time()
    pr = subprocess.Popen([SK, mode, ptsfile, str(maxhalf), str(capbits)],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = pr.communicate()
    dt = time.time() - t0
    signed = (mode == 'sweight')
    sols, raw = [], 0
    for line in out.split('\n'):
        f = line.split()
        if not f: continue
        if f[0] == 'WHIT':
            raw += 1
            q = int(f[1][2:]); t = int(f[2][2:])
            Sq = sum(decode(q, signed)); St = sum(decode(t, signed))
            for s in (1, -1):
                k = (Sq + s * St) % n
                if mul(k, G) == target:
                    sols.append(k)
        elif f[0] == 'WHITEXACT':
            raw += 1
            q = int(f[1][2:])
            k = sum(decode(q, signed)) % n
            if mul(k, G) == target:
                sols.append(k)
    return {'mode': mode, 'maxhalf': maxhalf, 'seconds': round(dt, 1),
            'raw_hits': raw, 'solutions': sorted(set(sols)), 'stderr': err.strip()}

if __name__ == '__main__':
    mode = sys.argv[1]; maxhalf = int(sys.argv[2]); capbits = int(sys.argv[3])
    signed = (mode == 'sweight')

    # ---- control: a planted scalar of exactly 2*maxhalf terms, over a REDUCED
    # generator set (first CN chain points) so the control run stays cheap while
    # exercising exactly the same C code path. ----
    CN = 48
    bits = [3, 9, 17, 22, 30, 35, 41, 46][:2 * maxhalf]
    assert all(b < CN for b in bits)
    if signed:
        ctrl = sum((-1) ** j * (1 << b) for j, b in enumerate(bits)) % n
    else:
        ctrl = sum(1 << b for b in bits)
    print(f"control scalar ({2*maxhalf} terms over {CN} bit positions): {ctrl}")
    cres = run(mode, maxhalf, min(capbits, 24), mul(ctrl, G), 'ctrl', npts=CN)
    ok = ctrl % n in cres['solutions']
    print(f"  control: {'RECOVERED' if ok else 'MISSED'}  raw_hits={cres['raw_hits']} "
          f"({cres['seconds']}s)")
    print("  " + cres['stderr'].replace('\n', '\n  '))
    if not ok:
        print("CONTROL FAILED -- results below are not trustworthy")

    # ---- the real instance ----
    rres = run(mode, maxhalf, capbits, T, 'real')
    print(f"  real: raw_hits={rres['raw_hits']} solutions={rres['solutions']} "
          f"({rres['seconds']}s)")
    print("  " + rres['stderr'].replace('\n', '\n  '))
    if rres['solutions']:
        k = rres['solutions'][0]
        print(f"\n*** SOLVED ***  k = {k}\n")
        json.dump({'solved': str(k)}, open(os.path.join(HERE, 'SOLVED.json'), 'w'))
    out = {'control_ok': ok, 'control': cres, 'real': rres,
           'covers': ('signed-digit' if signed else 'Hamming') + f' weight <= {2*maxhalf}'}
    json.dump(out, open(os.path.join(HERE, f'h3_{mode}_{maxhalf}.json'), 'w'), indent=1,
              default=str)
