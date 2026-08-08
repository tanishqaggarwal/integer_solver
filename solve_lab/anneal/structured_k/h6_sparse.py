#!/usr/bin/env python3
"""h6_sparse.py -- k is sparse in a larger radix: few NONZERO DIGITS base 2^w.

Hamming weight (h3) only covers digits in {0,1}.  Here a generator is
d * 2^(w*i) for every digit value d and position i, all generators at position i
sharing group id i, so a combination picks at most one digit per position.
A meet-in-the-middle with <= H generators per side therefore covers every k
with at most 2H nonzero base-2^w digits.

usage: h6_sparse.py <bytes|nibbles|sbytes|snibbles> <maxhalf> <capbits>
"""
import sys, os, subprocess, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from instance import p, n, G, T, PTS, add, mul, neg

SK = os.path.join(HERE, 'sk_new')

SPECS = {                      # name -> (digit-bits, digit values)
    'nibbles':  (4, list(range(1, 16))),
    'bytes':    (8, list(range(1, 256))),
    'snibbles': (4, [d for d in range(-15, 16) if d]),
    'sbytes':   (8, [d for d in range(-255, 256) if d]),
}

def make_generators(w, digits, ngroups):
    """returns (scalars, points, gids) for d*2^(w*i), computed by repeated addition"""
    scalars, points, gids = [], [], []
    for i in range(ngroups):
        base = PTS[w * i]
        pos = {}
        cur = None
        for d in range(1, max(digits) + 1):
            cur = base if cur is None else add(cur, base)
            pos[d] = cur
        for d in digits:
            P = pos[abs(d)]
            if d < 0: P = neg(P)
            scalars.append(d * (1 << (w * i)))
            points.append(P)
            gids.append(i)
    return scalars, points, gids

def run(name, maxhalf, capbits, target, scalars, points, gids, tag):
    path = os.path.join(HERE, f'gen_{tag}.txt')
    with open(path, 'w') as f:
        for gi, P in zip(gids, points):
            f.write(f"{gi} {P[0]:064x} {P[1]:064x}\n")
        f.write(f"T {target[0]:064x} {target[1]:064x}\n")
    t0 = time.time()
    pr = subprocess.run([SK, 'gweight', path, str(maxhalf), str(capbits)],
                        capture_output=True, text=True)
    dt = time.time() - t0
    def decode(v):
        cnt = (v >> 60) & 0xF
        return sum(scalars[(v >> (15 * z)) & 0x7FFF] for z in range(cnt))
    sols, raw = [], 0
    for line in pr.stdout.split('\n'):
        f2 = line.split()
        if not f2: continue
        if f2[0] == 'WHIT':
            raw += 1
            Sq = decode(int(f2[1][2:])); St = decode(int(f2[2][2:]))
            for s in (1, -1):
                k = (Sq + s * St) % n
                if mul(k, G) == target: sols.append(k)
        elif f2[0] == 'WHITEXACT':
            raw += 1
            k = decode(int(f2[1][2:])) % n
            if mul(k, G) == target: sols.append(k)
    return {'seconds': round(dt, 1), 'raw_hits': raw, 'solutions': sorted(set(sols)),
            'stderr': pr.stderr.strip(), 'ngen': len(scalars)}

if __name__ == '__main__':
    name = sys.argv[1]; maxhalf = int(sys.argv[2]); capbits = int(sys.argv[3])
    w, digits = SPECS[name]
    ngroups = 256 // w
    scalars, points, gids = make_generators(w, digits, ngroups)
    print(f"{name}: radix 2^{w}, {ngroups} digit positions, {len(digits)} digit values, "
          f"{len(scalars)} generators; MITM half={maxhalf} covers "
          f"<= {2*maxhalf} nonzero digits")

    # ---- control over a reduced set of digit positions ----
    CG = 6
    cs, cp, cg = make_generators(w, digits, CG)
    dv = digits[:4] if len(digits) >= 4 else digits
    ctrl = 0
    for j in range(2 * maxhalf):
        ctrl += dv[j % len(dv)] * (1 << (w * j))
    ctrl %= n
    if 2 * maxhalf <= CG:
        cres = run(name, maxhalf, min(capbits, 24), mul(ctrl, G), cs, cp, cg, 'ctrl')
        ok = ctrl in cres['solutions']
        print(f"  control {ctrl}: {'RECOVERED' if ok else 'MISSED'} ({cres['seconds']}s)")
        if not ok: print("  CONTROL FAILED -- results not trustworthy")
    else:
        ok = None; print("  control skipped (2*maxhalf > reduced positions)")

    rres = run(name, maxhalf, capbits, T, scalars, points, gids, 'real')
    print(f"  real: raw_hits={rres['raw_hits']} solutions={rres['solutions']} ({rres['seconds']}s)")
    print("  " + rres['stderr'].replace('\n', '\n  '))
    if rres['solutions']:
        print(f"\n*** SOLVED ***  k = {rres['solutions'][0]}\n")
        json.dump({'solved': str(rres['solutions'][0])},
                  open(os.path.join(HERE, 'SOLVED.json'), 'w'))
    json.dump({'spec': name, 'radix_bits': w, 'positions': ngroups,
               'digit_values': len(digits), 'generators': len(scalars),
               'maxhalf': maxhalf, 'covers_nonzero_digits': 2 * maxhalf,
               'control_ok': ok, 'real': rres},
              open(os.path.join(HERE, f'h6_{name}_{maxhalf}.json'), 'w'), indent=1, default=str)
