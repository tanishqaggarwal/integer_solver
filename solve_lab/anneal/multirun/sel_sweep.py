#!/usr/bin/env python3
"""sel_sweep.py -- fill in the one-hot look-up cost for every window width w = 1..11,
at s = 256 and at the REAL prime p = 2^256 - 2^32 - 977, and cross-check the
composite window cost against the independently measured window256_neq.json.

Merges into multirun/pieces256.json (key selW) and writes multirun/winreal.json.
"""
import json, os, sys, time

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, _HERE)

from pieces import piece_sel, S, P
from resources import marginal_window

PJ = os.path.join(_HERE, 'pieces256.json')
D = json.load(open(PJ))

for mode in ('binary', 'wallace'):
    for w in range(1, 12):
        k = f'sel{w}'
        if k in D[mode]:
            continue
        t0 = time.time()
        D[mode][k] = piece_sel(mode, w)
        print(f"{mode} sel w={w}: {D[mode][k]['vars']:,} qubits "
              f"({time.time()-t0:.0f}s)", flush=True)

json.dump(D, open(PJ, 'w'), indent=1)


# ---- composite window cost from the atoms, and the real-p correction --------
def V_scheme(mode, mu, w):
    """qubits for a comb over mu unknown scalar bits with window width w."""
    r = D[mode]
    M = -(-mu // w)                       # ceil
    if M < 1:
        return None
    return M * 2 * r[f'sel{w}']['vars'] + (M - 1) * r['add']['vars'] + r['final']['vars']


out = {'V': {}, 'atoms': {m: {k: v for k, v in D[m].items() if k != '_secs'}
                          for m in ('binary', 'wallace')}}
for mode in ('binary', 'wallace'):
    out['V'][mode] = {}
    for w in range(1, 12):
        out['V'][mode][str(w)] = {str(mu): V_scheme(mode, mu, w)
                                  for mu in list(range(1, 33)) + [40, 48, 64, 80, 96,
                                                                  112, 120, 128, 129,
                                                                  144, 160, 192, 224, 256]}

print("\n=== cross-check: full-instance total (mu = 256) ===")
print(f"{'mode':>8} {'w':>3} {'this file (real p)':>20} {'report.py (random p)':>22} {'ratio':>7}")
W = json.load(open(os.path.join(os.path.dirname(_HERE), 'window256_neq.json')))
chk = {}
for mode in ('binary', 'wallace'):
    for w in range(1, 12):
        k = f'{mode}_w{w}'
        if k not in W:
            continue
        M = -(-256 // w)
        rep = M * W[k]['vars']
        mine = V_scheme(mode, 256, w)
        chk[k] = dict(real_p=mine, random_p=rep, ratio=round(mine / rep, 4))
        print(f"{mode:>8} {w:3d} {mine:20,d} {rep:22,d} {mine/rep:7.3f}")
out['crosscheck_mu256'] = chk

# independent check that the atom decomposition reproduces marginal_window
print("\n=== atom decomposition vs. resources.marginal_window (random p, s=256) ===")
for mode in ('binary', 'wallace'):
    for w in (1, 8):
        t0 = time.time()
        v, c, jb = marginal_window(256, w, mode)
        r = D[mode]
        pred = 2 * r[f'sel{w}']['vars'] + r['add']['vars']
        print(f"{mode:>8} w={w}: direct(random p)={v:,}  atoms(real p)={pred:,}  "
              f"ratio={pred/v:.3f}  ({time.time()-t0:.0f}s)", flush=True)

json.dump(out, open(os.path.join(_HERE, 'winreal.json'), 'w'), indent=1)
print("\nwrote multirun/winreal.json")
