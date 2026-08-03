#!/usr/bin/env python3
"""Custom method: linear-response solve in bit-space over GF(P).

The circuit is multilinear in the 256 core bits. Around the all-bits-0 point,
compute for every main atom its first-order response D[atom][bit] to setting each
bit=1 (mod P, via propagation). Then solve D·x = -base (mod P) for x in GF(P)^256.
If the atoms are dominated by their bit-linear part, a 0/1 rounding of x is the
witness. Verify exactly. This uses no SAT/SMT — just the circuit + linear algebra."""
import json, time, sys
from collections import deque
import multiprocessing as mp
from propagate import load_atoms
from repair import boolean_vars
from modp import ModPEngine, P, NVARS

_G = {}
def init():
    atoms = load_atoms(); bset = boolean_vars(atoms)
    base = ModPEngine(atoms); base.propagate()
    mainv = set(json.load(open('main_comp.json'))['main_vars'])
    ma_idx = [i for i, poly in enumerate(atoms)
              if any(v in mainv for m in poly for v in m)]
    _G.update(atoms=atoms, bset=bset, bval=list(base.val), ma_idx=ma_idx)

def atomvals(ones):
    """Return dict {atom_index: value mod P} for main atoms, after setting `ones` bits=1."""
    atoms = _G['atoms']; bset = _G['bset']
    eng = ModPEngine(atoms); eng.val = list(_G['bval']); eng.wl = deque(); eng.inwl = [False]*len(atoms)
    for b in ones:
        if eng.val[b] is None: eng.assign(b, 1)
    eng.propagate()
    for v in [b for b in bset if eng.val[b] is None]:
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    for v in range(NVARS):
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    val = [x if x is not None else 0 for x in eng.val]
    out = {}
    for ai in _G['ma_idx']:
        s = 0
        for m, c in atoms[ai].items():
            t = c % P
            for x in m: t = (t*val[x]) % P
            s = (s+t) % P
        out[ai] = s
    return out

def worker(bit):
    return bit, atomvals([bit] if bit is not None else [])

def main():
    t0 = time.time()
    control = json.load(open('control_bits.json'))
    tasks = [None] + control  # None = baseline
    results = {}
    with mp.Pool(4, initializer=init) as pool:
        for k, (bit, av) in enumerate(pool.imap_unordered(worker, tasks)):
            results[bit] = av
            if (k+1) % 32 == 0: print(f"  computed {k+1}/{len(tasks)} responses ({time.time()-t0:.0f}s)", flush=True)
    base = results[None]
    # atoms nonzero at base = the ones we must fix; plus keep others at 0
    ma_idx = sorted(base.keys())
    # Build D (rows=atoms, cols=bits) mod P, and target t = -base
    import numpy as np
    # rows: only atoms that are nonzero at base OR change under some bit
    changed = set()
    for b in control:
        for ai in ma_idx:
            if results[b][ai] != base[ai]: changed.add(ai)
    rows = sorted(changed)
    print(f"main atoms: {len(ma_idx)}, atoms sensitive to some bit: {len(rows)} ({time.time()-t0:.0f}s)", flush=True)
    # GF(P) Gaussian elimination on D x = t
    ncol = len(control)
    D = [dict() for _ in rows]   # sparse row: col->coef
    tvec = []
    rowidx = {ai: r for r, ai in enumerate(rows)}
    for j, b in enumerate(control):
        for ai in rows:
            d = (results[b][ai] - base[ai]) % P
            if d: D[rowidx[ai]][j] = d
    for ai in rows:
        tvec.append((-base[ai]) % P)
    # solve
    def inv(a): return pow(a % P, P-2, P)
    pivcol = {}; used=[False]*len(rows)
    # forward elimination
    order = list(range(len(rows)))
    for j in range(ncol):
        pr = None
        for r in order:
            if not used[r] and D[r].get(j):
                pr = r; break
        if pr is None: continue
        used[pr]=True; pivcol[j]=pr
        f = inv(D[pr][j])
        D[pr] = {k:(v*f)%P for k,v in D[pr].items()}; tvec[pr]=(tvec[pr]*f)%P
        for r in order:
            if r!=pr and D[r].get(j):
                g=D[r][j]
                for k,v in D[pr].items():
                    D[r][k]=(D[r].get(k,0)-g*v)%P
                    if D[r][k]==0: del D[r][k]
                tvec[r]=(tvec[r]-g*tvec[pr])%P
    # check consistency and read solution
    incons=0
    for r in range(len(rows)):
        if not D[r] and tvec[r]!=0: incons+=1
    x=[0]*ncol
    for j,pr in pivcol.items():
        # x_j = tvec[pr] - sum_{k>j free} D[pr][k] x_k ; set free vars=0
        x[j]=tvec[pr]
    print(f"linear system: rank {len(pivcol)}, inconsistent rows {incons}")
    # round: bits should be 0/1
    ones=[control[j] for j in range(ncol) if x[j]==1]
    other=[control[j] for j in range(ncol) if x[j] not in (0,1)]
    print(f"solution x: {len(ones)} bits =1, {len(other)} bits not-0/1 (linear model imperfect if >0)")
    json.dump({"ones":ones, "nonbinary":len(other)}, open('delta_solution.json','w'))
    # verify the 0/1 part
    print("verifying (mod P) the bits set to 1...", flush=True)
    init(); av=atomvals(ones)
    viol=sum(1 for ai in av if av[ai]!=0)
    print(f"  bits-from-linear-solve -> {viol} violated main atoms mod P (baseline was {sum(1 for ai in base if base[ai]!=0)})")
    # save candidate for exact check
    from modp import ModPEngine
    json.dump(ones, open('delta_ones.json','w'))
    print(f"done {time.time()-t0:.0f}s; wrote delta_ones.json (run checker over Z if promising)")

if __name__ == '__main__':
    main()
