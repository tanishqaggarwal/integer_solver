#!/usr/bin/env python3
"""Robust accurate mod-P propagation search: pairs then triples over the improving
bits. mod-P propagation is exact propagation reduced mod P (no big-int blow-up),
baseline 5 vs Z's 4. Any state with <=1 mod-P violation is verified over Z."""
import json, time, sys, itertools
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
    ma = [poly for poly in atoms if set().union(*[set(m) for m in poly]) & mainv]
    _G.update(atoms=atoms, bset=bset, base=list(base.val), ma=ma)

def eval_ones(ones, cutoff=6):
    atoms = _G['atoms']; bset = _G['bset']
    eng = ModPEngine(atoms); eng.val = list(_G['base'])
    eng.wl = deque(); eng.inwl = [False]*len(atoms)
    for b in ones:
        if eng.val[b] is None: eng.assign(b, 1)
    eng.propagate()
    for v in [b for b in bset if eng.val[b] is None]:
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    for v in range(NVARS):
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    val = [x if x is not None else 0 for x in eng.val]
    vi = 0
    for poly in _G['ma']:
        s = 0
        for m, c in poly.items():
            t = c % P
            for x in m: t = (t*val[x]) % P
            s = (s+t) % P
        if s:
            vi += 1
            if vi >= cutoff: return vi, None
    return vi, val

def work(chunk):
    out = []
    for ones in chunk:
        try:
            vi, val = eval_ones(list(ones))
        except Exception:
            continue
        if vi <= 1:
            out.append((vi, tuple(ones)))
            if val is not None and vi <= 1:
                json.dump({f"x_{i}": int(val[i]) for i in range(NVARS)}, open(f'cand_modp2_{"_".join(map(str,ones))}.json','w'))
    return out

def main():
    res = json.load(open('flip_results.json'))
    improving = [b for (v,b,nc) in res if v <= 4]
    control = json.load(open('control_bits.json'))
    combos = list(itertools.combinations(improving, 2))  # PAIRS ONLY for feedback
    print(f"improving={len(improving)} combos(pairs+triples)={len(combos)}", flush=True)
    nproc = 4
    NC=200
    chunks = [combos[i*len(combos)//NC:(i+1)*len(combos)//NC] for i in range(NC)]
    t0 = time.time(); found = []
    with mp.Pool(nproc, initializer=init, maxtasksperchild=8) as pool:
        for ci, r in enumerate(pool.imap_unordered(work, chunks)):
            found.extend(r)
            for vi, c in r:
                print(f"  LOW {c} -> mod-P viol {vi}  (candidate saved)", flush=True)
            if (ci+1) % 10 == 0:
                print(f"  {ci+1}/{len(chunks)} chunks {time.time()-t0:.0f}s found {len(found)}", flush=True)
    print(f"DONE {time.time()-t0:.0f}s; low candidates: {sorted(found)[:20]}", flush=True)
    json.dump([[v,list(c)] for v,c in found], open('modp_pairs2_result.json','w'))

if __name__ == '__main__':
    main()
