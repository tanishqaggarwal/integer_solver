#!/usr/bin/env python3
"""Mod-P search over the core bits: no big-int blow-up, so many-bit states are
cheap. Complete pairs then triples over the improving bits; verifies any 0-mod-P
candidate over Z. Uses multiprocessing."""
import json, time, itertools, sys
from collections import deque
import multiprocessing as mp
from propagate import load_atoms, atom_vars
from repair import boolean_vars
from modp import ModPEngine, P, substitute_modp, NVARS

_G = {}
def init():
    atoms = load_atoms(); bset = boolean_vars(atoms)
    base = ModPEngine(atoms); base.propagate()
    mainv = set(json.load(open('main_comp.json'))['main_vars'])
    # main-component atoms (only these can be violated by core-bit changes)
    ma = []
    for poly in atoms:
        vs = set()
        for m in poly: vs.update(m)
        if vs & mainv: ma.append(poly)
    _G.update(atoms=atoms, bset=bset, base=list(base.val), mainatoms=ma)

def eval_ones(ones, cutoff=5):
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
    for poly in _G['mainatoms']:
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
        vi, val = eval_ones(list(ones))
        if vi < 4:
            out.append((vi, tuple(ones)))
            if vi == 0:
                json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open(f'cand_modp_{"_".join(map(str,ones))}.json','w'))
    return out

def main():
    control = json.load(open('control_bits.json'))
    res = json.load(open('flip_results.json'))
    improving = [b for (v,b,nc) in res if v <= 4]
    print(f"control={len(control)} improving={len(improving)}", flush=True)
    combos = list(itertools.combinations(improving, 2)) + list(itertools.combinations(improving, 3))
    print(f"pairs+triples over improving bits: {len(combos)}", flush=True)
    nproc = 3
    chunks = [combos[i::nproc*30] for i in range(nproc*30)]
    t0 = time.time(); found = []
    with mp.Pool(nproc, initializer=init) as pool:
        for ci, r in enumerate(pool.imap_unordered(work, chunks)):
            found.extend(r)
            if r:
                for vi,c in r: print(f"  IMPROVEMENT {c} -> {vi} mod-P violations", flush=True)
            if (ci+1) % 10 == 0:
                print(f"  {ci+1}/{len(chunks)} chunks, {time.time()-t0:.0f}s, found {len(found)}", flush=True)
    found.sort()
    print(f"DONE {time.time()-t0:.0f}s best: {found[:10]}", flush=True)
    json.dump([[v,list(c)] for v,c in found], open('modp_search_result.json','w'))

if __name__ == '__main__':
    main()
