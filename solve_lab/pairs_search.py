#!/usr/bin/env python3
"""Complete accurate pairs search over all 256 core bits (multiprocessed).
For each pair of control bits set to 1, propagate from the pin-baseline and
count violated atoms. Reports any pair with <4 violations (progress) or 0 (solve).
Also does singles for completeness."""
import json, time, sys, itertools
from collections import deque
import multiprocessing as mp
from propagate import load_atoms, atom_vars, substitute, solve_single
from repair import ProvEngine, boolean_vars

NVARS = 38748
_G = {}

def init():
    atoms = load_atoms()
    bset = boolean_vars(atoms)
    base = ProvEngine(atoms); base.propagate()
    _G['atoms'] = atoms; _G['bset'] = bset
    _G['bval'] = list(base.val); _G['bprov'] = list(base.prov); _G['bdom'] = dict(base.domain)

def evalset(ones):
    atoms = _G['atoms']; bset = _G['bset']
    eng = ProvEngine(atoms)
    eng.val = list(_G['bval']); eng.prov = list(_G['bprov']); eng.domain = dict(_G['bdom'])
    eng.wl = deque(); eng.inwl = [False] * len(atoms)
    for b in ones:
        if eng.val[b] is None: eng.assign(b, 1, ('f', ()))
    eng.propagate()
    for v in [b for b in bset if eng.val[b] is None]:
        if eng.val[v] is None: eng.assign(v, 0, ('f', ())); eng.propagate()
    for v in range(NVARS):
        if eng.val[v] is None: eng.assign(v, 0, ('f', ())); eng.propagate()
    val = [x if x is not None else 0 for x in eng.val]
    viol = 0
    for poly in atoms:
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s != 0: viol += 1
    return viol, val

def work(chunk):
    out = []
    for pair in chunk:
        viol, val = evalset(list(pair))
        if viol < 4:
            out.append((viol, pair))
            if viol == 0:
                json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open(f'cand_pairsolved_{pair[0]}_{pair[1]}.json', 'w'))
    return out

def main():
    control = json.load(open('control_bits.json'))
    pairs = list(itertools.combinations(control, 2))
    print(f"control bits: {len(control)}, pairs: {len(pairs)}", flush=True)
    # chunk
    nproc = 3
    chunks = [pairs[i::nproc*20] for i in range(nproc*20)]
    t0 = time.time()
    found = []
    with mp.Pool(nproc, initializer=init) as pool:
        for ci, res in enumerate(pool.imap_unordered(work, chunks)):
            found.extend(res)
            done = (ci + 1) / len(chunks)
            print(f"  chunk {ci+1}/{len(chunks)} done ({done*100:.0f}%), {time.time()-t0:.0f}s, improvements so far: {len(found)}", flush=True)
            for viol, pair in res:
                print(f"    PAIR x_{pair[0]},x_{pair[1]} -> {viol} violated", flush=True)
    found.sort()
    print(f"DONE pairs in {time.time()-t0:.0f}s. best: {found[:10]}", flush=True)
    json.dump([[v, list(p)] for v, p in found], open('pairs_result.json', 'w'))

if __name__ == '__main__':
    main()
