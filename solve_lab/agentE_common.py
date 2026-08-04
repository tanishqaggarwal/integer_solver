#!/usr/bin/env python3
"""Shared harness for agentE experiments: wire union-find, atoms, gates, best partial, checker."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS

p = 2**256 - 2**32 - 977
HERE = __file__.rsplit('/', 1)[0]
hc = json.load(open(HERE+'/huge_consts.json')); C1 = int(hc['C1']); C2 = int(hc['C2'])

def build_wire():
    """Signed union-find over 2-term identity atoms. Returns (wire dict v->sign, find2, par, sgn)."""
    A = load_atoms()
    par = list(range(NVARS)); sgn = [1]*NVARS
    def find2(x):
        s = 1; r = x
        while par[r] != r: s *= sgn[r]; r = par[r]
        return r, s
    def union(a, b, rel):
        ra, sa = find2(a); rb, sb = find2(b)
        if ra == rb: return
        par[ra] = rb; sgn[ra] = rel*sb*sa
    for pp in A:
        vs = atom_vars(pp)
        if len(vs) == 2 and pp.get((), 0) == 0:
            v1, v2 = sorted(vs); c1 = pp.get((v1,), 0); c2 = pp.get((v2,), 0)
            qok = all(pp.get(k, 0) == 0 for k in pp if isinstance(k, tuple) and len(k) == 2)
            if qok and c1 != 0 and c2 != 0 and abs(c1) == abs(c2):
                rel = (-c2)//c1
                if rel in (1, -1): union(v1, v2, rel)
    r0, _ = find2(26064)
    wire = {v: find2(v)[1] for v in range(NVARS) if find2(v)[0] == r0}
    return wire, find2, A

def load_gates():
    gates = []
    with open(HERE+'/atoms/gates.jsonl') as f:
        for line in f:
            d = json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
    return gates

def load_best():
    return {int(k[2:]): v for k, v in json.load(open(HERE+'/best/new_instance_partial_39013.json')).items()}

def load_lines():
    return [L for L in open(HERE+'/../EQUATIONS.txt').read().split('\n') if L.strip()]

VAR = re.compile(r'x_(\d+)')

CORE = [2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892]
NONCORE13 = [8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257]

if __name__ == '__main__':
    wire, find2, A = build_wire()
    print(f"wire size {len(wire)}")
    best = load_best()
    print(f"best nonzero {sum(1 for v in best.values() if v)}")
    print(f"x_26064 = {best.get(26064)} == p? {best.get(26064)==p}")
    for m in [5101,32017,26789]:
        print(f"wire member x_{m}: sign={wire.get(m)}, best_val_bits={best.get(m,0).bit_length()}")
