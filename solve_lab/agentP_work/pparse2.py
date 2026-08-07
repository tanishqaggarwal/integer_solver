#!/usr/bin/env python3
"""Agent P independent parse, v2: iterative peel of `scalar * L^k`."""
import sys, json, pickle
from collections import defaultdict
sys.setrecursionlimit(100000)
from pparse import parse, topoly, flat_add, flat_mul, is_const, constval, key, pmul, padd

EQ = '/home/user/integer_solver/EQUATIONS.txt'

def peel(ast):
    """Return (scalar, power, Lnode, distinct_factor_count_history)."""
    scal = 1
    pw = 1
    node = ast
    hist = []
    while True:
        summands = flat_add(node)
        parsed = []
        for sg, nd in summands:
            fac = flat_mul(nd)
            sc = sg; nc = []
            for f in fac:
                if f == ('NEG',): sc = -sc
                elif is_const(f): sc *= constval(f)
                else: nc.append(f)
            parsed.append((sc, nc))
        sigs = set(tuple(sorted(key(topoly(f)) for f in nc)) for sc, nc in parsed)
        if len(sigs) != 1:
            break
        nc0 = parsed[0][1]
        if len(nc0) == 0:
            # pure constant equation
            hist.append('CONST')
            scal *= sum(sc for sc, _ in parsed)
            return scal, 0, None, hist
        fk = set(key(topoly(f)) for f in nc0)
        if len(fk) != 1:
            hist.append(('MULTIFACTOR', len(fk)))
            break
        s = sum(sc for sc, _ in parsed)
        if len(parsed) == 1 and len(nc0) == 1 and s == 1 and nc0[0] is node:
            break  # no progress
        scal *= s
        pw *= len(nc0)
        node = nc0[0]
        if len(parsed) == 1 and len(nc0) == 1 and s == 1:
            # unwrapped a single parenthesis; keep going but guard
            pass
    return scal, pw, node, hist

def main():
    lines = [l.strip() for l in open(EQ) if l.strip()]
    atom_key_to_id = {}
    atom_polys = []
    eq_rows = []
    stats = defaultdict(int)
    for ei, line in enumerate(lines):
        lhs = line.rsplit('=', 1)[0]
        ast = parse(lhs)
        scal, pw, L, hist = peel(ast)
        stats[(scal != 1, pw, tuple(hist))] += 0
        stats[('pw', pw)] += 1
        if hist: stats[('hist', str(hist))] += 1
        row = []
        if L is not None:
            for sg, nd in flat_add(L):
                fac = flat_mul(nd)
                sc = sg; nc = []
                for f in fac:
                    if f == ('NEG',): sc = -sc
                    elif is_const(f): sc *= constval(f)
                    else: nc.append(f)
                if not nc:
                    ap = {(): 1}
                else:
                    ap = {(): 1}
                    for f in nc: ap = pmul(ap, topoly(f))
                k = key(ap)
                aid = atom_key_to_id.get(k)
                if aid is None:
                    aid = len(atom_polys); atom_key_to_id[k] = aid; atom_polys.append(ap)
                row.append((sc, aid))
        eq_rows.append({'scal': scal, 'pw': pw, 'row': row})
        stats[('natoms', len(row))] += 1
    print("power distribution:", {k[1]: v for k, v in stats.items() if k[0] == 'pw'})
    print("hist:", {k[1]: v for k, v in stats.items() if k[0] == 'hist'})
    nat = sorted((k[1], v) for k, v in stats.items() if k[0] == 'natoms')
    print("atoms-per-equation distribution:", nat)
    print("distinct atoms:", len(atom_polys))
    # atom degree profile
    degs = defaultdict(int)
    nv = defaultdict(int)
    for ap in atom_polys:
        d = max((len(m) for m in ap), default=0)
        degs[d] += 1
        vs = set()
        for m in ap: vs.update(m)
        nv[len(vs)] += 1
    print("atom degree:", dict(sorted(degs.items())))
    print("atom #vars:", dict(sorted(nv.items())))
    with open('/home/user/integer_solver/solve_lab/agentP_work/model2.pkl', 'wb') as f:
        pickle.dump({'eq_rows': eq_rows, 'atom_polys': atom_polys}, f)
    print("saved model2.pkl")

if __name__ == '__main__':
    main()
