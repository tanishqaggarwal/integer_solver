#!/usr/bin/env python3
"""Bit-blast the BV encoding to CNF and run the CDCL solvers bundled with pysat."""
import sys, time, json, os
from z3 import *
import sibling, z3enc, witness

def to_cnf(d, pin=None):
    C, bits = z3enc.build(d, 'bv', pin=pin)
    g = Goal()
    for c in C: g.add(c)
    t = Then(Tactic('simplify'), Tactic('bit-blast'), Tactic('simplify'), Tactic('tseitin-cnf'))
    sub = t(g)
    return sub[0]

def stats(name, cnf, tmo):
    from pysat.solvers import Solver as PS
    s = PS(name=name, bootstrap_with=cnf)
    t = time.time()
    r = s.solve()
    st = s.accum_stats()
    s.delete()
    return dict(solver=name, sat=r, time=round(time.time() - t, 2), stats=st)

if __name__ == '__main__':
    out = json.load(open('runs/cnf.json')) if os.path.exists('runs/cnf.json') else {}
    for m in (8, 10, 12, 16):
        d = sibling.instance(m); witness.witness(d, d['k'])
        for pin in (True, False):
            key = 'm%d_%s' % (m, 'pin' if pin else 'free')
            if key in out: continue
            t0 = time.time()
            gsub = to_cnf(d, pin=(d['k'] if pin else None))
            dim = gsub.dimacs()
            path = 'encodings/%s.cnf' % key
            open(path, 'w').write(dim)
            hdr = [l for l in dim.split('\n') if l.startswith('p cnf')]
            nv, nc = (int(hdr[0].split()[2]), int(hdr[0].split()[3])) if hdr else (0, 0)
            print('%s  bit-blast %.1fs  vars=%d clauses=%d' % (key, time.time() - t0, nv, nc), flush=True)
            cnf = [[int(x) for x in l.split()[:-1]] for l in dim.split('\n')
                   if l and not l[0] in 'pc']
            rec = {'vars': nv, 'clauses': nc, 'blast_time': round(time.time() - t0, 2), 'runs': []}
            for name in ('cadical195', 'kissat', 'cryptominisat', 'glucose42'):
                try:
                    r = stats(name, cnf, 300)
                except Exception as e:
                    r = {'solver': name, 'err': repr(e)[:120]}
                rec['runs'].append(r)
                print('   %-14s %s' % (name, r), flush=True)
            out[key] = rec
            json.dump(out, open('runs/cnf.json', 'w'), indent=1)
