#!/usr/bin/env python3
"""On-manifold experiment harness: definer DAG, free inputs, forward propagate."""
import os, pickle, json, sys, time
import jengine as E

HERE = os.path.dirname(os.path.abspath(__file__))
definer = E.build_definer()
order, cyc = E.topo(definer)
assert not cyc
ev = E.make_eval(definer)
FREE = set(range(E.NV)) - set(definer)
REF = E.load(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))
BASE = list(REF)
E.forward(BASE, order, ev, definer, FREE)
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663


def run(edits, base=None, verbose=True, tag=''):
    val = list(base if base is not None else BASE)
    for k, v in edits.items():
        val[k] = v
    bad = E.forward(val, order, ev, definer, FREE)
    s, fails, av = E.score(val)
    if verbose:
        nz = [i for i, x in enumerate(av) if x]
        print(f"  {tag} score={s} inexact={len(bad)} nz_atoms={len(nz)} "
              f"atoms={nz[:10]} nfail={len(fails)}")
    return s, fails, av, val


def save(val, s, name):
    p = os.path.join(HERE, f'J_{name}_{s}.json')
    E.save(val, p)
    return p


if __name__ == '__main__':
    s, f, av, v = run({}, tag='base')
    print('base fails', f)
