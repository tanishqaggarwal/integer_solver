#!/usr/bin/env python3
"""Experiment harness: load the fitted generative model, apply free-input edits,
propagate, score."""
import os, pickle, json, sys, time
import jengine as E

HERE = os.path.dirname(os.path.abspath(__file__))
F = pickle.load(open(os.path.join(HERE, 'jfit.pkl'), 'rb'))
definer, order, free = F['definer'], F['order'], F['free']
ev = E.make_eval(definer)
REF = E.load(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))


def run(edits, base=None, verbose=True):
    val = list(base if base is not None else REF)
    for k, v in edits.items():
        val[k] = v
    bad = E.forward(val, order, ev, definer, free)
    s, fails, av = E.score(val)
    if verbose:
        nz = [i for i, x in enumerate(av) if x]
        print(f"  score={s} inexact={len(bad)} nonzero_atoms={len(nz)} fails={fails[:15]}")
    return s, fails, av, val


def free_report(names):
    for n in names:
        print(f"x_{n}: free={n in free} val={str(REF[n])[:30]}")


if __name__ == '__main__':
    import ast
    print("free residual vars:")
    free_report([4432, 7068, 642, 9413, 17325, 29854, 31864, 1329, 10903,
                 9118, 8731, 2081, 4287, 31861, 14865, 6418, 12553, 19964, 2099])
    print("\nbaseline:")
    run({})
