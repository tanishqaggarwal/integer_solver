#!/usr/bin/env python3
"""Deterministic endgame finisher.

Given a candidate free-input assignment `sol` (dict {varid:int}) that:
  (1) satisfies all 39022 non-gap equations exactly over Z, AND
  (2) has x_7068 ≡ M1 (mod p) and x_4432 ≡ M2 (mod p)
      where M1 = sol's x_2099 (=x_6418) value, M2 = sol's x_19964 (=x_12553),
this applies the four verified side-effect-free absorbers to zero G1,G2 exactly,
producing a full 39033 solution. Returns the completed dict or None if the
mod-p / mod-3 carry conditions are not met.

Absorbers (all verified side-effect-free; multipliers x_38744=x_22972=p):
  x_2099  = M1 + 15804267*p*x_3387     (x_3387 free; x_6418 follows the pin)
  x_19964 = M2 + p*x_5081              (x_5081 free; x_12553 follows the pin)
  G1 slack 7376877*p*x_17325 ; G2 slack p*x_9413.
Closure: G1=0 needs (x_7068-M1) = p*(15804267*x_3387 + 7376877*x_17325);
  since gcd(15804267,7376877)=3, requires 3 | (x_7068-M1)/p.  G2=0 needs
  (x_4432-M2) = p*(x_5081 + x_9413), any carry (coeff 1).
"""
import sys, os
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
os.chdir('/home/user/integer_solver/solve_lab')
import heal_harness as H
from math import gcd

p = H.p

def egcd(a, b):
    if b == 0: return (a, 1, 0)
    d, x, y = egcd(b, a % b); return (d, y, x - (a // b) * y)

def finish(sol, verbose=True):
    """sol: dict {varid:int} over free inputs. Returns completed dict or None."""
    for v in H.freeinp: H.val[v] = sol.get(v, 0)
    H.forward()
    M1 = H.val[2099]; M2 = H.val[19964]
    d1 = H.val[7068] - M1
    d2 = H.val[4432] - M2
    if d1 % p != 0:
        if verbose: print('FAIL condition (A): x_7068 != M1 mod p'); return None
    if d2 % p != 0:
        if verbose: print('FAIL condition (C): x_4432 != M2 mod p'); return None
    k1 = d1 // p; k2 = d2 // p
    a, b = 15804267, 7376877
    g = gcd(a, b)
    if k1 % g != 0:
        if verbose: print(f'FAIL condition (B\'): 3 does not divide (x_7068-M1)/p (k1 mod 3={k1%g})'); return None
    dd, x0, y0 = egcd(a, b); mult = k1 // g
    x3387, x17325 = x0 * mult, y0 * mult
    # apply G1 absorbers
    H.val[3387] = x3387; H.forward(); x26777 = H.val[26777]
    H.val[6418] = M1 + 15804267 * x26777
    H.val[17325] = x17325
    # apply G2 absorbers
    H.val[5081] = k2; H.forward(); x13458 = H.val[13458]
    H.val[12553] = M2 + x13458
    H.val[9413] = 0
    H.forward()
    nf = len(H.fails())
    if verbose: print(f'after finisher: fails = {nf}')
    if nf != 0:
        if verbose: print('  remaining fails:', sorted(H.fails())[:20]); return None
    return {v: H.val[v] for v in H.freeinp}

if __name__ == '__main__':
    import json
    path = sys.argv[1] if len(sys.argv) > 1 else 'best_agentA_39022.json'
    sol = H.loadd(path)
    out = finish(sol)
    if out is not None:
        outpath = 'SOLVED_finished.json'
        json.dump({('x_%d' % v): out[v] for v in H.freeinp}, open(outpath, 'w'))
        print('WROTE', outpath, '-> verify: python3 ../checker.py solve_lab/'+outpath)
    else:
        print('candidate does not satisfy the mod-p preconditions (expected for best_agentA).')
