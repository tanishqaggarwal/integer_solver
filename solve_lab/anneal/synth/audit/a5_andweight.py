#!/usr/bin/env python3
"""ATTACK 2: AND-weight sufficiency.

Two questions:
  (i)  Is E=0 <=> solution independent of W?  YES if W>=1: every square>=0 and
       every Rosenberg AND penalty g=ij-2iz-2jz+3z >=0 (=0 iff z=ij), so E=0
       forces every gate correct for ANY W>=1.  We verify W>=1 always.
  (ii) finalize()'s LOCAL load bound (W = 1 + max_z sum|coef in pen containing z|)
       claims rigidity: it OMITS any contribution from z feeding ANOTHER AND gate.
       That omission is only safe if AND OUTPUTS never feed AND inputs.  We check
       that structurally on a real modmul and a real ladder, and also brute-check
       that no wrong-gate state can undercut a right one at the chosen W.
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
SQ = os.path.join(HERE, '..', '..', 'squeeze')
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
sys.path.insert(0, os.path.abspath(SQ))
from mmqb import MMQB
from mm import build_modmul
from ladder2 import build_win2
from ecsmall import curve, find


def and_of_and(Q):
    """count AND gates whose input i or j is itself an AND output."""
    andouts = set(Q.andcache.values())
    bad = []
    for (i, j), z in Q.andcache.items():
        if i in andouts or j in andouts:
            bad.append((i, j, z))
    return andouts, bad


def check_modmul(p, mult='schoolbook', red='naf', mode='wallace', square=False):
    s = p.bit_length()
    Q = MMQB(chunk=16, mode=mode)
    A = Q.mkword('A', s, lambda wv: wv['_a'])
    B = A if square else Q.mkword('B', s, lambda wv: wv['_b'])
    C = Q.mkword('C', s, lambda wv: wv['_c'])
    build_modmul(Q, p, A, B, C, mult=mult, leaf=4, red=red)
    Q.finalize()
    andouts, bad = and_of_and(Q)
    return Q.W, len(andouts), bad, Q.and_weight_ok()


def check_ladder(p, B, m, w, mode='wallace', mult='schoolbook'):
    add, mul = curve(p, B)
    G, order = find(p, B)
    M = (m + w - 1) // w
    table = [[mul(((t + 1) << (w * j)) % order, G) for t in range(1 << w)] for j in range(M)]
    off = sum(1 << (w * j) for j in range(M))
    Tp = add(mul(1 % order, G), mul(off % order, G))
    Q, U = build_win2(p, table, Tp, w, mode=mode, mult=mult, leaf=4, red='naf')
    Q.finalize()
    andouts, bad = and_of_and(Q)
    return Q.W, len(andouts), bad, Q.and_weight_ok()


if __name__ == '__main__':
    print("=== (i) W>=1 and (ii) no AND-of-AND on real modmuls ===")
    okall = True
    for p in (13, 61, 251):
        for mult in ('schoolbook', 'karatsuba', 'toom3'):
            for red in ('naf', 'quotient', 'fold'):
                for sq in (False, True):
                    W, na, bad, awok = check_modmul(p, mult=mult, red=red, mode='wallace', square=sq)
                    okall &= (W >= 1 and not bad and awok)
                    if bad or W < 1 or not awok:
                        print(f"  p={p} {mult} {red} sq={sq}: W={W} andvars={na} "
                              f"AND-of-AND={len(bad)} and_weight_ok={awok}  <-- ISSUE")
    print(f"  modmuls: W>=1 & no-AND-of-AND & and_weight_ok all hold = {okall}")

    print("\n=== no AND-of-AND on real ladders (squeezed comb) ===")
    lok = True
    for (p, B, m, w) in [(127,3,3,1),(251,1,4,2),(1021,3,3,1)]:
        W, na, bad, awok = check_ladder(p, B, m, w)
        lok &= (W >= 1 and not bad and awok)
        print(f"  p={p} m={m} w={w}: W={W} andvars={na} AND-of-AND={len(bad)} and_weight_ok={awok}"
              + ("  <-- ISSUE" if bad or W < 1 or not awok else ""))
    print(f"  ladders OK = {lok}")

    print("\n=== brute: at W=1 (minimum), does E=0 still force all gates correct? (p=3 modmul) ===")
    # exhaustive: force W_and=1 and confirm zero-energy set unchanged vs W=default
    import verify
    for red in ('naf','quotient','fold'):
        Q1, z1, t1 = None, None, None
        s = 2  # p=3
        # default
        Qd, zd, td = verify.L0X(3, mult='schoolbook', leaf=3, red=red, mode='binary')
        # forced W=1
        from mmqb import MMQB as MM
        p=3
        Q = MM(chunk=16, mode='binary'); Q.W_and = 1
        A=Q.mkword('A',2,lambda wv:wv['_a']); Bw=Q.mkword('B',2,lambda wv:wv['_b']); C=Q.mkword('C',2,lambda wv:wv['_c'])
        build_modmul(Q,p,A,Bw,C,mult='schoolbook',leaf=3,red=red); Q.finalize()
        zeros=set()
        for x in verify.zero_states(Q):
            assert Q.energy(x)==0
            a=sum(x[v]<<t for t,v in enumerate(A.bits)); b=sum(x[v]<<t for t,v in enumerate(Bw.bits)); c=sum(x[v]<<t for t,v in enumerate(C.bits))
            zeros.add((a,b,c))
        print(f"  red={red}: W_default={Qd.W} zeros(default)={len(zd)} ; W_forced=1 zeros={len(zeros)} equal={zeros==zd} == truth={zeros==td}")
    print("\nDONE")
