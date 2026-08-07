#!/usr/bin/env python3
"""agent V -- base loader.

Loads L's calibrated engine out of a PRIVATE mirror under agentV_work/mirror/L so that nothing
outside agentV_work is ever written.  The mirror's *.pkl caches were copied from agentT_work's
mirror (which rebuilt them from cold after the container restart); v_check.py below re-derives
L's own census numbers from them so the copy is validated rather than trusted.

Exposes: E, SL, SHIFT, p, NV, M, relift, vars_of, atomvalvars, influences, nzcount,
         assignment, ORIENT, T1, T2, factor, crt_list, fit, probe, peval, CGT2.
"""
import os, sys, json, time, itertools, random, collections
from math import gcd

V = '/home/user/integer_solver/solve_lab/agentV_work'
L = os.path.join(V, 'mirror', 'L')
_cwd0 = os.getcwd()
os.chdir(L)
sys.path.insert(0, L)
_g = {'__name__': 'v_drv'}
exec(compile(open(os.path.join(L, 'closeS4.py')).read().split("if __name__")[0], 'c4', 'exec'), _g)
os.chdir(_cwd0)

E = _g['E']; SL = _g['SL']; SHIFT = _g['SHIFT']; p = _g['p']; NV = _g['NV']; M = _g['M']
relift = _g['relift']; vars_of = _g['vars_of']; atomvalvars = _g['atomvalvars']
influences = _g['influences']; nzcount = _g['nzcount']; assignment = _g['assignment']
ORIENT = _g['ORIENT']; T1 = _g['T1']; T2 = _g['T2']
factor = _g['factor']; crt_list = _g['crt_list']
fit = _g['fit']; probe = _g['probe']; peval = _g['peval']
CGT2 = _g['CGT2']
GL = _g

TGTW = ('x24468', 'x18956')      # the two target-congruence wires


def wires_of(a):
    """the SHIFT wires an atom a can be moved on (L's own candidate set)."""
    return (set(q for q in vars_of(E.atoms[a]) if q in SHIFT) |
            set(q for q in atomvalvars[a] if q in SHIFT))


def greedy_init(S):
    """L's construction: build the assignment for ON-set S, then run the greedy lift to fixpoint."""
    v, isl, valn = assignment(set(S), ORIENT)
    v[24468] = T1; v[18956] = T2
    vv = [0]*NV
    for k, x in v.items():
        vv[k] = x
    for rd in range(60):
        bad = relift(vv)
        if not bad:
            break
        r = E.run(vv); fx = 0
        for a in bad:
            i = E.residx[a]; cur = r[i]; sm = abs(SL[a])
            if cur % p:
                continue
            imm = [q for q in vars_of(E.atoms[a]) if q in SHIFT]
            for w in imm + [q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
                old = vv[w]; vv[w] = old + p
                d = E.run(vv)[i] - cur
                vv[w] = old
                if d == 0:
                    continue
                g = gcd(d, sm)
                if cur % g:
                    continue
                mm = sm//g
                t = (-(cur//g))*pow((d//g) % mm, -1, mm) % mm if mm > 1 else 0
                vv[w] = old + p*t; fx += 1
                break
        if fx == 0:
            break
    return vv


def violated(vv, r=None):
    """the c>1 conditions currently NOT discharged (residual nonzero and not divisible by c*p)."""
    if r is None:
        r = E.run(vv)
    return [a for a in SL if r[E.residx[a]] != 0 and SL[a] and r[E.residx[a]] % abs(SL[a]) != 0]


def onset(n, seed=7):
    """L's own Random(7) ON-set convention, so sizes are comparable across agents."""
    if n == 2:
        return [24601, 2081]
    return random.Random(seed).sample(M['live'], n)


if __name__ == '__main__':
    print('NV=%d  atoms=%d  SHIFT wires=%d  c>1 atoms(CGT2)=%d' % (
        NV, len(E.res), len(SHIFT), len(CGT2)))
    print('p bits = %d' % p.bit_length())
