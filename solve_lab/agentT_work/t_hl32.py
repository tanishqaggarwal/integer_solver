#!/usr/bin/env python3
"""AUDIT T35a -- diagnose why the |S|=32 handle-less atom survives, and check the mirror.

Loads close_M32.json (the |S|=32 end state left by t_close2wj.py), confirms the mirror reproduces
the published 3-nonzero-atom state, then interrogates the surviving HANDLE-LESS atom:
  * which wires influence it, R(t) degree and integer root on each;
  * what the forced root breaks (the collateral, by name);
  * whether the collateral shares a wire with it (the mixed-condition pair the tasking names);
  * and, for that pair, exactly which branch of joint_pair()'s mixed path bails out.
"""
import os, sys, json, collections, itertools
T = '/home/user/integer_solver/solve_lab/agentT_work'
sys.path.insert(0, T)
import t_close2w as C
import t_close2wj as J
E = C.E; SL = C.SL; SHIFT = C.SHIFT; p = C.p; NV = C.NV
relift = C.relift; vars_of = C.vars_of; atomvalvars = C.atomvalvars
influences = C.influences; nzcount = C.nzcount

STATE = sys.argv[1] if len(sys.argv) > 1 else 'close_M32.json'
vv = [0]*NV
for k, val in json.load(open(os.path.join(T, STATE))).items():
    vv[int(k[2:])] = int(val)
relift(vv)
r = E.run(vv)
nz = [E.res[i] for i, x in enumerate(r) if x]
print('%s reloaded in the MIRROR: %d nonzero atoms of %d' % (STATE, len(nz), len(E.res)))
for a in nz:
    kind = 'EXACT(handle-less)' if J.is_exact(a) else 'c=%d' % (abs(SL[a])//p)
    print('   %-70s %s   R%%p==0: %s' % (a[:70], kind, r[E.residx[a]] % p == 0))
base = nzcount(vv)
print('baseline global nonzero = %d' % base)

HL = [a for a in nz if a not in SL]
print('\nhandle-less nonzero atoms: %d' % len(HL))

def wires_of(a):
    ws = set(x for x in vars_of(E.atoms[a]) if x in SHIFT) | \
         set(x for x in atomvalvars[a] if x in SHIFT)
    return sorted(ws)

def newton(ys):
    d = [list(ys)]
    for k in range(len(ys)-1):
        d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
    return [d[k][0] for k in range(len(ys))]

ZERO0 = set(i for i, x in enumerate(r) if x == 0)

for a in HL:
    i = E.residx[a]
    ws = wires_of(a)
    inf = [w for w in ws if influences(vv, a, w)]
    print('\n=== %s' % a)
    print('    residual = %d   (mod p == 0: %s)' % (r[i], r[i] % p == 0))
    print('    syntactic shift wires %d, influencing %d: %s' % (len(ws), len(inf), inf))
    for w in inf:
        ys = []
        for t in range(8):
            ow = vv[w]; vv[w] = ow + p*t; ys.append(E.run(vv)[i]); vv[w] = ow
        cf = newton(ys)
        top = max([k for k in range(8) if cf[k]], default=0)
        rs = J.newton_int_roots(cf[:6])
        print('    wire x%-6d deg=%d  newton=%s  int roots=%s'
              % (w, top, [str(x)[:14] for x in cf[:top+2]], rs))
        if isinstance(rs, list) and rs:
            for t in rs:
                snap = vv[:]
                vv[w] += p*t
                y = E.run(vv)[i]
                n = nzcount(vv)
                rn = E.run(vv)
                broke = [E.res[i2] for i2 in ZERO0 if rn[i2]]
                print('      t=%d -> atom = %d ; global %d -> %d ; breaks %d: %s'
                      % (t, y, base, n, len(broke), [b[:60] for b in broke]))
                for b in broke:
                    bw = [x for x in wires_of(b) if influences(vv, b, x)]
                    print('         broken %-58s %s  wires=%s  shares w? %s'
                          % (b[:58],
                             'EXACT' if J.is_exact(b) else 'c=%d' % (abs(SL[b])//p),
                             bw, w in bw))
                vv[:] = snap
