#!/usr/bin/env python3
"""AUDIT T33c -- can the two handle-less atoms at |S|=32 be zeroed by their one admitted wire?
They must go to EXACTLY zero over Z (no cofactor can absorb them), so the question is whether the
fitted integer polynomial R(t)/p has an integer root.  Fit exactly, then test every rational-root
candidate by DIRECT RECOMPUTATION."""
import os, sys, json
T = '/home/user/integer_solver/solve_lab/agentT_work'
sys.path.insert(0, T)
import t_close2w as C
E = C.E; SL = C.SL; p = C.p; NV = C.NV; relift = C.relift; SHIFT = C.SHIFT
vars_of = C.vars_of; atomvalvars = C.atomvalvars; influences = C.influences

vv = [0]*NV
for k, val in json.load(open(os.path.join(T, 'close_T32g.json'))).items():
    vv[int(k[2:])] = int(val)
relift(vv)
r = E.run(vv)
HL = [a for i, a in ((E.residx[a], a) for a in E.res) if r[E.residx[a]] and a not in SL]
print('handle-less nonzero atoms: %d' % len(HL))
for a in HL:
    i = E.residx[a]
    ws = sorted(set(x for x in vars_of(E.atoms[a]) if x in SHIFT) |
                set(x for x in atomvalvars[a] if x in SHIFT))
    ws = [w for w in ws if influences(vv, a, w)]
    print('\n%s   wires %s' % (a[:60], ws))
    for w in ws:
        ys = []
        okp = True
        for t in range(8):
            y = C.probe2(vv, i, w, t, w if False else ws[0], 0) if False else None
            ow = vv[w]; vv[w] = ow + p*t; y = E.run(vv)[i]; vv[w] = ow
            if y % p:
                okp = False
            ys.append(y)
        # exact Newton differences of R(t) itself (not R/p): degree and monotonicity
        d = [ys[:]]
        for k in range(7):
            d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
        cf = [d[k][0] for k in range(8)]
        top = max([k for k in range(8) if cf[k]], default=0)
        print('   wire x%-6d  R(t) degree %d   R(0) has %d digits   all R(t)%%p==0: %s'
              % (w, top, len(str(abs(ys[0]))), okp))
        print('      Newton coeffs (first 4): %s' % [str(c)[:26] for c in cf[:4]])
        # a linear R(t) = a + b*t has the integer root -a/b iff b | a
        if top == 1:
            a0, b0 = cf[0], cf[1]
            print('      linear: root exists over Z ? %s   (a=%d digits, b=%d digits)'
                  % (a0 % b0 == 0 if b0 else False, len(str(abs(a0))), len(str(abs(b0)))))
            if b0 and a0 % b0 == 0:
                t = -a0//b0
                ow = vv[w]; vv[w] = ow + p*t
                print('      DIRECT RECOMPUTATION at t=%d : R = %s' % (t, E.run(vv)[i]))
                vv[w] = ow
