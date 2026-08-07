#!/usr/bin/env python3
"""AUDIT T33b -- what the two extra nonzero atoms at |S|=32 are.
They are NOT c>1 divisibility conditions (those all discharged); classify them: is the residual
divisible by p at all (a lift/divisibility problem) or not (a genuine mod-p failure, a different
obstruction class entirely)?"""
import os, sys, json
T = '/home/user/integer_solver/solve_lab/agentT_work'
sys.path.insert(0, T)
import t_close2w as C
E = C.E; SL = C.SL; p = C.p; NV = C.NV; relift = C.relift; SHIFT = C.SHIFT
vars_of = C.vars_of; atomvalvars = C.atomvalvars; influences = C.influences
for tag in sys.argv[1:] or ['T32g']:
    vv = [0]*NV
    for k, val in json.load(open(os.path.join(T, 'close_%s.json' % tag))).items():
        vv[int(k[2:])] = int(val)
    relift(vv)
    r = E.run(vv)
    nz = [E.res[i] for i, x in enumerate(r) if x]
    print('\n== %s : %d nonzero atoms ==' % (tag, len(nz)))
    for a in nz:
        val = r[E.residx[a]]
        s = SL.get(a)
        cls = ('c>1 divisibility, c=%d' % (abs(s)//p) if s and s % p == 0 and abs(s)//p > 1
               else 'c==1 (slope +-p)' if s and s % p == 0
               else 'slope 0' if s == 0
               else 'NO single handle' if s is None else 'slope not p-divisible')
        ws = [w for w in (set(x for x in vars_of(E.atoms[a]) if x in SHIFT) |
                          set(x for x in atomvalvars[a] if x in SHIFT))]
        print('  %-52s' % a[:52])
        print('     class            : %s' % cls)
        print('     residual %% p == 0 : %s   (residual has %d digits)'
              % (val % p == 0, len(str(abs(val)))))
        print('     shift wires      : %d total, %d admitted by influences()'
              % (len(ws), sum(1 for w in ws if influences(vv, a, w))))
