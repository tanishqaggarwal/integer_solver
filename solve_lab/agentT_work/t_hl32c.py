#!/usr/bin/env python3
"""AUDIT T35c -- forced-exact + freeze + guarded repair, tested on the SAVED |S|=32 end state.
Runs the full outer loop of the patched closer starting from close_M32.json rather than from
cold, so the mechanism can be checked in ~1 min instead of ~8.
"""
import os, sys, json, collections
from math import gcd
T = '/home/user/integer_solver/solve_lab/agentT_work'
sys.path.insert(0, T)
import t_close2w as C
import t_close2wj as J
E = C.E; SL = C.SL; SHIFT = C.SHIFT; p = C.p; NV = C.NV
relift = C.relift; vars_of = C.vars_of; atomvalvars = C.atomvalvars
influences = C.influences; nzcount = C.nzcount; solve_group3 = C.solve_group3

TAG = sys.argv[1] if len(sys.argv) > 1 else 'T32f'
vv = [0]*NV
for k, val in json.load(open(os.path.join(T, 'close_M32.json'))).items():
    vv[int(k[2:])] = int(val)
relift(vv)


def log(s):
    print(s, flush=True)

TGT = ('x24468', 'x18956')
gen = 0
for outer in range(20):
    base = nzcount(vv); r = E.run(vv); gen += 1
    viol = [a for a in SL if r[E.residx[a]] != 0 and SL[a] and r[E.residx[a]] % abs(SL[a]) != 0
            and not any(t in a for t in TGT)]
    hl0 = [a for a in E.res if r[E.residx[a]] and a not in SL]
    log('outer %d: global nonzero %d, violated c>1 (non-target) %d, nonzero handle-less %d'
        % (outer, base, len(viol), len(hl0)))
    if not viol and not hl0:
        log('   CLOSED'); break
    prog = 0
    if viol:
        wires = collections.defaultdict(list)
        for a in viol:
            for w in (set(q for q in vars_of(E.atoms[a]) if q in SHIFT) |
                      set(q for q in atomvalvars[a] if q in SHIFT)):
                wires[w].append(a)
        for w, ats in sorted(wires.items(), key=lambda kv: -len(kv[1])):
            V = [a for a in ats if influences(vv, a, w)]
            if not V:
                continue
            t = solve_group3(vv, V, w, gen, base)
            if t:
                log('   single-wire: x%d += p*%d accepted' % (w, t))
                prog += 1; base = nzcount(vv); gen += 1
    if prog:
        continue
    if J.handleless_pass(vv, nzcount(vv), log):
        continue
    if any(J.joint_pair(vv, a, base, log) for a in viol+hl0):
        gen += 1; continue
    if hl0 and J.forced_exact_pass(vv, hl0, log):
        gen += 1; continue
    log('   nothing moves -> stop'); break

relift(vv); r = E.run(vv)
nz = [E.res[i] for i, x in enumerate(r) if x]
print('\nFROZEN wires: %s' % sorted(J.FROZEN))
print('final nonzero atoms: %d' % len(nz))
for a in nz:
    print('   ', a[:100])
json.dump({'x_%d' % i: vv[i] for i in range(NV) if vv[i]},
          open(os.path.join(T, 'close_%s.json' % TAG), 'w'))
print('dumped close_%s.json' % TAG)
