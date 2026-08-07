#!/usr/bin/env python3
"""K26: drive the real equations from the 256 leaf selector bits alone (the only true
boolean inputs), close mod p, and read the root's two input pairs."""
import sys, os, json, time
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
from cascadep import CascadeP, NV, P

C = CascadeP()
vc = json.load(open(K + '/varclass2.json'))
handles, leafsel, otherbools, wires = vc['handles'], vc['leafsel'], vc['otherbools'], vc['wires']
defvars = [u for u in range(NV) if u not in set(C.E.free)]
ORDER = handles + leafsel + otherbools + defvars + wires


TX = '91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002'
TY = '125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626'
FORBID = tuple(i for i, n in enumerate(C.names) if TX in n or TY in n)


def drive(on, forward_only=True):
    seed = {u: 0 for u in handles}
    for u in leafsel: seed[u] = 1 if u in on else 0
    for u in otherbools: seed[u] = 0
    v, _ = C.close(seed, ORDER, forbid=FORBID if forward_only else ())
    return v


S = FD.SHIFT
def rootpair(v):
    return (((v[12186] + S) % P, v[16742]), ((v[14853] + S) % P, v[24908]))

def oncurve(pt):
    return (pt[1] ** 2 - pow(pt[0], 3, P) - FD.B) % P == 0


if __name__ == '__main__':
    D = FD.points()
    ch = json.load(open(K + '/chain.json'))
    bypow = {}
    for i_s, e in ch['exp'].items():
        bypow[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
    sel2exp = {ch['sel'][str(i)]: ch['exp'][str(i)] for i in range(256)}
    exp2sel = {e: s for s, e in sel2exp.items()}
    rs = json.load(open(K + '/rootsupport.json'))
    IA = set(sel2exp[s] for s in rs['A.x']) | set(sel2exp[s] for s in rs['A.y'])
    IB = set(sel2exp[s] for s in rs['B.x']) | set(sel2exp[s] for s in rs['B.y'])

    def foldexp(es):
        R = FD.INF
        for e in es: R = FD.add(R, bypow[e])
        return R

    t0 = time.time()
    v = drive(set(leafsel))
    A, B = rootpair(v)
    print('ALL ON: closed %.1fs' % (time.time() - t0))
    print('  A =', A, 'on curve', oncurve(A))
    print('  B =', B, 'on curve', oncurve(B))
    for lab, ia, ib in [('163->A', IA | {163}, IB), ('163->B', IA, IB | {163})]:
        print('  %-8s A match %s   B match %s' % (lab, foldexp(sorted(ia)) == A, foldexp(sorted(ib)) == B))

    # deliverable ON-set
    v2 = drive({2081, 24601})
    A2, B2 = rootpair(v2)
    print('ON={2081,24601}:  A =', A2)
    print('                  B =', B2)
    print('  leaf2081 =', bypow[sel2exp[2081]], ' leaf24601 =', bypow[sel2exp[24601]])
