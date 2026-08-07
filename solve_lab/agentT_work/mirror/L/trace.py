import sys, os, json, re, collections, pickle
F = '/home/user/integer_solver/solve_lab/agentT_work/mirror/F'
sys.path.insert(0, F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E = Engine()
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
defvar = {}   # var -> atom string
for a in E.order:
    c = E.cls[a]
    defvar[c[1]] = a
defrhs = {c[1]: c[2] for c in (E.cls[a] for a in E.order)}
resby = collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)

def dtree(v, depth=0, maxd=6, seen=None):
    pad = '  '*depth
    if v in defrhs:
        s = node_str(defrhs[v])
        print('%sx%d := %s' % (pad, v, s if len(s)<160 else s[:160]+'...'))
        if depth < maxd:
            for u in sorted(vars_of(defrhs[v])):
                dtree(u, depth+1, maxd)
    else:
        print('%sx%d  FREE   res-atoms: %s' % (pad, v, [a[:90] for a in resby.get(v,[])][:4]))

if __name__ == '__main__':
    for v in map(int, sys.argv[1:]):
        dtree(v, 0, int(os.environ.get('MAXD','3')))
        print('---')
