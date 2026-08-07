"""Agent B: backward cone trace of a variable through the gate graph."""
import pickle, collections, json, sys

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model5.pkl','rb'))
facs, atoms, eqs = M['facs'], M['atoms'], M['eqs']
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663

occ = collections.defaultdict(list)
for i, p in enumerate(facs):
    vs = set()
    for m in p: vs.update(m)
    for v in vs: occ[v].append(i)

def load(path='/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'):
    val = [0]*38748
    for k, v in json.load(open(path)).items(): val[int(k.split('_')[1])] = int(v)
    return val

def s(p, lim=160):
    return ''.join('%s*%s' % (('%+d'%c) if abs(c) < 10**8 else ('+K' if c > 0 else '-K'),
                               '*'.join('x%d'%v for v in m) if m else '1')
                   for m, c in sorted(p.items()))[:lim]

def fv(p, val):
    t = 0
    for m, c in p.items():
        z = c
        for v in m:
            z *= val[v]
            if z == 0: break
        t += z
    return t

def desc(v, val):
    x = val[v]
    tag = ''
    if x == P: tag = ' =p'
    elif x == 0: tag = ' =0'
    elif x == 1: tag = ' =1'
    return 'x%d(%db%s)' % (v, x.bit_length(), tag)

def trace(roots, val, depth=3):
    seen = set(); frontier = set(roots)
    for d in range(depth):
        print("--- depth %d, %d vars" % (d, len(frontier)))
        nxt = set()
        for v in sorted(frontier):
            if v in seen: continue
            seen.add(v)
            print(" %s in %d gates:" % (desc(v, val), len(occ[v])))
            if len(occ[v]) > 8:
                print("     (hub, skipped)"); continue
            for f in occ[v]:
                vs = set()
                for m in facs[f]: vs.update(m)
                print("    f%-6d [%s] %s  ->%d" % (f, s(facs[f]), ' '.join(desc(u, val) for u in sorted(vs)), fv(facs[f], val)))
                nxt |= vs
        frontier = nxt - seen

if __name__ == '__main__':
    val = load(sys.argv[1] if len(sys.argv) > 1 else
               '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
    roots = [int(a) for a in sys.argv[2].split(',')]
    trace(roots, val, int(sys.argv[3]) if len(sys.argv) > 3 else 3)
