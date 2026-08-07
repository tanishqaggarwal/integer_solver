"""S11 step 71: DECOMPILE the seven residual constraints.

Printed as polynomials the residual atoms are tiny:

    a22229  x7068 - x2099 - 7376877*x642          a22230  x28730 - x9413*x17499
    a35758  x29854 - x1329*x22665                 a35759  5113045*x7075*x9118 - x29854
    a35760  x31864 - x10903*x28961                a35761  x7075*x8731 + x31864
    a35762  x642 - x17325*x28599

Each is one multiplication or one sum, and they come in PAIRS that compute the same
wire two ways -- classic R1CS redundancy.  So the residual is not an exotic object;
it is a handful of products that must agree.  Whether that is reachable depends on
what the operands are, which means walking the definition chain back to free inputs.

Usage: decompile.py [DEPTH]
"""
import os, sys, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
DEPTH = int(sys.argv[1]) if len(sys.argv) > 1 else 6
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
FREE = set(t for t in range(L.NVARS) if t not in L.definer)


def poly(a):
    ts = L.polys[a]
    it = ts.items() if isinstance(ts, dict) else ts
    out = []
    for mono, c in it:
        if isinstance(mono, tuple) and mono:
            out.append('%+d*%s' % (c, '*'.join('x%d' % m for m in mono)))
        else:
            out.append('%+d' % c)
    return ' '.join(out)


def kind(t):
    if t in FREE:
        return 'FREE'
    return 'a%d' % L.definer[t]


SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
seed = set()
for a in SEVEN:
    seed |= set(L.avars[a])
print('=== residual atoms ===')
for a in SEVEN:
    print('  a%-6d = %s' % (a, poly(a)))
print()
print('=== operand provenance, %d levels ===' % DEPTH)
seen, front = set(), sorted(seed)
for lv in range(DEPTH):
    nxt = set()
    print('--- level %d (%d vars) ---' % (lv, len(front)))
    nfree = 0
    for t in front:
        if t in seen:
            continue
        seen.add(t)
        if t in FREE:
            nfree += 1
            if lv < 3:
                print('   x%-6d FREE   value bits %d' % (t, v[t].bit_length()))
            continue
        a = L.definer[t]
        if lv < 3:
            print('   x%-6d <- a%-6d : %s' % (t, a, poly(a)[:150]))
        for w in L.avars[a]:
            if w != t and w not in seen:
                nxt.add(w)
    if lv >= 3:
        print('   ... %d vars, %d of them free' % (len(front), nfree))
    front = sorted(nxt)
    if not front:
        break
print()
print('total vars reached: %d, free among them: %d'
      % (len(seen), sum(1 for t in seen if t in FREE)))
