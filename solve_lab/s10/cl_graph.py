"""CL: the global congruence graph.  Every gadget asserts A == B (mod p); pins assert
A == CONST (mod p).  Build the graph, propagate, and locate the contradiction."""
import os, sys, json, re, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0, nz0, S0, bad0 = E.stats(v0)

# gadget shapes
re_cc  = re.compile(r'^(\d+) \* \(x_(\d+) - x_(\d+)\) - x_(\d+)$')          # c*(A-B) - H
re_sel = re.compile(r'^x_(\d+) \* \(x_(\d+) - (\d+)\)(?: - (?:(\d+) \* )?x_(\d+))?$')  # sel*(A-C) - H
re_scc = re.compile(r'^x_(\d+) \* \(x_(\d+) - x_(\d+)\)(?: - (?:(\d+) \* )?x_(\d+))?$')

gad, pin, other = [], [], []
for a in range(L.NA):
    if a in atom_out: continue
    s = L.atom_src[a].strip()
    m = re_cc.match(s)
    if m:
        gad.append((a, int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)))); continue
    m = re_sel.match(s)
    if m:
        pin.append((a, int(m.group(1)), int(m.group(2)), int(m.group(3)), m.group(5))); continue
    m = re_scc.match(s)
    if m:
        gad.append((a, ('sel', int(m.group(1))), int(m.group(2)), int(m.group(3)), m.group(5))); continue

print(f'gadgets  c*(A-B)-H : {len([g for g in gad if not isinstance(g[1],tuple)])}')
print(f'gadgets  s*(A-B)-H : {len([g for g in gad if isinstance(g[1],tuple)])}')
print(f'pins     s*(A-C)-H : {len(pin)}')

# ---- resolve each var to its residue mod p
def r(u): return v0[u] % P

print('\n=== PINS (variable forced to a constant mod p) ===')
consts = {}
for a, sel, u, C, h in pin:
    ok = av0[a] == 0
    live = v0[sel] % P != 0
    print(f'  a{a:<6} eqs={len(L.atom2eq.get(a,{})):<3} sel=x_{sel}({v0[sel]}) x_{u} == {str(C)[:24]}..  '
          f'({C.bit_length()}b)  sat={ok} live={live}  x_{u} mod p == C mod p ? {r(u)==C%P}')
    if live: consts.setdefault(u, []).append((C % P, a))

print('\n=== GADGETS: A == B (mod p) ===')
uf = {}
def find(x):
    while uf.get(x, x) != x: uf[x] = uf.get(uf[x], uf[x]); x = uf[x]
    return x
def union(x, y):
    rx, ry = find(x), find(y)
    if rx != ry: uf[rx] = ry
edges = []
nsat = nfail = 0
for g in gad:
    a, c, A, B, h = g
    sel = None
    if isinstance(c, tuple):
        sel = c[1]
        if v0[sel] % P == 0:
            continue                      # gadget is switched off
    ok = av0[a] == 0
    cong = r(A) == r(B)
    edges.append((a, A, B, ok, cong))
    if ok: nsat += 1
    else: nfail += 1
    if cong: union(A, B)
print(f'  live gadgets: {len(edges)}  satisfied {nsat}  failing {nfail}')
for a, A, B, ok, cong in edges:
    if not ok or not cong:
        print(f'  a{a:<6} x_{A} == x_{B} ?  sat={ok} congruent={cong}  eqs={len(L.atom2eq.get(a,{}))}')

# classes
cls = collections.defaultdict(list)
for a, A, B, ok, cong in edges:
    cls[find(A)].append(A); cls[find(B)].append(B)
print(f'\ncongruence classes among gadget endpoints: {len(cls)}')
big = sorted(cls.items(), key=lambda kv: -len(set(kv[1])))
for k, mem in big[:12]:
    ms = sorted(set(mem))
    print(f'  class rep x_{k}: {len(ms)} members, residue {r(ms[0])==r(ms[-1])}, '
          f'value mod p = {str(r(ms[0]))[:26]}..  members {ms[:12]}')

print('\n=== the failing gadgets and their class neighbourhood ===')
for a, A, B, ok, cong in edges:
    if ok: continue
    print(f'\na{a}: x_{A} (mod p {r(A)}) vs x_{B} (mod p {r(B)})')
    for X in (A, B):
        cl = sorted(set(cls.get(find(X), [])))
        print(f'   x_{X} class: {len(cl)} members {cl[:14]}')
        for u in cl:
            if u in consts:
                print(f'      PINNED x_{u} == {consts[u][0][0]} (a{consts[u][0][1]})')
json.dump({'pins': [[a, sel, u, str(C)] for a, sel, u, C, h in pin],
           'gadgets': [[e[0], e[1], e[2], e[3], e[4]] for e in edges]},
          open(os.path.join(HERE,'cl_graph.json'),'w'))
