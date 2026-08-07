#!/usr/bin/env python3
"""Decode every stage input slot: find, for each stage input wire w, the residual atom that ties w to a
source expression, and unfold that source into its gated `selector * value` terms.

Handles are recognised structurally: Z = wires that are == 0 (mod p) for EVERY assignment (closure over
the definition DAG seeded by the literal p).  Any term containing a Z-wire drops out mod p, so what is
left in the atom is  c1*w + c2*z  and z is the slot's source.
"""
import sys, os, json, re, collections, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E = Engine()
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
defmap = {E.cls[a][1]: a for a in E.order}
defrhs = {v: E.cls[a][2] for a, v in ((a, E.cls[a][1]) for a in E.order)}

# ---- assignment-independent Z closure ----
Z = set()
def isZ(n):
    o = n[0]
    if o == 'c': return n[1] % p == 0
    if o == 'v': return n[1] in Z
    if o == 'neg': return isZ(n[1])
    if o == '*': return isZ(n[1]) or isZ(n[2])
    return isZ(n[1]) and isZ(n[2])
ch = True
while ch:
    ch = False
    for v, r in defrhs.items():
        if v not in Z and isZ(r): Z.add(v); ch = True

def lf(n):
    """(dict var->coef, const) mod p, or None if genuinely nonlinear."""
    o = n[0]
    if o == 'v':
        return ({}, 0) if n[1] in Z else ({n[1]: 1}, 0)
    if o == 'c': return ({}, n[1] % p)
    if o == 'neg':
        r = lf(n[1]); return None if r is None else ({k: (-x) % p for k, x in r[0].items()}, (-r[1]) % p)
    if o in '+-':
        a = lf(n[1]); b = lf(n[2])
        if a is None or b is None: return None
        s = 1 if o == '+' else -1
        d = dict(a[0])
        for k, x in b[0].items():
            d[k] = (d.get(k, 0) + s * x) % p
            if d[k] == 0: del d[k]
        return (d, (a[1] + s * b[1]) % p)
    a = lf(n[1]); b = lf(n[2])
    if a is None or b is None: return None
    if not a[0]: return ({k: a[1] * x % p for k, x in b[0].items()} if a[1] else {}, a[1] * b[1] % p)
    if not b[0]: return ({k: b[1] * x % p for k, x in a[0].items()} if b[1] else {}, b[1] * a[1] % p)
    return None

resby = collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)


def source_of(w):
    """Return (z, coef) with  w == coef*z (mod p) forced by a residual atom, or a tag."""
    cands = []
    for a in resby.get(w, []):
        f = lf(E.atoms[a])
        if f is None: continue
        d, c = f
        if w not in d: continue
        if len(d) == 2 and c == 0:
            z = [k for k in d if k != w][0]
            cands.append((z, (-d[z]) % p * pow(d[w], p - 2, p) % p))
        elif len(d) == 1 and c:
            cands.append(('CONST', (-c) % p * pow(d[w], p - 2, p) % p))
    return cands


def mux_terms(z):
    """Unfold z's definition into gated terms."""
    terms = []
    def walk(n):
        if n[0] in ('+', '-'): walk(n[1]); walk(n[2]); return
        if n[0] == 'v':
            u = n[1]
            if u in Z: return
            r = defrhs.get(u)
            if r is None: terms.append(('free', u)); return
            if r[0] == '*': terms.append(('gated', node_str(r))); return
            if r[0] in ('+', '-'): walk(r); return
            terms.append(('other', node_str(r)[:50])); return
        if n[0] == 'c' and n[1] % p: terms.append(('const', n[1]))
    r = defrhs.get(z)
    if r is None: return [('free', z)]
    walk(r)
    return terms


if __name__ == '__main__':
    roles = json.load(open(os.path.join(HERE, 'stage_roles.json')))
    print('|Z| (wires == 0 mod p for every assignment):', len(Z))
    stat = collections.Counter(); shapes = collections.Counter(); wired = {}
    for g, rs in roles.items():
        r = rs[0]
        for slot in ('inA', 'inB'):
            for w in r[slot]:
                c = source_of(w)
                if not c:
                    stat['no_source'] += 1
                    for a in resby.get(w, []):
                        shapes[re.sub(r'\d{2,}', 'C', re.sub(r'x\d+', 'X', a))[:70]] += 1
                    continue
                z = c[0][0]
                if z == 'CONST':
                    stat['literal'] += 1; wired.setdefault(g, {}).setdefault(slot, []).append(('CONST', str(c[0][1])))
                    continue
                t = mux_terms(z)
                stat['terms=%d' % len(t)] += 1
                wired.setdefault(g, {}).setdefault(slot, []).append((z, t))
    print('slot-wire decode:', stat.most_common())
    print('unresolved atom shapes:')
    for k, v in shapes.most_common(8): print('   %4d  %s' % (v, k))
    json.dump({k: {s: [[str(z), [list(map(str, x)) for x in t]] for z, t in v] for s, v in d.items()}
               for k, d in wired.items()}, open(os.path.join(HERE, 'mux_wiring.json'), 'w'))
    print('wrote mux_wiring.json for %d stages' % len(wired))
