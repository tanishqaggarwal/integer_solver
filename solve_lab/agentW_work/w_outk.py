"""W stage 19:
 (A) THE FOUR OUT-OF-K BLOCKS: price an injection at 3227/4429/30886/31606 exactly.
 (B) A GAP IN MY OWN ROUND-1 KNOB SET: K was 'free inputs reaching a NONZERO atom'.  To repair
     a broken equation you may instead move a DIFFERENT atom in that equation to compensate --
     and those movers were treated as constants.  Measure K+ = all free inputs touching any
     atom of the 7 failing equations."""
import sys, os, json, itertools, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB
PVAL = 115792089237316195423570985008687907853269984665640564039457584007908834671663
A = frameB.atom_src; EQT = frameB.eq_terms
B = {b['E']: b for b in json.load(open('w_blocks4.json'))}
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v0 = [0]*frameB.NV
for k, val in W.items(): v0[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fr = frameB.Frame([642, 28730, 29854, 31864])
fv = {u: v0[u] for u in fr.free if v0[u] != 0}
st0 = frameB.State(fr, fv)
assert st0.score() == 39026
V = 9470933787112530972301743374550701798831        # the deliverable's own injected magnitude

print('=== (A) injecting at the four equation-disjoint minimum-incidence blocks ===')
print('    control  E=7181: the deliverable already injects there; site damage = 5 of its 9')
for E in (3227, 4429, 30886, 31606):
    b = B[E]
    free_slots = [b['i%d' % k] for k in (5, 6) if b['i%d' % k] in fr.free]
    best = None
    for comb in ([free_slots[0]] if free_slots else []) + \
                ([free_slots[1]] if len(free_slots) > 1 else []) + \
                ([free_slots] if len(free_slots) > 1 else []):
        comb = comb if isinstance(comb, list) else [comb]
        s = st0.clone().set_free({u: V for u in comb})
        new = sorted(set(s.fails) - set(st0.fails))
        if best is None or len(new) < best[0]: best = (len(new), comb, s.score(), new)
    print('  E=%-6d free output slots %s -> cheapest injection breaks %d NEW equations %s (score %d)'
          % (E, free_slots, best[0], best[3], best[2]))

print()
print('=== (B) K  vs  K+ : who can move the atoms of the 7 failing equations? ===')
K = set(json.load(open('w_K.json'))['K'])
FAILEQ = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
Kplus = set(); atoms = set()
for e in FAILEQ:
    m, sq, tl = EQT[e]
    for c, a in tl:
        atoms.add(a)
        Kplus |= set(fr.SUPV.get(a, []))
print('  atoms in the 7 failing equations: %d' % len(atoms))
print('  |K| (round 1)      = %d' % len(K))
print('  |K+| (all movers)  = %d' % len(Kplus))
NEW = sorted(Kplus - K)
print('  free inputs OUTSIDE round-1 K that can move an atom of a failing equation: %d' % len(NEW))
print('  ', NEW[:40])
print()
print('  screening each new mover singly (delta = +1, +2, and the injected magnitude):')
hits = []
for u in NEW:
    o = fv.get(u, 0)
    for dv in (o + 1, o + 2, o - 1, V, 0):
        s = st0.clone().set_free({u: dv})
        if s.score() > 39026:
            hits.append((u, dv, s.score(), sorted(s.fails))); print('   *** IMPROVE', u, s.score())
        elif s.score() == 39026 and sorted(s.fails) != FAILEQ:
            hits.append((u, dv, s.score(), sorted(s.fails)))
print('  single-knob improvements:', [h for h in hits if h[2] > 39026])
print('  single-knob lateral moves (same score, different failing set):', len([h for h in hits if h[2] == 39026]))
json.dump({'K': sorted(K), 'Kplus': sorted(Kplus), 'NEW': NEW}, open('w_outk.json', 'w'))
