"""Is the mod-p barrier ambient — i.e. independent of collateral and of the lattice entirely?

`fr.csup` is the EXACT syntactic free-input support of an atom, so the 68 candidates are provably
every free input of Frame(POOL) that can change any atom of the region: no other free input can move
a region equation at all.  So if a region row's exact polynomial over those 68 knobs has every
non-constant coefficient divisible by p while its constant is not, that row is nonzero mod p for
EVERY integer point of Z^68 — the row cannot be zeroed by any assignment whatsoever, with or without
collateral, and independently of the lattice.

That is a claim of the form "nothing can move X", so the knob set and the configuration are stated
with it and verified by direct recomputation at random points.
"""
import os, sys, json, random, pickle
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import ev, optN
from optN import make, build, WIT, fr, FREE, FR0, atom_eqs, _bits, inner
from polyexact import P
from polyfull import exact_polys, evalP

Pp = 115792089237316195423570985008687907853269984665640564039457584007908834671663

st = make(WIT)
b0 = build(st)
Rl = b0['R']
atoms_R = set()
for e in Rl:
    for c, a in ev.eq_terms[e][2]:
        atoms_R.add(a)
cands = set()
for q in atoms_R:
    if q in fr.csup:
        cands.update(FR0[bb] for bb in _bits(fr.csup[q]))
cands = sorted(y for y in cands if y in FREE)
k = len(cands)

print('CONFIGURATION', flush=True)
print('  frame        : Frame(POOL), 65 pool variables detached', flush=True)
print('  detach set D : %s (the witness set)' % WIT, flush=True)
print('  selectors    : best/new_instance_partial_39026.json', flush=True)
print('  knob set     : ALL %d free inputs of the frame that syntactically support any atom of any'
      % k, flush=True)
print('                 region equation — i.e. every free input that can move the region at all',
      flush=True)
print('  region       : |R| = %d  %s' % (len(Rl), Rl), flush=True)

polys = exact_polys(st, Rl, cands)

print('\nrow-by-row reduction of the EXACT polynomial mod p:', flush=True)
print('%-8s %-5s %-7s %-14s %-14s %s' %
      ('eq', 'deg', 'terms', 'const mod p', 'nonconst mod p', 'verdict'), flush=True)
stuck = []
for e in Rl:
    pol = polys[e]
    c0 = pol.c.get((0,) * k, 0)
    nonc = [v for m, v in pol.c.items() if sum(m)]
    allp = all(v % Pp == 0 for v in nonc)
    verdict = ''
    if allp and c0 % Pp:
        verdict = 'UNZEROABLE mod p — no integer knob vector can zero it'
        stuck.append(e)
    elif allp and c0 % Pp == 0:
        verdict = 'entire row == 0 mod p'
    print('%-8d %-5d %-7d %-14s %-14s %s'
          % (e, pol.deg(), pol.nterms(), '0' if c0 % Pp == 0 else 'nonzero',
             'all 0' if allp else '%d/%d nonzero' % (sum(1 for v in nonc if v % Pp), len(nonc)),
             verdict), flush=True)

print('\nrows unzeroable mod p over the COMPLETE knob set: %d -> %s' % (len(stuck), stuck),
      flush=True)

# --- direct recomputation, never trust the expansion ---------------------------------------
print('\ndirect recomputation at random integer points (all %d knobs moved at once):' % k,
      flush=True)
rnd = random.Random(7)
bad = 0
for tr in range(12):
    t = [rnd.randint(-10 ** 6, 10 ** 6) for _ in cands]
    h = st.clone().set_free({Y: st.fv.get(Y, 0) + t[j] for j, Y in enumerate(cands)})
    for e in Rl:
        a = inner(h, e)
        b = evalP(polys[e], t)
        if a != b:
            bad += 1
    for e in stuck:
        if inner(h, e) % Pp == 0:
            print('   COUNTEREXAMPLE: eq %d hit 0 mod p' % e, flush=True)
            bad += 1
print('   %d evaluations, %d mismatches/counterexamples' % (12 * (len(Rl) + len(stuck)), bad),
      flush=True)

# --- how many equations does the mod-p barrier account for? ---------------------------------
print('\nblast radius: the deliverable fails on %d equations; this barrier accounts for %d of them'
      % (7, len(stuck)), flush=True)
print('   failing set: [12231, 12270, 12350, 14584, 18673, 22044, 29125]', flush=True)
print('   stuck set  : %s' % stuck, flush=True)

json.dump(dict(knobs=k, region=Rl, stuck=stuck, p=str(Pp), mismatches=bad),
          open(os.path.join(HERE, 'runs', 'pambient.json'), 'w'), indent=1)
print('\nwrote runs/pambient.json', flush=True)
