"""bl_cong: the two binding congruences, and which load pins feed them.

C0 = x_7068 - x_2099   (atom a22229: C0 - 7376877*x_642)
a7930 = 9367949*(x_24548 - x_25442) - x_7927
For each of the 727 conditional constant loads, decide whether its pinned wire is
an ANCESTOR of C0 / of a7930 / of a21617,a29539.
"""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, CANON, F2, pot, FORBID
P = 2**256-2**32-977

w = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json')); F2.fwd(w)
v0 = L.load(os.path.join(HERE,'mod9118_0.json')); CANON.fwd(v0)
LD = json.load(open(os.path.join(HERE,'bl_pins3.json')))['loads']

def vcone(F, seeds):
    c, st = set(), list(seeds)
    while st:
        t = st.pop()
        if t in c: continue
        c.add(t)
        d = F.definer.get(t)
        if d is None: continue
        for x in L.avars[d]:
            if x != t: st.append(x)
    return c

CONES = {
 'C0_canon'  : vcone(CANON, [7068, 2099]),
 'C0_f2'     : vcone(F2,    [7068, 2099]),
 'a7930'     : vcone(CANON, sorted(L.avars[7930])),
 'a21617'    : vcone(CANON, sorted(L.avars[21617])),
 'a29539'    : vcone(CANON, sorted(L.avars[29539])),
 'seven_f2'  : vcone(F2, sorted(set().union(*[L.avars[a] for a in
                 (22229,22230,35758,35759,35760,35761,35762)]))),
 'a37662'    : vcone(CANON, sorted(L.avars[37662])),
 'a40826'    : vcone(CANON, sorted(L.avars[40826])),
}
for k, c in CONES.items():
    print(f'cone {k:<10}: {len(c)} vars, {len(c & BOOL)} booleans, '
          f'{len(c & BOOL & CANON.FREE)} free booleans')

print('\nfor each cone: which load-pin GATES sit inside it (flipping the gate '
      'unpins / pins a constant that the cone depends on)')
for k, c in CONES.items():
    hits = collections.defaultdict(set)
    for a, b, K, pr, n in LD:
        for x in pr:
            if x in c: hits[b].add((a, x))
    print(f'\n  {k}: {len(hits)} gate booleans')
    for b in sorted(hits):
        st = 'ON ' if w[b] == 1 else 'off'
        print(f'    x_{b:<6} [{st}] -> {sorted(hits[b])[:6]}')

print('\n--- current congruence residues ---')
print(f'  C0 = x_7068-x_2099  (F2 witness) mod p = {(w[7068]-w[2099])%P}')
print(f'  x_642 mod p = {w[642]%P}   7376877*x_642 mod p = {7376877*w[642]%P}')
print(f'  a22229 mod p = {L.evalpoly(L.polys[22229], w)%P}')
print(f'  a7930 value (canon) = {L.evalpoly(L.polys[7930], v0)}')
print(f'  a7930 value (F2 witness) = {L.evalpoly(L.polys[7930], w)}')
for t in (24548, 25442, 7927):
    print(f'   x_{t}: canonFREE={t in CANON.FREE} definer={CANON.definer.get(t)} '
          f'src={L.atom_src[CANON.definer[t]][:80] if t in CANON.definer else "FREE"}')
