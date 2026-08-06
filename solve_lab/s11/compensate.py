"""Equation-space compensation: an atom a may be NONZERO if some other atom b appears in
   exactly the same equations with a proportional coefficient vector (then b = -k*a cancels it).
   Search for such 'shadow' atoms for each candidate absorber."""
import sys, os, json
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
CAND = [14445, 27139, 34580, 33796, 26719, 26721, 26723, 21050, 7881, 26839, 40065]

for a in CAND:
    ea = L.atom2eq.get(a, {})
    S = set(ea)
    shadows = []
    for b in range(L.NA):
        if b == a:
            continue
        eb = L.atom2eq.get(b, {})
        if not eb:
            continue
        if not S.issubset(set(eb)):
            continue
        # proportional on S ?
        rat = None
        ok = True
        for e in S:
            r = Fraction(eb[e], ea[e])
            if rat is None:
                rat = r
            elif r != rat:
                ok = False
                break
        if ok:
            extra = len(set(eb) - S)
            shadows.append((b, rat, extra, L.atom_out.get(b) is None))
    shadows.sort(key=lambda t: t[2])
    print(f"a{a} [{len(S)} eqs]: {len(shadows)} proportional shadows"
          f" -> {[(b, str(r), x, 'chk' if c else 'gate') for b, r, x, c in shadows[:6]]}")
