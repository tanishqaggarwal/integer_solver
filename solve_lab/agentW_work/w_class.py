"""W stage 10: THE CLASSIFICATION.  Solve  N1 = E*A^2 - B^2 = 0  and  N2 = A*H - B*J = 0
(a) as polynomial identities (sympy),
(b) EXHAUSTIVELY over small fields F_p, every Q, every (i1..i6)  -- validates that the case
    analysis has no missing case, mechanically,
(c) over Z, where the extra integrality A | B appears."""
import sys, itertools, json
sys.set_int_max_str_digits(20_000_000)
import sympy as sp
from collections import Counter

i1,i2,i3,i4,i5,i6,Q = sp.symbols('i1 i2 i3 i4 i5 i6 Q')
Av = i1-i2; Bv = i4-i3; Ev = i1+i2+i5+Q; Hv = i3+i6; Jv = i2-i5
N1 = sp.expand(Ev*Av**2 - Bv**2); N2 = sp.expand(Av*Hv - Bv*Jv)

print('=== (a) polynomial identities ===')
print('N1 at A=0 :', sp.simplify(N1.subs(i1, i2)), '   (= -B^2)')
print('N2 at A=0,B=0 :', sp.simplify(N2.subs({i1: i2, i4: i3})))
Jac = sp.Matrix([[sp.diff(N1,i5), sp.diff(N1,i6)],[sp.diff(N2,i5), sp.diff(N2,i6)]])
print('Jacobian d(N1,N2)/d(i5,i6) =', Jac.tolist(), ' det =', sp.factor(Jac.det()))
lam = sp.Symbol('lam')
i5s = lam**2 - i1 - i2 - Q
i6s = lam*(i2 - i5s) - i3
chk1 = sp.simplify(N1.subs({i5: i5s, i4: i3 + lam*(i1-i2)}))
chk2 = sp.simplify(N2.subs({i5: i5s, i6: i6s, i4: i3 + lam*(i1-i2)}))
print('chord substitution B=lam*A, i5,i6 as claimed:  N1 ->', chk1, ' N2 ->', chk2)
# rank of the coordinate map (i1..i6) -> (A,B,E,H,J)
M = sp.Matrix([[sp.diff(f, v) for v in (i1,i2,i3,i4,i5,i6)] for f in (Av,Bv,Ev,Hv,Jv)])
print('coordinate map matrix rank over Q:', M.rank(), ' Smith-ish dets of 5x5 minors:',
      sorted({abs(M[:, c].det()) for c in itertools.combinations(range(6),5)}))

print()
print('=== (b) EXHAUSTIVE over small fields: every p, every Q, every (i1..i6) ===')
rows = []
for p in (5, 7, 11, 13):
    tot = Counter()
    for q in range(p):
        for a1 in range(p):
         for a2 in range(p):
          A_ = (a1-a2) % p
          for a3 in range(p):
           for a4 in range(p):
            B_ = (a4-a3) % p
            for a5 in range(p):
             E_ = (a1+a2+a5+q) % p
             n1 = (E_*A_*A_ - B_*B_) % p
             if n1: continue
             for a6 in range(p):
              n2 = (A_*(a3+a6) - B_*(a2-a5)) % p
              if n2: continue
              if A_ == 0 and B_ == 0: tot['degeneracy A=B=0'] += 1
              elif A_ == 0: tot['A=0,B!=0  (IMPOSSIBLE?)'] += 1
              else:
                  lam_ = B_ * pow(A_, p-2, p) % p
                  okE = (E_ - lam_*lam_) % p == 0
                  okH = ((a3+a6) - lam_*(a2-a5)) % p == 0
                  tot['chord (lam=B/A)' if (okE and okH) else 'THIRD FAMILY'] += 1
    rows.append((p, dict(tot)))
    print('  p=%2d  solutions over all Q and all (i1..i6):' % p, dict(tot), flush=True)

print()
print('=== (c) over Z: the integrality refinement ===')
import math
R = 24; cnt = Counter(); ex = []
for A_ in range(-R, R+1):
    for B_ in range(-R, R+1):
        # E*A^2 = B^2 solvable in Z ?
        if A_ == 0:
            cnt['A=0 -> needs B=0: %s' % ('B=0' if B_ == 0 else 'no solution')] += 1
            continue
        if (B_*B_) % (A_*A_) == 0:
            E_ = (B_*B_)//(A_*A_)
            lam2 = E_
            div = (B_ % A_ == 0)
            cnt['A|B' if div else 'A^2|B^2 but A does NOT divide B'] += 1
            if not div: ex.append((A_, B_))
        else:
            cnt['no integer E'] += 1
print('  box |A|,|B| <= %d :' % R, dict(cnt))
print('  counterexamples to  A^2|B^2 => A|B :', ex[:5], '(count %d)' % len(ex))
json.dump({'small_field': rows}, open('w_class.json', 'w'))
