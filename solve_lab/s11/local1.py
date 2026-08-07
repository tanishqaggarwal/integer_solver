"""EXACT reduction of the checkpoint's neighbourhood.

Nine variables are FULLY LOCAL (they occur in no atom outside the 15 equations that touch
the broken atoms):  x1844 x21574 x642 x9413 x1329 x29854 x10903 x31864 x17325.
Holding everything else fixed, the nine broken/at-risk atom values are

    z0 = a35756 = x1844 - p*x21574        z1 = a35757 = S - x1844          (S = x1956*x17065)
    y3 = a35758 = x29854 - p*x1329        y4 = a35759 = T - x29854         (T = 5113045*x7075*x9118)
    y5 = a35760 = x31864 - p*x10903       y6 = a35761 = U + x31864         (U = x7075*x8731)
    y7 = a35762 = x642 - p*x17325         y1 = a22229 = A - 7376877*x642   (A = x7068 - x2099)
    y2 = a22230 = B - p*x9413                                              (B = x28730)

Eliminating the nine knobs, the REACHABLE SET of (z0,z1,y3..y7,y1,y2) is exactly the coset

    z0 + z1 == S            (mod p)
    y3 + y4 == T            (mod p)
    y6 - y5 == U            (mod p)
    y2      == B            (mod p)
    y1 + 7376877*y7 == A    (mod 7376877*p)

and nothing else -- every other integer combination is free.  The 15 equations are linear
forms in these nine values.  So "how many equations must fail here" is an exact, tiny
integer feasibility question, and we can enumerate the answer.
"""
import sys, os, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
from zsolve import solve_int
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)

src = sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
v = load_raw(src)
AV=[atomval(a,v) for a in range(L.NA)]
KEY=[35756,35757,35758,35759,35760,35761,35762,22229,22230]   # -> z0 z1 y3 y4 y5 y6 y7 y1 y2
NK=len(KEY)
EQS=sorted(set().union(*[set(L.atom2eq.get(a,{})) for a in KEY]))
print("equations touched:",len(EQS))
# constant part = contribution of every OTHER atom (must be unchanged)
rows=[]; const=[]
for e in EQS:
    m,sq,co=L.eq_atoms[e]
    assert not sq
    rows.append([co.get(a,0) for a in KEY])
    const.append(sum(c*AV[a] for a,c in co.items() if a not in KEY))
print("constant (non-KEY) contributions all zero:", all(c==0 for c in const))
cur=[AV[a] for a in KEY]
now=[i for i,e in enumerate(EQS) if sum(rows[i][j]*cur[j] for j in range(NK))+const[i]!=0]
print("currently failing among them:",len(now),[EQS[i] for i in now])

S = v[1956]*v[17065]; T = 5113045*v[7075]*v[9118]; U = v[7075]*v[8731]
A = v[7068]-v[2099];  B = v[28730]; Q = 7376877
# congruence rows:  coeff vector on KEY, extra p-multiple unknown, rhs
CONG=[([1,1,0,0,0,0,0,0,0], P,       S),
      ([0,0,1,1,0,0,0,0,0], P,       T),
      ([0,0,0,0,-1,1,0,0,0],P,       U),
      ([0,0,0,0,0,0,0,0,1], P,       B),
      ([0,0,0,0,0,0,Q,1,0], Q*P,     A)]
# sanity: the current point satisfies them
for cv,mod,rhs in CONG:
    assert (sum(cv[j]*cur[j] for j in range(NK))-rhs)%mod==0, cv
print("current point lies in the reachable coset: OK")

NT=len(CONG)
def feasible(keep):
    """integer system: kept rows =0 ; congruence rows with slack unknowns"""
    M=[]; r=[]
    for i in keep:
        M.append(rows[i]+[0]*NT); r.append(-const[i])
    for t,(cv,mod,rhs) in enumerate(CONG):
        M.append(cv+[(-mod if k==t else 0) for k in range(NT)]); r.append(rhs)
    return solve_int(M,r)

t0=time.time()
best=None
for nd in range(0,len(EQS)+1):
    hit=None
    for D in itertools.combinations(range(len(EQS)), nd):
        keep=[i for i in range(len(EQS)) if i not in D]
        x=feasible(keep)
        if x is not None:
            hit=(D,x); break
    print(f"  drop {nd}: {'FEASIBLE  drop-set '+str([EQS[i] for i in hit[0]]) if hit else 'infeasible'}"
          f"  ({time.time()-t0:.0f}s)", flush=True)
    if hit:
        best=hit; break
D,x=best
print()
print("MINIMUM failures in this neighbourhood =",len(D)," (currently 7)")
print("y* =",x[:NK])
import json
json.dump({'EQS':EQS,'drop':[EQS[i] for i in D],'y':[int(t) for t in x[:NK]],'KEY':KEY},
          open(os.path.join(HERE,'data','local1.json'),'w'))
