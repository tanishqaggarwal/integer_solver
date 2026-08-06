"""Exact integer optimisation over the 9-atom / 8-knob lattice with the two extra
congruence-breaking knobs x_9118 and x_8731.  Maximise the number of the 13 equations that vanish."""
import pickle, itertools, sys, os
from fractions import Fraction
HERE=os.path.dirname(os.path.abspath(__file__)); S9=os.path.dirname(HERE)
sys.path.insert(0,S9); os.chdir(S9)
import harness as H
P = 2**256-2**32-977
d = pickle.load(open('atoms.pkl','rb')); eq_terms = d['eq_terms']
v = H.load_assignment('../best/new_instance_partial_39022.json')
D1 = v[7068]-v[2099]-7376877*v[642]
D2 = v[4432]-v[19964]
ATOMS = [22229,35762,22230,22231,35758,35759,35760,35761,37887]   # y1..y9
S13 = [2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125,9123,18673]

# y = M z + c ,  z = (t, h, s, k, y5, y6, y7, y8)
#   y1 = D1 - 7376877 t ; y2 = t - p h ; y3 = s - p k ; y4 = D2 - s ; y5..y8 free ; y9 = y4^2
NZ = 8
def ycol(j):        # coefficient vector of z_j in (y1..y8)
    m=[0]*8
    if j==0: m[0]=-7376877; m[1]=1
    if j==1: m[1]=-P
    if j==2: m[2]=1; m[3]=-1
    if j==3: m[2]=-P
    if j>=4: m[j]=1
    return m
CONST=[D1,0,0,D2,0,0,0,0]

rows=[]
for i in S13:
    m,sq,tl=eq_terms[i]
    co={a:0 for a in ATOMS}
    for c,a in tl:
        if a in co: co[a]+=c
    rows.append((i,sq,co))
    others=[a for c,a in tl if a not in co]
print('13 equations, atom coefficients over the 9-atom set:')
for i,sq,co in rows: print(f'  eq {i} sq={sq}: ' + ', '.join(f'{a}:{c}' for a,c in co.items() if c))

from snf import solve_int
def _unused(A, b):
    """Integer solvability of A z = b (A: list of rows, ints). Returns a solution or None."""
    import copy
    n = len(A[0]); m = len(A)
    M = [row[:] + [b[k]] for k, row in enumerate(A)]
    piv = []
    r = 0
    for c in range(n):
        # find pivot with nonzero entry, use gcd elimination
        while True:
            nz = [k for k in range(r, m) if M[k][c]]
            if len(nz) <= 1: break
            nz.sort(key=lambda k: abs(M[k][c]))
            k0 = nz[0]
            for k in nz[1:]:
                q = M[k][c] // M[k0][c]
                for j in range(c, n+1): M[k][j] -= q*M[k0][j]
        nz = [k for k in range(r, m) if M[k][c]]
        if not nz: continue
        k0 = nz[0]
        M[r], M[k0] = M[k0], M[r]
        piv.append((r, c)); r += 1
        if r == m: break
    for k in range(r, m):
        if M[k][n] != 0: return None
    # back-substitute over the pivot rows (upper-triangular in pivot columns)
    z = [0]*n
    for (rr, cc) in reversed(piv):
        s = M[rr][n] - sum(M[rr][j]*z[j] for j in range(cc+1, n))
        if s % M[rr][cc]: return None
        z[cc] = s // M[rr][cc]
    for k in range(m):
        if sum(A[k][j]*z[j] for j in range(n)) != b[k]: return None
    return z

def try_subset(T):
    """Can all equations in T vanish simultaneously?  Returns z or None."""
    A=[]; b=[]
    for i in T:
        sq = dict((x,(y,z)) for x,y,z in [])  # placeholder
    for i in T:
        _,sq,co = next(r for r in rows if r[0]==i)
        if co.get(37887,0):
            # equation value is m*(sum)^2 with sum containing y9=y4^2 -> needs y4 = 0
            A.append([1 if j==2 else 0 for j in range(NZ)]); b.append(D2)   # s = D2 => y4 = 0
            co = {a:c for a,c in co.items() if a != 37887}
            if not any(co.values()): continue
        row=[0]*NZ; rhs=0
        for a,c in co.items():
            if not c: continue
            k=ATOMS.index(a)
            if k>=8: continue
            for j in range(NZ): row[j]+=c*ycol(j)[k]
            rhs-=c*CONST[k]
        if any(row): A.append(row); b.append(rhs)
    if not A: return [0]*NZ
    return solve_int(A,b)

# validate: the known 39,024 witness zeroes {2554,6816,8124,8680}
known=[2554,6816,8124,8680]
print('\nvalidation: known 4-subset solvable?', try_subset(known) is not None)
print('validation: 5-subsets containing it:',
      sum(1 for e in S13 if e not in known and try_subset(known+[e]) is not None))

best=None
for size in range(len(S13),0,-1):
    found=None
    for T in itertools.combinations(S13,size):
        z=try_subset(list(T))
        if z is not None: found=(T,z); break
    if found:
        best=found; print(f'\nMAX simultaneously-zeroable equations = {size}: {found[0]}')
        print(f'   z = (t,h,s,k,y5,y6,y7,y8) = {found[1]}')
        print(f'   => failing equations = {len(S13)} - {size} = {len(S13)-size}  -> score {39033-(len(S13)-size)}')
        break
pickle.dump(best, open('kernel/opt6.pkl','wb'))
