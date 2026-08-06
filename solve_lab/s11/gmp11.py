"""The three coupled rows: a7930, a21617, a33796.

Each is 'free mirror minus computed value' with its own free input (x24548, x14623, x31339),
so naively they are three independent knobs.  Every obstruction certificate contains all three,
which says the responses are linearly dependent.  Look at the actual matrix.
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
from gmp7 import allchk, failing, resp_at
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
bd=allchk(base)
R=[7930,21617,33796]
K=[24548,14623,31339]
print("row residues:", [('FAIL' if bd[a] else 'ok') for a in R])
for u in K:
    d=resp_at(base,bd,u)
    print(f"x{u}: touches {sorted(d)}")
    print(f"    on R: {[d.get(a,0)%P for a in R]}")
M=[[resp_at(base,bd,u).get(a,0)%P for u in K] for a in R]
def rank(M,n):
    A=[r[:] for r in M]; m=len(A); r=0
    for c in range(n):
        pr=next((i for i in range(r,m) if A[i][c]),None)
        if pr is None: continue
        A[r],A[pr]=A[pr],A[r]
        inv=pow(A[r][c],-1,P); A[r]=[x*inv%P for x in A[r]]
        for i in range(m):
            if i!=r and A[i][c]:
                f=A[i][c]; A[i]=[(A[i][k]-f*A[r][k])%P for k in range(n)]
        r+=1
    return r
print("3x3 rank:",rank(M,3))
for i,a in enumerate(R): print(f"  row a{a}: {[str(x)[:14]+'..' for x in M[i]]}")
