import json, sympy
from sympy import Matrix, Rational, zeros
p=2**256-2**32-977
D=json.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/dioph.json'))
A=[[int(x) for x in row] for row in D['A']]
b0=[int(x) for x in D['b0']]
AFF=D['AFF']; F=D['F']
m=len(F); n=len(AFF)
print(f"system {m} eqs x {n} unknowns")
Am=Matrix(A); bm=Matrix([-b0[k] for k in range(m)])  # want A delta = -b0
# rank over Q
r=Am.rank(); print("rank_Q(A)=",r)
Aug=Am.row_join(bm); ra=Aug.rank(); print("rank_Q([A|b])=",ra)
if r==ra:
    print("=> CONSISTENT over Q")
    sol,params=Am.gauss_jordan_solve(bm)
    print("particular rational solution delta:")
    dvec=[sol[i] for i in range(n)]
    for i in range(n):
        val=dvec[i]
        print(f"  d(x_{AFF[i]}) = {val}  {'INTEGER' if val.is_integer else 'FRACTION den='+str(val.q)}")
    print("num free params:",params.shape)
else:
    print("=> INCONSISTENT over Q (overdetermined). Checking mod p.")
# mod p rank/consistency
Amp=Matrix([[A[k][j]%p for j in range(n)] for k in range(m)])
bmp=Matrix([(-b0[k])%p for k in range(m)])
def rank_modp(M):
    M=M.copy(); rows,cols=M.shape; r=0
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if M[i,c]%p!=0: piv=i;break
        if piv is None: continue
        M.row_swap(r,piv); inv=pow(int(M[r,c]),p-2,p)
        for j in range(cols): M[r,j]=(M[r,j]*inv)%p
        for i in range(rows):
            if i!=r and M[i,c]%p!=0:
                f=M[i,c]
                for j in range(cols): M[i,j]=(M[i,j]-f*M[r,j])%p
        r+=1
    return r,M
rp,_=rank_modp(Amp)
rap,_=rank_modp(Amp.row_join(bmp))
print(f"\nmod p: rank(A)={rp}, rank([A|b])={rap} -> {'CONSISTENT mod p' if rp==rap else 'INCONSISTENT mod p'}")
