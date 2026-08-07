"""Exact integer solve of A d = b, plus integer kernel basis."""
from flint import fmpz_mat
def hnf_t(A):
    M=fmpz_mat(A)
    try:
        H,U=M.hnf(transform=True)
    except TypeError:
        H=M.hnf(); U=None
    return H.tolist(), (U.tolist() if U is not None else None)

def solve_int(A, b):
    """A: m x n list of ints; b: m ints.  Return (d, kernel_basis) or (None, kernel_basis)."""
    m=len(A); n=len(A[0]) if m else 0
    AT=[[A[i][j] for i in range(m)] for j in range(n)]   # n x m
    H,U=hnf_t(AT)                                        # U (n x n) * AT = H (n x m)
    # pivots
    piv=[]
    for j in range(n):
        row=H[j]
        p=None
        for i in range(m):
            if row[i]!=0: p=i; break
        piv.append(p)
    z=[0]*n
    rhs=list(b)
    # forward substitution: process rows j in order of pivot
    order=sorted([j for j in range(n) if piv[j] is not None], key=lambda j: piv[j])
    for j in order:
        i=piv[j]
        cur=rhs[i]
        if cur % H[j][i]!=0:
            return None, [[int(y) for y in U[j]] for j in range(n) if piv[j] is None]
        c=cur//H[j][i]
        z[j]=c
        if c:
            for k in range(m):
                rhs[k]-=c*H[j][k]
    if any(rhs): return None, [[int(y) for y in U[j]] for j in range(n) if piv[j] is None]
    # d = U^T z  =>  d_k = sum_j z_j * U[j][k]
    d=[0]*n
    for j in range(n):
        if z[j]:
            uj=U[j]
            for k in range(n): d[k]+=z[j]*uj[k]
    ker=[U[j] for j in range(n) if piv[j] is None]
    return [int(x) for x in d], [[int(y) for y in r] for r in ker]
if __name__=='__main__':
    A=[[2,1,0],[0,3,1]]; b=[4,6]
    d,k=solve_int(A,b); print(d,k)
    print([sum(A[i][j]*d[j] for j in range(3)) for i in range(2)])
