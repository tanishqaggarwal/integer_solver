"""Exact integer linear solve M x = b via column-style Hermite reduction."""
def solve_int(M, b):
    """M: list of rows (lists), b: list. Return integer x or None."""
    m=len(M); n=len(M[0]) if m else 0
    A=[row[:] for row in M]
    U=[[1 if i==j else 0 for j in range(n)] for i in range(n)]   # columns ops
    rhs=b[:]
    piv=[]
    r=0
    for c in range(m):            # each row of A in turn -> pivot column
        # find column >= r with nonzero entry in row c
        while True:
            nz=[j for j in range(r,n) if A[c][j]!=0]
            if not nz: break
            if len(nz)==1: break
            # euclidean reduce among nz columns
            nz.sort(key=lambda j: abs(A[c][j]))
            j0=nz[0]
            done=True
            for j in nz[1:]:
                q=A[c][j]//A[c][j0]
                if q:
                    for i in range(m): A[i][j]-=q*A[i][j0]
                    for i in range(n): U[i][j]-=q*U[i][j0]
                    done=False
            if done: break
        nz=[j for j in range(r,n) if A[c][j]!=0]
        if not nz: 
            piv.append(None); continue
        j0=nz[0]
        if j0!=r:
            for i in range(m): A[i][r],A[i][j0]=A[i][j0],A[i][r]
            for i in range(n): U[i][r],U[i][j0]=U[i][j0],U[i][r]
        piv.append(r); r+=1
    # now A is lower-triangular-ish: row c has pivot at column piv[c]
    y=[0]*n
    for c in range(m):
        s=rhs[c]-sum(A[c][j]*y[j] for j in range(n))
        j=piv[c]
        if j is None:
            if s!=0: return None
            continue
        if A[c][j]==0:
            if s!=0: return None
            continue
        if s % A[c][j]: return None
        y[j]=s//A[c][j]
    # verify
    x=[sum(U[i][j]*y[j] for j in range(n)) for i in range(n)]
    for c in range(m):
        if sum(M[c][i]*x[i] for i in range(n))!=b[c]: return None
    return x
if __name__=='__main__':
    print(solve_int([[2,3],[4,5]],[8,14]))
    print(solve_int([[2,4]],[7]))
    print(solve_int([[2,4]],[6]))
