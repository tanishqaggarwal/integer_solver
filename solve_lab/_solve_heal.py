import pickle
p=2**256-2**32-977
J=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/jacmat.pkl','rb'))
M=J['M']; rhs=J['rhs']; frees=J['frees']; nfrag=J['nfrag']
nrow=len(M); ncol=len(frees)
# augmented [M|rhs]; track which original row is S (nfrag) and T (nfrag+1)
A=[row[:]+[rhs[r]] for r,row in enumerate(M)]
# Gaussian elimination mod p
def inv(x): return pow(x%p,p-2,p)
pivcol=[]
r=0
rowsrc=list(range(nrow))  # track provenance
for c in range(ncol):
    # find pivot in column c at or below r
    piv=-1
    for rr in range(r,nrow):
        if A[rr][c]%p!=0: piv=rr; break
    if piv<0: continue
    A[r],A[piv]=A[piv],A[r]
    rowsrc[r],rowsrc[piv]=rowsrc[piv],rowsrc[r]
    iv=inv(A[r][c])
    A[r]=[(x*iv)%p for x in A[r]]
    for rr in range(nrow):
        if rr!=r and A[rr][c]%p!=0:
            f=A[rr][c]
            A[rr]=[(A[rr][k]-f*A[r][k])%p for k in range(ncol+1)]
    pivcol.append(c); r+=1
    if r==nrow: break
rank=r
# check inconsistency: any row all-zero in cols but nonzero rhs
inconsistent=[]
for rr in range(nrow):
    if all(A[rr][c]%p==0 for c in range(ncol)) and A[rr][ncol]%p!=0:
        inconsistent.append(rr)
print(f"rank={rank}, nrow={nrow}, ncol={ncol}")
print(f"inconsistent rows: {len(inconsistent)}")
# Is the FULL system consistent?
print("HEAL SYSTEM CONSISTENT (partial guards):", len(inconsistent)==0)
# Now the key: even ignoring RHS, are S,T rows independent of guard rows?
# rank of guards only vs guards+S+T
import copy
def rankof(rows):
    B=[row[:] for row in rows]; nc=len(B[0]); rr=0
    for c in range(nc):
        piv=-1
        for k in range(rr,len(B)):
            if B[k][c]%p!=0: piv=k;break
        if piv<0: continue
        B[rr],B[piv]=B[piv],B[rr]
        ivv=inv(B[rr][c]); B[rr]=[(x*ivv)%p for x in B[rr]]
        for k in range(len(B)):
            if k!=rr and B[k][c]%p!=0:
                f=B[k][c]; B[k]=[(B[k][j]-f*B[rr][j])%p for j in range(nc)]
        rr+=1
        if rr==len(B): break
    return rr
guards=[M[r] for r in range(nfrag)]
rg=rankof(guards)
rgS=rankof(guards+[M[nfrag]])
rgT=rankof(guards+[M[nfrag+1]])
rgST=rankof(guards+[M[nfrag],M[nfrag+1]])
print(f"\nrank(guards)={rg}")
print(f"rank(guards+S)={rgS}  -> S {'INDEPENDENT (movable given these guards)' if rgS>rg else 'in rowspace (S LOCKED by these guards)'}")
print(f"rank(guards+T)={rgT}  -> T {'INDEPENDENT (movable)' if rgT>rg else 'in rowspace (T LOCKED)'}")
print(f"rank(guards+S+T)={rgST}  (guards rank + {rgST-rg} core dims free)")
