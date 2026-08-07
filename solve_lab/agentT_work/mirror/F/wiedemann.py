#!/usr/bin/env python3
"""Independent cross-check of singularity of M by Wiedemann over GF(q)."""
import numpy as np, scipy.sparse as sp, os, sys, time
HERE=os.path.dirname(os.path.abspath(__file__))
M=sp.load_npz(os.path.join(HERE,'M.npz')).tocsr()
n=M.shape[0]

def bm(seq,q):
    """Berlekamp-Massey over GF(q). Returns C (numpy int64 array, C[0]=1)."""
    C=np.zeros(len(seq)+1,dtype=np.int64); C[0]=1
    B=np.zeros(len(seq)+1,dtype=np.int64); B[0]=1
    L=0; mm=1; b=1
    for N in range(len(seq)):
        # discrepancy d = sum_{i=0..L} C[i]*seq[N-i]
        d=int((C[:L+1]*seq[N-L:N+1][::-1] % q).sum() % q)
        if d==0:
            mm+=1
        elif 2*L<=N:
            T=C.copy()
            coef=d*pow(int(b),q-2,q)%q
            C[mm:mm+len(B)-mm]=(C[mm:]-coef*B[:len(B)-mm])%q
            L=N+1-L; B=T; b=d; mm=1
        else:
            coef=d*pow(int(b),q-2,q)%q
            C[mm:]=(C[mm:]-coef*B[:len(B)-mm])%q
            mm+=1
    return C[:L+1],L

def run(q,seed):
    rng=np.random.default_rng(seed)
    v=rng.integers(1,q,size=n).astype(np.int64)
    u=rng.integers(1,q,size=n).astype(np.int64)
    seq=np.zeros(2*n,dtype=np.int64)
    x=v.copy(); t0=time.time()
    for i in range(2*n):
        seq[i]=int((u*x % q).sum() % q)
        x=M.dot(x)%q
        if i%20000==0: print('   krylov',i,'%.0fs'%(time.time()-t0),flush=True)
    print('  sequence done %.0fs'%(time.time()-t0),flush=True)
    t1=time.time()
    C,L=bm(seq,q)
    print('  BM done %.0fs  minpoly degree L=%d'%(time.time()-t1,L),flush=True)
    # C[0]=1 is the leading coeff of the reversed poly; the annihilator is
    # sum_i C[i] a_{k-i}; A is singular iff x | minpoly, i.e. C[L]==0
    print('  q=%d  deg=%d (n=%d)  trailing coeff C[L]=%d  =>  M %s over GF(q)'%(
        q,L,n,int(C[L]), 'SINGULAR' if int(C[L])%q==0 else 'NONSINGULAR'),flush=True)
    return L,int(C[L])

if __name__=='__main__':
    for q,seed in [(2147483647,1),(2147483629,2)]:
        print('=== q =',q,flush=True)
        run(q,seed)
