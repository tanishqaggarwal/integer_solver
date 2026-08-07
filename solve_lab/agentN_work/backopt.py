"""BACKWARDS step 4: exact integer optimum over the realizable knob lattice, and the price of
   the atoms the lattice cannot reach."""
import frameB as FB, ev, json, time
from frameB import Frame, State
from collections import defaultdict
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
D=json.load(open('backreal.json')); R=D['R']; S=D['S']
KN=sorted(int(k) for k in D['knobs'])
fr=Frame([642,28730,29854,31864])
W=json.load(open('../best/new_instance_partial_39026.json'))
v=[0]*38748
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
base=State(fr,{u:v[u] for u in fr.free if v[u]!=0})
# inner sums of the region's equations, as linear functions of the 9 knob offsets
def inner(st,e):
    m,sq,tl=ev.eq_terms[e]; t=0
    for c,a in tl:
        x=st.av.get(a)
        if x: t+=c*x
    return t
b=[inner(base,e) for e in R]
A=[]
for X in KN:
    g=base.clone().set_free({X:base.fv.get(X,0)+1})
    A.append([inner(g,e)-b[i] for i,e in enumerate(R)])
# verify linearity of the whole map
g=base.clone()
for X in KN: g.set_free({X:base.fv.get(X,0)+1})
lin=all(inner(g,e)==b[i]+sum(A[j][i] for j in range(len(KN))) for i,e in enumerate(R))
print('region rows %d, knobs %d, jointly linear: %s'%(len(R),len(KN),lin))
print('rows currently zero:',sum(1 for x in b if x==0),'of',len(R))

def solve_int(rows):
    """solve  sum_j t_j A[j][i] = -b[i] for i in rows, over Z.  Return a solution or None."""
    n=len(KN); M=[[A[j][i] for j in range(n)]+[-b[i]] for i in rows]
    # integer Gaussian elimination (fraction-free), then back-substitute requiring integrality
    piv=[]; r=0
    for c in range(n):
        k=None
        for i in range(r,len(M)):
            if M[i][c]: k=i; break
        if k is None: continue
        M[r],M[k]=M[k],M[r]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f1=M[r][c]; f2=M[i][c]
                M[i]=[f1*M[i][j]-f2*M[r][j] for j in range(n+1)]
        piv.append(c); r+=1
        if r==len(M): break
    for i in range(r,len(M)):
        if all(M[i][j]==0 for j in range(n)) and M[i][n]!=0: return None
    t=[0]*n
    for i,c in enumerate(piv):
        num=M[i][n]; den=M[i][c]
        if den==0: return None
        if num%den: return None
        t[c]=num//den
    # verify
    for i in rows:
        if b[i]+sum(t[j]*A[j][i] for j in range(n))!=0: return None
    return t

best=(sum(1 for x in b if x==0),None)
t0=time.time(); tested=0
import itertools
nrows=len(R)
for size in range(min(9,nrows),0,-1):
    if size<=best[0]: break
    for rows in itertools.combinations(range(nrows),size):
        tested+=1
        t=solve_int(list(rows))
        if t is None: continue
        z=sum(1 for i in range(nrows) if b[i]+sum(t[j]*A[j][i] for j in range(len(KN)))==0)
        if z>best[0]: best=(z,t); print('  size %d -> %d rows zeroed'%(size,z),flush=True)
    print(' size %d done, tested %d, best %d, %.0fs'%(size,tested,best[0],time.time()-t0),flush=True)
print('EXACT integer optimum over the realizable lattice: %d of %d rows zeroed -> failing %d in region'%(best[0],nrows,nrows-best[0]))
print('total failing would be', nrows-best[0], '(region) + 0 outside  =>  score', 39033-(nrows-best[0]))
if best[1]:
    st=base.clone()
    for j,X in enumerate(KN):
        if best[1][j]: st.set_free({X:base.fv.get(X,0)+best[1][j]})
    print('CONSTRUCTED state score:',st.score(),'failing',len(st.fails))
    if st.score()>39026:
        out='H_%d_backwards.json'%st.score()
        json.dump({('x_%d'%i):st.v[i] for i in range(38748) if st.v[i]!=0},open(out,'w'))
        print('WROTE',out)
