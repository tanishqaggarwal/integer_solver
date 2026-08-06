"""Two affine knobs, three congruences -- is the target in their span mod p?"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','fix7_29539_7930.json'))
T=[11150,25739,37758]
c=[v[t]%P for t in T]
K=[22162,30213]
R=[]
for u in K:
    v2=list(v); L.ripple(v2,{u:v[u]+1})
    R.append([(v2[t]-v[t])%P for t in T])
print("target  c =",[str(x)[:20]+'..' for x in c])
for u,col in zip(K,R): print(f"col x{u} =",[str(x)[:20]+'..' for x in col])
# solve R^T . delta = -c  (3 eqs, 2 unknowns) over GF(p)
A=[[R[0][i],R[1][i],(-c[i])%P] for i in range(3)]
r=0
for col in range(2):
    pv=next((i for i in range(r,3) if A[i][col]),None)
    if pv is None: continue
    A[r],A[pv]=A[pv],A[r]
    inv=pow(A[r][col],-1,P); A[r]=[x*inv%P for x in A[r]]
    for i in range(3):
        if i!=r and A[i][col]:
            f=A[i][col]; A[i]=[(A[i][k]-f*A[r][k])%P for k in range(3)]
    r+=1
bad=[i for i in range(3) if not A[i][0] and not A[i][1] and A[i][2]]
print("rank",r,"  consistent:",not bad)
if not bad:
    d=[0,0]
    # read off pivots
    rr=0
    for col in range(2):
        if rr<r and A[rr][col]==1:
            d[col]=A[rr][2]; rr+=1
    print("delta =",d)
    L.ripple(v,{K[0]:v[K[0]]+d[0], K[1]:v[K[1]]+d[1]})
    print("  x11150%p==0:",v[11150]%P==0," x25739%p==0:",v[25739]%P==0," x37758%p==0:",v[37758]%P==0)
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    B=[a for a in range(L.NA) if AV[a]!=0]
    print("  broken:",B," failing:",len(L.failing_eqs(AV)))
    json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','modp3_out.json'),'w'))
