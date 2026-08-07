import sys, json, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close import closure, evalat, CHK
from ort import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BASE={542:1, 91:1, 22162:K2, 30213:K1}
def st(extra,verbose=False):
    s=dict(BASE); s.update(extra)
    return closure(s,verbose=verbose)
seeds,v,sc,nz=st({})
x1=v[12186]%P; y2=v[24908]%P; x3=v[22162]%P; y3=v[30213]%P
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
c=(x1-x3)%P; T=(x3+2*x1+K)%P; M=(y3+y2)%P
# w^3 - 3M w^2 + (3M^2 + c^3 - c^2 T) w - M^3 + c^2 T M = 0
co=[1,(-3*M)%P,(3*M*M+pow(c,3,P)-c*c%P*T)%P,(-pow(M,3,P)+c*c%P*T%P*M)%P]
print('cubic coeffs',co)
# roots mod p: gcd(f, x^p-x)
def pmul(a,b,f):
    r=[0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        if ai:
            for j,bj in enumerate(b): r[i+j]=(r[i+j]+ai*bj)%P
    return pmod(r,f)
def pmod(a,f):
    a=a[:]; df=len(f)-1
    inv=pow(f[0],P-2,P)
    while len(a)-1>=df and len(a)>1:
        if a[0]==0: a.pop(0); continue
        k=len(a)-1-df; co2=a[0]*inv%P
        for i in range(len(f)): a[i+0]=(a[i]-co2*f[i])%P
        a.pop(0)
    while len(a)>1 and a[0]==0: a.pop(0)
    return a
def ppow(base,e,f):
    r=[1]; b=base[:]
    while e:
        if e&1: r=pmul(r,b,f)
        b=pmul(b,b,f); e>>=1
    return r
def pgcd(a,b):
    while len(b)>1 or (len(b)==1 and b[0]!=0):
        a=pmod(a,b) if len(a)>=len(b) else a
        if len(a)>=len(b) and not(len(a)==1 and a[0]==0):
            a,b=b,a
        else: a,b=b,a
        if len(b)==1 and b[0]==0: break
    return a
f=co[:]
xp=ppow([1,0],P,f)
g=[(xp+[0]*(2-len(xp)))[i] if False else 0 for i in range(0)]
# compute xp - x
d=xp[:]
while len(d)<2: d=[0]+d
d[-2]=(d[-2]-1)%P
while len(d)>1 and d[0]==0: d.pop(0)
# gcd(f, d)
def gcdp(a,b):
    a=a[:]; b=b[:]
    while not(len(b)==1 and b[0]==0):
        r=pmod(a,b); a,b=b,r
    inv=pow(a[0],P-2,P) if a[0] else 1
    return [x*inv%P for x in a]
g=gcdp(f,d)
print('gcd degree',len(g)-1)
def evalp(f,x):
    s=0
    for co_ in f: s=(s*x+co_)%P
    return s
roots=[]
if len(g)-1==1: roots=[(-g[1])%P]
elif len(g)-1>=2:
    # brute: try random splitting
    for _ in range(200):
        a=random.randrange(P)
        h=ppow([1,a],(P-1)//2,g)
        h=h[:]; 
        while len(h)<2: h=[0]+h
        h[-1]=(h[-1]-1)%P
        while len(h)>1 and h[0]==0: h.pop(0)
        gg=gcdp(g,h)
        if 1<=len(gg)-1<len(g)-1:
            # recurse crudely
            for part in [gg,pmod(g,gg) if False else None]:
                pass
            # just collect linear factors by repeated splitting
            stack=[g]; out=[]
            while stack:
                q=stack.pop()
                if len(q)-1==1: out.append((-q[1]*pow(q[0],P-2,P))%P); continue
                for _ in range(100):
                    a2=random.randrange(P)
                    h2=ppow([1,a2],(P-1)//2,q)
                    while len(h2)<2: h2=[0]+h2
                    h2[-1]=(h2[-1]-1)%P
                    while len(h2)>1 and h2[0]==0: h2.pop(0)
                    q1=gcdp(q,h2)
                    if 1<=len(q1)-1<len(q)-1:
                        # divide
                        num=q[:]; q2=[]
                        # polynomial division q / q1
                        aa=q[:]; bb=q1[:]; quot=[]
                        while len(aa)>=len(bb) and not(len(aa)==1 and aa[0]==0):
                            co3=aa[0]*pow(bb[0],P-2,P)%P
                            quot.append(co3)
                            for i in range(len(bb)): aa[i]=(aa[i]-co3*bb[i])%P
                            aa.pop(0)
                            while len(aa)>1 and aa[0]==0: aa.pop(0); quot.append(0)
                        stack.append(q1); stack.append(quot if quot else [1])
                        break
                else: break
            roots=out; break
print('roots of cubic:',len(roots))
cands=[]
# degenerate w=0
cands.append(('w=0 (P1=P2)',y2,x1))
for w in roots:
    s=(M-w)%P
    if s==0: continue
    y1=(y2-w)%P; x2=(x1+c*w%P*pow(s,P-2,P))%P
    cands.append(('cubic w=%d'%(w%10**8),y1,x2))
print('candidates',len(cands))
res=[]
for tag,y1,x2 in cands:
    seeds,v,sc,nz=st({16742:y1,14853:x2})
    print('[%s] score=%d nz=%d  A=%d B=%d'%(tag,sc,len(nz),v[35389]%P,v[6671]%P))
    res.append((sc,tag,y1,x2))
    if sc>=39013:
        json.dump({f'x_{i}':v[i] for i in range(L.NVARS) if v[i]!=0},
                  open('/home/user/integer_solver/solve_lab/agentC_work/cand_%d.json'%sc,'w'))
print(sorted(res,reverse=True)[:3])
