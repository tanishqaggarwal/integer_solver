import sys, math, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
N=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
x1=v[12186]%P; y1=v[16742]%P; b1=(y1*y1-pow(x1,3,P))%P
def add(A,B):
    if A is None: return B
    if B is None: return A
    ax,ay=A; bx,by=B
    if ax==bx and (ay+by)%P==0: return None
    lam = 3*ax*ax%P*pow(2*ay,P-2,P)%P if A==B else (by-ay)*pow(bx-ax,P-2,P)%P
    cx=(lam*lam-ax-bx)%P; return (cx,(lam*(ax-cx)-ay)%P)
def mul(k,A):
    R=None;Q=A
    while k:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    return R
t=P+1-N
f2=(4*P-t*t)
assert f2%3==0
f=math.isqrt(f2//3); assert f*f==f2//3
print('t_secp=%d  f=%d  (t^2-4p=-3f^2 check %s)'%(t,f,t*t-4*P==-3*f*f))
traces=set()
for s in (1,-1):
    traces.add(s*t)
    if (s*t+3*f)%2==0: traces.add((s*t+3*f)//2)
    if (s*t-3*f)%2==0: traces.add((s*t-3*f)//2)
G=(x1,y1)
import sympy
print('the six j=0 orders and whether (x1,y1) has that order:')
for tt in sorted(traces):
    m=P+1-tt
    r=mul(m,G)
    print('  t=%-80d order=%-80d  m*G=O? %-6s prime? %s'%(tt,m,r is None,sympy.isprime(m)))
# pin-constant test against b1
lits=set()
for a in range(L.NA):
    for m_,c in L.polys[a].items():
        if abs(c)>=10**40: lits.add(abs(c)%P)
lits=sorted(lits)
print('\n%d distinct large-literal residues'%len(lits))
onx=[X for X in lits if pow((pow(X,3,P)+b1)%P,(P-1)//2,P)==1]
print('literals X with X^3+b1 a QR (valid x-coord on y^2=x^3+b1): %d  (chance ~%.0f)'%(len(onx),len(lits)/2))
S=set(lits)
pairs=[(X,Y) for X in lits for Y in (None,) if False]
cnt=0; ex=[]
sq={}
for Y in lits: sq.setdefault(pow(Y,2,P),[]).append(Y)
for X in lits:
    key=(pow(X,3,P)+b1)%P
    if key in sq:
        cnt+=1; ex.append((X,sq[key][0]))
print('literal pairs (X,Y) with Y^2 = X^3 + b1 : %d   examples %s'%(cnt,ex[:3]))
# also: do the literals lie on the SAME curve as (x3,y3)?
x3=v[22162]%P; y3=v[30213]%P; b3=(y3*y3-pow(x3,3,P))%P
print('b3 = %d ; b3==b1 ? %s'%(b3,b3==b1))
for name,bb in [('b1',b1),('b3',b3),('7',7)]:
    c2=0
    for X in lits:
        key=(pow(X,3,P)+bb)%P
        if key in sq: c2+=1
    print('   pairs on y^2=x^3+%s : %d'%(name,c2))
