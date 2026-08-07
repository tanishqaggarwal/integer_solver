import sys, json, math
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
P=2**256-2**32-977
C=json.load(open('/home/user/integer_solver/solve_lab/agentC_work/curve.json'))
a2=int(C['KA']); a4=int(C['a4']); a6=int(C['a6'])
# short Weierstrass: X = x + a2/3
inv3=pow(3,P-2,P)
A=(a4-a2*a2%P*inv3)%P
B=(2*pow(a2,3,P)*pow(inv3,3,P)-a4*a2%P*inv3+a6)%P
print('short form: A =',A,' B =',B)
assert A==0, 'not j=0'
n_secp=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
t=P+1-n_secp
print('t_secp =',t)
# 4p = L^2+27M^2 ; for secp256k1 t = L (up to sign)
L=t
M2=(4*P-L*L)
print('(4p - t^2) % 27 ==0 ?', M2%27==0)
M=math.isqrt(M2//27)
print('M ok?',M*M*27+L*L==4*P, 'M =',M)
cands=set()
for s in (1,-1):
    cands.add(s*L)
    cands.add(s*(L+9*M)//2 if (L+9*M)%2==0 else None)
    cands.add(s*(L-9*M)//2 if (L-9*M)%2==0 else None)
cands={c for c in cands if c is not None}
# also with 3M
for s in (1,-1):
    for e in [(L+9*M),(L-9*M),(-L+9*M),(-L-9*M)]:
        if e%2==0: cands.add(s*e//2)
orders=sorted(set(P+1-c for c in cands))
print('candidate orders:',len(orders))
def add(Pt,Qt):
    if Pt is None: return Qt
    if Qt is None: return Pt
    x1,y1=Pt; x2,y2=Qt
    if (x1-x2)%P==0:
        if (y1+y2)%P==0: return None
        lam=(3*x1*x1+2*a2*x1+a4)*pow(2*y1,P-2,P)%P
    else:
        lam=(y2-y1)*pow(x2-x1,P-2,P)%P
    x3=(lam*lam-a2-x1-x2)%P; y3=(lam*(x1-x3)-y1)%P
    return (x3,y3)
def mul(k,Pt):
    if k==0: return None
    neg=k<0; k=abs(k)
    R=None; Q=Pt
    while k:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    if neg and R: R=(R[0],(-R[1])%P)
    return R
G=(int(C['P1'][0]),int(C['P1'][1]))
Q=(int(C['Q'][0]),int(C['Q'][1]))
good=[]
for N in orders:
    if N<=0: continue
    if mul(N,G) is None: good.append(N); print('ORDER FOUND:',N)
print('n_secp works?', mul(n_secp,G) is None)
json.dump({'A':str(A),'B':str(B),'orders':[str(o) for o in good]},open('/home/user/integer_solver/solve_lab/agentC_work/order.json','w'))
