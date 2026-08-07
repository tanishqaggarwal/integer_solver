import json
P=2**256-2**32-977
C=json.load(open('/home/user/integer_solver/solve_lab/agentC_work/curve.json'))
a2=int(C['KA']); a4=int(C['a4']); a6=int(C['a6'])
N=115792089237316195423570985008687907852837564279074904382605163141518161494337
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
def neg(Pt): return None if Pt is None else (Pt[0],(-Pt[1])%P)
def mul(k,Pt):
    if k==0 or Pt is None: return None
    ng=k<0; k=abs(k)%N
    R=None; Q=Pt
    while k:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    return neg(R) if ng else R
def leafpoints():
    d=json.load(open('/home/user/integer_solver/solve_lab/agentC_work/leafpts2.json'))
    return {int(k):(int(v[0]),int(v[1])) for k,v in d.items()}
