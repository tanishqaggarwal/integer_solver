import sys, json, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
GX=0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY=0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
c=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/curve_final.json'))
sh=int(c['shift']); a2=int(c['a2']); a4=int(c['a4']); a6=int(c['a6'])
d=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/weier.json'))
raw=[(int(x)%P,int(y)%P) for x,y in d['points']]
u=12830242018875522506555146473674089970775060590290859819641972374662130570109
def add(A_,B_):
    if A_ is None: return B_
    if B_ is None: return A_
    ax,ay=A_; bx,by=B_
    if ax==bx and (ay+by)%P==0: return None
    lam=(3*ax*ax%P)*pow(2*ay,-1,P)%P if A_==B_ else ((by-ay)*pow(bx-ax,-1,P))%P
    cx=(lam*lam-ax-bx)%P; return (cx,(lam*(ax-cx)-ay)%P)
idx={}; Q=(pow(u,2,P)*GX%P, pow(u,3,P)*GY%P)
for i in range(260):
    idx.setdefault(Q[0],i); Q=add(Q,Q)
got={}
for k,(X,Y) in enumerate(raw):
    for (xx,yy) in ((X,Y),(Y,X)):
        if (pow(yy,2,P)-pow(xx,3,P)-a2*pow(xx,2,P)-a4*xx-a6)%P==0:
            sx=(xx+sh)%P
            if sx in idx: got[k]=idx[sx]
            break
print('leaf pins identified as 2^i*G (allowing either pin order): %d of %d'%(len(got),len(raw)))
ex=sorted(got.values())
print('exponents: %d distinct, min %d max %d, complete 0..255 ? %s'%(
      len(set(ex)),min(ex),max(ex),sorted(set(ex))==list(range(256))))
