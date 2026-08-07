"""Reconstruct the doubling chain and identify each leaf point as 2^i * G."""
import sys, json, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L, sympy
P=env.P
NSEC=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
GX=0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY=0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
c=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/curve_final.json'))
B=int(c['B']); sh=int(c['shift']); a2=int(c['a2']); a4=int(c['a4']); a6=int(c['a6'])
d=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/weier.json'))
raw=[(int(x)%P,int(y)%P) for x,y in d['points']]
onmask=[(pow(Y,2,P)-pow(X,3,P)-a2*pow(X,2,P)-a4*X-a6)%P==0 for X,Y in raw]
on=[((X+sh)%P,Y) for (X,Y),m in zip(raw,onmask) if m]
u=12830242018875522506555146473674089970775060590290859819641972374662130570109
def sec2ours(pt): return (pow(u,2,P)*pt[0]%P, pow(u,3,P)*pt[1]%P)
def add(A_,B_):
    if A_ is None: return B_
    if B_ is None: return A_
    ax,ay=A_; bx,by=B_
    if ax==bx and (ay+by)%P==0: return None
    lam=(3*ax*ax%P)*pow(2*ay,-1,P)%P if A_==B_ else ((by-ay)*pow(bx-ax,-1,P))%P
    cx=(lam*lam-ax-bx)%P; return (cx,(lam*(ax-cx)-ay)%P)
# map 2^i * G for i=0..300 into our model, index by x-coordinate
idx={}
Q=sec2ours((GX,GY))
for i in range(300):
    idx.setdefault(Q[0], i)
    Q=add(Q,Q)
found=collections.Counter(); which={}
for k,(X,Y) in enumerate(on):
    if X in idx:
        found[idx[X]]+=1; which[idx[X]]=(X,Y)
print('leaf points identified as 2^i * G : %d of %d'%(len(found),len(on)))
ii=sorted(found)
print('exponents i present: min=%d max=%d count=%d'%(min(ii),max(ii),len(ii)))
missing=[i for i in range(min(ii),max(ii)+1) if i not in found]
print('gaps in 2^i coverage:',missing[:20],'(%d total)'%len(missing))
# the 34 off-curve leaf points -- what are they?
off=[(X,Y) for (X,Y),m in zip(raw,onmask) if not m]
print('\noff-curve leaf pins: %d'%len(off))
# maybe they are on the curve with x and y swapped, or are y-only / scalar constants
sw=sum(1 for X,Y in off if (pow(X,2,P)-pow(Y,3,P)-a2*pow(Y,2,P)-a4*Y-a6)%P==0)
print('   on the curve with the pair SWAPPED: %d'%sw)
