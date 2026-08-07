"""Is the 256-point gated set a double-and-add ladder?  Test 2*P_i membership, and
compare against the secp256k1 generator."""
import sys, json, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
NSEC=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
GX=0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY=0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
c=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/curve_final.json'))
A=int(c['A']); B=int(c['B']); sh=int(c['shift']); a2=int(c['a2']); a4=int(c['a4']); a6=int(c['a6'])
d=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/weier.json'))
pts=[(int(x)%P,int(y)%P) for x,y in d['points']]
on=[((X+sh)%P,Y) for X,Y in pts if (pow(Y,2,P)-pow(X,3,P)-a2*pow(X,2,P)-a4*X-a6)%P==0]
print('%d points on the curve (short form Y^2 = X^3 + %d)'%(len(on),B))
def add(Pt,Qt):
    if Pt is None: return Qt
    if Qt is None: return Pt
    ax,ay=Pt; bx,by=Qt
    if ax==bx and (ay+by)%P==0: return None
    lam=(3*ax*ax%P)*pow(2*ay,-1,P)%P if Pt==Qt else ((by-ay)*pow(bx-ax,-1,P))%P
    cx=(lam*lam-ax-bx)%P; return (cx,(lam*(ax-cx)-ay)%P)
def mul(k,Pt):
    R=None;Q=Pt
    while k:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    return R
S={p_[0] for p_ in on}          # by x-coordinate (sign-insensitive)
dbl=sum(1 for p_ in on if add(p_,p_) is not None and add(p_,p_)[0] in S)
print('points whose DOUBLE is also in the set (x-coord match): %d of %d'%(dbl,len(on)))
# sums of pairs
Sx=set(S); hits=0; ex=[]
for i in range(min(40,len(on))):
    for j in range(i+1,min(40,len(on))):
        r=add(on[i],on[j])
        if r is not None and r[0] in Sx:
            hits+=1
            if len(ex)<3: ex.append((i,j))
print('P_i + P_j landing back in the set (first 40x40 pairs): %d  %s'%(hits,ex))
# is the secp256k1 generator among them?  (need the iso: our curve is B, secp is 7)
import math
u6=B*pow(7,-1,P)%P
# find u with u^6 = B/7 : u = (B/7)^((2*(p-1)/3+1)/6)? do it by solving via a 6th root
def sixth_root(a):
    # p = 2 mod 3? p mod 3:
    for t in range(1,2):
        pass
    # brute: use that x -> x^6 is 6-to-1; solve by taking cube root then square root
    # cube root: p mod 9 handling; use sympy nthroot_mod
    import sympy
    rs=sympy.ntheory.residue_ntheory.nthroot_mod(a,6,P,all_roots=True)
    return rs
try:
    rs=sixth_root(u6)
    print('sixth roots of B/7 : %d found'%(len(rs) if rs else 0))
    if rs:
        u=rs[0]
        # iso: secp point (x,y) -> our curve (u^2 x, u^3 y)
        gx=(pow(u,2,P)*GX)%P; gy=(pow(u,3,P)*GY)%P
        print('   secp G mapped onto our curve: x=%d'%gx)
        print('   is it a leaf point ? %s'%(gx in Sx))
        # try all six roots
        for u in rs:
            gx=(pow(u,2,P)*GX)%P
            if gx in Sx: print('   *** MATCH with u=%d'%u)
except Exception as e:
    print('sixth-root failed:',e)
