"""Exact determination of the curve from the 256 gated leaf points, its short form,
its isomorphism class relative to secp256k1, and its group order."""
import sys, json, math, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L, sympy
P=env.P
NSEC=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
d=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/weier.json'))
pts=[(int(x)%P,int(y)%P) for x,y in d['points']]
a2,a4,a6=[int(x)%P for x in d['a']]
on=[(X,Y) for X,Y in pts if (pow(Y,2,P)-pow(X,3,P)-a2*pow(X,2,P)-a4*X-a6)%P==0]
print('leaf points: %d ; on the fitted curve: %d'%(len(pts),len(on)))
# refit using ALL inliers, exactly, and confirm the fit is over-determined & consistent
M=[[pow(X,2,P),X%P,1,(pow(Y,2,P)-pow(X,3,P))%P] for X,Y in on]
piv=[];r=0
for c in range(3):
    pr=None
    for i in range(r,len(M)):
        if M[i][c]%P: pr=i;break
    M[r],M[pr]=M[pr],M[r]
    inv=pow(M[r][c],-1,P); M[r]=[x*inv%P for x in M[r]]
    for i in range(len(M)):
        if i!=r and M[i][c]%P:
            f=M[i][c]; M[i]=[(u-f*w)%P for u,w in zip(M[i],M[r])]
    piv.append(c); r+=1
inc=sum(1 for i in range(r,len(M)) if M[i][3]%P)
print('refit over %d inliers: rank=%d inconsistent=%d  (0 => all lie on ONE curve)'%(len(on),r,inc))
print('a2 = %d'%a2); print('a4 = %d'%a4); print('a6 = %d'%a6)
i3=pow(3,-1,P)
A=(a4-a2*a2%P*i3)%P
B=(a6 + 2*pow(a2,3,P)%P*pow(27,-1,P) - a2*a4%P*i3)%P
print('\nshort form after X = x + a2/3 :  Y^2 = X^3 + A X + B')
print('A = %d   %s'%(A,'<== A == 0' if A==0 else ''))
print('B = %d'%B)
print('B == 7 ? %s'%(B==7))
for k in (2,3,6):
    e=(P-1)//math.gcd(P-1,k)
    print('B/7 is a %d-th power mod p : %s'%(k, pow(B*pow(7,-1,P)%P,e,P)==1))
# curve group order on the short form
def add(Pt,Qt,a=A):
    if Pt is None: return Qt
    if Qt is None: return Pt
    ax,ay=Pt; bx,by=Qt
    if ax==bx and (ay+by)%P==0: return None
    lam=((3*ax*ax+a)%P)*pow(2*ay,-1,P)%P if Pt==Qt else ((by-ay)*pow(bx-ax,-1,P))%P
    cx=(lam*lam-ax-bx)%P; return (cx,(lam*(ax-cx)-ay)%P)
def mul(k,Pt):
    R=None;Q=Pt
    while k:
        if k&1: R=add(R,Q)
        Q=add(Q,Q); k>>=1
    return R
sh=(a2*i3)%P
G=((on[0][0]+sh)%P, on[0][1]%P)
print('\nsample short-form point on Y^2=X^3+A X+B :', (pow(G[1],2,P)-pow(G[0],3,P)-A*G[0]-B)%P==0)
print('n_secp * G = %s   %s'%(mul(NSEC,G), '<== ORDER IS n_secp (prime)' if mul(NSEC,G) is None else ''))
print('n_secp prime? %s'%sympy.isprime(NSEC))
# also test all inliers have order dividing n_secp
bad=sum(1 for X,Y in on if mul(NSEC,((X+sh)%P,Y))is not None)
print('inlier points with n_secp*Pt != O : %d of %d'%(bad,len(on)))
# where does the pinned point (x12186,x16742) sit?
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
x1,y1=v[12186]%P,v[16742]%P
print('\npinned (x12186,x16742) on the GENERAL curve ? %s'%(
      (pow(y1,2,P)-pow(x1,3,P)-a2*pow(x1,2,P)-a4*x1-a6)%P==0))
G1=((x1+sh)%P,y1)
print('   on the short form ? %s ; n_secp*G1 = %s'%(
      (pow(G1[1],2,P)-pow(G1[0],3,P)-A*G1[0]-B)%P==0, mul(NSEC,G1)))
json.dump({'a2':str(a2),'a4':str(a4),'a6':str(a6),'A':str(A),'B':str(B),'shift':str(sh),
           'n_inliers':len(on)},open('/home/user/integer_solver/solve_lab/agentA_work/curve_final.json','w'))
