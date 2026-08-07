"""Extract the two load-pin constants per selector bit as a POINT (X,Y) and fit the
general Weierstrass curve  y^2 = x^3 + a2 x^2 + a4 x + a6  by EXACT linear algebra over F_p.
Each point contributes a LINEAR equation  y^2 - x^3 = a2 x^2 + a4 x + a6  in (a2,a4,a6)."""
import sys, collections, json, itertools; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
# parse pins:  bit*(H + <deg-2 partner>) - s*x_T   (exactly one huge coeff, on a deg-1 monomial)
bybit=collections.defaultdict(list)
for a in range(L.NA):
    Pp=L.polys[a]
    hs=[(m,c) for m,c in Pp.items() if abs(c)>=10**40]
    if len(hs)!=1: continue
    m,H=hs[0]
    if len(m)!=1: continue
    bit=m[0]
    rest=[(mm,cc) for mm,cc in Pp.items() if mm!=m]
    deg2=[(mm,cc) for mm,cc in rest if len(mm)==2 and bit in mm]
    deg1=[(mm,cc) for mm,cc in rest if len(mm)==1]
    if len(deg2)!=1 or len(deg1)!=1: continue
    (m2,c2),(m1,s)=deg2[0],deg1[0]
    partner=[u for u in m2 if u!=bit]
    partner=partner[0] if partner else bit
    tgt=m1[0]
    bybit[bit].append({'atom':a,'H':H,'c2':c2,'partner':partner,'s':s,'tgt':tgt})
two={b:x for b,x in bybit.items() if len(x)==2}
print('gating bits with exactly two clean load pins: %d  (all bits: %d)'%(len(two),len(bybit)))
def loaded(pin):
    # a = H*bit + c2*bit*partner + s*tgt = 0 with bit=1 -> tgt = -(H + c2*partner)/s
    val=(-(pin['H'] + pin['c2']*v[pin['partner']]) * pow(pin['s'],-1,P))%P
    return val
pts=[]
for b,x in sorted(two.items()):
    x=sorted(x,key=lambda d:d['atom'])
    pts.append((b,loaded(x[0]),loaded(x[1]),x[0]['tgt'],x[1]['tgt']))
print('extracted %d candidate leaf points'%len(pts))
for b,X,Y,t1,t2 in pts[:3]:
    print('   bit x%-6d X=%d  Y=%d  (targets x%d,x%d)'%(b,X,Y,t1,t2))

def fit(points):
    """solve a2*x^2 + a4*x + a6 = y^2 - x^3 over F_p, exactly, overdetermined."""
    rows=[[pow(X,2,P)%P, X%P, 1, (pow(Y,2,P)-pow(X,3,P))%P] for X,Y in points]
    M=[r[:] for r in rows]; piv=[]; r=0
    for c in range(3):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]%P: pr=i; break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],-1,P); M[r]=[x*inv%P for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%P:
                f=M[i][c]; M[i]=[(a-f*bb)%P for a,bb in zip(M[i],M[r])]
        piv.append(c); r+=1
    inc=sum(1 for i in range(r,len(M)) if M[i][3]%P)
    sol=[0,0,0]
    for i,c in enumerate(piv): sol[c]=M[i][3]%P
    return sol,r,inc,len(M)

for order,name in [((1,2),'(pin0=X, pin1=Y)'),((2,1),'(pin0=Y, pin1=X)')]:
    pl=[(p[order[0]],p[order[1]]) for p in pts]
    sol,r,inc,tot=fit(pl)
    print('\nfit %s over %d points: rank=%d inconsistent rows=%d'%(name,tot,r,inc))
    if inc==0:
        a2,a4,a6=sol
        print('   *** ALL %d POINTS ON ONE CURVE ***'%tot)
        print('   a2=%d\n   a4=%d\n   a6=%d'%(a2,a4,a6))
        json.dump({'a2':str(a2),'a4':str(a4),'a6':str(a6),
                   'points':[[str(x),str(y)] for x,y in pl]},
                  open('/home/user/integer_solver/solve_lab/agentA_work/weier.json','w'))
    else:
        # how many points fit the curve determined by the first 3?
        a2,a4,a6=sol
        good=sum(1 for X,Y in pl if (pow(Y,2,P)-pow(X,3,P)-a2*pow(X,2,P)-a4*X-a6)%P==0)
        print('   points satisfying the fitted curve: %d of %d'%(good,tot))
