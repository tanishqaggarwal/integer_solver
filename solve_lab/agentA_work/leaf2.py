"""RANSAC over every plausible extraction of the leaf-point coordinates.
A general Weierstrass fit y^2 = x^3+a2 x^2+a4 x+a6 absorbs ANY x-shift (so the
depression question cannot cause a false negative here); what it does NOT absorb is a
per-pin scaling, so we sweep the scaling conventions explicitly."""
import sys, collections, json, random, itertools; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
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
    pr=[u for u in m2 if u!=bit]; pr=pr[0] if pr else bit
    bybit[bit].append({'atom':a,'H':H,'c2':c2,'partner':pr,'s':s,'tgt':m1[0]})
two={b:sorted(x,key=lambda d:d['atom']) for b,x in bybit.items() if len(x)==2}
print('bits with 2 clean pins: %d'%len(two))
pv=collections.Counter()
for b,x in two.items():
    for pin in x: pv[('partner=0' if v[pin['partner']]==0 else ('partner=p' if v[pin['partner']]==P else 'partner=other'))]+=1
print('partner values:',dict(pv))
sv=collections.Counter('s=+-1' if abs(x[0]['s'])==1 else 's=multi' for b,x in two.items())
print('s(pin0):',dict(sv))
VARIANTS={
 'neg(H+c2*prt)/s': lambda q:(-(q['H']+q['c2']*v[q['partner']])*pow(q['s'],-1,P))%P,
 'pos(H+c2*prt)/s': lambda q:(( q['H']+q['c2']*v[q['partner']])*pow(q['s'],-1,P))%P,
 'negH/s'         : lambda q:((-q['H'])*pow(q['s'],-1,P))%P,
 'posH/s'         : lambda q:(( q['H'])*pow(q['s'],-1,P))%P,
 'negH'           : lambda q:(-q['H'])%P,
 'posH'           : lambda q:( q['H'])%P,
 'neg(H+c2*prt)'  : lambda q:(-(q['H']+q['c2']*v[q['partner']]))%P,
 'negH*s'         : lambda q:((-q['H'])*q['s'])%P,
}
def fit3(pts3):
    M=[[pow(X,2,P), X%P, 1, (pow(Y,2,P)-pow(X,3,P))%P] for X,Y in pts3]
    piv=[]; r=0
    for c in range(3):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]%P: pr=i;break
        if pr is None: return None
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],-1,P); M[r]=[x*inv%P for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%P:
                f=M[i][c]; M[i]=[(a-f*bb)%P for a,bb in zip(M[i],M[r])]
        piv.append(c); r+=1
    sol=[0,0,0]
    for i,c in enumerate(piv): sol[c]=M[i][3]%P
    return sol
random.seed(1)
best=None
for name,f in VARIANTS.items():
    for order in [(0,1),(1,0)]:
        pl=[]
        for b,x in sorted(two.items()):
            pl.append((f(x[order[0]]), f(x[order[1]])))
        bi=0; ba=None
        for _ in range(200):
            s3=random.sample(pl,3)
            sol=fit3(s3)
            if sol is None: continue
            a2,a4,a6=sol
            cnt=sum(1 for X,Y in pl if (pow(Y,2,P)-pow(X,3,P)-a2*pow(X,2,P)-a4*X-a6)%P==0)
            if cnt>bi: bi,ba=cnt,sol
        print('%-18s order=%s  best inliers %d / %d'%(name,order,bi,len(pl)),flush=True)
        if best is None or bi>best[0]: best=(bi,name,order,ba,pl)
print('\nBEST: %d inliers with variant %s order %s'%(best[0],best[1],best[2]))
if best[0]>=10:
    json.dump({'variant':best[1],'order':list(best[2]),'a':[str(x) for x in best[3]],
               'points':[[str(x),str(y)] for x,y in best[4]]},
              open('/home/user/integer_solver/solve_lab/agentA_work/weier.json','w'))
    print('saved weier.json')
