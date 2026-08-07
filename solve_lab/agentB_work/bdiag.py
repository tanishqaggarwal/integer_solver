import os,sys,json,collections
os.environ['ORIENT']=os.environ.get('ORIENT','orient4.pkl'); sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentB_work')
import beval as E, bfix as F
P=F.P; Q=F.Q
base='/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
v0=E.load(base)
fv0=[E.fval(f,v0) for f in range(len(E.facs))]
def s(p,lim=130):
    return ''.join('%s*%s'%(('%+d'%c) if abs(c)<10**8 else ('+K' if c>0 else '-K'),'*'.join('x%d'%v for v in m) if m else '1') for m,c in sorted(p.items()))[:lim]
def go(knobs,iters=5):
    freeval={v:v0[v] for v in E.free}
    val,nd,_=E.forward(freeval,default=v0)
    for it in range(iters):
        if 4432 in knobs: freeval[4432]=freeval[4432]-(val[28730]%P)
        if 3349 in knobs: freeval[3349]=freeval[3349]-(val[8731]%P)
        if 9118 in knobs: freeval[9118]=freeval[9118]-(val[9118]%P)
        if 7068 in knobs: freeval[7068]=freeval[7068]-((val[7068]-val[2099])%Q)
        val,nd,_=E.forward(freeval,default=v0)
    val=F.set_handles(val)
    ok,fail,fv=E.score(val)
    av=[]
    for a in E.atoms:
        t=1
        for f in a:
            t*=fv[f]
            if t==0: break
        av.append(t)
    nz=[i for i,x in enumerate(av) if x]
    print('knobs',knobs,'score',ok,'nfail',len(fail),'nonzero atoms',len(nz),'nd',nd)
    for a in nz:
        for f in E.atoms[a]:
            if fv[f]!=0:
                st = 'WAS-BAD' if fv0[f]!=0 else 'NEW-BAD'
                print('   atom%-6d f%-6d %-70s %s def=%s'%(a,f,s(E.facs[f],70),st,E.matchF[f]))
for k in [[3349,9118],[4432,3349,9118,7068],[4432],[7068]]:
    go(k); print()
