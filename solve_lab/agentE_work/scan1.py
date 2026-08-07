import sys, pickle, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import harness as H
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
base={18956:C}
v=H.forward(base); f,av=H.eqfails(v)
print("base fails",len(f),"nz",sorted(av.items())[:5])
BAD=sorted(av)
free=pickle.load(open('conefree.pkl','rb'))
res={}
t0=time.time()
for k,fv in enumerate(free):
    s=dict(base); s[fv]=base.get(fv,0)+1
    vv=H.forward(s); ff,aav=H.eqfails(vv)
    res[fv]=(len(ff), sorted(aav), {a:aav[a] for a in BAD if a in aav})
    if k%50==0: print(k,len(free),f"{time.time()-t0:.0f}s",file=sys.stderr)
pickle.dump(res, open('scan1.pkl','wb'))
# report movers
for fv,(nf,nz,d) in sorted(res.items(), key=lambda kv: kv[1][0]):
    if nz!=BAD or nf!=len(f):
        print(f"x_{fv}: fails={nf} nzatoms={nz}")
