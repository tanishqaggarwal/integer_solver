import sys, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, bitfeas2 as B
R=pickle.load(open('scan_B.pkl','rb'))
bad=[b for b,d in R.items() if not d['feas']]
out={}
for b in sorted(bad):
    t0=time.time()
    base={18956:B.C, b:1}
    bad0,res,S,cols,nonlin,rounds,v0=B.analyse(base,{18956,b})
    maxr,nv,nr,nd,msg,sol=res[-1]
    ex=None
    if sol is not None:
        ns=dict(base)
        for f,d in sol.items():
            if d: ns[f]=v0[f]+d
        v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
        ex=(len(ff),sorted(set(av)-B.CORE))
    print(f"x_{b}\tpins={sorted(set(bad0)-B.CORE)}\tFEAS={sol is not None}\tmsg={msg[:60]}\texact={ex}\t{time.time()-t0:.0f}s",flush=True)
    out[b]=(sol is not None,msg,ex)
    pickle.dump(out,open('recheckB.pkl','wb'))
