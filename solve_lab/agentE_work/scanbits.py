import sys, time, pickle, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, bitfeas2 as B
import signal
class TO(Exception): pass
def _h(a,b): raise TO()
signal.signal(signal.SIGALRM,_h)
LIMIT=int(__import__('os').environ.get('LIMIT','75'))
L=pickle.load(open('bitlists.pkl','rb'))
which=sys.argv[1]
bits=L['A'] if which=='A' else L['B']
lo=int(sys.argv[2]) if len(sys.argv)>2 else 0
hi=int(sys.argv[3]) if len(sys.argv)>3 else len(bits)
res={}
for b in bits[lo:hi]:
    t0=time.time()
    base={18956:B.C, b:1}
    signal.alarm(LIMIT)
    try:
        bad0,out,S,cols,nonlin,rounds,v0=B.analyse(base,{18956,b})
        signal.alarm(0)
    except Exception as e:
        signal.alarm(0)
        print(f"x_{b}\tTIMEOUT/ERROR {type(e).__name__}",flush=True); res[b]=dict(feas=None); pickle.dump(res,open(f'scan_{which}.pkl','wb')); continue
    maxr,nv,nr,nd,msg,sol=out[-1]
    pins=sorted(set(bad0)-B.CORE)
    feas=sol is not None
    exact=None
    if feas:
        ns=dict(base)
        for f,d in sol.items():
            if d: ns[f]=v0[f]+d
        v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
        exact=(len(ff), sorted(set(av)-B.CORE))
    print(f"x_{b}\tpins={pins}\tvars={nv}\trows={nr}\tnl={nd}\tFEAS={feas}\tmsg={msg[:70]}\texact={exact}\t{time.time()-t0:.0f}s",flush=True)
    res[b]=dict(pins=pins,feas=feas,msg=msg,vars=nv,rows=nr,exact=exact,
                sol={int(k):int(x) for k,x in (sol or {}).items() if x})
    pickle.dump(res,open(f'scan_{which}.pkl','wb'))
