"""Fast beam search over pin-cascade repairs (incremental evaluator)."""
import ev, fast, json, os, sys, time, math, pickle
from fast import St, chk, csup, inv
HERE=os.path.dirname(os.path.abspath(__file__))
S=pickle.load(open(os.path.join(HERE,'support.pkl'),'rb'))
def bl(x):
    o=[]
    while x: b=x&-x; o.append(inv[b.bit_length()-1]); x^=b
    return o
SUPV={a:bl(s) for a,s in csup.items()}

def isqrt_exact(n):
    if n<0: return None
    x=math.isqrt(n); return x if x*x==n else None

def solve_u(st,a,u,deg=3):
    base=st.fv.get(u,0)
    ys=[]
    tmp=st.clone()
    for k in range(deg+1):
        tmp.set_free({u:base+k}); ys.append(tmp.av[a])
    d=[ys[:]]
    for k in range(deg):
        p=d[-1]; d.append([p[i+1]-p[i] for i in range(len(p)-1)])
    c=[d[k][0] for k in range(deg+1)]
    if all(x==0 for x in c): return []
    order=max(k for k in range(deg+1) if c[k]!=0)
    if order==0: return []
    if order==1:
        if c[0]%c[1]: return []
        return [base-c[0]//c[1]]
    if order==2:
        aa=c[2]; bb=2*c[1]-c[2]; cc=2*c[0]
        disc=bb*bb-4*aa*cc
        r=isqrt_exact(disc)
        if r is None: return []
        out=[]
        for s in (-bb+r,-bb-r):
            if s%(2*aa)==0: out.append(base+s//(2*aa))
        return out
    return []

def gen_moves(st):
    out=[]
    for a in st.nz():
        for u in SUPV[a]:
            for rt in solve_u(st,a,u):
                if st.fv.get(u,0)!=rt: out.append((u,rt))
    return out

def beam(fv0,width=30,depth=80,verbose=True,log=None):
    st0=St(fv0)
    frontier=[st0]
    best=st0.clone(); bestk=(len(st0.fails),len(st0.nz()))
    seen=set()
    for d in range(depth):
        cand=[]
        for st in frontier:
            for (u,rt) in gen_moves(st):
                g=st.clone().set_free({u:rt})
                sig=(tuple(sorted(g.fails)),)
                k=(len(g.fails),len(g.nz()))
                cand.append((k,g))
        if not cand: break
        cand.sort(key=lambda t:t[0])
        # dedupe by fail-set signature
        fr=[]; sg=set()
        for k,g in cand:
            s=tuple(sorted(g.fails))
            if s in sg: continue
            sg.add(s); fr.append(g)
            if len(fr)>=width: break
        frontier=fr
        k=(len(frontier[0].fails),len(frontier[0].nz()))
        if k<bestk: bestk=k; best=frontier[0].clone()
        m='d%02d best=(%d fail,%d nz) score=%d cands=%d'%(d,k[0],k[1],39033-k[0],len(cand))
        if verbose: print(m,flush=True)
        if log: log.write(m+'\n'); log.flush()
        if k[0]==0: break
    return best,bestk

if __name__=='__main__':
    C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    ub=int(sys.argv[1]); wb=int(sys.argv[2])
    W=int(os.environ.get('W','30')); D=int(os.environ.get('D','60'))
    fv={5096:C2,21589:C1,16742:C1,12186:C2,18956:C1,24468:C2,ub:1,wb:1}
    t0=time.time()
    log=open(os.path.join(HERE,'runs_beam_%d_%d.log'%(ub,wb)),'w')
    best,k=beam(fv,W,D,log=log)
    print('BEST score',39033-k[0],'nz',k[1],'%.1fs'%(time.time()-t0))
    out=os.path.join(HERE,'B_%d_%d_%d.json'%(39033-k[0],ub,wb))
    json.dump({('x_%d'%i):best.v[i] for i in range(38748) if best.v[i]!=0},open(out,'w'))
    json.dump({str(a):b for a,b in best.fv.items()},open(out.replace('.json','_fv.json'),'w'))
    print('wrote',out)
