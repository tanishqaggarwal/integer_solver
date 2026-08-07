"""Repair engine in agent-H frame: greedy exact zeroing of nonzero check atoms via free inputs."""
import ev, pickle, json, os, sys, time
from collections import defaultdict
HERE=os.path.dirname(os.path.abspath(__file__))
S=pickle.load(open(os.path.join(HERE,'support.pkl'),'rb'))
csup=S['csup']; fidx=S['fidx']; inv={i:v for v,i in fidx.items()}
def bl(x):
    o=[]
    while x: b=x&-x; o.append(inv[b.bit_length()-1]); x^=b
    return o
SUPV={a:bl(s) for a,s in csup.items()}
eq_of=defaultdict(list)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: eq_of[a].append(i)
# how many check atoms each free input touches (cheapness heuristic)
touch=defaultdict(int)
for a,vs in SUPV.items():
    for u in vs: touch[u]+=1

def state(fv):
    v=ev.forward(fv); e,av=ev.eqvals(v)
    fails=[i for i,x in enumerate(e) if x!=0]
    nz={a:x for a,x in av.items() if x!=0}
    return v,fails,nz

def atomval(fv,a):
    v=ev.forward(fv)
    ns={'v':v,'__builtins__':{}}
    return eval(ev.CHECKCODE[a],ns)

def solve_u(fv,a,u,deg=3):
    """integer roots of atom a as a polynomial in free input u (degree<=deg)."""
    base=fv.get(u,0)
    ys=[]
    g=dict(fv)
    for k in range(deg+1):
        g[u]=base+k; ys.append(atomval(g,a))
    # finite differences -> Newton form
    d=[ys[:]]
    for k in range(deg):
        prev=d[-1]; d.append([prev[i+1]-prev[i] for i in range(len(prev)-1)])
    c=[d[k][0] for k in range(deg+1)]
    if all(x==0 for x in c): return 'ANY'
    # highest nonzero order
    order=max(k for k in range(deg+1) if c[k]!=0)
    if order==0: return []
    if order==1:
        # y(t) = c0 + c1*t  (t = u-base)
        if c[0] % c[1]: return []
        return [base - c[0]//c[1]]
    # order>=2: brute rational root not attempted; try small window + exact quadratic
    if order==2:
        # y(t)=c0 + c1*t + c2*t(t-1)/1  (Newton with step 1): c2 = second difference
        # expand: c2/2*? -- use Newton basis binom
        A=c[2]//2 if c[2]%2==0 else None
        # generic: y(t)=c0 + c1*t + c2*t*(t-1)/2
        # multiply by 2: 2c0 + 2c1 t + c2 t^2 - c2 t = 0
        aa=c[2]; bb=2*c[1]-c[2]; cc=2*c[0]
        disc=bb*bb-4*aa*cc
        if disc<0: return []
        r=isqrt_exact(disc)
        if r is None: return []
        out=[]
        for s in (-bb+r,-bb-r):
            if aa and s%(2*aa)==0: out.append(base+s//(2*aa))
        return out
    return []

def isqrt_exact(n):
    if n<0: return None
    import math
    x=math.isqrt(n)
    return x if x*x==n else None

def greedy(fv, rounds=12, verbose=True, forbid=()):
    fv=dict(fv)
    v,fails,nz=state(fv)
    best=len(fails); bestfv=dict(fv)
    for r in range(rounds):
        improved=False
        order=sorted(nz, key=lambda a: len(SUPV[a]))
        for a in order:
            if atomval(fv,a)==0: continue
            cands=sorted(SUPV[a], key=lambda u: touch[u])
            for u in cands:
                if u in forbid: continue
                roots=solve_u(fv,a,u)
                if roots=='ANY' or not roots: continue
                for rt in roots:
                    g=dict(fv); g[u]=rt
                    _,f2,nz2=state(g)
                    if len(f2)<len(fails):
                        fv=g; fails=f2; nz=nz2; improved=True
                        if verbose: print('  fix a%d via x_%d -> failing %d'%(a,u,len(fails)))
                        break
                if improved: break
            if improved: break
        if len(fails)<best:
            best=len(fails); bestfv=dict(fv)
        if not improved: break
    return bestfv, best

if __name__=='__main__':
    C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    fv={5096:C2,21589:C1,16742:C1,12186:C2,18956:C1,24468:C2,542:1,1438:1}
    t0=time.time()
    v,fails,nz=state(fv)
    print('start score',39033-len(fails),'nz',sorted(nz))
    bfv,best=greedy(fv,rounds=30)
    print('after greedy: score',39033-best,'time %.1f'%(time.time()-t0))
    v=ev.forward(bfv)
    json.dump({('x_%d'%i):v[i] for i in range(38748) if v[i]!=0},open(os.path.join(HERE,'g1.json'),'w'))
    json.dump({str(k):v2 for k,v2 in bfv.items()},open(os.path.join(HERE,'g1_fv.json'),'w'))
