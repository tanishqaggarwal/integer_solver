"""Wide oracle: classify leaves by their FULL bad-atom delta signature, not E's 5-row projection.
   Then compare the induced partition against tree96 subtree supports."""
import sys,json,collections,time,pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentM_work')
import mcore as M, xcompare as X
import engine as E, fast

def wide(seed, mode='support'):
    """mode 'support': signature = frozenset of atoms whose residual changed.
       mode 'full'   : signature = tuple of (atom, delta) pairs."""
    v0=E.forward(seed); bad0=E.badatoms(v0)
    sig={}
    for f in M.bools():
        if v0[f]!=0: sig[f]='ON'; continue
        b1,_=fast.resid_delta(v0,bad0,{f:1})
        keys=set(b1)|set(bad0)
        d={a:(b1.get(a,0)-bad0.get(a,0)) for a in keys}
        d={a:x for a,x in d.items() if x}
        sig[f]=frozenset(d) if mode=='support' else tuple(sorted(d.items()))
    return v0,bad0,sig

def part(sig):
    cls=collections.defaultdict(list)
    for f,s in sig.items(): cls[s].append(f)
    return cls

def report(tag,seed,mode='full'):
    t0=time.time()
    v0,bad0,sig=wide(seed,mode)
    cls=part(sig)
    live=[(s,v) for s,v in cls.items() if s!='ON' and s]
    inert=[v for s,v in cls.items() if s!='ON' and not s]
    on=[v for s,v in cls.items() if s=='ON']
    print('=== %s  bad0=%d  classes=%d  (%.0fs)'%(tag,len(bad0),len(live),time.time()-t0))
    parts=[('ch',set(v)) for _,v in sorted(live,key=lambda kv:-len(kv[1]))]
    if inert: parts.append(('INERT',set(inert[0])))
    if on: parts.append(('ON',set(on[0])))
    for lab,s in parts:
        nd,left=X.decompose(s)
        ex=[k for k,g in X.SUB if g==s]
        print('  %-6s n=%3d exact=%s cover=%s loose=%d'%(lab,len(s),ex or '-',
            ','.join('%s(%d)'%(a,b) for a,b in nd) or '-',len(left)))
    onset=parts[-1][1] if on else set()
    ps=[s for lab,s in parts if lab!='ON']
    cross=[]
    for k,g in X.SUB:
        g2=g-onset
        if not g2: continue
        hit=[len(g2&s) for s in ps if g2&s]
        if len(hit)>1: cross.append((k,len(g2),sorted(hit,reverse=True)))
    print('  crossings(ON removed): %d  %s'%(len(cross),cross[:8]))
    return parts,cross

if __name__=='__main__':
    s0=M.load_seed(); BASE=dict(s0); BASE[1530]=0; BASE[1603]=0
    report('alloff (wide)',dict(BASE))
    c=dict(BASE); c[2081]=1; report('one b-side leaf 2081 (wide)',c)
    c=dict(BASE); c[24601]=1; report('one a-side leaf 24601 (wide)',c)
    sd={int(k):int(v) for k,v in json.load(open('/home/user/integer_solver/solve_lab/agentM_work/deliv_seed.json')).items()}
    report('deliverable cfg (wide)',sd)
