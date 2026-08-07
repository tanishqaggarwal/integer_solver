"""Exact cascade repair on a full assignment: zero a target atom by solving it for one of its
   variables, then recursively repair the atoms that break.  Scored at EQUATION level."""
import sys, json, collections, math, time, random, heapq
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E
P=C.P
NA=len(H.atoms)

def loadv(path):
    d=json.load(open(path)); v=[0]*E.NV
    for k,val in d.items(): v[int(k.split('_')[1])]=int(val)
    return v

def aval(v,i):
    return eval(H.acodes[i],{'v':v,'__builtins__':{}})

def solve_for(v,i,u):
    """All integer values of v[u] making atom i zero (atom is <=deg2 in u). Returns list."""
    o=v[u]
    v[u]=0; c0=aval(v,i); v[u]=1; c1=aval(v,i); v[u]=2; c2=aval(v,i); v[u]=o
    A2=c2-2*c1+c0
    if A2==0:
        sl=c1-c0
        if sl==0: return [] if c0 else 'ANY'
        return [-c0//sl] if c0%sl==0 else []
    A=A2//2; B=c1-c0-A; Cc=c0
    disc=B*B-4*A*Cc
    if disc<0: return []
    r=math.isqrt(disc)
    if r*r!=disc: return []
    return sorted({(-B+s)//(2*A) for s in (r,-r) if (-B+s)%(2*A)==0})

def full_bad(v):
    ns={'v':v,'__builtins__':{}}
    return {i:eval(H.acodes[i],ns) for i in range(NA) if eval(H.acodes[i],ns)}

def eqscore(bad):
    return 39033-len(E.eqfails(bad))

class Rep:
    def __init__(self, v0, bad0):
        self.v0=v0; self.bad0=dict(bad0)
    def evaluate(self, changes):
        v=list(self.v0)
        for u,x in changes.items(): v[u]=x
        dirty=set(self.bad0)
        for u in changes: dirty|=set(H.occ[u])
        ns={'v':v,'__builtins__':{}}
        bad={}
        for i in dirty:
            r=eval(H.acodes[i],ns)
            if r: bad[i]=r
        for i,r in self.bad0.items():
            if i not in dirty: bad[i]=r
        return v,bad

def dfs(rep, changes, open_bad, depth, budget, best, seen, order_rand):
    key=tuple(sorted(changes.items()))
    if key in seen: return
    seen.add(key)
    v,bad=rep.evaluate(changes)
    ff=len(E.eqfails(bad))
    if ff<best[0]:
        best[0]=ff; best[1]=dict(changes); best[2]=sorted(bad)
        print(f"    depth{depth} fails={ff} SCORE={39033-ff} atoms={sorted(bad)} changes={len(changes)}",flush=True)
    if depth>=budget or not bad: return
    # choose a bad atom with the fewest variables/occurrences
    cands=sorted(bad, key=lambda a:(sum(len(H.occ[u]) for u in H.avars[a]), a))
    for a in cands[:3]:
        vs=sorted(H.avars[a], key=lambda u:len(H.occ[u]))
        for u in vs[:4]:
            sols=solve_for(v,a,u)
            if sols=='ANY': sols=[0]
            for s in sols[:2]:
                if s==v[u]: continue
                if abs(s).bit_length()>60000: continue
                nc=dict(changes); nc[u]=s
                dfs(rep,nc,None,depth+1,budget,best,seen,order_rand)

if __name__=='__main__':
    src=sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
    bud=int(sys.argv[2]) if len(sys.argv)>2 else 4
    v0=loadv(src); bad0=full_bad(v0)
    print("start atoms",sorted(bad0),"fails",len(E.eqfails(bad0)))
    rep=Rep(v0,bad0)
    best=[len(E.eqfails(bad0)),{},sorted(bad0)]
    t0=time.time()
    dfs(rep,{},None,0,bud,best,set(),random.Random(1))
    print("BEST fails=%d score=%d changes=%d atoms=%s (%.0fs)"%(best[0],39033-best[0],len(best[1]),best[2],time.time()-t0))
    if best[0]<len(E.eqfails(bad0)):
        v,bad=rep.evaluate(best[1])
        json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('S_cascade_%d.json'%(39033-best[0]),'w'))
        print("wrote S_cascade_%d.json"%(39033-best[0]))
