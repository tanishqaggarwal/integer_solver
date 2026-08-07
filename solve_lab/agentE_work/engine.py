"""Repair engine: cone-restricted forward evaluation + exact single-free-input solve."""
import sys, pickle, math, json, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import harness as H

atoms=H.atoms; acodes=H.acodes; avars=H.avars; definer=H.definer; NV=H.NV
SEQ=H.SEQ
pos={u:k for k,u in enumerate(SEQ)}
FREE=set(H.FREE)|{u for u in range(NV) if definer[u] is None}

_conecache={}
def cone(aid):
    c=_conecache.get(aid)
    if c is not None: return c
    seen=set(); st=list(avars[aid]); fr=set()
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        dv=definer[u]
        if dv is None: fr.add(u); continue
        for w in avars[dv[0]]:
            if w!=u: st.append(w)
    order=sorted([u for u in seen if definer[u] is not None], key=lambda u: pos[u])
    c=(order, sorted(fr), seen)
    _conecache[aid]=c
    return c

def _solvevar(v, ns, u, i, kd):
    c=acodes[i]
    v[u]=0; c0=eval(c,ns)
    v[u]=1; c1=eval(c,ns)
    if kd=='lin':
        sl=c1-c0
        v[u]=-c0//sl if sl and c0%sl==0 else 0
    else:
        v[u]=2; c2=eval(c,ns)
        A2=c2-2*c1+c0; A=A2//2; B=c1-c0-A; C=c0
        disc=B*B-4*A*C
        if disc<0 or A==0: v[u]=0; return
        r=math.isqrt(disc)
        if r*r!=disc: v[u]=0; return
        rts={(-B+s)//(2*A) for s in (r,-r) if (-B+s)%(2*A)==0}
        v[u]=rts.pop() if len(rts)==1 else 0

def forward(seed):
    v=[0]*NV
    for k,val in seed.items(): v[k]=val
    ns={'v':v,'__builtins__':{}}
    for u in SEQ:
        i,kind=definer[u]
        _solvevar(v,ns,u,i,kind[0])
    return v

def resid_cone(aid, seed, base_v=None):
    """Evaluate only the cone of atom aid; returns residual."""
    order,fr,seen=cone(aid)
    v=list(base_v) if base_v is not None else [0]*NV
    for k,val in seed.items(): v[k]=val
    ns={'v':v,'__builtins__':{}}
    for u in order:
        i,kind=definer[u]
        _solvevar(v,ns,u,i,kind[0])
    return eval(acodes[aid],ns), v

def badatoms(v):
    ns={'v':v,'__builtins__':{}}
    out={}
    for i in range(len(atoms)):
        r=eval(acodes[i],ns)
        if r: out[i]=r
    return out

def eqfails(av):
    fails=[]
    for e,(issq,outer,terms) in enumerate(H.eqt):
        s=0
        for c,a in terms:
            if a<0: s+=c
            elif a in av: s+=c*av[a]
        if s: fails.append(e)
    return fails

def score(seed):
    v=forward(seed); av=badatoms(v); f=eqfails(av)
    return len(f), av, v

def solve_for(aid, f, seed):
    """Return candidate value of free var f that zeroes atom aid, or None."""
    s=dict(seed); c=s.get(f,0)
    vals=[]
    for d in (0,1,2):
        s[f]=c+d
        r,_=resid_cone(aid,s)
        vals.append(r)
    c0,c1,c2=vals
    A2=c2-2*c1+c0
    if A2==0:
        sl=c1-c0
        if sl==0: return None
        if c0%sl: return None
        return c - c0//sl
    A=A2//2; B=c1-c0-A; C=c0
    disc=B*B-4*A*C
    if disc<0: return None
    r=math.isqrt(disc)
    if r*r!=disc: return None
    rts=sorted({(-B+s2)//(2*A) for s2 in (r,-r) if (-B+s2)%(2*A)==0})
    if not rts: return None
    return c + rts[0]
