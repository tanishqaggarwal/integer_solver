"""Clean forward evaluator + check reporter, with cycle handling."""
import sys, os, collections, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L

NV = L.NVARS
P = L.P

# ---- build topo over defined vars, and SCCs for the cyclic remainder ----
definer = L.definer
avars = L.avars
deps = {}
for v,a in definer.items():
    deps[v] = [u for u in avars[a] if u!=v and u in definer]

def tarjan(nodes, dep):
    index={}; low={}; onstk={}; stk=[]; out=[]; counter=[0]
    for root in nodes:
        if root in index: continue
        work=[(root,0)]
        while work:
            v,pi = work[-1]
            if pi==0:
                index[v]=low[v]=counter[0]; counter[0]+=1
                stk.append(v); onstk[v]=True
            rec=False
            for i in range(pi, len(dep[v])):
                w=dep[v][i]
                if w not in index:
                    work[-1]=(v,i+1); work.append((w,0)); rec=True; break
                elif onstk.get(w):
                    low[v]=min(low[v], index[w])
            if rec: continue
            if low[v]==index[v]:
                comp=[]
                while True:
                    w=stk.pop(); onstk[w]=False; comp.append(w)
                    if w==v: break
                out.append(comp)
            work.pop()
            if work:
                u=work[-1][0]; low[u]=min(low[u], low[v])
    return out

SCCS = tarjan(list(definer), deps)   # reverse topological order
ORDER = SCCS  # each comp: list of vars

def evalpoly(Pp, v):
    s=0
    for m,c in Pp.items():
        t=c
        for u in m: t*=v[u]
        s+=t
    return s

def solve_lin(a, t, v):
    """solve atom a == 0 for var t (must be linear in t)."""
    c=0
    for m,cc in L.polys[a].items():
        k=m.count(t)
        if k==0: continue
        if k>1: return None
        term=cc
        for u in m:
            if u!=t: term*=v[u]
        c+=term
    if c==0: return None
    old=v[t]; v[t]=0
    rest=evalpoly(L.polys[a], v)
    v[t]=old
    if rest % c: return None
    return -rest//c

def forward(v, iters=40):
    for comp in ORDER:
        if len(comp)==1:
            u=comp[0]
            x=solve_lin(definer[u], u, v)
            if x is not None: v[u]=x
        else:
            for _ in range(iters):
                ch=False
                for u in comp:
                    x=solve_lin(definer[u], u, v)
                    if x is not None and x!=v[u]:
                        v[u]=x; ch=True
                if not ch: break
    return v

CHECKS = [a for a in range(L.NA) if a not in L.atom_out or L.atom_out.get(a) is None]
CHECKS = [a for a in range(L.NA) if L.atom_out.get(a) is None]

def bad_checks(v):
    return [a for a in CHECKS if evalpoly(L.polys[a], v)!=0]

if __name__=='__main__':
    import json
    v=[0]*NV
    sel = {}
    if len(sys.argv)>1 and sys.argv[1]!='-':
        for kv in sys.argv[1].split(','):
            k,val=kv.split('='); sel[int(k)]=int(val)
    for k,val in sel.items(): v[k]=val
    forward(v)
    for k,val in sel.items(): assert v[k]==val, (k,v[k],val)
    b=bad_checks(v)
    av=L.all_atom_values(v)
    f=L.failing_eqs(av)
    print(f"sel={sel}  bad_checks={len(b)}  failing_eqs={len(f)}  score={L.NEQ-len(f)}")
    print("first bad checks:", b[:40])
