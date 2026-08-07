"""bn_defic: exact minimum of |E(S)| - |S| over nonempty sets S of boolean atoms,
via maximum-closure / min-cut.  Also greedy growth for small blocks."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad

cen = json.load(open(os.path.join(HERE,'bn_census.json')))
bools = {int(a):tuple(t) for a,t in cen['bools'].items()}
BA = sorted(bools)
FREESET = set(ad.FREE)

aeq = {a: set(L.atom2eq.get(a,{})) for a in BA}
allE = set()
for a in BA: allE |= aeq[a]
print(f'boolean atoms {len(BA)}, union of their equations {len(allE)}, '
      f'incidences {sum(len(aeq[a]) for a in BA)}')
print(f'global deficiency |E(ALL)|-|ALL| = {len(allE)-len(BA)}')

# how many boolean atoms per equation
per = collections.Counter()
for a in BA:
    for e in aeq[a]: per[e]+=1
h = collections.Counter(per.values())
print('boolean-atoms-per-equation histogram:', sorted(h.items())[:25])

# ------------------------------------------------------------------
# Dinic max-flow for maximum closure: max (|S| - |E(S)|)
# ------------------------------------------------------------------
class Dinic:
    def __init__(s,n):
        s.n=n; s.g=[[] for _ in range(n)]
    def add(s,u,v,c):
        s.g[u].append([v,c,len(s.g[v])])
        s.g[v].append([u,0,len(s.g[u])-1])
    def bfs(s,S,T):
        s.lvl=[-1]*s.n; s.lvl[S]=0; q=[S]
        for u in q:
            for e in s.g[u]:
                if e[1]>0 and s.lvl[e[0]]<0:
                    s.lvl[e[0]]=s.lvl[u]+1; q.append(e[0])
        return s.lvl[T]>=0
    def dfs(s,u,T,f):
        if u==T: return f
        while s.it[u]<len(s.g[u]):
            e=s.g[u][s.it[u]]
            if e[1]>0 and s.lvl[e[0]]==s.lvl[u]+1:
                d=s.dfs(e[0],T,min(f,e[1]))
                if d>0:
                    e[1]-=d; s.g[e[0]][e[2]][1]+=d; return d
            s.it[u]+=1
        return 0
    def flow(s,S,T):
        fl=0
        while s.bfs(S,T):
            s.it=[0]*s.n
            while True:
                f=s.dfs(S,T,10**18)
                if f==0: break
                fl+=f
        return fl

def maxclosure(atoms, label=''):
    """max over nonempty S subset of `atoms` of (|S| - |E(S)|)."""
    E = sorted(set().union(*[aeq[a] for a in atoms])) if atoms else []
    eidx = {e:i for i,e in enumerate(E)}
    n = len(atoms)+len(E)+2
    S=n-2; T=n-1
    d=Dinic(n)
    INF=10**9
    ai={a:i for i,a in enumerate(atoms)}
    for a in atoms: d.add(S, ai[a], 1)
    for e in E: d.add(len(atoms)+eidx[e], T, 1)
    for a in atoms:
        for e in aeq[a]:
            d.add(ai[a], len(atoms)+eidx[e], INF)
    sys.setrecursionlimit(100000)
    f = d.flow(S,T)
    best = len(atoms)-f
    # recover S: nodes reachable from source in residual
    seen=[False]*n; st=[S]; seen[S]=True
    while st:
        u=st.pop()
        for e in d.g[u]:
            if e[1]>0 and not seen[e[0]]:
                seen[e[0]]=True; st.append(e[0])
    Ssel=[a for a in atoms if seen[ai[a]]]
    Esel=set().union(*[aeq[a] for a in Ssel]) if Ssel else set()
    print(f'[{label}] maxflow={f}  max(|S|-|E(S)|) = {best}  '
          f'recovered |S|={len(Ssel)} |E(S)|={len(Esel)} defic={len(Esel)-len(Ssel)}')
    return best, Ssel, Esel

print()
print('=== ALL boolean atoms ===')
b1,S1,E1 = maxclosure(BA,'all')
print()
print('=== boolean atoms on FREE vars only ===')
BF = [a for a in BA if bools[a][0] in FREESET]
b2,S2,E2 = maxclosure(BF,'free')

json.dump({'all':{'best':b1,'S':S1,'E':sorted(E1)},
           'free':{'best':b2,'S':S2,'E':sorted(E2)}},
          open(os.path.join(HERE,'bn_defic.json'),'w'))
