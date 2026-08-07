"""The complete set of atoms whose corruption could possibly affect the target equations.
This is more fundamental than a site list: it bounds the entire space of useful corruptions."""
import sys, json, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentT_work/mirror/F')
from fwd import Engine,NV
from circ2 import vars_of
E=Engine()
EQ={a:set(v) for a,v in pickle.load(open('eqmap_exact.pkl','rb')).items()}
M25=set([2554,5324,6816,8124,8680,9041,9123,9421,11226,12231,12270,12350,14584,15558,18673,
         21000,22044,22534,22997,28929,29125,29330,32026,35512,38051])
MINE13=set(json.load(open('filterset.json')))
FS=M25 | MINE13          # UNION: never discard something that could help either baseline
print('M baseline 25; my baseline %d (strict subset of M25: %s); union used for filtering: %d'%(
    len(MINE13), MINE13<=M25, len(FS)))
M=pickle.load(open('full_model.pkl','rb')); NODE=M['NODE']; OUT=M['OUT']; ROOT=M['ROOT']
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
sys.setrecursionlimit(100000)
Hd=pickle.load(open('handles.pkl','rb')); handle=set(Hd['handle'])
freeall={}
def fa(v):
    if v in freeall: return freeall[v]
    if v not in defrhs: freeall[v]={v}; return freeall[v]
    freeall[v]=set(); s=set()
    for u in vars_of(defrhs[v]): s|=fa(u)
    freeall[v]=s; return s
atomu={}
for a in E.res:
    s=set()
    for v in vars_of(E.atoms[a]): s|=fa(v)
    hs=[v for v in s if v in handle]
    if len(hs)==1: atomu[a]=hs[0]
inc={a:len(EQ.get(a,set())&FS) for a in atomu}
inc13={a:len(EQ.get(a,set())&MINE13) for a in atomu}
hot=sorted([(n,a) for a,n in inc.items() if n>0],reverse=True)
print('FILTER SET: %d equations'%len(FS))
print('atoms with a unique cofactor: %d'%len(atomu))
print('ATOMS INCIDENT TO THE FILTER SET: %d   (all others cannot change any of them)'%len(hot))
print('incidence histogram:',sorted(collections.Counter(n for n,_ in hot).items(),reverse=True))
# locate each hot atom in the tree
wire2role={}
for n in NODE:
    for i,d in enumerate(OUT[n]):
        for k in ('va','vb','vab','out'): wire2role.setdefault(d[k],[]).append('x%d.%s[%d]'%(n,k,i))
print('\nALL incident atoms (h = the P-multiple to corrupt, u = its free cofactor):')
rows=[]
for n,a in hot:
    u=atomu[a]
    hm=[v for v in set(vars_of(E.atoms[a])) if v in defrhs and fa(v)=={u}]
    h=min(hm) if hm else None
    roles=[]
    for v in set(vars_of(E.atoms[a])):
        roles+=wire2role.get(v,[])
    rows.append(dict(atom=a,rows_target_union=n,rows_target_mine13=inc13[a],h=h,u=u,roles=sorted(set(roles))[:3]))
    print('  rt(union) %-3d rt(mine13) %-3d h=x%-6s u=x%-6d %-52s %s'%(n,inc13[a],h,u,a[:52],','.join(sorted(set(roles))[:2])))
json.dump(rows,open('incident_atoms.json','w'),indent=0)
