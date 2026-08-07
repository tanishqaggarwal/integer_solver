"""Emit candidate handle-corruption SETS, one per cut site, in the deliverable's exact shape:
   2 slot-link handles (parent's slot wires for child c) + 2 guard handles (c's own vab wires).
Calibration: the deliverable is the site c=x27994 / parent x4971.va -> {642,28730,29854,31864}."""
import sys, pickle, collections, re, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentF_work')
from fwd import Engine,NV
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
M=pickle.load(open('full_model.pkl','rb'))
NODE=M['NODE']; OUT=M['OUT']; tree=M['tree']; ROOT=M['ROOT']; sub=M['sub']; liveset=set(M['live'])
PR=pickle.load(open('price.pkl','rb')); eqs_of={a:set(v) for a,v in PR['eqs_of'].items()}
Hd=pickle.load(open('handles.pkl','rb')); handle=set(Hd['handle'])
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
import sys as _s; _s.setrecursionlimit(100000)
freeall={}
def fa(v):
    if v in freeall: return freeall[v]
    if v not in defrhs: freeall[v]={v}; return freeall[v]
    freeall[v]=set(); s=set()
    for u in vars_of(defrhs[v]): s|=fa(u)
    freeall[v]=s; return s
atomh={}
for a in E.res:
    s=set()
    for u in vars_of(E.atoms[a]): s|=fa(u)
    hs=[u for u in s if u in handle]
    if len(hs)==1: atomh[a]=hs[0]
byvar=collections.defaultdict(list)
for a in E.res:
    for u in set(vars_of(E.atoms[a])): byvar[u].append(a)
parent={}; side_of={}
for n in NODE:
    for s,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])): parent[ch]=n; side_of[ch]=s
depth={ROOT:0}
def sd(n,d):
    depth[n]=d
    if tree[n] is not None:
        for c in tree[n]: sd(c,d+1)
sd(ROOT,0)
LINK=PR['atom_of']              # slot wire -> its link atom
rows=[]
for c in NODE:
    if c==ROOT: continue
    n=parent[c]; s=side_of[c]
    slots=[d[s] for d in OUT[n]]
    vabs=[d['vab'] for d in OUT[c]]
    la=[LINK.get(w) for w in slots]
    ga=[]
    for w in vabs:
        cands=[a for a in byvar[w] if a in atomh and ('*x%d'%w in a or 'x%d*'%w in a)]
        ga.append(cands[0] if len(cands)==1 else (cands[0] if cands else None))
    atoms=[x for x in la+ga if x]
    if len(atoms)!=4: continue
    hs=sorted({atomh[a] for a in atoms if a in atomh})
    if len(hs)!=4: continue
    # M corrupts the DEFINED P-multiple h (the var appearing in the atom), not the free cofactor u
    hm=[]
    for a in atoms:
        u=atomh.get(a)
        cand=[v for v in set(vars_of(E.atoms[a])) if v in defrhs and fa(v)=={u}]
        hm.append(min(cand) if cand else None)
    if any(x is None for x in hm) or len(set(hm))!=4: continue
    hm=sorted(set(hm))
    u=set()
    for a in atoms: u|=eqs_of.get(a,set())
    nlive=len([x for x in sub[c] if x in liveset])
    rows.append(dict(site_child=c,parent=n,side=s,slot_wires=slots,vab_wires=vabs,
                     handles_h_Pmultiple=hm,handles_u_cofactor=hs,incidence=len(u),depth=depth[c],live_leaves_under=nlive))
import itertools
rows.sort(key=lambda r:(r['incidence'],-r['depth']))
D=[r for r in rows if r['site_child']==27994]
print('CALIBRATION row (the deliverable):',json.dumps(D[0]) if D else 'NOT FOUND')
print('total candidate sites of the deliverable shape:',len(rows))
json.dump(rows,open('candidates.json','w'),indent=0)
print('\nTop 25 by incidence (inflated - ordinal only):')
for r in rows[:25]:
    print('  handles %-28s site c=x%-6d parent x%-6d.%s  incid %-3d depth %d  liveleaves %d'%(
        ','.join('x%d'%h for h in r['handles_h_Pmultiple']),r['site_child'],r['parent'],r['side'],
        r['incidence'],r['depth'],r['live_leaves_under']))
