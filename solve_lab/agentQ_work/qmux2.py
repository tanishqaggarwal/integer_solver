#!/usr/bin/env python3
"""Q-11c: close the mux layer at all 383 slots, with a cross-check.

A slot matches only if BOTH coordinate muxes are found with the IDENTICAL coefficient wires:
    Xout = cA*Xa + cB*Xb + cC*Xout_chord
    Yout = cA*Ya + cB*Yb + cC*Yout_chord
    cC = s1*s2 ,  cA = s1*(1-s2) ,  cB = s2*(1-s1)
Requiring the same cA,cB,cC on both coordinates rules out accidental structural matches.
"""
import pickle,collections,re,json
terms=pickle.load(open('qterms.pkl','rb'))
PROD=re.compile(r'^x_(\d+) - x_(\d+) \* x_(\d+)$')
SUM =re.compile(r'^x_(\d+) - \(x_(\d+) \+ x_(\d+)\)$')
NOT =re.compile(r'^x_(\d+) - \(1 - x_(\d+)\)$')
EQ  =re.compile(r'^x_(\d+) - x_(\d+)$')
BOOL=[re.compile(r'^x_(\d+) \* \(x_(\d+) - 1\)$'),re.compile(r'^x_(\d+) \* x_(\d+) - x_(\d+)$')]
prod={}; summ={}; notof={}; alias={}; boolw=set()
pbf=collections.defaultdict(list); sba=collections.defaultdict(list)
for s,vs in terms:
    m=PROD.match(s)
    if m:
        t,u,v=map(int,m.groups()); prod[t]=(u,v); pbf[u].append(t); pbf[v].append(t); continue
    m=SUM.match(s)
    if m:
        t,u,v=map(int,m.groups()); summ[t]=(u,v); sba[u].append(t); sba[v].append(t); continue
    m=NOT.match(s)
    if m: t,u=map(int,m.groups()); notof[t]=u; continue
    m=EQ.match(s)
    if m: t,u=map(int,m.groups()); alias[t]=u; continue
    for B in BOOL:
        m=B.match(s)
        if m and len(set(m.groups()))==1: boolw.add(int(m.group(1)))
def root(w):
    seen=set()
    while w in alias and w not in seen: seen.add(w); w=alias[w]
    return w
def expand(w,d=0):
    if d>5 or w not in summ: return [w]
    a,b=summ[w]; return expand(a,d+1)+expand(b,d+1)
def mux(pa,pb,po):
    """find (cA,cB,cC) with  out = cA*pa + cB*pb + cC*po ."""
    out=[]
    for t in pbf.get(po,[]):
        u,v=prod[t]; cC = v if u==po else u
        cands=set(); fr=[t]
        for _ in range(3):
            nxt=[]
            for x in fr:
                for y in sba.get(x,[]): cands.add(y); nxt.append(y)
            fr=nxt
        for X in cands:
            lv=expand(X)
            if t not in lv or len(lv)!=3: continue
            oth=[x for x in lv if x!=t]
            if any(x not in prod for x in oth): continue
            pair={}
            for rr in oth:
                aa,bb=prod[rr]
                if aa in (pa,pb): pair[aa]=bb
                elif bb in (pa,pb): pair[bb]=aa
            if set(pair)!={pa,pb}: continue
            out.append((pair[pa],pair[pb],cC,X))
    return out
def isnot(n,s):  return n in notof and root(notof[n])==root(s)
def quad(cA,cB,cC):
    if cC not in prod or cA not in prod or cB not in prod: return None
    s1,s2=[root(x) for x in prod[cC]]
    for (S1,S2) in ((s1,s2),(s2,s1)):
        a1,a2=prod[cA]; b1,b2=prod[cB]
        okA=(root(a1)==S1 and isnot(a2,S2)) or (root(a2)==S1 and isnot(a1,S2))
        okB=(root(b1)==S2 and isnot(b2,S1)) or (root(b2)==S2 and isnot(b1,S1))
        if okA and okB: return (S1,S2)
    return None
ST=[x for x in json.load(open('qstages.json'))['stages'] if 'u3' in x]
res=collections.Counter(); bad=[]
for g in ST:
    U,V=(g['ua'],g['ub'],g['u3']),(g['ya'],g['yb'],g['y3'])
    W,Z=(g['ua'],g['ub'],g['y3']),(g['ya'],g['yb'],g['u3'])
    hit=None
    for P,Q in ((U,V),(V,U),(W,Z),(Z,W)):
        for cA,cB,cC,_ in mux(*P):
            q=quad(cA,cB,cC)
            if not q: continue
            if any((cA,cB,cC)==(a2,b2,c2) for a2,b2,c2,_ in mux(*Q)):
                hit=(q,cA,cB,cC); break
        if hit: break
    if hit:
        s1,s2=hit[0]
        res['BOTH coordinate muxes, same coefficients'+(' [selectors boolean-pinned]' if (s1 in boolw and s2 in boolw) else '')]+=1
    else: res['unmatched']+=1; bad.append(g['u3'])
print('slots: %d'%len(ST))
for k,v in res.most_common(): print('   %-58s %d'%(k,v))
if bad: print('unmatched slot outputs:',bad[:20])
