#!/usr/bin/env python3
"""Q-11: solve the mux layer SYMBOLICALLY (do not propagate through it).

Read off one slot completely, then check the same shape at all 383.  The claim under test:
    Xout = cA*Xa + cB*Xb + cC*u3 ,  Yout = cA*Ya + cB*Yb + cC*y3
    cA = s1*(1-s2)   cB = s2*(1-s1)   cC = s1*s2      s1,s2 boolean
so the four quadrants are  (0,0)->(0,0) identity   (1,0)->A   (0,1)->B   (1,1)->chord(A,B)=A+B.
"""
import pickle,collections,re,json
terms=pickle.load(open('qterms.pkl','rb'))
PROD=re.compile(r'^x_(\d+) - x_(\d+) \* x_(\d+)$')
SUM =re.compile(r'^x_(\d+) - \(x_(\d+) \+ x_(\d+)\)$')
NOT =re.compile(r'^x_(\d+) - \(1 - x_(\d+)\)$')
EQ  =re.compile(r'^x_(\d+) - x_(\d+)$')
BOOL=[re.compile(r'^x_(\d+) \* \(x_(\d+) - 1\)$'),re.compile(r'^x_(\d+) \* x_(\d+) - x_(\d+)$')]
prod={}; summ={}; notof={}; alias={}; boolw=set()
prod_by_factor=collections.defaultdict(list); sum_by_arg=collections.defaultdict(list)
for s,vs in terms:
    m=PROD.match(s)
    if m:
        t,u,v=map(int,m.groups()); prod[t]=(u,v)
        prod_by_factor[u].append(t); prod_by_factor[v].append(t); continue
    m=SUM.match(s)
    if m:
        t,u,v=map(int,m.groups()); summ[t]=(u,v)
        sum_by_arg[u].append(t); sum_by_arg[v].append(t); continue
    m=NOT.match(s)
    if m: t,u=map(int,m.groups()); notof[t]=u; continue
    m=EQ.match(s)
    if m: t,u=map(int,m.groups()); alias[t]=u; continue
    for B in BOOL:
        m=B.match(s)
        if m and len(set(m.groups()))==1: boolw.add(int(m.group(1)))
def root(w,seen=None):
    seen=seen or set()
    while w in alias and w not in seen: seen.add(w); w=alias[w]
    return w
ST=[x for x in json.load(open('qstages.json'))['stages'] if 'u3' in x]
res=collections.Counter(); detail=[]
for g in ST:
    ua,ub,u3=g['ua'],g['ub'],g['u3']
    ok=False; why='no cC*u3 product'
    for t in prod_by_factor.get(u3,[]):
        u,v=prod[t]; cC = v if u==u3 else u
        # Xout = (cA*Xa + cB*Xb) + cC*u3
        cands=set()
        fr=[t]
        for _ in range(3):
            nxt=[]
            for x in fr:
                for y in sum_by_arg.get(x,[]): cands.add(y); nxt.append(y)
            fr=nxt
        for Xout in cands:
            def expand(w,d=0):
                if d>4 or w not in summ: return [w]
                a,b=summ[w]; return expand(a,d+1)+expand(b,d+1)
            leaves=expand(Xout)
            if t not in leaves or len(leaves)!=3: continue
            others=[x for x in leaves if x!=t]
            if any(x not in prod for x in others): continue
            pair={}
            for rr in others:
                aa,bb=prod[rr]
                if aa in (ua,ub): pair[aa]=bb
                elif bb in (ua,ub): pair[bb]=aa
            if set(pair)!={ua,ub}: continue
            cA,cB=pair[ua],pair[ub]
            if cC not in prod: why='cC not a product'; continue
            s1,s2=[root(x) for x in prod[cC]]
            def isand_not(c,s,other):
                if c not in prod: return False
                x,y=[root(z) for z in prod[c]]
                for A,B in ((x,y),(y,x)):
                    if A==s and B in [root(k) for k in notof if root(notof[k])==other and root(k)==B]: return True
                    if A==s and (B in notof and root(notof[B])==other): return True
                return False
            okA=isand_not(cA,s1,s2) or isand_not(cA,s2,s1)
            okB=isand_not(cB,s1,s2) or isand_not(cB,s2,s1)
            bl=(s1 in boolw and s2 in boolw)
            if okA and okB:
                ok=True; why='ok' if bl else 'ok (selectors not both boolean-pinned)'
                res['quadrant mux confirmed'+('' if bl else ' [sel not boolean]')]+=1
                break
        if ok: break
    if not ok: res['NOT matched: '+why]+=1; detail.append(g['u3'])
print('slots: %d'%len(ST))
for k,v in res.most_common(): print('   %-55s %d'%(k,v))
print('unmatched output wires (first 10):',detail[:10])
