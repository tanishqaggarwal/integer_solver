#!/usr/bin/env python3
"""Q-9c: find EVERY chord-law stage in the circuit and verify each one on RANDOM curve points.

Gadget (division-free chord law; wires hold raw coords u = X - K/3, and 3*(K/3) = K):
   dx = ua - ub                dy = ya - yb
   R1 = S*dx^2 - dy^2   with  S = u3 + ua + ub + K        <=> lambda^2 = u3+ua+ub+K
   R2 = A*dx  - B*dy    with  A = y3 + yb,  B = ub - u3    <=> y3+yb = lambda*(ub-u3)
Both residuals are driven to 0 by the equations.  Verification is a Schwartz-Zippel identity test:
put RANDOM points of the cubic on the two inputs, set (u3,y3) from the group law, evaluate the
actual sub-DAG taken from EQUATIONS.txt, and require R1 = R2 = 0 mod p.
"""
import json,collections,random,sys,os,re
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from qgrp import add,neg,oncur,p,cs,A_,B_
G=[json.loads(l) for l in open('atoms/gates.jsonl')]
defs={}; uses=collections.defaultdict(list)
for d in G:
    defs.setdefault(d['t'],d)
    for v in d['vids']: uses[v].append(d)
SQ=re.compile(r'^x_(\d+) \* x_\1$'); SUB=re.compile(r'^x_(\d+) - x_(\d+)$')
MUL=re.compile(r'^x_(\d+) \* x_(\d+)$')
leaf={int(g):v for g,v in json.load(open('qleaf.json')).items()}
lad=json.load(open('qladder.json')); sel2exp={int(k):v for k,v in lad['sel2exp'].items()}
X2e={v[2]:sel2exp[int(g)] for g,v in leaf.items()}; Y2e={v[3]:sel2exp[int(g)] for g,v in leaf.items()}
LEAFW=set(X2e)|set(Y2e)
sqof={}                      # w -> wire holding w*w
for d in G:
    m=SQ.match(d['rhs'])
    if m: sqof.setdefault(int(m.group(1)),d['t'])
issq={t:w for w,t in sqof.items()}
def sub(w):
    d=defs.get(w); m=SUB.match(d['rhs']) if d else None
    return (int(m.group(1)),int(m.group(2))) if m else None
stages=[]
for d in G:                                     # R1 = m - sq(dy)
    m=SUB.match(d['rhs'])
    if not m: continue
    mw,s2=int(m.group(1)),int(m.group(2))
    if s2 not in issq: continue
    dmw=defs.get(mw); mm=MUL.match(dmw['rhs']) if dmw else None
    if not mm: continue
    a,b=int(mm.group(1)),int(mm.group(2))
    S,s1=(a,b) if b in issq else ((b,a) if a in issq else (None,None))
    if S is None: continue
    dx,dy=issq[s1],issq[s2]
    pa,pb=sub(dx),sub(dy)
    if not pa or not pb: continue
    stages.append({'R1':d['t'],'S':S,'dx':dx,'dy':dy,'ua':pa[0],'ub':pa[1],'ya':pb[0],'yb':pb[1]})
print('chord-law stage gadgets found:',len(stages))
def findR2(st):
    for d in uses[st['dx']]:
        mm=MUL.match(d['rhs'])
        if not mm: continue
        for e in uses[d['t']]:
            me=SUB.match(e['rhs'])
            if not me: continue
            o=int(me.group(2)) if int(me.group(1))==d['t'] else int(me.group(1))
            do=defs.get(o); mo=MUL.match(do['rhs']) if do else None
            if mo and st['dy'] in (int(mo.group(1)),int(mo.group(2))): return e['t']
    return None
def cone(root,cut):
    seen=set(); stk=[root]
    while stk:
        w=stk.pop()
        if w in seen or w in cut: continue
        seen.add(w)
        d=defs.get(w)
        if d:
            for v in d['vids']: stk.append(v)
    done=set(cut); out=[]; prog=True
    while prog:
        prog=False
        for w in seen:
            if w in done or w not in defs: continue
            if all(v in done for v in defs[w]['vids']): out.append(w); done.add(w); prog=True
    return out,seen
def ev(order,val):
    for w in order:
        d=defs[w]; e=d['rhs']
        for v in sorted(set(d['vids']),reverse=True): e=e.replace('x_%d'%v,'(%d)'%val[v])
        val[w]=eval(e,{'__builtins__':{}})%p
random.seed(11)
def randpt():
    while True:
        x=random.randrange(p); r=(pow(x,3,p)+A_*x+B_)%p
        y=pow(r,(p+1)//4,p)
        if y*y%p==r: return (x,y)
kinds=collections.Counter(); vk=collections.Counter(); orients=collections.Counter(); bad=[]
for st in stages:
    kind=('leaf-adjacent' if (st['ua'] in LEAFW and st['ub'] in LEAFW) else
          'mixed' if (st['ua'] in LEAFW or st['ub'] in LEAFW) else 'internal')
    kinds[kind]+=1
    cut={st['ua'],st['ub'],st['ya'],st['yb']}
    R2=findR2(st)
    o1,c1=cone(st['R1'],cut); o2,c2=cone(R2,cut) if R2 else ([],set())
    free=sorted(w for w in (c1|c2) if w not in defs and w not in cut)
    if len(free)!=2: bad.append((st['R1'],kind,'free=%d'%len(free))); continue
    hit=None
    for trial in range(2):
        Pa=randpt(); Pb=randpt()
        for sa in (1,-1):
            for sb in (1,-1):
                P3=add(Pa if sa>0 else neg(Pa),Pb if sb>0 else neg(Pb))
                if P3 is None: continue
                for u3,y3 in ((free[0],free[1]),(free[1],free[0])):
                    val=collections.defaultdict(int)
                    val[st['ua']]=(Pa[0]-cs)%p; val[st['ya']]=Pa[1]%p
                    val[st['ub']]=(Pb[0]-cs)%p; val[st['yb']]=Pb[1]%p
                    val[u3]=(P3[0]-cs)%p; val[y3]=P3[1]%p
                    ev(o1,val)
                    if R2: ev(o2,val)
                    if val[st['R1']]%p==0 and (R2 is None or val[R2]%p==0):
                        hit=(sa,sb); break
                if hit: break
            if hit: break
        if hit: break
    if hit: vk[kind]+=1; orients[hit]+=1; st['u3'],st['y3']=u3,y3; st['kind']=kind
    else: bad.append((st['R1'],kind,'no orientation'))
print('classification :',dict(kinds))
print('VERIFIED chord law on random points :',dict(vk))
print('input orientations (sign of each input point):',dict(orients))
print('unverified:',len(bad)); 
for b in bad[:8]: print('   ',b)
json.dump({'stages':stages},open('qstages.json','w'))
