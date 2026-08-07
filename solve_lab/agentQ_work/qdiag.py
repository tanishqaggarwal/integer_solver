import json,collections,random,sys
sys.path.insert(0,'.')
import qsolve as Q
from qgrp import cs
leaf={int(g):v for g,v in json.load(open('qleaf.json')).items()}
LP={g:(int(v[0]),int(v[1])) for g,v in leaf.items()}; LW={g:(v[2],v[3]) for g,v in leaf.items()}
random.seed(3); S=set(random.sample(sorted(LP),64))
V=[None]*Q.NV
for g in LP:
    V[g]=1 if g in S else 0
    wx,wy=LW[g]; V[wx]=(LP[g][0]-cs)%Q.p; V[wy]=LP[g][1]%Q.p
s,c=Q.propagate(V)
gate=set()
for i,(st,vs) in enumerate(Q.terms): pass
G=[json.loads(l) for l in open('atoms/gates.jsonl')]
tg={d['t'] for d in G}
free=[v for v in range(Q.NV) if v not in tg]
print('solved %d  contradictions %d'%(s,c))
print('wires known: %d / %d'%(sum(1 for v in V if v is not None),Q.NV))
print('free inputs known: %d / %d'%(sum(1 for v in free if V[v] is not None),len(free)))
ST=[x for x in json.load(open('qstages.json'))['stages'] if 'u3' in x]
print('gadget outputs known: %d / %d'%(sum(1 for x in ST if V[x['u3']] is not None),len(ST)))
# terms blocked with exactly 2 unknowns
c2=collections.Counter()
for st,vs in Q.terms:
    u=[v for v in vs if V[v] is None]
    if len(u)==2: c2['2-unknown']+=1
    elif len(u)>2: c2['3+']+=1
    elif len(u)==1: c2['1']+=1
print('terms by unknown count:',dict(c2))
