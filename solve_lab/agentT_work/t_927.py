#!/usr/bin/env python3
"""AUDIT T17 -- item 3: is the 927 decomposition-dependent?

Agent L: of its 3,681 residual atoms, the measured multiplier c is 1 for 2,747 and >1 for 927.
c>1 means the integrality condition c*p | R is STRICTLY STRONGER than R == 0 mod p, so the lift
to Z is not free.  L's 927 matches P's 927 from an unshared decomposition.  Two agents agreeing
is evidence, but per my B1 finding a count matching across models is not the same as a count
being intrinsic -- both could be reading the same literal off the same file text.

Rebuild it in agent F's 39,033-atom decomposition, the one certified faithful in audit T2:
  u (free cofactor) -> its unique atom, of shape  (h - (P * u))   [P is a p-valued wire]
  h -> the OTHER atom containing h, the GUARD
  c := the literal multiplying h INSIDE THE GUARD (bare h => c = 1)
Then count c > 1.  Independent of L's and P's atom algebra; only the cofactor list is borrowed."""
import os,sys,pickle,collections,re,json
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
from circ2 import vars_of
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
v2a=collections.defaultdict(list)
for i,a in enumerate(names):
    for u in vars_of(atoms[a]): v2a[u].append(i)
H=pickle.load(open(os.path.join(LAB,'agentL_work','handles.pkl'),'rb'))
U=sorted(set(H['handle']))
print("cofactors u borrowed from L: %d"%len(U))
defpat=re.compile(r'^\(x_(\d+)-\(x_(\d+)\*x_(\d+)\)\)$')
stats=collections.Counter(); cvals=collections.Counter(); rows=[]
noguard=[]; oddshape=[]
for u in U:
    ai=v2a[u]
    if len(ai)!=1: stats['u not in exactly 1 atom']+=1; continue
    s=names[ai[0]].replace(' ','')
    m=defpat.match(s)
    if not m: oddshape.append((u,s)); stats['definition not (h-(P*u))']+=1; continue
    h=int(m.group(1)); Pv=int(m.group(2)); uu=int(m.group(3))
    if uu!=u: h,Pv=int(m.group(1)),int(m.group(2))
    others=[j for j in v2a[h] if j!=ai[0]]
    if len(others)!=1: noguard.append((u,h,len(others))); stats['h not in exactly 2 atoms']+=1; continue
    g=names[others[0]].replace(' ','')
    mm=re.search(r'\((\d+)\*x_%d\)'%h, g)
    if mm: c=int(mm.group(1))
    elif re.search(r'(?<![0-9_])x_%d(?![0-9])'%h, g): c=1
    else: stats['h not found in guard']+=1; continue
    cvals[c]+=1; stats['ok']+=1
    rows.append((u,h,Pv,c,g))
print('\nparse outcome:', dict(stats))
if oddshape: print('  sample odd definition shapes:', oddshape[:3])
if noguard: print('  sample h-without-unique-guard:', noguard[:3])
c1=sum(n for c,n in cvals.items() if c==1)
cg=sum(n for c,n in cvals.items() if c>1)
print('\n== RESULT, rebuilt in F\'s 39,033-atom parse ==')
print('   c == 1 (condition collapses to the mod-p congruence) : %d'%c1)
print('   c >  1 (a genuine extra integer condition)           : %d'%cg)
print('   total parsed                                          : %d'%(c1+cg))
print('\n   agent L reports 2,747 / 927 ; agent P reports 927.')
print('   MATCH on c>1 ?  %s'%('YES' if cg==927 else 'NO  (F-parse says %d)'%cg))
print('   MATCH on c==1 ? %s'%('YES' if c1==2747 else 'NO  (F-parse says %d)'%c1))
print('\n   distinct multiplier values: %d ; 10 commonest: %s'%(len(cvals),cvals.most_common(10)))
json.dump({'c1':c1,'cgt1':cg,'multipliers':{str(k):v for k,v in cvals.items()}},
          open(os.path.join(T,'t_927.json'),'w'))
