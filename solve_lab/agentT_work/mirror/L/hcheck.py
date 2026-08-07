"""Distinguish the free cofactor u from the defined P-multiple h. Check P's 927."""
import sys, pickle, collections, re
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentT_work/mirror/F')
from fwd import Engine,NV
from circ2 import vars_of
from parse import node_str
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
SL=pickle.load(open('slopes.pkl','rb'))
Hd=pickle.load(open('handles.pkl','rb')); handle=set(Hd['handle'])
print('my "handle" class size: %d  -- are these FREE vars? %s'%(len(handle), all(h not in defrhs for h in handle)))
# occurrence counts
occ_atoms=collections.Counter(); occ_defs=collections.Counter()
for a in E.res:
    for u in set(vars_of(E.atoms[a])): occ_atoms[u]+=1
for w,r in defrhs.items():
    for u in set(vars_of(r)): occ_defs[u]+=1
print('free cofactors u (my handle class): atom-occurrence histogram',
      collections.Counter(occ_atoms[h] for h in handle).most_common())
print('   ... plus definition-occurrences histogram',
      collections.Counter(occ_defs[h] for h in handle).most_common())
# the DEFINED p-multiples: vars defined as p*u  (or c*p*u)
pm=[w for w,r in defrhs.items() if r[0]=='*' and (('c',p) in (r[1],r[2]))]
print('defined vars of the form p*<var>: %d'%len(pm))
print('   their TOTAL occurrence histogram (defs+atoms)',
      collections.Counter(occ_atoms[w]+occ_defs[w] for w in pm).most_common())
# multiplier c per atom, from the measured slope
cs=collections.Counter()
for a,s in SL.items():
    if s==0: cs['slope0']+=1; continue
    if s% p: cs['not_p_mult']+=1; continue
    cs[abs(s)//p]+=1
c1=cs[1]; cgt=sum(v for k,v in cs.items() if isinstance(k,int) and k>1)
print('atoms with a single free cofactor: %d'%len(SL))
print('   c == 1 : %d   (integer condition collapses to the mod-p congruence)'%c1)
print('   c  > 1 : %d   (real extra condition  c*p | R)'%cgt)
print('   slope 0: %d   other: %d'%(cs['slope0'],cs['not_p_mult']))
