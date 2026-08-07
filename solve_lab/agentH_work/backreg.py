"""BACKWARDS step 1: enumerate target regions R (<=12 equations) and score them by
   un-cancellable rows and by the balance count |R| - |S(R)| + c."""
import ev, json, time, pickle
from collections import defaultdict, Counter
eq_atoms=[]            # equation -> set of atoms
atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    s=set(a for c,a in tl)
    eq_atoms.append(s)
    for a in s: atom_eqs[a].add(i)
NA=len(ev.atom_src)
print('atoms',NA,'equations',len(eq_atoms))
small=[a for a in atom_eqs if 1<=len(atom_eqs[a])<=12]
print('atoms living in <=12 equations:',len(small))

def score_region(R):
    """R: frozenset of equations.  S = atoms whose whole eq-footprint is inside R."""
    cand=set()
    for e in R: cand|=eq_atoms[e]
    S=[a for a in cand if atom_eqs[a]<=R]
    if not S: return None
    Sset=set(S)
    unc=sum(1 for e in R if len(eq_atoms[e]&Sset)==1)
    return (len(R),len(S),unc,S)

t0=time.time()
best=[]
seen=set()
for a in small:
    R=frozenset(atom_eqs[a])
    if R in seen: continue
    seen.add(R)
    r=score_region(R)
    if r is None: continue
    nR,nS,unc,S=r
    best.append((nR-nS,unc,nR,nS,tuple(sorted(S))))
best.sort()
print('distinct regions from single-atom footprints:',len(seen),'%.0fs'%(time.time()-t0))
print('\ntop 25 by  |R|-|S|  (lower is better; failing >= |R|-|S|+c):')
for b in best[:25]:
    print('  |R|-|S|=%3d  uncancellable=%2d  |R|=%2d |S|=%2d  atoms=%s'%(b[0],b[1],b[2],b[3],list(b[4])[:9]))
print('\ntop 25 by uncancellable rows:')
for b in sorted(best,key=lambda t:(t[1],t[0]))[:25]:
    print('  unc=%2d  |R|-|S|=%3d  |R|=%2d |S|=%2d  atoms=%s'%(b[1],b[0],b[2],b[3],list(b[4])[:9]))
# where does the deliverable's region sit?
DEL={22229,22230,35758,35759,35760,35761,35762}
RD=set()
for a in DEL: RD|=atom_eqs[a]
r=score_region(frozenset(RD))
print('\nDELIVERABLE region: |R|=%d |S|=%d uncancellable=%d  S=%s'%(r[0],r[1],r[2],sorted(r[3])))
pickle.dump(best,open('backreg.pkl','wb'))
