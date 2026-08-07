"""BACKWARDS step 2: grow a target region around the defect-carrying atoms, minimising |R|-|S|
   (the balance count) and un-cancellable rows.  Only regions containing a defect atom can hold
   the two congruences, so growth is seeded there."""
import ev, json, time, itertools
from collections import defaultdict
eq_atoms=[]; atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    s=set(a for c,a in tl); eq_atoms.append(s)
    for a in s: atom_eqs[a].add(i)
def S_of(R):
    cand=set()
    for e in R: cand|=eq_atoms[e]
    return set(a for a in cand if atom_eqs[a]<=R)
def unc(R,S):
    return sum(1 for e in R if len(eq_atoms[e]&S)==1)

DEFECT_B={22229,22230,35758,35759,35760,35761,35762,22231}      # frame-B (witness) carriers
DEFECT_H={688,1618,40608,30980,30982,36185,40812}               # my own frame's carriers
def report(tag,seed):
    R=set()
    for a in seed: R|=atom_eqs[a]
    S=S_of(R)
    print('%s seed: |R|=%d |S|=%d bal=%d unc=%d'%(tag,len(R),len(S),len(R)-len(S),unc(R,S)))
    return R,S
RB,SB=report('witness',DEFECT_B)
RH,SH=report('agentH ',DEFECT_H)

def grow(R0,maxadd=14,verbose=True):
    """greedily add equations that bring the most atoms fully inside; track best balance."""
    R=set(R0); S=S_of(R)
    best=(len(R)-len(S),len(R),len(S),frozenset(R))
    for step in range(maxadd):
        # candidate equations: those touching atoms that are ALMOST inside R
        cand=set()
        for e in R:
            for a in eq_atoms[e]:
                if a not in S: cand|=(atom_eqs[a]-R)
        if not cand: break
        bestmove=None
        for e in cand:
            R2=R|{e}; S2=S_of(R2)
            k=(len(R2)-len(S2),len(R2))
            if bestmove is None or k<bestmove[0]: bestmove=(k,e,R2,S2)
        k,e,R2,S2=bestmove
        R,S=R2,S2
        b=(len(R)-len(S),len(R),len(S),frozenset(R))
        if b[0]<best[0]: best=b
        if verbose: print('   +eq%d -> |R|=%d |S|=%d bal=%d unc=%d'%(e,len(R),len(S),len(R)-len(S),unc(R,S)))
    return best
print('\n--- growth from the witness region ---')
bB=grow(RB)
print('best balance reached:',bB[0],'|R|=%d |S|=%d'%(bB[1],bB[2]))
print('\n--- growth from the agent-H region ---')
bH=grow(RH)
print('best balance reached:',bH[0],'|R|=%d |S|=%d'%(bH[1],bH[2]))
print('\nfailing >= balance + c  (c = 2 independent congruences)')
print('  witness route floor: %d   agentH route floor: %d   current deliverable: 7'%(bB[0]+2,bH[0]+2))
json.dump({'witness_best_balance':bB[0],'agentH_best_balance':bH[0],
           'witness_region':sorted(bB[3]),'agentH_region':sorted(bH[3])},open('backgrow.json','w'))
