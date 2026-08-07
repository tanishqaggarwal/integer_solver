#!/usr/bin/env python3
"""Agent X: validate the 128-rotation splitting system BEFORE trusting any negative from it.
A_j = {(j+t) mod 256 : t = 0..127}.  Claim: every 10-subset S has some j in [0,128) with
|S n A_j| = 5 (hence |S n B_j| = 5)."""
import random, itertools
def f(S,j): return sum(1 for i in S if ((i-j) % 256) < 128)
def bal(S): return [j for j in range(128) if f(S,j)==5]
random.seed(3)
bad=0; nb=[]
for _ in range(20000):
    S=random.sample(range(256),10); b=bal(S)
    if not b: bad+=1
    nb.append(len(b))
print('random 10-sets tested: 20000   with NO balanced rotation: %d'%bad)
print('balanced rotations per set: min %d  mean %.1f  max %d'%(min(nb),sum(nb)/len(nb),max(nb)))
# adversarial families
tight=list(range(10))
print()
print('TIGHT set S = {0..9}  -> balanced rotations:',bal(tight),'  (unique; +-1 in the rotation index MISSES it)')
print('  f(j) for j=0..11:',[f(tight,j) for j in range(12)])
anti=[i for k in range(5) for i in (k,k+128)]
print('ANTIPODAL set S =',sorted(anti),' -> #balanced rotations:',len(bal(anti)),'(all of them: S is +128-invariant, f == 5)')
print('  f(j) for j=0..5:',[f(anti,j) for j in range(6)])
# is 4|6 ever unavailable?  (justifies choosing the 5|5 split)
print()
print('ANTIPODAL set: any rotation with |S n A_j| = 4 ?', any(f(anti,j)==4 for j in range(256)),
      '  -> the 4|6 split does NOT always exist, so 5|5 is forced')
# exhaustive-ish worst case: sets inside a short arc
worst=[]
for w in range(10,40):
    S=list(range(0,w,max(1,w//10)))[:10]
    if len(S)==10: worst.append((w,len(bal(S))))
print('sets inside a short arc: (arc width, #balanced rotations) =',worst[:8])
print()
print('=> the covering claim holds on every case tested; a UNIQUE-rotation set exists ({0..9}, j=5),')
print('   which is exactly the planted test an off-by-one in the rotation index would fail.')
