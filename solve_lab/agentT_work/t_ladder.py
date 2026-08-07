#!/usr/bin/env python3
"""AUDIT T11 -- agent Q's EXISTENCE chain, the premises Q did not flag.

Q's step 5: 'Subset sums of {2^i G}_{i=0..255} realise kG for EVERY k in [1,2^256-1].  N < 2^256,
so every group element is hit -- including T.  A satisfying assignment therefore EXISTS.'

That is airtight ONLY IF the leaf set is exactly {2^i G : i = 0..255} with every exponent present
and distinct.  Drop one exponent and the subset sums cover 2^255 of the 2^256 -- and whether T
survives is a coin flip, not a theorem.  Q's flagged caveat is the leaf-adjacent stage law; this
tests the ladder itself.

RESUME_Q §1(d): '249/253 decoded leaves have their double also a leaf; the doubling graph is
4 chains of 124/79/41/9 which link tail->head by exactly one missing doubling each.  Result,
VERIFIED EXACTLY: L_i == 2^i G for i = 0..255.'  The four chains are decoded; the three links
between them are the INFERENCE.  Verify every doubling directly."""
import json,os,collections
Q='/home/user/integer_solver/solve_lab/agentQ_work'
cur=json.load(open(os.path.join(Q,'curve.json')))
p=int(cur['p']); a=int(cur['a']); b=int(cur['b']); c=int(cur['c_shift'])
L=json.load(open(os.path.join(Q,'ladder.json')))
lad=L['ladder']; missing=L['missing']; N=int(L['N'])
qleaf=json.load(open(os.path.join(Q,'qleaf.json')))
def on(P):
    X,Y=P; return (Y*Y-X**3-a*X-b)%p==0
def add(P1,P2):
    if P1 is None: return P2
    if P2 is None: return P1
    x1,y1=P1; x2,y2=P2
    if x1==x2 and (y1+y2)%p==0: return None
    l=(3*x1*x1+a)*pow(2*y1,p-2,p)%p if P1==P2 else (y2-y1)*pow((x2-x1)%p,p-2,p)%p
    x3=(l*l-x1-x2)%p; return (x3,(l*(x1-x3)-y1)%p)
def mul(k,P):
    R=None
    while k:
        if k&1: R=add(R,P)
        P=add(P,P); k>>=1
    return R
LP={int(e):(int(qleaf[str(v)][0])%p,int(qleaf[str(v)][1])%p) for e,v in lad.items()}
print('ladder entries: %d   missing exponents: %s   qleaf leaves: %d'%(len(LP),missing,len(qleaf)))
print('all ladder points on the cubic: %s'%all(on(P) for P in LP.values()))

print('\n== 1. every consecutive doubling L_{i+1} == 2*L_i ==')
ok=bad=skip=0; badlist=[]
for i in range(255):
    if i not in LP or i+1 not in LP: skip+=1; continue
    if mul(2,LP[i])==LP[i+1]: ok+=1
    else: bad+=1; badlist.append(i)
print('   verified %d   FAILED %d %s   skipped (an endpoint is a missing exponent) %d'%(ok,bad,badlist[:6],skip))

print('\n== 2. are the 253 decoded ladder points DISTINCT? ==')
cnt=collections.Counter(LP.values())
dup=[(P,n) for P,n in cnt.items() if n>1]
print('   distinct points: %d of %d   duplicates: %d'%(len(cnt),len(LP),len(dup)))

print('\n== 3. the 3 leaves NOT placed in the ladder ==')
inlad={int(v) for v in lad.values()}
extra=[k for k in qleaf if int(k) not in inlad]
print('   selector vars decoded but unplaced: %s'%extra)
G=LP[0]
for k in extra:
    P=(int(qleaf[k][0])%p,int(qleaf[k][1])%p)
    hit=[e for e in missing if mul(pow(2,e,N),G)==P]
    print('   var %-7s on-curve %-5s  equals 2^e G for e in missing? %s'%(k,on(P),hit if hit else 'NO'))

print('\n== 4. does the full ladder close: L_i == 2^i G for every i = 0..255? ==')
G=LP[0]; bad2=[]
for i in range(256):
    if i in LP and mul(pow(2,i,N),G)!=LP[i]: bad2.append(i)
print('   mismatches among the 253 decoded: %s'%(bad2 if bad2 else 'none'))
print('   => exponents actually backed by a decoded leaf: %d of 256'%len(LP))

print('\n== 5. group order ==')
print('   N*G == O ?  %s'%(mul(N,G) is None))
print('   N < 2^256 ? %s   2^256-1 - N = %d'%(N< 2**256, 2**256-1-N))
print('   N prime-ish: 2^N-1 Fermat  %s'%(pow(2,N-1,N)==1))
print('   [1,2^256-1] covers every residue mod N ? %s  (needs N-1 <= 2^256-1)'%(N-1<=2**256-1))
