#!/usr/bin/env python3
"""Agent X standing SIGNED validation (design supplied by agent Z; pass criterion from agent Y).

Replaces the vacuous 1-term plant in srep_c.txt.  That test was worthless: a 1-term k makes
EVERY scan point a genuine 2-term table hit, so it passed independently of sign handling.

Here the plant has m=5, so a b=2 scan forces the table to supply exactly the other 3 digits.
PASS requires the EXACT set of C(5,2)=10 splits to appear -- not merely that some hit appeared.
Three sign patterns are run so that sign bookkeeping is actually exercised:
  (a) lowest digit NEGATIVE   (b) ALL digits negative   (c) all positive (control)
"""
import json,os,subprocess,sys,itertools
d=json.load(open('xdata.json')); p=int(d['p']); A_=int(d['a']); N=int(d['N'])
lad=[(int(a),int(b)) for a,b in d['ladder']]
def inv(z):return pow(z,p-2,p)
def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x1,y1=P;x2,y2=Q
    if (x1-x2)%p==0:
        if (y1+y2)%p==0:return None
        l=(3*x1*x1+A_)%p*inv(2*y1%p)%p
    else: l=(y2-y1)%p*inv((x2-x1)%p)%p
    x3=(l*l-x1-x2)%p;return (x3,(l*(x1-x3)-y1)%p)
def mul(k,P):
    k%=N; R=None;Q=P
    while k>0:
        if k&1:R=add(R,Q)
        Q=add(Q,Q);k>>=1
    return R
E=[11,63,128,190,240]
CASES={'lowest_negative':[-1,1,1,1,1],'all_negative':[-1,-1,-1,-1,-1],'all_positive_control':[1,1,1,1,1]}
allpass=True
for name,eps in CASES.items():
    terms=list(zip(E,eps)); k=sum(s*(1<<e) for e,s in terms)
    T=mul(k%N,lad[0])
    fn='splant_z_%s.txt'%name; rep='srep_z_%s.txt'%name
    open(fn,'w').write('%d %d\n'%T + ''.join('%d %d\n'%(x,y) for x,y in lad))
    if os.path.exists(rep): os.remove(rep)
    subprocess.run(['./xsigned','scan',fn,'2','stbls.bin','sbm.bin',rep],
                   stderr=subprocess.DEVNULL,check=True)
    hits=[l.split() for l in open(rep) if l.startswith('HIT')]
    got=set()
    for h in hits:
        code=int(h[2]); sc=[(code>>(16*i))&0xFFFF for i in range(2)]
        got.add(tuple(sorted((s>>1, -1 if (s&1) else 1) for s in sc)))
    # expected: for each 2-subset of the planted digits, sigma = -eps
    exp=set()
    for pair in itertools.combinations(terms,2):
        exp.add(tuple(sorted((e,-s) for e,s in pair)))
    uw=bin(k%N).count('1')
    ok = (len(hits)==10) and (got==exp)
    allpass &= ok
    print('%-22s  unsigned wt(k mod N)=%3d   HIT lines=%2d (want 10)   exact splits %2d/%2d   %s'%(
          name,uw,len(hits),len(got&exp),len(exp),'PASS' if ok else 'FAIL'))
    if not ok:
        print('   missing:',sorted(exp-got)); print('   extra  :',sorted(got-exp))
print()
print('OVERALL SIGNED VALIDATION:','PASS' if allpass else 'FAIL')

# --- why fixing the table's leading sign is LOSSLESS (Z's point; verify, do not assert) ---
import random
random.seed(11); M=(1<<64)-1; bad=0
for _ in range(200):
    a=sorted(random.sample(range(256),3))
    sg=[random.choice((1,-1)) for _ in a]
    P=None
    for e,s in zip(a,sg):
        Q=lad[e] if s>0 else (lad[e][0],(-lad[e][1])%p)
        P=add(P,Q)
    Pn=(P[0],(-P[1])%p)          # negation flips EVERY sign, incl. the leading one
    if (P[0]&M)!=(Pn[0]&M): bad+=1
print('leading-sign restriction lossless: x(P)==x(-P) on 200 random signed 3-term sums, mismatches=%d'%bad)
print('  => every leading-NEGATIVE signed sum has the same 64-bit key as its leading-POSITIVE negation,')
print('     so the table storing only leading-positive representatives loses no key.')
