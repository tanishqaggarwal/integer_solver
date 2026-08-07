#!/usr/bin/env python3
"""Solve mod many primes, save the vectors, CRT them into a solution modulo a large composite Q,
verify all 39,033 equations vanish mod Q, then test the balanced integer representative."""
import sys,os,json,time,pickle
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
from modm import solve_modm, ModM
E=Engine()
PRIMES=[1000003,1000033,1000037,1000039,1000081,1000099,1000117,1000121,
        1000133,1000151,1000159,1000171,1000183,1000187,1000193,1000199,
        1000211,1000213,1000231,1000249]
out=os.path.join(HERE,'crt_sols'); os.makedirs(out,exist_ok=True)
vecs=[]
for q in PRIMES:
    f=os.path.join(out,'v_%d.json'%q)
    if os.path.exists(f):
        vecs.append((q,json.load(open(f)))); print('cached',q,flush=True); continue
    t0=time.time()
    v,M,ok=solve_modm(E,q)
    r=M.run(v); bad=M.score(r)
    assert ok and not bad, (q,ok,len(bad))
    json.dump(v,open(f,'w'))
    vecs.append((q,v))
    print('solved mod %d (0 failing) t=%.1f'%(q,time.time()-t0),flush=True)
Q=1
for q,_ in vecs: Q*=q
print('Q has %d bits'%Q.bit_length(),flush=True)
# CRT
def crt_pair(a1,m1,a2,m2):
    from math import gcd
    def egcd(a,b):
        if b==0: return (a,1,0)
        g,x,y=egcd(b,a%b); return (g,y,x-(a//b)*y)
    g,s,t=egcd(m1,m2)
    return (a1 + m1*(((a2-a1)//g*s) % (m2//g))) % (m1*m2//g)
V=[0]*NV
mods=[q for q,_ in vecs]
t0=time.time()
for i in range(NV):
    a=vecs[0][1][i]; m=mods[0]
    for k in range(1,len(vecs)):
        a=crt_pair(a,m,vecs[k][1][i],mods[k]); m*=mods[k]
    V[i]=a
print('CRT done t=%.1f'%(time.time()-t0),flush=True)
MQ=ModM(E,Q)
r=MQ.run(list(V)); bad=MQ.score(r)
print('mod Q: nonzero residual atoms=%d  failing equations=%d  (of 39033)'%(sum(1 for x in r if x),len(bad)),flush=True)
json.dump({'Q_bits':Q.bit_length(),'nonzero_atoms_modQ':sum(1 for x in r if x),'failing_eqs_modQ':len(bad)},
          open(os.path.join(HERE,'crt_report.json'),'w'))
# balanced representative over Z
B=[a-Q if a>Q//2 else a for a in V]
json.dump({'x_%d'%i:B[i] for i in range(NV) if B[i]},open(os.path.join(HERE,'crt_balanced.json'),'w'))
print('wrote crt_balanced.json (check with checker.py)',flush=True)
