"""Remaining hard facts: the OR constraint, minimum weight, and what each bit loads."""
import sys, os, json, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
real=[r[1] for r in json.load(open(os.path.join(HERE,'data','gmp16.json')))]
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
D=json.load(open(os.path.join(HERE,'data','bits1.json')))
PIN={int(k):v for k,v in D['pins'].items()}
off=list(base)
for b in real: off[b]=0
forwardp(off)
print("a23000 (the OR gate) with ALL bits off:", evalp(L.polys[23000],off), " -- nonzero means violated")
print("  a23000 =", L.polys[23000], " vars:", {u:off[u]%P for u in L.avars(23000)} if False else "")
for b in (2081,24601):
    v=list(off); v[b]=1; forwardp(v)
    print(f"  only x{b} on:  U={v[7715]} V={v[34554]}  a23000={evalp(L.polys[23000],v)}  "
          f"channel x15298={v[15298]} x5647={v[5647]} x34606={v[34606]}")
# what does each bit load, and is it private?
loadvars=collections.Counter()
for b in real:
    for a,x,C in PIN[b]: loadvars[x]+=1
print(f"\ndistinct loaded variables across all 512 pins: {len(loadvars)}")
print(f"   loaded by exactly one bit: {sum(1 for x,c in loadvars.items() if c==1)}")
print(f"   shared: {sum(1 for x,c in loadvars.items() if c>1)}")
# the small multiplier m in each pin
mult=collections.Counter()
for b in real:
    for a,x,C in PIN[b]:
        Pp=L.polys[a]
        ms=[abs(c) for m,c in Pp.items() if len(m)==1 and m[0]!=b]
        if ms: mult[ms[0]]+=1
print(f"\npin multipliers m (count): {len(mult)} distinct; sample {sorted(mult.items())[:6]}")
print(f"   m == 1: {mult.get(1,0)} of 512")
