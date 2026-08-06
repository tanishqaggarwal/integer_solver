"""The obstruction as an INVARIANT of the message.

A certificate y satisfies y.J = 0, so the weighted combination  sum_a y_a * r_a  of check
residues is CONSTANT under every continuous move.  A solution needs all r_a = 0, hence needs the
invariant to vanish.  It is computable from one 0.08 s forward evaluation, which turns the search
from 'count failures' into 'hit zero in GF(p)' -- and lets us see whether the invariants have
structure or are effectively random in the message.
"""
import sys, os, json, time, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
print(f"{len(CERT)} certificates, supports {[len(c) for c in CERT]}")
def inv(v):
    return [sum(y*evalp(L.polys[a],v) for a,y in c.items())%P for c in CERT]
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
forwardp(base)
print("invariants at the checkpoint base:", [str(x)[:16]+'..' if x else 'ZERO' for x in inv(base)])
for S in [{2081,24601},{4287,24601},{13195,24601},{12054,24601},{2081},{24601}]:
    v=msg(S)
    I=inv(v)
    print(f"  msg {sorted(S)}: nonzero invariants {sum(1 for x in I if x)}/{len(I)}  "
          f"{[('ZERO' if not x else str(x)[:10]+'..') for x in I]}")
# is the invariant invariant under moves?  perturb a free knob and re-check
v=msg({2081,24601}); I0=inv(v)
import random
rnd=random.Random(1)
FREEV=[u for u in range(L.NVARS) if u not in L.definer]
same=0; tot=0
for u in rnd.sample(FREEV,40):
    w=list(v); w[u]=(w[u]+rnd.randrange(1,10**6))%P; forwardp(w)
    tot+=1
    if inv(w)==I0: same+=1
print(f"invariants unchanged under {same}/{tot} random single-knob perturbations")
