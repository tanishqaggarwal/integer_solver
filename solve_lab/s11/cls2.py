"""Class sweep with the BEST bit from each tree, and all pairs of best bits."""
import sys, os, json, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
w1=json.load(open(os.path.join(HERE,'data','bits5.json')))   # [count, bit, tree]
best=collections.defaultdict(list)
for n,b,tk in sorted(w1): best[tk].append((n,b))
print("best 4 bits per tree (weight-1 score):")
for tk in 'ABCD': print(f"   {tk}: {best[tk][:4]}")
TOP={tk:[b for n,b in best[tk][:3]] for tk in 'ABCD'}
rows=[]
for us in [(), ('A',), ('B',), ('A','B')]:
    for vs in [(), ('C',), ('D',), ('C','D')]:
        if not us and not vs: continue
        combos=[dict(zip(us+vs,c)) for c in itertools.product(*[TOP[t] for t in us+vs])]
        bestn=None
        for c in combos:
            S=set(c.values())
            v=msg(S); F=fails(v)
            if bestn is None or len(F)<bestn[0]: bestn=(len(F), sorted(S), v[38170], v[3896], F[:6])
        rows.append((bestn[0], us, vs, bestn[1], bestn[2], bestn[3], bestn[4]))
rows.sort()
print(f"\n{'U-side':10s} {'V-side':10s} m2 m1  best-of-top3 message           failing")
for n,us,vs,S,m2,m1,F in rows:
    print(f"{str(us):10s} {str(vs):10s} {m2}  {m1}   {str(S):30s} {n:3d}  {F}")
