"""Re-score the 256 real bits, repairing the handles a flip switches on.

Turning a message bit on activates its load pins.  A pin is  bit*(x - C) = m*handle, and both
x and the handle are free inputs -- so the check it creates is meant to be SATISFIED by setting
them, not counted as damage.  The previous scan never did that, which made every flip look bad.

'Private' knobs are free inputs whose base response touches exactly one check (485 of the 1,726);
those are precisely the handles.  After each flip, repair with them to a fixpoint, then score.
"""
import sys, os, json, time, pickle, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
D=pickle.load(open(os.path.join(HERE,'data','resp_modp.pkl'),'rb'))
priv=collections.defaultdict(list)
for u,d in D['cols'].items():
    if len(d)==1: priv[next(iter(d))].append(u)
print("checks with a private knob:",len(priv))
# also allow free inputs syntactically inside the atom
def cands(a):
    s=list(priv.get(a,()))
    s+= [u for u in sorted(L.avars[a]) if u in FREE and u not in s]
    return s
def fails(v): return [a for a in CHK if evalp(L.polys[a],v)]
def repair(v, rounds=8):
    for _ in range(rounds):
        F=fails(v)
        if not F: return v,F
        did=False
        for a in F:
            r0=evalp(L.polys[a],v)
            for u in cands(a):
                w=list(v); w[u]=(w[u]+1)%P; forwardp(w)
                c=(evalp(L.polys[a],w)-r0)%P
                if not c: continue
                # does it disturb any currently-satisfied check?
                ok=True
                for b2 in CHK:
                    if b2==a: continue
                    if evalp(L.polys[b2],w)!=evalp(L.polys[b2],v): ok=False; break
                if not ok: continue
                v[u]=(v[u]-r0*pow(c,-1,P))%P; forwardp(v); did=True; break
        if not did: break
    return v, fails(v)
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
forwardp(base)
print("base failing:",fails(base))
real=[r[1] for r in json.load(open(os.path.join(HERE,'data','gmp16.json')))]
t0=time.time(); res=[]
for i,u in enumerate(real):
    v=list(base); v[u]=(1-base[u])%P; forwardp(v)
    v,F=repair(v)
    res.append((len(F),u,F[:8]))
    if i%25==0: print(f"   {i}/{len(real)} best={min(r[0] for r in res)} ({time.time()-t0:.0f}s)", flush=True)
res.sort()
print(f"done ({time.time()-t0:.0f}s).  distribution:", dict(sorted(collections.Counter(r[0] for r in res).items())))
for n,u,F in res[:15]: print(f"   x{u} -> failing {n}  {F}")
json.dump([[r[0],r[1]] for r in res], open(os.path.join(HERE,'data','gmp17.json'),'w'))
