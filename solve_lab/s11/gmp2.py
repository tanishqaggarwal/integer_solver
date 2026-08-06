"""The barrier, in its smallest form.

A global mod-p forward evaluation from the checkpoint satisfies EVERY gate and all but SIX
checks:  a7930, a29539, a35759, a35760, a40826, a41512.  Zero those six residues in GF(p) and
every atom is 0 mod p; each equation is then p*r and the p-quantised handles absorb it exactly.

So: for each free input, perturb it, re-run the mod-p forward evaluation, and record
  (a) how many checks fail mod p afterwards  -- perturbing an input can break the checks that
      constrain it (booleanity, load pins), and
  (b) the response of the six targets.
Inputs that keep the failing count at six are the usable knobs.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
TGT=[7930,29539,35759,35760,40826,41512]
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
FREE=[u for u in range(L.NVARS) if u not in L.definer]
print("free inputs:",len(FREE))
t0=time.time()
v=list(base); forwardp(v)
print(f"one forward eval: {time.time()-t0:.2f}s")
def probe(u, step):
    v=list(base); v[u]=(v[u]+step)%P
    forwardp(v)
    res=[evalp(L.polys[a],v) for a in TGT]
    nf=sum(1 for a in CHK if evalp(L.polys[a],v))
    return res, nf
b0=[evalp(L.polys[a],base) for a in TGT]
nb=sum(1 for a in CHK if evalp(L.polys[a],base))
print("base residues nonzero:",[r!=0 for r in b0]," base failing checks:",nb)
rows=[]
t0=time.time()
for i,u in enumerate(FREE):
    r1,n1=probe(u,1)
    if r1==b0 and n1==nb: continue           # inert
    rows.append((u,n1,[ (r1[j]-b0[j])%P for j in range(6)]))
    if len(rows)%50==0:
        print(f"   {i}/{len(FREE)} live={len(rows)} ({time.time()-t0:.0f}s)", flush=True)
print(f"live free inputs: {len(rows)}  ({time.time()-t0:.0f}s)")
clean=[r for r in rows if r[1]<=nb]
print(f"  of which keep failing checks <= {nb}: {len(clean)}")
print("  examples:", [(r[0],r[1]) for r in clean[:25]])
json.dump([[r[0],r[1]]+[str(x) for x in r[2]] for r in rows],
          open(os.path.join(HERE,'data','gmp2.json'),'w'))
