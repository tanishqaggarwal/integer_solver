"""Cache the full mod-p response matrix: free input -> (check -> delta).

One global mod-p forward evaluation costs 0.08s, so all 7,273 free inputs can be probed exactly.
This is the object every earlier attempt lacked: the true, global, exact linearisation of the
instance in the layer where the obstruction lives.
"""
import sys, os, json, time, pickle, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=[u for u in range(L.NVARS) if u not in L.definer]
src=sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'data','gmp1_state.json')
out=sys.argv[2] if len(sys.argv)>2 else os.path.join(HERE,'data','resp_modp.pkl')
base=[int(x) for x in json.load(open(src))]
forwardp(base)
bd={a:evalp(L.polys[a],base) for a in CHK}
print("base failing checks:",sorted(a for a,x in bd.items() if x), flush=True)
cols={}
t0=time.time()
for i,u in enumerate(FREE):
    v=list(base); v[u]=(v[u]+1)%P; forwardp(v)
    d={}
    for a in CHK:
        x=evalp(L.polys[a],v)
        if x!=bd[a]: d[a]=(x-bd[a])%P
    if d: cols[u]=d
    if i%400==0:
        print(f"  {i}/{len(FREE)} live={len(cols)} ({time.time()-t0:.0f}s)", flush=True)
print(f"live knobs: {len(cols)}  ({time.time()-t0:.0f}s)")
sz=collections.Counter()
for u,d in cols.items(): sz[min(len(d),50)]+=1
print("  response widths (capped at 50):", dict(sorted(sz.items())[:20]))
pickle.dump({'base':base,'bd':bd,'cols':cols}, open(out,'wb'))
print("saved", out)
