import sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
def cone_free(v, maxn=200000):
    seen=set(); out=set(); st=[v]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        d=L.definer.get(u)
        if d is None: out.add(u); continue
        for w in L.avars[d]:
            if w!=u: st.append(w)
    return out, seen
for tgt in [8599,21839,7304,25956]:
    fr,seen=cone_free(tgt)
    print(f"x{tgt}: cone={len(seen)} vars, free inputs={len(fr)}: {sorted(fr)}")
