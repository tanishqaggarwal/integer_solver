import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
bd={a:evalp(L.polys[a],base) for a in CHK}
print("base failing:",sorted(a for a,x in bd.items() if x))
for u in [14853,24548,9118,8731]:
    v=list(base); v[u]=(v[u]+1)%P; forwardp(v)
    d={a:(evalp(L.polys[a],v)-bd[a])%P for a in CHK}
    d={a:x for a,x in d.items() if x}
    print(f"\nx{u}: response touches {len(d)} checks: {sorted(d)}")
    for a in sorted(d):
        fv=[t for t in sorted(L.avars[a]) if t in FREE]
        print(f"    a{a} (in {len(L.atom2eq.get(a,{}))} eqs) currently {'FAIL' if bd[a] else 'ok'}"
              f"  free vars in atom: {fv[:8]}")
