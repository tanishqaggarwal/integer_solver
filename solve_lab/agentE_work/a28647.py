"""Enumerate EVERY knob reaching a28647 (affine and boolean) at several selector
   configurations, and ask whether a20215 is simultaneously reachable with content coprime to p."""
import sys, json, math, re, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, harness as H
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
TGT=[28647,20215]
base={int(k):int(v) for k,v in json.load(open('triple8_seed.json')).items()}
def isb(f):
    for i in H.occ[f]:
        t=re.sub(r'x_%d\b'%f,'X',H.atoms[i])
        if t in ('X - X * X','X * X - X','X * (X - 1)','2 * X * (1 - X)'): return True
    return False
CAND=sorted(set(E.cone(28647)[1])|set(E.cone(20215)[1]))
print(f"union cone of a28647,a20215: {len(CAND)} free vars ({sum(1 for f in CAND if not isb(f))} non-boolean)",flush=True)
def sweep(seed, tag):
    v0=E.forward(seed); bad0=E.badatoms(v0)
    r={a:bad0.get(a,0) for a in TGT}
    rows=collections.defaultdict(list)
    for f in CAND:
        o=v0[f]
        if isb(f):
            probes=[(1 if o!=1 else 0,'bool')]
        else:
            probes=[(o+1,'aff')]
        for val,kind in probes:
            b1,_=fast.resid_delta(v0,bad0,{f:val})
            for a in TGT:
                d=b1.get(a,0)-bad0.get(a,0)
                if d and math.gcd(abs(d),P)==1: rows[a].append((f,kind,d))
    print(f"--- {tag}: bad={sorted(bad0)}  R(28647)%p={r[28647]%P if r[28647] else 0}  R(20215)%p={r[20215]%P if r[20215] else 0}",flush=True)
    for a in TGT:
        ks=rows[a]
        print(f"    knobs on a{a} with delta COPRIME to p: {len(ks)}  "
              f"(non-bool {sum(1 for x in ks if x[1]=='aff')})  e.g. {[x[0] for x in ks[:8]]}",flush=True)
    both=set(f for f,k,d in rows[28647]) & set(f for f,k,d in rows[20215])
    print(f"    knobs coprime-to-p on BOTH rows: {len(both)} {sorted(both)[:10]}",flush=True)
    return rows,bad0
CFG=[({}, "cfg0: selectors 1530,1603 (baseline)"),
     ({490:1}, "cfg1: +x_490 (b-tree)"),
     ({2081:1}, "cfg2: +x_2081 (b-tree)"),
     ({47:1},   "cfg3: +x_47 (a-tree)"),
     ({22106:1},"cfg4: +x_22106 (a-tree)"),
     ({1530:0}, "cfg5: drop x_1530 (single selector 1603)"),
     ({1603:0}, "cfg6: drop x_1603 (single selector 1530)"),
     ({1530:0,1603:0}, "cfg7: no selectors")]
for extra,tag in CFG:
    s=dict(base); s.update(extra)
    try: sweep(s,tag)
    except Exception as e: print(f"--- {tag}: ERR {type(e).__name__}",flush=True)
