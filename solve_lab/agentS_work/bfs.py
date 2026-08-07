"""BFS over selector configurations; state key = (5 cluster rows mod p, plus which atoms are bad).
   Records the full reachable image rather than sampling it."""
import sys, json, collections, pickle, time, random
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
base=dict(C.BASE); v0=E.forward(base); bad0=E.badatoms(v0)
BOOLS=[f for f in C.cluster_cone() if C.isbool(f)]
MOVES=BOOLS+[30163]
ROWS=C.ROWS
def key_of(assign):
    v,_=fast.apply_delta(v0,assign)
    ns={'v':v,'__builtins__':{}}
    return tuple(eval(H.acodes[a],ns)%P for a in ROWS)
start={}
K0=key_of(start)
seen={K0:dict(start)}
frontier=[dict(start)]
t0=time.time(); nev=0; LIMIT=float(sys.argv[1]) if len(sys.argv)>1 else 3000
gen=0
while frontier and time.time()-t0<LIMIT:
    gen+=1
    nf=[]
    for st in frontier:
        for f in MOVES:
            cur=st.get(f, v0[f])
            new=0 if cur else 1
            a=dict(st); a[f]=new
            try: k=key_of(a)
            except Exception: continue
            nev+=1
            if k not in seen:
                seen[k]=dict(a); nf.append(dict(a))
                if k[3]==0 or k[4]==0:
                    print("*** ZERO ROW HIT",k,flush=True)
                    json.dump({str(x):str(y) for x,y in a.items()},open('bfs_hit_%d.json'%len(seen),'w'))
            if time.time()-t0>LIMIT: break
        if time.time()-t0>LIMIT: break
    print(f"gen{gen}: frontier {len(frontier)} -> {len(nf)} new; total distinct keys {len(seen)}; evals {nev}; {time.time()-t0:.0f}s",flush=True)
    frontier=nf[:40]
    if not nf: break
print("DONE distinct mod-p 5-tuples:",len(seen))
vals=[collections.Counter() for _ in ROWS]
for k in seen:
    for i,x in enumerate(k): vals[i][x]+=1
for i,a in enumerate(ROWS):
    print(f"  a{a}: {len(vals[i])} distinct values mod p; contains 0? {0 in vals[i]}")
pickle.dump(seen,open('bfs_image.pkl','wb'))
