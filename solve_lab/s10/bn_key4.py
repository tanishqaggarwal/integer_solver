"""bn_key4: exhaustive sweep of the only boolean vars that can reach the residual.

Chunked + checkpointed to JSONL (absolute paths).
usage: bn_key4.py START END   (index into the candidate list)
"""
import os, sys, json, itertools, time
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

OUT=os.path.join(HERE,'bn_key4_%s.jsonl'%(sys.argv[1] if len(sys.argv)>1 else '0'))
inf=json.load(open(os.path.join(HERE,'bn_infl.json')))
KEY=sorted(set(inf['freebool_anc']) | set(inf['bool_anc2']))
print('key boolean vars:',KEY, flush=True)

v0=L.load(B.BEST); BG=B.broken_gates(v0)
base=B.score(v0)[0]
print('baseline',base, flush=True)

# candidate assignments: each key var over a range of integers
VALS=list(range(-6,9))
cands=[]
for u in KEY:
    for x in VALS:
        if x!=v0[u]: cands.append(((u,),(x,)))
for u1,u2 in itertools.combinations(KEY,2):
    for x1 in VALS:
        for x2 in VALS:
            if x1==v0[u1] and x2==v0[u2]: continue
            cands.append(((u1,u2),(x1,x2)))
if len(KEY)>=3:
    for tri in itertools.combinations(KEY,3):
        for xs in itertools.product(range(-3,6),repeat=3):
            cands.append((tri,xs))
if len(KEY)>=4:
    for xs in itertools.product(range(-2,5),repeat=4):
        cands.append((tuple(KEY),xs))
print('candidates',len(cands), flush=True)

done=set()
if os.path.exists(OUT):
    for ln in open(OUT):
        try: r=json.loads(ln); done.add(r['key'])
        except Exception: pass
S=int(sys.argv[1]); E=min(int(sys.argv[2]),len(cands))
t0=time.time(); best=base
f=open(OUT,'a')
for i in range(S,E):
    us,xs=cands[i]
    key='|'.join(f'{u}={x}' for u,x in zip(us,xs))
    if key in done: continue
    w=list(v0)
    for u,x in zip(us,xs): w[u]=x
    B.fwdb(w,BG,1)
    s,fl,av=B.score(w)
    nz=sum(1 for z in av if z)
    f.write(json.dumps({'key':key,'score':s,'nz':nz,'fail':fl if s>=base else None})+'\n'); f.flush()
    if s>best: best=s; print(f'  NEW BEST {s} at {key}',flush=True)
f.close()
print(f'range [{S},{E}) done in {time.time()-t0:.0f}s best={best}',flush=True)
