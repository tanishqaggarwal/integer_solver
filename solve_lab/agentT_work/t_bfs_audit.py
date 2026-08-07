#!/usr/bin/env python3
"""AUDIT T8 -- agent S's central negative claim.

RESUME_S §3: the BFS over selector configurations 'terminated by exhaustion (gen6 produced 0
new), giving exactly 48 distinct mod-p 5-tuples', with a20215 never 0 -> 'a20215 mod p = 0 is
the one thing nothing reached', which contradicts agent Q's existence proof.

bfs.py dedups and expands on `key_of(assign)` = the 5 ROW VALUES MOD P, not on the assignment.
Any newly generated assignment whose 5-tuple was already seen is DISCARDED, never expanded.
That is a valid closure only if the 5-tuple determines the set of reachable 5-tuples, i.e. only
if the quotient is compatible with the moves.  It is not: the state is 258 booleans.
This script tests the quotient directly, then re-runs the BFS deduping on the ASSIGNMENT.
Read-only w.r.t. agentS_work / agentE_work (PYTHONDONTWRITEBYTECODE=1 in the runner)."""
import sys, os, json, collections, time, pickle
S='/home/user/integer_solver/solve_lab/agentS_work'
os.chdir(S); sys.path.insert(0,S)
import common as C
import harness as H, engine as E, fast
P=C.P; ROWS=C.ROWS
base=dict(C.BASE); v0=E.forward(base)
BOOLS=[f for f in C.cluster_cone() if C.isbool(f)]
MOVES=BOOLS+[30163]
print('move set: %d booleans + switch 30163 = %d moves'%(len(BOOLS),len(MOVES)),flush=True)
codes=[H.acodes[a] for a in ROWS]
def key_of(assign):
    v,_=fast.apply_delta(v0,assign)
    ns={'v':v,'__builtins__':{}}
    return tuple(eval(c,ns)%P for c in codes)
def akey(a):  # canonical key of the ASSIGNMENT itself
    return tuple(sorted((f,val) for f,val in a.items() if val!=v0[f]))

MODE=sys.argv[1] if len(sys.argv)>1 else 'quotient'
LIMIT=float(sys.argv[2]) if len(sys.argv)>2 else 1200

if MODE=='quotient':
    # ---- Test 1: does the 5-tuple determine the successor 5-tuples?  Collect several
    # DISTINCT assignments sharing one key, then compare their single-flip successor sets.
    print('\n== TEST 1: is bfs.py\'s dedup key a valid quotient? ==',flush=True)
    bykey=collections.defaultdict(list)
    st={}
    for f in MOVES:
        a=dict(st); a[f]=0 if st.get(f,v0[f]) else 1
        try: bykey[key_of(a)].append(a)
        except Exception: pass
    cands=[(k,v) for k,v in bykey.items() if len(v)>=2]
    cands.sort(key=lambda kv:-len(kv[1]))
    print('depth-1 states: %d assignments -> %d distinct 5-tuples'%(len(MOVES),len(bykey)),flush=True)
    print('   (bfs.py keeps ONE assignment per 5-tuple and discards the rest)',flush=True)
    k,group=cands[0]
    print('largest collision class: %d distinct assignments share one 5-tuple'%len(group),flush=True)
    succ=[]
    for a in group[:6]:
        s=set()
        for f in MOVES:
            b=dict(a); b[f]=0 if a.get(f,v0[f]) else 1
            try: s.add(key_of(b))
            except Exception: pass
        succ.append(s)
        print('   rep %-28s -> %d distinct successor 5-tuples'%(str(akey(a))[:28],len(s)),flush=True)
    same=all(s==succ[0] for s in succ)
    print('ALL successor sets identical? %s'%same,flush=True)
    if not same:
        u=set().union(*succ)
        print('   union %d vs individual %s  -> QUOTIENT IS INVALID:'%(u and len(u),[len(s) for s in succ]),flush=True)
        print('   bfs.py explores ONE branch of a %d-way fork and calls the result a closure.'%len(group),flush=True)
        a20=set(t[3] for s in succ for t in s)
        print('   distinct a20215 values across the union of these successor sets: %d'%len(a20),flush=True)
else:
    # ---- Test 2: re-run the BFS deduping on the ASSIGNMENT, no frontier truncation.
    print('\n== TEST 2: BFS deduped on the ASSIGNMENT (bfs.py used nf[:40] and key dedup) ==',flush=True)
    t0=time.time()
    start={}
    seenA={akey(start)}
    K=key_of(start)
    keys={K}; a20={K[3]}
    frontier=[start]; gen=0; nev=0; hit=None
    while frontier and time.time()-t0<LIMIT:
        gen+=1; nf=[]
        for stt in frontier:
            for f in MOVES:
                b=dict(stt); b[f]=0 if stt.get(f,v0[f]) else 1
                ak=akey(b)
                if ak in seenA: continue
                seenA.add(ak)
                try: k=key_of(b)
                except Exception: continue
                nev+=1
                keys.add(k); a20.add(k[3]); nf.append(b)
                if k[3]==0:
                    hit=b; print('*** a20215 == 0 REACHED:',ak,flush=True); break
                if time.time()-t0>LIMIT: break
            if hit or time.time()-t0>LIMIT: break
        print('gen%d: %d -> %d new states; distinct 5-tuples %d; distinct a20215 %d; evals %d; %.0fs'%(
              gen,len(frontier),len(nf),len(keys),len(a20),nev,time.time()-t0),flush=True)
        if hit: break
        frontier=nf
        if not nf: break
    print('\nRESULT distinct 5-tuples reached: %d   (bfs.py reported 48 as the CLOSURE)'%len(keys),flush=True)
    print('RESULT distinct a20215 values : %d   (RESUME_S §3 reports 2)'%len(a20),flush=True)
    print('a20215 == 0 reached? %s'%(hit is not None),flush=True)
    pickle.dump(sorted(a20),open('/home/user/integer_solver/solve_lab/agentT_work/t_a20215_values.pkl','wb'))
