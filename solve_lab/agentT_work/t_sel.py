#!/usr/bin/env python3
"""AUDIT T12 -- the premise agent Q did NOT flag.

Q's step 5 needs: for ANY k in [1,2^256-1], setting the 256 leaf selectors to the binary
expansion of k makes the circuit fold to kG.  That presupposes the 256 selectors INDEPENDENTLY
control which leaves are live.  Q verified the leaf VALUES (the ladder) and the LAW's algebra,
but never that the selector vector actually drives the fold.

Test: at several bases, flip each of the 256 selectors and measure (i) how many wires change,
(ii) whether that leaf's own x-coordinate appears on any wire.  A selector that changes nothing
cannot contribute a bit to k.  Multiple bases -- a single base is exactly the trap agent S fell
into.  Read-only w.r.t. other agents' dirs."""
import sys,os,json,collections,time
Q='/home/user/integer_solver/solve_lab/agentQ_work'; S='/home/user/integer_solver/solve_lab/agentS_work'
cur=json.load(open(os.path.join(Q,'curve.json')))
p=int(cur['p']); c=int(cur['c_shift'])
lad=json.load(open(os.path.join(Q,'ladder.json')))['ladder']
qleaf=json.load(open(os.path.join(Q,'qleaf.json')))
SELof={int(e):int(v) for e,v in lad.items()}
RAWX={int(e):(int(qleaf[str(v)][0])-c)%p for e,v in lad.items()}
os.chdir(S); sys.path.insert(0,S); sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E
import common as C
ALLSEL={int(k) for k in qleaf}
asg=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
D={int(k[2:]):int(v) for k,v in asg.items()}
allon={s:1 for s in ALLSEL}
BASES=[('deliverable',dict(D)),
       ('deliverable, all selectors OFF',{**D,**{s:0 for s in ALLSEL}}),
       ('triple8_seed (agent S base)',dict(C.BASE)),
       ('all-zero',{})]
for bname,bseed in BASES:
    try: v0=E.forward(bseed)
    except Exception as e:
        print('%-32s forward FAILED %s'%(bname,e),flush=True); continue
    nz0=sum(1 for x in v0 if x)
    vs0=collections.Counter(x%p for x in v0 if x)
    on0=sorted(e for e,s in SELof.items() if bseed.get(s,0))
    live0=[e for e in SELof if vs0[RAWX[e]]]
    t0=time.time(); inert=0; moved=0; became=0; ch=[]
    for e in sorted(SELof):
        s=SELof[e]
        b=dict(bseed); b[s]=0 if bseed.get(s,0) else 1
        try: v=E.forward(b)
        except Exception: continue
        d=sum(1 for i in range(len(v0)) if v[i]!=v0[i])
        ch.append(d)
        if d==0: inert+=1
        else: moved+=1
        vsn=collections.Counter(x%p for x in v if x)
        if vsn[RAWX[e]] and not vs0[RAWX[e]]: became+=1
    ch.sort()
    print('%-32s selectors ON=%-2d  live leaves=%-3d | flipping each of 256: inert %3d, moves wires %3d, makes its OWN leaf appear %3d | wires changed med=%d max=%d  (%.0fs)'%(
        bname,len(on0),len(live0),inert,moved,became,ch[len(ch)//2] if ch else -1,ch[-1] if ch else -1,time.time()-t0),flush=True)
