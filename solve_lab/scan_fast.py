#!/usr/bin/env python3
"""Fast 2^12 scan: incremental forward (only recompute gates downstream of varied inputs)
+ evaluate only the equations those inputs can affect. Find any combo below agentA's 16."""
import heal_harness as H, json, itertools, sys
p=H.p
pins=json.load(open('pinrec.json'))
sels=[4287,5910,11368,12054,13195,16586,17406,18022,22562,23751,24365,28005]
pinmap={}
for pr in pins:
    ai,sel,tgt,cst,coef,hnd=pr[:6]
    if sel in sels: pinmap.setdefault(sel,{}).setdefault(tgt,cst%p)
pinmap={s:list(d.items()) for s,d in pinmap.items()}
tgts=set(t for s in sels for t,_ in pinmap.get(s,[]))
varied=set(sels)|tgts|{7068,4432,17325,9413,2081}

# affected gates: those with a varied free-ancestor; recompute only these in topo order
aff_gate_idx=[k for k,t in enumerate(H.order) if (H.anc.get(t,set()) & varied)]
aff_targets=set(H.order[k] for k in aff_gate_idx)
# affected equations: those touching a varied input or affected gate (via ancestry ∩ varied)
aff_eqs=[i for i in range(len(H.eqcode)) if any((H.anc.get(v,{v})&varied) for v in H.eqvars[i])]
print(f'affected gates {len(aff_gate_idx)}/{len(H.order)}, affected eqs {len(aff_eqs)}/{len(H.eqcode)}',flush=True)

vA=H.loadd('best_agentA_39022.json')
base=[0]*H.NVARS
for v in H.freeinp: base[v]=vA.get(v,0)
val=H.val
def inc_forward():
    ns={'v':val,'__builtins__':{}}
    for k in aff_gate_idx: val[H.order[k]]=eval(H.gcode[k],ns)
def aff_fail_count():
    ns={'v':val,'__builtins__':{}}
    return sum(1 for i in aff_eqs if eval(H.eqcode[i],ns)!=0)

# baseline sanity: agentA all-zero combo
for v in range(H.NVARS): val[v]=base[v]
val[2081]=1
H.forward()  # full once
val[17325]=0; val[9413]=0
H.forward()
val[7068]=val[2099]; val[4432]=val[19964]
H.forward()
print('agentA-combo affected-fail count',aff_fail_count(),' (should ~16)',flush=True)

best=99; bestcombo=None; hits=[]
for cnt,combo in enumerate(itertools.product([0,1],repeat=len(sels))):
    # reset varied inputs to base (only those we touch)
    for v in varied: val[v]=base[v]
    val[2081]=1; val[17325]=0; val[9413]=0
    for s,bit in zip(sels,combo):
        val[s]=bit
        if bit:
            for t,c in pinmap.get(s,[]): val[t]=c
    inc_forward()
    val[7068]=val[2099]; val[4432]=val[19964]
    inc_forward()
    nf=aff_fail_count()
    if nf<best: best=nf; bestcombo=combo
    if nf<=15: hits.append((nf,combo))
    if (cnt+1)%1024==0: print(f'  {cnt+1}/4096 best={best}',flush=True)
print(f'DONE best={best} combo={dict(zip(sels,bestcombo))}',flush=True)
print(f'combos <=15: {len(hits)}',flush=True)
for nf,combo in sorted(hits)[:25]:
    print(f'  fails={nf}: ON={[s for s,b in zip(sels,combo) if b]}',flush=True)
# save best combo state for exact follow-up
if bestcombo is not None:
    for v in varied: val[v]=base[v]
    val[2081]=1; val[17325]=0; val[9413]=0
    for s,bit in zip(sels,bestcombo):
        val[s]=bit
        if bit:
            for t,c in pinmap.get(s,[]): val[t]=c
    H.forward(); val[7068]=val[2099]; val[4432]=val[19964]; H.forward()
    json.dump({str(v):val[v] for v in range(H.NVARS)}, open('scan_best.json','w'))
    print('saved scan_best.json, full fails:',len(H.fails()),flush=True)
