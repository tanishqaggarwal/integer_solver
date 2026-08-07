#!/usr/bin/env python3
"""Q-11d: close the induction.  For every slot recover (s1,s2,cA,cB,cC,Xout,Yout); then check
  (i)  every live bit is either a boolean-pinned leaf selector or the OR wire of a lower slot,
  (ii) every slot input is either a leaf coordinate wire or a lower slot's mux output,
  (iii) the whole thing is one tree over the 256 leaves whose top mux output is the ROOT PIN.
"""
import pickle,collections,re,json
exec(open('qmux2.py').read().split('ST=[x for x')[0])
ST=[x for x in json.load(open('qstages.json'))['stages'] if 'u3' in x]
leaf={int(g):v for g,v in json.load(open('qleaf.json')).items()}
LEAFSEL={int(g) for g in leaf}; LEAFW={v[2] for v in leaf.values()}|{v[3] for v in leaf.values()}
ROOTX,ROOTY=24468,18956
# OR wires:  x_t - (x_u - x_v)  with summ[u]==prod[v]=={a,b}
SUBP=re.compile(r'^x_(\d+) - \(x_(\d+) - x_(\d+)\)$')
ORof={}
for s,vs in terms:
    m=SUBP.match(s)
    if not m: continue
    t,u,v=map(int,m.groups())
    if u in summ and v in prod and set(summ[u])==set(prod[v]): ORof[t]=tuple(sorted(summ[u]))
info={}
for g in ST:
    U,V=(g['ua'],g['ub'],g['u3']),(g['ya'],g['yb'],g['y3'])
    W,Z=(g['ua'],g['ub'],g['y3']),(g['ya'],g['yb'],g['u3'])
    for P,Q in ((U,V),(V,U),(W,Z),(Z,W)):
        got=None
        for cA,cB,cC,X in mux(*P):
            q=quad(cA,cB,cC)
            if not q: continue
            for a2,b2,c2,Y in mux(*Q):
                if (cA,cB,cC)==(a2,b2,c2): got=(q,cA,cB,cC,X,Y,P,Q); break
            if got: break
        if got: info[id(g)]=got; break
print('slots characterised: %d / %d'%(len(info),len(ST)))
orpair={tuple(sorted((root(a),root(b)))):t for t,(a,b) in ORof.items()}
livesrc=collections.Counter(); insrc=collections.Counter()
muxout={}
for g in ST:
    (s1,s2),cA,cB,cC,X,Y,P,Q=info[id(g)]
    muxout[X]=g; muxout[Y]=g
for g in ST:
    (s1,s2),cA,cB,cC,X,Y,P,Q=info[id(g)]
    for s in (s1,s2):
        livesrc['boolean-pinned leaf selector' if root(s) in LEAFSEL else
                ('OR wire of a lower slot' if root(s) in ORof or s in ORof else 'other')]+=1
    for w in (P[0],P[1],Q[0],Q[1]):
        insrc['leaf coordinate wire' if w in LEAFW else
              ('lower slot mux output' if w in muxout else
               ('hard zero / dummy' if w in [k for k in (P[1],Q[1])] and w not in muxout and w not in LEAFW else 'other'))]+=1
print('live-bit sources (2 per slot, %d total):'%(2*len(ST)))
for k,v in livesrc.most_common(): print('    %-38s %d'%(k,v))
print('slot input sources (4 per slot, %d total):'%(4*len(ST)))
for k,v in insrc.most_common(): print('    %-38s %d'%(k,v))
print('ROOT: x_%d is a slot mux output: %s ; x_%d: %s'%(ROOTX,ROOTX in muxout,ROOTY,ROOTY in muxout))
tops=[g for g in ST if info[id(g)][4] not in [i for gg in ST for i in (info[id(gg)][6][0],info[id(gg)][6][1],info[id(gg)][7][0],info[id(gg)][7][1])]]
print('slots whose mux output feeds no other slot (tree tops): %d'%len(tops))
for g in tops[:6]:
    (s1,s2),cA,cB,cC,X,Y,P,Q=info[id(g)]
    print('    top slot mux outputs: x_%d , x_%d   %s'%(X,Y,'== ROOT PIN' if {X,Y}=={ROOTX,ROOTY} else ''))
