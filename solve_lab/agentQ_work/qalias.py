#!/usr/bin/env python3
"""Q-12: characterise layer (e) -- what sits between a slot's mux output and its parent's input.

Observed shapes, read verbatim:
    x_17675 - x_20820 - x_36780              ->  parent_in = muxout + (x_4116*x_22163)
    6910381*(x_15439 - x_18440) - x_11630    ->  parent_in = muxout + (x_1962*x_10858)/6910381
    x_24468 - x_13682 - 12354891*x_34243     ->  ROOT PIN  = muxout + 12354891*x_34243
i.e. an AFFINE ALIAS: parent input = mux output + (a multiple of) one further wire.
"""
import pickle,collections,re,json
exec(open('qmux2.py').read().split('ST=[x for x')[0])
ST=[x for x in json.load(open('qstages.json'))['stages'] if 'u3' in x]
byvar=collections.defaultdict(list)
for i,(s,vs) in enumerate(terms):
    for v in vs: byvar[v].append(i)
SLOTIN=set()
for g in ST: SLOTIN|={g['ua'],g['ub'],g['ya'],g['yb']}
ROOT={24468,18956}
A1=re.compile(r'^x_(\d+) - x_(\d+) - x_(\d+)$')
A2=re.compile(r'^(\d+) \* \(x_(\d+) - x_(\d+)\) - x_(\d+)$')
A3=re.compile(r'^x_(\d+) - x_(\d+) - (\d+) \* x_(\d+)$')
A4=re.compile(r'^x_(\d+) - \(x_(\d+) \+ x_(\d+)\)$')
def alias_of(M):
    """return list of (parent_wire, slack_wire, form)"""
    out=[]
    for i in byvar[M]:
        s=terms[i][0]
        m=A1.match(s)
        if m:
            P,X,Q=map(int,m.groups())
            if X==M: out.append((P,Q,'P = M + Q'))
            elif Q==M: out.append((P,X,'P = Q + M'))
            continue
        m=A2.match(s)
        if m:
            k,P,X,Q=map(int,m.groups())
            if X==M: out.append((P,Q,'k*(P-M) = Q'))
            continue
        m=A3.match(s)
        if m:
            P,X,k,Q=map(int,m.groups())
            if X==M: out.append((P,Q,'P = M + k*Q'))
            continue
        m=A4.match(s)
        if m:
            P,u,v=map(int,m.groups())
            if u==M: out.append((P,v,'P = M + Q (sum)'))
            elif v==M: out.append((P,u,'P = Q + M (sum)'))
    return out
# recover both mux outputs per slot
outs=[]
for g in ST:
    for P,Q in (((g['ua'],g['ub'],g['u3']),(g['ya'],g['yb'],g['y3'])),
                ((g['ya'],g['yb'],g['y3']),(g['ua'],g['ub'],g['u3'])),
                ((g['ua'],g['ub'],g['y3']),(g['ya'],g['yb'],g['u3'])),
                ((g['ya'],g['yb'],g['u3']),(g['ua'],g['ub'],g['y3']))):
        got=None
        for cA,cB,cC,X in mux(*P):
            if not quad(cA,cB,cC): continue
            for a2,b2,c2,Y in mux(*Q):
                if (cA,cB,cC)==(a2,b2,c2): got=(X,Y); break
            if got: break
        if got: break
    outs.append(got)
res=collections.Counter(); slack=[]; forms=collections.Counter(); unres=[]
for (X,Y) in outs:
    for M in (X,Y):
        al=alias_of(M)
        tgt=[(P,Q,f) for (P,Q,f) in al if P in SLOTIN or P in ROOT]
        if tgt:
            P,Q,f=tgt[0]
            res['ROOT PIN' if P in ROOT else 'parent slot input']+=1
            forms[f]+=1; slack.append(Q)
        else:
            res['no aliased parent found']+=1; unres.append(M)
print('mux outputs examined: %d (2 per slot)'%(2*len(ST)))
for k,v in res.most_common(): print('   %-32s %d'%(k,v))
print('alias forms:',dict(forms))
print('distinct slack wires:',len(set(slack)))
# what are the slack wires?
PROD2=re.compile(r'^x_(\d+) - x_(\d+) \* x_(\d+)$')
kinds=collections.Counter()
for q in set(slack):
    ss=[terms[i][0] for i in byvar[q] if len(terms[i][0])<45]
    kinds['product of two wires' if any(PROD2.match(x) and int(PROD2.match(x).group(1))==q for x in ss) else 'other']+=1
print('slack wire shapes:',dict(kinds))
if unres: print('unresolved mux outputs:',unres[:10])
