#!/usr/bin/env python3
"""AUDIT T14 -- the routing layer: do the leaf pins actually pin the leaves?
Shape found for exp16:  sel*(w1 - X_leaf) - 4949965*x_5923   and   sel*(w2 - Y_leaf) - x_33102.
So atom=0 gives w1 = X_leaf + 4949965*x_5923, NOT w1 = X_leaf.  Census this over all 256 leaves
and report whether the slack variables are free."""
import sys,os,json,re,collections
S='/home/user/integer_solver/solve_lab/agentS_work'
os.chdir(S); sys.path.insert(0,S); sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import harness as H, engine as E
q=json.load(open('../agentQ_work/qleaf.json'))
lad=json.load(open('../agentQ_work/ladder.json'))['ladder']
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
pat=re.compile(r'^x_(\d+) \* \(x_(\d+) - (\d+)\)(.*)$')
kinds=collections.Counter(); slackvars=set(); nopin=[]
exact=0; withslack=0
for e,selv in sorted(lad.items(), key=lambda kv:int(kv[0])):
    sel=int(selv); ent=q[str(sel)]
    found=0
    for i in H.occ[sel]:
        m=pat.match(H.atoms[i].strip())
        if not m: continue
        if int(m.group(1))!=sel: continue
        found+=1
        tail=m.group(4).strip()
        if tail=='':
            exact+=1; kinds['sel*(w - CONST)  [exact pin]']+=1
        else:
            withslack+=1
            kinds['sel*(w - CONST) - <slack>']+=1
            for v in re.findall(r'x_(\d+)',tail): slackvars.add(int(v))
    if found==0: nopin.append(e)
print('leaf pin atoms of the form  sel * (wire - CONST) ...')
for k,v in kinds.most_common(): print('   %-34s %d'%(k,v))
print('   leaves with NO pin atom of this shape: %d %s'%(len(nopin),nopin[:6]))
print()
print('distinct slack variables appearing on the right of a leaf pin: %d'%len(slackvars))
fr=[v for v in slackvars if E.definer[v] is None]
print('   of which FREE variables (independently settable): %d'%len(fr))
print('   sample: %s'%sorted(slackvars)[:10])
print()
print('CONSEQUENCE: setting atom=0 gives  wire = CONST + <slack>.  The leaf coordinate is')
print('recovered only when every slack var is 0.  Agent Q\'s ladder is the slack=0 section.')
