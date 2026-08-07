#!/usr/bin/env python3
"""AUDIT T13 -- 'liveness is fully determined by the selectors': measured, not assumed.
qleaf[selvar] = [X, Y, w1, w2]; w1,w2 are the two wire vars agent Q associates with the leaf.
Turn each selector ON alone and ask whether the leaf's coordinates actually arrive on w1/w2."""
import sys,os,json,collections,random
Q='/home/user/integer_solver/solve_lab/agentQ_work'; S='/home/user/integer_solver/solve_lab/agentS_work'
p=int(json.load(open(os.path.join(Q,'curve.json')))['p']); c=int(json.load(open(os.path.join(Q,'curve.json')))['c_shift'])
lad=json.load(open(os.path.join(Q,'ladder.json')))['ladder']
qleaf=json.load(open(os.path.join(Q,'qleaf.json')))
os.chdir(S); sys.path.insert(0,S); sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E
ALLSEL={int(k) for k in qleaf}
asg=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
D={int(k[2:]):int(v) for k,v in asg.items()}
ent={int(e):qleaf[str(v)] for e,v in lad.items()}
SELof={int(e):int(v) for e,v in lad.items()}
print('are the leaf wire vars w1,w2 FREE or DEFINED in E\'s parse?')
w1f=sum(1 for e in ent if E.definer[int(ent[e][2])] is None)
w2f=sum(1 for e in ent if E.definer[int(ent[e][3])] is None)
print('   w1 free in %d/253, w2 free in %d/253'%(w1f,w2f))
random.seed(3)
SAMPLE=sorted(random.sample(sorted(ent),12))
for bname,bseed in [('all selectors OFF',{**D,**{s:0 for s in ALLSEL}}),('all-zero',{})]:
    print('\nbase = %s'%bname)
    hit=0
    for e in SAMPLE:
        b=dict(bseed); b[SELof[e]]=1
        try: v=E.forward(b)
        except Exception as ex: print('   exp %-4d forward FAILED'%e); continue
        w1,w2=int(ent[e][2]),int(ent[e][3])
        X,Y=int(ent[e][0])%p,int(ent[e][1])%p
        got1,got2=v[w1]%p,v[w2]%p
        want={'X':X,'Y':Y,'X-c':(X-c)%p,'Y-c':(Y-c)%p}
        m1=[k for k,val in want.items() if got1==val]; m2=[k for k,val in want.items() if got2==val]
        anywhere = collections.Counter(x%p for x in v if x)
        appears = sum(anywhere[val] for val in want.values())
        if m1 or m2 or appears: hit+=1
        print('   exp %-4d sel x%-6d w1=x%-6d -> %-6s  w2=x%-6d -> %-6s | leaf coords anywhere in circuit: %d wires'%(
            e,SELof[e],w1,(m1 or ['no'])[0],w2,(m2 or ['no'])[0],appears))
    print('   leaves that arrived anywhere: %d/%d'%(hit,len(SAMPLE)))
