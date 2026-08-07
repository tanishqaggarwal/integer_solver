#!/usr/bin/env python3
"""AUDIT T10b -- same test as t_fold.py, but perturbing the KNOWN-LIVE deliverable state
(whose leaf values demonstrably reach wires) instead of building a seed from scratch."""
import sys,os,json,collections,random,time
Q='/home/user/integer_solver/solve_lab/agentQ_work'; S='/home/user/integer_solver/solve_lab/agentS_work'
cur=json.load(open(os.path.join(Q,'curve.json')))
p=int(cur['p']); a=int(cur['a']); b=int(cur['b']); c=int(cur['c_shift'])
lad=json.load(open(os.path.join(Q,'ladder.json')))['ladder']
qleaf=json.load(open(os.path.join(Q,'qleaf.json')))
LP={int(e):(int(qleaf[str(v)][0])%p,int(qleaf[str(v)][1])%p) for e,v in lad.items()}
SELof={int(e):int(v) for e,v in lad.items()}
def add(P1,P2):
    if P1 is None: return P2
    if P2 is None: return P1
    x1,y1=P1; x2,y2=P2
    if x1==x2 and (y1+y2)%p==0: return None
    l=(3*x1*x1+a)*pow(2*y1,p-2,p)%p if P1==P2 else (y2-y1)*pow((x2-x1)%p,p-2,p)%p
    x3=(l*l-x1-x2)%p; return (x3,(l*(x1-x3)-y1)%p)
def rawx(P1): return None if P1 is None else (P1[0]-c)%p
os.chdir(S); sys.path.insert(0,S); sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E
ALLSEL={int(k) for k in qleaf}
asg=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
D={int(k[2:]):int(v) for k,v in asg.items()}
on0=sorted(s for s in ALLSEL if D.get(s,0))
e_on0=sorted(e for e,s in SELof.items() if D.get(s,0))
print('deliverable: selectors ON = %s  -> ladder exponents %s'%(on0,e_on0),flush=True)
def probe(exps,label):
    seed=dict(D)
    for s in ALLSEL: seed[s]=0
    for e in exps: seed[SELof[e]]=1
    try: v=E.forward(seed)
    except Exception as ex: print('%-30s forward FAILED %s'%(label,ex),flush=True); return
    vs=collections.Counter(x%p for x in v if x)
    G=None
    for e in exps: G=add(G,LP[e])
    nsum=vs[rawx(G)] if G is not None else -1
    ind={e:vs[rawx(LP[e])] for e in exps}
    tot=sum(ind.values())
    print('%-30s indiv-leaf wires %-22s sum-wires %-4d  %s'%(
        label,str(ind)[:22],nsum,'GROUP SUM PRESENT' if nsum else ('leaves only' if tot else 'NO leaf reaches a wire')),flush=True)
# --- are the 256 selectors even independently assignable free variables? ---
free=[s for s in ALLSEL if E.definer[s] is None]
print('\n256 leaf selectors that are FREE variables in agent E\'s parse: %d of %d'%(len(free),len(ALLSEL)),flush=True)
notfree=[s for s in ALLSEL if E.definer[s] is not None]
if notfree: print('   DEFINED (not freely assignable): %s'%notfree[:20],flush=True)
random.seed(7)
print('\n-- control: reproduce the deliverable\'s own ON-set through this code path --',flush=True)
probe(e_on0,'deliverable ON-set %s'%e_on0)
print('\n-- fresh ON-sets --',flush=True)
for T in [[0],[1],[72],[235],[0,1],[72,235],[5,92],[3,72,152],sorted(random.sample(sorted(LP),5))]:
    probe(T,'ON = %s'%T)
