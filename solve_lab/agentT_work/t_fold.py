#!/usr/bin/env python3
"""AUDIT T10 -- agent Q's premise NOBODY TESTED.

Q's chain: (1) stage law has a cubic invariant; (2) the law is an abelian group op; (3) THEREFORE
'the fold of a leaf subset is its GROUP SUM'; (4) leaves are the ladder 2^i G; (5) subset sums hit
every k, N < 2^256, so T is hit -> A SATISFYING ASSIGNMENT EXISTS.

Steps 1,2,4 are measured in agentQ_work.  Step 3 is an INFERENCE from the law's algebra -- it
assumes every stage of the real circuit implements that law, including the 24 leaf-adjacent
stages agent F explicitly SKIPPED ('fewer than 6 free inputs') and Q flagged as its one caveat.
Nothing in agentQ_work runs the actual circuit on a chosen ON-set and compares the result to the
predicted group sum.  Q's own §5b in fact reports the group sum appearing on ZERO wires for the
deliverable's ON-set {72,235}, and explains it away as the deliverable being degenerate.

This runs the test on FRESH ON-sets, through agent E's independent forward engine.
Read-only w.r.t. other agents' dirs."""
import sys,os,json,collections,random,time
Q='/home/user/integer_solver/solve_lab/agentQ_work'
S='/home/user/integer_solver/solve_lab/agentS_work'
cur=json.load(open(os.path.join(Q,'curve.json')))
p=int(cur['p']); a=int(cur['a']); b=int(cur['b']); c=int(cur['c_shift'])
lad=json.load(open(os.path.join(Q,'ladder.json')))['ladder']
qleaf=json.load(open(os.path.join(Q,'qleaf.json')))
def oncurve(X,Y): return (Y*Y-(X*X*X+a*X+b))%p==0
def add(P1,P2):
    if P1 is None: return P2
    if P2 is None: return P1
    x1,y1=P1; x2,y2=P2
    if x1==x2 and (y1+y2)%p==0: return None
    if P1==P2: l=(3*x1*x1+a)*pow(2*y1,p-2,p)%p
    else:      l=(y2-y1)*pow((x2-x1)%p,p-2,p)%p
    x3=(l*l-x1-x2)%p; return (x3,(l*(x1-x3)-y1)%p)
# leaf points in the shifted (X,Y) chart
LP={}
for e,var in lad.items():
    xs,ys=qleaf[str(var)][0],qleaf[str(var)][1]
    LP[int(e)]=(int(xs)%p,int(ys)%p)   # qleaf already stores the SHIFTED (X,Y) chart
bad=[e for e,P1 in LP.items() if not oncurve(*P1)]
print('ladder leaves on the cubic: %d/%d  (off: %s)'%(len(LP)-len(bad),len(LP),bad[:5]),flush=True)
os.chdir(S); sys.path.insert(0,S)
import common as C
import harness as H, engine as E
SEL={int(v) for v in lad.values()}                     # 253 decoded ladder selectors
ALLSEL={int(k) for k in qleaf}                         # all 256 leaf selectors
base=dict(C.BASE)
exps=sorted(LP)
random.seed(11)
TESTS=[[exps[0]],[exps[1]],[exps[0],exps[1]],[exps[5],exps[90]],
       [exps[3],exps[70],exps[150]],random.sample(exps,5),random.sample(exps,8)]
print('\nON-set -> does the predicted GROUP SUM appear on any wire?\n',flush=True)
hdr='%-26s %8s %8s %10s %10s'%('ON-set (ladder exps)','#wires','#wires','#wires','verdict')
print('%-26s %-9s %-9s %-9s %s'%('ON-set (ladder exps)','sum','indiv','chord2','verdict'),flush=True)
for T in TESTS:
    seed=dict(base)
    for s in ALLSEL: seed[s]=0
    for e in T: seed[int(lad[str(e)])]=1
    t0=time.time()
    try: v=E.forward(seed)
    except Exception as ex:
        print('%-26s forward FAILED: %s'%(str(T)[:26],ex),flush=True); continue
    G=None
    for e in T: G=add(G,LP[e])
    want=set(); 
    if G is not None: want.add((G[0]-c)%p)
    indiv={(LP[e][0]-c)%p for e in T}
    vs=collections.Counter(x%p for x in v if x)
    nsum=sum(vs[w] for w in want)
    nind=sum(vs[w] for w in indiv)
    verdict='GROUP SUM PRESENT' if nsum else ('only individual leaves' if nind else 'neither')
    print('%-26s %-9d %-9d %-9s %s   (%.0fs)'%(str(T)[:26],nsum,nind,'-',verdict,time.time()-t0),flush=True)

# ---- CONTROL: reproduce Q's own §5b numbers on the deliverable, to validate the wire scan.
print('\n== CONTROL: agent Q §5b says at the deliverable, L_72 x-coord is on 92 wires,',flush=True)
print('   L_235 on 5, the group sum L_72+L_235 on 0, target C1 on 4. ==',flush=True)
asg=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
seed={int(k[2:]):int(v) for k,v in asg.items()}
v=E.forward(seed)
vs=collections.Counter(x%p for x in v if x)
def rawx(P1): return (P1[0]-c)%p
s72,s235=LP[72],LP[235]
print('   L_72   on %d wires'%vs[rawx(s72)],flush=True)
print('   L_235  on %d wires'%vs[rawx(s235)],flush=True)
print('   L_72+L_235 on %d wires'%vs[rawx(add(s72,s235))],flush=True)
C1=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002%p
print('   target C1 on %d wires'%vs[C1],flush=True)
