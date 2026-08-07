"""Tree-native oracle: for each leaf, which of the 96 stages' wires change when it is switched on.
   Prediction from the circuit model: that set is exactly the leaf's ancestor path."""
import sys,json,collections,time,pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentM_work')
import mcore as M, xcompare as X
import engine as E, fast

T,NODES=M.tree()
STAGE={k:set(v['six'])|{int(k)} for k,v in T.items()}
ALLW=sorted(set().union(*STAGE.values()))

def anc(seed):
    v0=E.forward(seed)
    out={}
    for f in M.bools():
        if v0[f]!=0: out[f]='ON'; continue
        v1,aff=fast.apply_delta(v0,{f:1})
        ch={w for w in ALLW if v1[w]!=v0[w]}
        out[f]=frozenset(k for k,ws in STAGE.items() if ws&ch)
    return v0,out

if __name__=='__main__':
    s0=M.load_seed(); BASE=dict(s0); BASE[1530]=0; BASE[1603]=0
    t0=time.time()
    v0,A=anc(dict(BASE))
    print('probe time %.0fs'%(time.time()-t0))
    sizes=collections.Counter(len(s) for s in A.values() if s!='ON')
    print('ancestor-set sizes:',dict(sorted(sizes.items())))
    # test: is A(f) exactly {stages X : f in gsup(X)} ?
    ok=bad=0; mism=[]
    for f,s in A.items():
        if s=='ON': continue
        truth=frozenset(k for k in NODES if f in NODES[k])
        if s==truth: ok+=1
        else: bad+=1; mism.append((f,sorted(s-truth),sorted(truth-s)))
    print('leaves whose changed-stage set == its tree96 ancestor set: %d ; mismatched: %d'%(ok,bad))
    for m in mism[:15]: print('   leaf %d  extra=%s  missing=%s'%m)
    pickle.dump({f:(sorted(s) if s!='ON' else 'ON') for f,s in A.items()},open('/home/user/integer_solver/solve_lab/agentM_work/anc_alloff.pkl','wb'))
