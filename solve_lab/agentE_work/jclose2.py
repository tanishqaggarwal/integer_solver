"""Fast closure of the exact linear map free-vars -> atom residuals (incremental probing)."""
import sys, json, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, fast
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
def mux11(s,n=5):
    v=E.forward(s)
    for _ in range(n):
        s[22162]=v[13682]; s[30213]=v[18956]-v[32237]; v=E.forward(s)
    return v,s
def cone_free(a): return set(E.cone(a)[1])

def run(abit,bbit,maxvars=20000,log=sys.stdout):
    s={18956:C,abit:1,bbit:1}
    v0,s=mux11(s)
    bad0=E.badatoms(v0)
    print("base bad",sorted(bad0),file=log,flush=True)
    FROZEN={abit,bbit}
    S=[]; cols={}; nonlin=set(); processed=set(); rounds={}
    pending=set(bad0)
    t0=time.time()
    for rnd in range(20):
        newS=set()
        for a in pending: newS|=cone_free(a)
        newS-=set(S)|FROZEN
        newS=sorted(newS)
        if not newS: break
        print(f"round {rnd}: +{len(newS)} vars (total {len(S)+len(newS)})",file=log,flush=True)
        for f in newS:
            b1,_=fast.resid_delta(v0,bad0,{f:v0[f]+1})
            b2,_=fast.resid_delta(v0,bad0,{f:v0[f]+2})
            col={}
            for a in set(b1)|set(bad0):
                d=b1.get(a,0)-bad0.get(a,0)
                if d: col[a]=d
            for a in set(b2)|set(bad0)|set(col):
                if b2.get(a,0)-bad0.get(a,0)!=2*col.get(a,0): nonlin.add((f,a))
            cols[f]=col; S.append(f); rounds.setdefault(rnd,[]).append(f)
        aff=set()
        for f in newS: aff|=set(cols[f])
        processed|=pending
        pending=(aff|set(bad0))-processed
        print(f"  affected atoms so far {len(processed|pending)}, nonlin {len(nonlin)}, {time.time()-t0:.0f}s",file=log,flush=True)
        if len(S)>maxvars: print("  var cap hit",file=log,flush=True); break
    allat=set(bad0)
    for f in cols: allat|=set(cols[f])
    print(f"FINAL: vars {len(S)} atoms {len(allat)} nonlin {len(nonlin)} {time.time()-t0:.0f}s",file=log,flush=True)
    pickle.dump({'r0':bad0,'cols':cols,'nonlin':nonlin,'S':S,'base':s,'atoms':sorted(allat),'rounds':rounds},
                open(f'jac_{abit}_{bbit}.pkl','wb'))
if __name__=='__main__':
    run(int(sys.argv[1]),int(sys.argv[2]))
