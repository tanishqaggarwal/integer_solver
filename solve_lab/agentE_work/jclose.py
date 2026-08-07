import sys, json, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, jsolve as J
C=J.C
import ast
base=json.load(open(sys.argv[1])) if len(sys.argv)>1 else {"18956":C,"4279":1,"26005":1}
base={int(k):int(v) for k,v in base.items()}
r0,_=J.resid(base)
print("base bad", sorted(r0), flush=True)
S=set(); done=set(); cols={}; nonlin=set()
pending=set(r0)
for rnd in range(12):
    newS=set()
    for a in pending: newS|=J.cone_free(a)
    newS -= set(base)|S|{22162,30213}
    newS=sorted(newS)
    if not newS: break
    print(f"round {rnd}: adding {len(newS)} free vars", flush=True)
    _,c2,nl=J.build(base,newS)
    cols.update(c2); nonlin|=nl; S|=set(newS)
    aff=set()
    for f in newS: aff|=set(cols[f])
    pending=aff-done; done|=aff
    print(f"  affected now {len(done)} atoms; nonlin {len(nonlin)}", flush=True)
    if len(S)>1200: break
allat=set(r0)
for f in cols: allat|=set(cols[f])
print("FINAL: vars",len(S),"atoms",len(allat),"nonlin",len(nonlin))
pickle.dump({'r0':r0,'cols':cols,'nonlin':nonlin,'S':sorted(S),'base':base,'atoms':sorted(allat)},open(sys.argv[2] if len(sys.argv)>2 else 'jacC.pkl','wb'))
