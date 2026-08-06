"""Reusable engine for the 16-ripple heal on the core-solved / residue-G1G2 branch."""
import heal_harness as H
from math import isqrt
import time
p=H.p
RIP=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
CORE=[2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892]
G1G2=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]

def setup():
    """core-solved (39022) with x_7068,x_4432 at residues -> exactly the 16 ripple fail."""
    d=H.loadd('best_agentA_39022.json')
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.val[7068]=None; H.val[4432]=None  # placeholder
    # residues are x_2099,x_19964 which are gates; compute by forward with compensators first
    H.val[7068]=d.get(7068,0); H.val[4432]=d.get(4432,0)
    H.forward()
    r7=H.val[2099]; r4=H.val[19964]
    H.val[7068]=r7; H.val[4432]=r4
    H.forward()
    return r7,r4

_ns={'v':H.val,'__builtins__':{}}
def LHS(idxs): return {i:eval(H.eqcode[i],_ns) for i in idxs}

# eq type: 'lin' -> residual=LHS (affine); 'sq' -> residual=isqrt(LHS)
_TYPE={}
def classify():
    """determine lin vs sq for each RIP eq by bumping x_7068 from compensator."""
    d=H.loadd('best_agentA_39022.json')
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.forward()
    base=LHS(RIP)
    H.val[7068]+=1; H.forward(); v1=LHS(RIP)
    H.val[7068]+=1; H.forward(); v2=LHS(RIP)
    H.val[7068]-=2
    for i in RIP:
        a=v1[i]-base[i]; b=v2[i]-base[i]
        if a!=0 and b==2*a: _TYPE[i]='lin'
        elif a!=0 and b==4*a: _TYPE[i]='sq'
        else:
            # depends on x_4432 not x_7068; test with x_4432
            _TYPE[i]='sq'  # default; refine below
    # refine unknowns via x_4432
    H.val[4432]+=1; H.forward(); w1=LHS(RIP)
    H.val[4432]+=1; H.forward(); w2=LHS(RIP)
    H.val[4432]-=2; H.forward()
    for i in RIP:
        if _TYPE[i]=='sq' and (v1[i]-base[i])==0:
            a=w1[i]-base[i]; b=w2[i]-base[i]
            if a!=0 and b==2*a: _TYPE[i]='lin'
            elif a!=0 and b==4*a: _TYPE[i]='sq'
    return dict(_TYPE)

def inner(idxs):
    """signed inner residual: lin -> LHS ; sq -> isqrt(LHS) (>=0)."""
    L=LHS(idxs); out={}
    for i in idxs:
        if _TYPE.get(i)=='lin': out[i]=L[i]
        else:
            v=L[i]; out[i]=isqrt(v) if v>=0 else -isqrt(-v)
    return out

if __name__=='__main__':
    t=time.time(); classify(); print("classify",{i:_TYPE[i] for i in RIP})
    r7,r4=setup()
    t0=time.time()
    for _ in range(50): H.forward()
    print(f"forward x50: {time.time()-t0:.2f}s  -> {(time.time()-t0)/50*1000:.1f}ms each")
    F=H.fails(); print("fails now",len(F))
    ri=inner(RIP)
    print("inner residuals nonzero:",sum(1 for v in ri.values() if v!=0))
    # VALIDATION: build inner-Jacobian over knobs {7068,4432}, solve, expect compensator gap
    d=H.loadd('best_agentA_39022.json'); C2=d[7068]; C1=d[4432]
    g1=C2-r7; g2=C1-r4   # moving residue->compensator should zero all 16
    base=inner(RIP)
    H.val[7068]+=1; H.forward(); j7=inner(RIP)
    H.val[7068]-=1; H.val[4432]+=1; H.forward(); j4=inner(RIP); H.val[4432]-=1; H.forward()
    # check linear model: base + a*g1 + b*g2 == 0 ?
    ok=0
    for i in RIP:
        a=j7[i]-base[i]; b=j4[i]-base[i]
        pred=base[i]+a*g1+b*g2
        if pred==0: ok+=1
    print(f"linear model recovers compensator (inner->0) for {ok}/16 eqs")
