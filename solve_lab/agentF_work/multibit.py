#!/usr/bin/env python3
"""Multi-bit ladder consistency test: 2+ bits on, exact integer Jacobian, exact component solve."""
import sys,os,json,pickle,collections,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
from jac import build_jac
from intsolve import solve_int
E=Engine(); NR=len(E.res)
p=115792089237316195423570985008687907853269984665640564039457584007908834671663

def analyze(bits, tag):
    v=[0]*NV
    for b in bits: v[b]=1
    t0=time.time()
    F0,cols,nonaff=build_jac(v)
    rows=collections.defaultdict(dict)
    for f,c in cols.items():
        for i,val in c.items(): rows[i][f]=val
    colrows=collections.defaultdict(set)
    for i,r in rows.items():
        for f in r: colrows[f].add(i)
    nz=[i for i in range(NR) if F0[i]]
    print(tag,'bits',bits,'broken atoms',len(nz),'jac t=%.0fs'%(time.time()-t0),flush=True)
    seen=set()
    for s in nz:
        if s in seen or s not in rows: 
            if s not in rows: print('   atom',s,E.res[s][:80],'has NO affine knob',flush=True)
            continue
        stack=[s]; comp=set()
        while stack:
            i=stack.pop()
            if i in seen: continue
            seen.add(i); comp.add(i)
            for f in rows[i]:
                for j in colrows[f]:
                    if j not in seen: stack.append(j)
        comp=sorted(comp)
        knobs=sorted(set(f for i in comp for f in rows[i]))
        A=[[rows[i].get(f,0) for f in knobs] for i in comp]
        b=[-F0[i] for i in comp]
        x=solve_int(A,b)
        print('   component rows=%d knobs=%d  integrally solvable=%s'%(len(comp),len(knobs),x is not None),flush=True)
        if x is None:
            for i in comp: print('      %5d %-90s rhs=%s'%(i,E.res[i][:90],str(-F0[i])[:26]),flush=True)
    return F0,cols

if __name__=='__main__':
    d=json.load(open(os.path.join(HERE,'curve.json')))
    # reconstruct chain order
    pts={int(k):(int(v[0]),int(v[1])) for k,v in d['pts'].items()}
    sup=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))
    A=set(sup['7715'])
    def add(P,Q):
        x1,y1=P;x2,y2=Q
        if x1==x2:
            if (y1+y2)%p==0: return None
            l=3*x1*x1%p*pow(2*y1,p-2,p)%p
        else: l=(y2-y1)*pow(x2-x1,p-2,p)%p
        x3=(l*l-x1-x2)%p; return (x3,(l*(x1-x3)-y1)%p)
    inv={v:k for k,v in pts.items()}
    succ={}
    for bit,P in pts.items():
        D=add(P,P)
        if D in inv: succ[bit]=inv[D]
    chain=[2779]
    while chain[-1] in succ: chain.append(succ[chain[-1]])
    print('chain from 2779 length',len(chain),flush=True)
    aa=[b for b in chain if b in A][:2]
    print('two tree-A chain bits:',aa,flush=True)
    analyze(aa,'TWO-A')
