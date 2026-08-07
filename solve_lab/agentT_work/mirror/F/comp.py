#!/usr/bin/env python3
import sys,os,pickle,collections,json,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
from jac import build_jac
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663

def analyze(v0,label,save=None):
    t0=time.time()
    F0,cols,nonaff=build_jac(v0)
    rows=collections.defaultdict(dict)
    for f,c in cols.items():
        for i,val in c.items(): rows[i][f]=val
    colrows=collections.defaultdict(set)
    for i,r in rows.items():
        for f in r: colrows[f].add(i)
    nz=[i for i in range(len(F0)) if F0[i]]
    print(label,'nonzero atoms',len(nz),'affine knobs',len(cols),'nonaff',len(nonaff),'t',time.time()-t0)
    seen=set(); out=[]
    for s in nz:
        if s in seen: continue
        stack=[s]; comp=set()
        while stack:
            i=stack.pop()
            if i in seen: continue
            seen.add(i); comp.add(i)
            for f in rows.get(i,{}):
                for j in colrows[f]:
                    if j not in seen: stack.append(j)
        out.append(sorted(comp))
    for comp in out:
        ks=sorted(set(f for i in comp for f in rows[i]))
        print('  COMPONENT rows=%d knobs=%d'%(len(comp),len(ks)))
        for i in comp:
            print('    %5d %-95s rhs=%s'%(i,E.res[i][:95],str(-F0[i])[:34]))
            print('          row:',{f:(str(vv//p)+'*p' if vv%p==0 else str(vv)) for f,vv in rows[i].items()})
    if save: pickle.dump({'v0':v0,'F0':F0,'cols':cols,'nonaff':nonaff},open(os.path.join(HERE,save),'wb'))
    return F0,cols,nonaff,rows

if __name__=='__main__':
    import json
    d=json.load(open(os.path.join(HERE,'theirs_closed.json')))
    v=[0]*NV
    for k,x in d.items(): v[int(k[2:])]=int(x)
    v0=[0]*NV
    for f in E.free: v0[f]=v[f]
    analyze(v0,'THEIRS','jacT.pkl')
