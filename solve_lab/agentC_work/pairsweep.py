"""Sweep over leaf-bit pairs (u on s1 side, w on s2 side); force P_u = P_w by overriding
one bit's pinned coordinate free-variables, close, and score."""
import sys, json, re, time, collections, os
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close2 import *
from ort import leaves
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
W='/home/user/integer_solver/solve_lab/agentC_work/'
lp={int(k):v for k,v in json.load(open(W+'leafpts2.json')).items()}   # bit -> [x,y,varx,vary]
side={}
for r,tag in [(8599,'s1'),(21839,'s1'),(25956,'s2'),(7304,'s2')]:
    for x in leaves(r):
        if x in lp: side[x]=tag
S1=[b for b in lp if side[b]=='s1']; S2=[b for b in lp if side[b]=='s2']
def build(u,w,direction=0):
    ctrl={22162:K2, 30213:K1, u:1, w:1}
    if direction==0: src,dst=w,u        # overwrite u's coords with w's point
    else: src,dst=u,w
    xs,ys,vx_s,vy_s=lp[src]; xd,yd,vx_d,vy_d=lp[dst]
    ctrl[int(vx_d)]=int(xs); ctrl[int(vy_d)]=int(ys)
    return ctrl
if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='one':
        u,w=int(sys.argv[2]),int(sys.argv[3])
        for d in (0,1):
            sc,v,nz=closure2(build(u,w,d),rounds=12)
            print('pair',u,w,'dir',d,'score',sc,'nz',len(nz),[ (a,len(L.atom2eq.get(a,{}))) for a in nz])
            if sc>=39020:
                json.dump({f'x_{i}':v[i] for i in range(L.NVARS) if v[i]!=0},open(W+'PS_%d_%d_%d_%d.json'%(sc,u,w,d),'w'))
    else:
        shard=int(sys.argv[2]); nsh=int(sys.argv[3]); lim=int(sys.argv[4]) if len(sys.argv)>4 else 10**9
        res=[]; best=0; t0=time.time(); k=0
        pairs=[(u,w) for u in sorted(S1) for w in sorted(S2)]
        for i,(u,w) in enumerate(pairs):
            if i%nsh!=shard: continue
            k+=1
            if k>lim: break
            for d in (0,1):
                try: sc,v,nz=closure2(build(u,w,d),rounds=10)
                except Exception as e: continue
                res.append((sc,u,w,d,len(nz)))
                if sc>best:
                    best=sc; print('NEW BEST',sc,u,w,d,'t=%.0f'%(time.time()-t0),flush=True)
                    if sc>=39020:
                        json.dump({f'x_{i2}':v[i2] for i2 in range(L.NVARS) if v[i2]!=0},open(W+'PS_%d_%d_%d_%d.json'%(sc,u,w,d),'w'))
            if k%25==0:
                print('  %d pairs, best %d, %.0fs'%(k,best,time.time()-t0),flush=True)
                json.dump(res,open(W+'runs/pairsweep_%d.json'%shard,'w'))
        json.dump(res,open(W+'runs/pairsweep_%d.json'%shard,'w'))
        print('DONE',best)
