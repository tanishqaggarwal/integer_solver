import sys, pickle, random, json, os
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from lib2 import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
base=forward([0]*L.NVARS)
base_av=L.all_atom_values(base)
BASESC=L.NEQ-len(L.failing_eqs(base_av))
FREE=[u for u in range(L.NVARS) if u not in outs]
def probe(u,val):
    v=list(base); v[u]=val
    forward(v)
    av=L.all_atom_values(v)
    sc=L.NEQ-len(L.failing_eqs(av))
    return (sc,(v[18956]-K1)%P,(v[24468]-K2)%P,v[7715],v[34554])
if __name__=='__main__':
    shard=int(sys.argv[1]); nsh=int(sys.argv[2]); mode=sys.argv[3]
    random.seed(12345)
    R=random.randrange(1,P)
    out=[]
    for i,u in enumerate(FREE):
        if i%nsh!=shard: continue
        val=1 if mode=='one' else R
        try: r=probe(u,val)
        except Exception as e: r=('ERR',str(e)[:40],0,0,0)
        out.append((u,)+tuple(r))
    json.dump(out,open('/home/user/integer_solver/solve_lab/agentC_work/runs/probe1_%s_%d.json'%(mode,shard),'w'))
    print('done',shard,mode,len(out))
