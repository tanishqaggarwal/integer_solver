import os, sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
import suppfree
src = sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
v = L.load(src); ad.fwd(v, rounds=6)
vm=[x%P for x in v]
TARG=[688,1618,19297,19299,30984,36185,40608,40812]
for modp in (True,False,None):
    idx, freelist, vs = suppfree.build(vm, modp=modp)
    tot=set()
    for a in TARG:
        m=0
        for w in L.avars[a]:
            m |= vs[w] if w<len(vs) else 0
        s={freelist[i] for i in range(len(freelist)) if (m>>i)&1}
        # also w itself if free
        s |= {w for w in L.avars[a] if w not in L.definer}
        tot |= s
        print('modp=%s a%-6d supp=%d' % (modp,a,len(s)))
    print('modp=%s UNION over the 8 = %d free inputs' % (modp,len(tot)))
    print(sorted(tot)[:80])
    print()
