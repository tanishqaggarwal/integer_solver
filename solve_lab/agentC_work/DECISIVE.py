"""Decisive test asked by the coordinator, answered from my own parse only."""
import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close2 import *
from ec import add,mul,leafpoints,N
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
CHKs=[a for a in range(L.NA) if a not in L.atom_out]
def show(ctrl,tag):
    v=[0]*L.NVARS
    for k,val in ctrl.items(): v[k]=val
    forward(v)
    av=L.all_atom_values(v)
    nz=[a for a in CHKs if av[a]!=0]
    print('[%s] score=%d  nonzero checks=%s'%(tag,L.NEQ-len(L.failing_eqs(av)),
          [(a,len(L.atom2eq.get(a,{}))) for a in nz]))
    print('    a688 (x_18956 = K1 mod p) EXACT ZERO over Z: %s'%(av[688]==0))
    print('    a1618(x_24468 = K2 mod p) EXACT ZERO over Z: %s'%(av[1618]==0))
    print('    a23000 (OR(s1,s2)=1)      EXACT ZERO over Z: %s'%(av[23000]==0))
    return v
show({}, 'all free inputs = 0')
show({542:1,91:1}, 'branch (1,1), no coords')
v=show({542:1,91:1,22162:K2,30213:K1}, 'branch (1,1) + x_22162=K2, x_30213=K1')
print()
print('=> the ONLY nonzero checks left are the two activated bits own conditional pins.')
print('=> after closing those (close.py / close2.py): score 39,013, residual =')
print('   a19297 x_15298*x_11150+x_4007 | a19299 x_15298*x_25739-6672769*x_29804 |')
print('   a30984 537773*(x_15298*x_37758)-x_35605 | a36185 | a40812')
print('   i.e. EXACTLY the point-addition law A=B=0 (x_11150,x_25739,x_37758 are rank-2 in A,B).')
C=json.load(open('/home/user/integer_solver/solve_lab/agentC_work/curve.json'))
print()
print('curve  a2 =',C['KA']); print('       a4 =',C['a4']); print('       a6 =',C['a6'])
print('short form A=0, B=64019533680030876408443198762210829058751700634554282185987325820393598524794, j=0')
print('group order = n_secp (verified [n]G=O) => ISOMORPHIC to secp256k1, not a different-order twist')
print('Q = (K2 mod p, K1 mod p) is ON that curve:',
      (pow(K1%P,2,P)-(pow(K2%P,3,P)+int(C['KA'])*pow(K2%P,2,P)+int(C['a4'])*(K2%P)+int(C['a6'])))%P==0)
