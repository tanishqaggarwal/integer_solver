import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H, full11 as F
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*E.NV
for k,x in d.items(): vd[int(k.split('_')[1])]=int(x)
r=F.solve_pair(24601,2081,verbose=False); n,ns,av=r
vp=E.forward(ns)
VS=[7075,8731,9118,9413,1329,10903,17325,22665,28961,28599,17499,29854,31864,642,2099,7068,28730,11150,21889,32405,25739,25156,3023,37758,2287,15663,15298]
print(f'{"var":>10} {"free?":>6} {"deliv_bits":>11} {"pair_bits":>10}  equal')
for u in VS:
    fr = E.definer[u] is None
    print(f'x_{u:<8} {"FREE" if fr else "def":>6} {abs(vd[u]).bit_length():>11} {abs(vp[u]).bit_length():>10}  {vd[u]==vp[u]}')
