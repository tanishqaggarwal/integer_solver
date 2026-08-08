import os,sys,time,json
sys.path.insert(0,'.')
from lp_core import build_modmul_instance
# secp256k1 prime
p = 2**256 - 2**32 - 977
for red in ['naf']:
    t=time.time()
    inst = build_modmul_instance(p, mult='schoolbook', red=red)
    Q=inst['Q']; n=Q.n; s=inst['s']
    nsq=len(Q.squares); nand=len(Q.andcache)
    ncarry=sum(1 for v in range(n) if Q.kind[v] in ('adder','carry','chunk'))
    ninp=sum(1 for v in range(n) if Q.kind[v]=='input')
    print(json.dumps(dict(p='secp256k1',s=s,red=red,n=n,n_eq=nsq,n_and=nand,
        n_carry=ncarry,n_inputs=ninp,
        linear_reduction_frac=round(nsq/n,4),
        eff_dim=n-nsq, eff_dim_minus_inputs=n-nsq-ninp,
        secs=round(time.time()-t,1))))
